use std::fs::File;
use std::io::{self, BufReader, Read, Write};
use std::path::Path;

use crate::ffi;

pub struct ByteWriter {
    pub buf: Vec<u8>,
}

impl ByteWriter {
    pub fn with_capacity(cap: usize) -> Self {
        Self { buf: Vec::with_capacity(cap) }
    }

    pub fn write_u8(&mut self, v: u8) {
        self.buf.push(v);
    }

    pub fn write_bytes(&mut self, s: &[u8]) {
        self.buf.extend_from_slice(s);
    }

    pub fn write_u32_le(&mut self, v: u32) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    pub fn write_u64_le(&mut self, v: u64) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    pub fn write_var_u64(&mut self, mut v: u64) {
        while v >= 0x80 {
            self.buf.push(((v & 0x7F) as u8) | 0x80);
            v >>= 7;
        }
        self.buf.push(v as u8);
    }
}

pub struct ByteReader<'a> {
    pub data: &'a [u8],
    pub pos: usize,
}

impl<'a> ByteReader<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0 }
    }

    pub fn has_remaining(&self) -> bool {
        self.pos < self.data.len()
    }

    pub fn read_u8(&mut self) -> Result<u8, String> {
        if self.pos >= self.data.len() {
            return Err("Unexpected EOF".into());
        }
        let b = self.data[self.pos];
        self.pos += 1;
        Ok(b)
    }

    pub fn read_bytes(&mut self, len: usize) -> Result<&'a [u8], String> {
        if self.pos + len > self.data.len() {
            return Err(format!("Unexpected EOF reading {len} bytes at {}", self.pos));
        }
        let sl = &self.data[self.pos..self.pos + len];
        self.pos += len;
        Ok(sl)
    }

    pub fn read_u32_le(&mut self) -> Result<u32, String> {
        let s = self.read_bytes(4)?;
        Ok(u32::from_le_bytes([s[0], s[1], s[2], s[3]]))
    }

    pub fn read_u64_le(&mut self) -> Result<u64, String> {
        let s = self.read_bytes(8)?;
        Ok(u64::from_le_bytes([s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]]))
    }

    pub fn read_var_u64(&mut self) -> Result<u64, String> {
        let mut res = 0u64;
        let mut shift = 0u32;
        loop {
            let b = self.read_u8()?;
            res |= ((b & 0x7F) as u64) << shift;
            if (b & 0x80) == 0 {
                return Ok(res);
            }
            shift += 7;
            if shift > 63 {
                return Err("Malformed varlong".into());
            }
        }
    }
}

// Java VarIntCodec Mode 3 & Mode 6 RangeCoder Decoder
fn range_aux_decode(src: &[u8], raw_len: usize, count_model: bool) -> Vec<u8> {
    const TOP: u32 = 1 << 24;
    const MOVE_BITS: u32 = 5;
    let mut probs = vec![1024u16; 1 << 16];
    let mut c0 = vec![0u16; if count_model { 1 << 16 } else { 0 }];
    let mut c1 = vec![0u16; if count_model { 1 << 16 } else { 0 }];

    let mut pos = 0usize;
    let mut next_b = || -> u32 {
        if pos < src.len() {
            let b = src[pos] as u32;
            pos += 1;
            b
        } else {
            0
        }
    };
    let mut code = 0u32;
    for _ in 0..5 {
        code = (code << 8) | next_b();
    }
    let mut range = 0xFFFF_FFFFu32;

    let mut decode_bit = |idx: usize| -> u32 {
        let prob = if !count_model {
            probs[idx] as u32
        } else {
            let z = c0[idx] as u64;
            let o = c1[idx] as u64;
            (((z + 1) << 11) / (z + o + 2)).clamp(1, 2047) as u32
        };
        let bound = (range >> 11) * prob;
        let bit = if code < bound {
            range = bound;
            0u32
        } else {
            code = code.wrapping_sub(bound);
            range = range.wrapping_sub(bound);
            1u32
        };
        if !count_model {
            let mut p = probs[idx] as i32;
            p += if bit == 0 { (2048 - p) >> MOVE_BITS } else { -(p >> MOVE_BITS) };
            probs[idx] = p as u16;
        } else {
            let mut z = c0[idx];
            let mut o = c1[idx];
            if bit == 0 { z += 1; } else { o += 1; }
            if z + o > 16384 {
                z = (z + 1) >> 1;
                o = (o + 1) >> 1;
            }
            c0[idx] = z;
            c1[idx] = o;
        }
        while range < TOP {
            range <<= 8;
            code = (code << 8) | next_b();
        }
        bit
    };

    let mut out = vec![0u8; raw_len];
    let mut prev = 0usize;
    for i in 0..raw_len {
        let base = prev << 8;
        let mut node = 1usize;
        for _ in 0..8 {
            let b = decode_bit(base | node) as usize;
            node = (node << 1) | b;
        }
        let sym = (node & 0xFF) as u8;
        out[i] = sym;
        prev = sym as usize;
    }
    out
}

// Java VarIntCodec Mode 5 LzRangeCoder Decoder
fn lz_range_decode(src: &[u8], raw_len: usize) -> Result<Vec<u8>, String> {
    const TOP: u32 = 1 << 24;
    const MOVE_BITS: u32 = 5;
    const MIN_MATCH: usize = 4;
    const OFF_IS_MATCH: usize = 0;
    const OFF_LITERAL: usize = 2;
    const OFF_LEN: usize = OFF_LITERAL + (256 << 8);
    const OFF_SLOT: usize = OFF_LEN + 1024;
    const OFF_MANT: usize = OFF_SLOT + (32 * 32);
    const MODEL_SIZE: usize = OFF_MANT + (32 * 32);

    let mut p = vec![1024i32; MODEL_SIZE];
    let mut pos = 0usize;
    let mut next_b = || -> u32 {
        if pos < src.len() {
            let b = src[pos] as u32;
            pos += 1;
            b
        } else {
            0
        }
    };
    let mut code = 0u32;
    for _ in 0..5 {
        code = (code << 8) | next_b();
    }
    let mut range = 0xFFFF_FFFFu32;

    let mut dec_bit = |idx: usize| -> u32 {
        let prob = p[idx] as u32;
        let bound = (range >> 11) * prob;
        let bit = if code < bound {
            range = bound;
            0u32
        } else {
            code = code.wrapping_sub(bound);
            range = range.wrapping_sub(bound);
            1u32
        };
        let v = p[idx];
        p[idx] = if bit == 0 { v + ((2048 - v) >> MOVE_BITS) } else { v - (v >> MOVE_BITS) };
        while range < TOP {
            range <<= 8;
            code = (code << 8) | next_b();
        }
        bit
    };

    let mut out = vec![0u8; raw_len];
    let mut i = 0usize;
    let mut prev_byte = 0usize;
    let mut prev_slot = 0usize;
    let mut last_match = false;
    while i < raw_len {
        let is_m = dec_bit(OFF_IS_MATCH + usize::from(last_match));
        if is_m == 0 {
            let base = OFF_LITERAL + (prev_byte << 8);
            let mut node = 1usize;
            for _ in 0..8 {
                let b = dec_bit(base + node) as usize;
                node = (node << 1) | b;
            }
            let sym = (node & 0xFF) as u8;
            out[i] = sym;
            i += 1;
            prev_byte = sym as usize;
            last_match = false;
        } else {
            let mut node = 1usize;
            for _ in 0..9 {
                let b = dec_bit(OFF_LEN + node) as usize;
                node = (node << 1) | b;
            }
            let mut len = (node & 0x1FF) + MIN_MATCH;
            let sctx = prev_slot * 32;
            node = 1;
            for _ in 0..5 {
                let b = dec_bit(OFF_SLOT + sctx + node) as usize;
                node = (node << 1) | b;
            }
            let slot = node & 0x1F;
            let mut mant = 0usize;
            for bit in (0..slot).rev() {
                let b = dec_bit(OFF_MANT + slot * 32 + (slot - 1 - bit)) as usize;
                mant |= b << bit;
            }
            let dist = (1usize << slot) + mant;
            if dist == 0 || dist > i {
                return Err(format!("Invalid LZ distance {dist} at {i}"));
            }
            if i + len > raw_len {
                len = raw_len - i;
            }
            for _ in 0..len {
                out[i] = out[i - dist];
                i += 1;
            }
            prev_byte = out[i - 1] as usize;
            prev_slot = slot;
            last_match = true;
        }
    }
    Ok(out)
}

pub fn write_compressed_aux(out: &mut ByteWriter, raw: &[u8]) {
    if raw.is_empty() {
        out.write_u8(0);
        return;
    }
    if raw.len() <= 12 {
        out.write_u8(1);
        out.write_var_u64(raw.len() as u64);
        out.write_bytes(raw);
        return;
    }
    if let Ok(comp) = ffi::deflate_raw(raw) {
        if comp.len() + 4 < raw.len() {
            out.write_u8(2);
            out.write_var_u64(raw.len() as u64);
            out.write_var_u64(comp.len() as u64);
            out.write_bytes(&comp);
            return;
        }
    }
    out.write_u8(1);
    out.write_var_u64(raw.len() as u64);
    out.write_bytes(raw);
}

pub fn read_compressed_aux(inp: &mut ByteReader<'_>) -> Result<Vec<u8>, String> {
    let mode = inp.read_u8()?;
    if mode == 0 {
        return Ok(Vec::new());
    }
    let raw_len = inp.read_var_u64()? as usize;
    match mode {
        1 => Ok(inp.read_bytes(raw_len)?.to_vec()),
        2 => {
            let clen = inp.read_var_u64()? as usize;
            let comp = inp.read_bytes(clen)?;
            let dec = ffi::inflate_raw(comp, raw_len).map_err(|e| e.to_string())?;
            if dec.len() != raw_len {
                return Err(format!("Deflate length mismatch: expected {raw_len}, got {}", dec.len()));
            }
            Ok(dec)
        }
        3 => {
            let clen = inp.read_var_u64()? as usize;
            let comp = inp.read_bytes(clen)?;
            Ok(range_aux_decode(comp, raw_len, false))
        }
        6 => {
            let clen = inp.read_var_u64()? as usize;
            let comp = inp.read_bytes(clen)?;
            Ok(range_aux_decode(comp, raw_len, true))
        }
        5 => {
            let clen = inp.read_var_u64()? as usize;
            let comp = inp.read_bytes(clen)?;
            lz_range_decode(comp, raw_len)
        }
        4 => {
            let values = inp.read_var_u64()? as usize;
            let planes = inp.read_u8()? as usize;
            let mut digits = Vec::with_capacity(planes);
            for k in 0..planes {
                let len = inp.read_var_u64()? as usize;
                let d = read_compressed_aux(inp)?;
                if d.len() != len {
                    return Err(format!("Plane {k} length mismatch"));
                }
                digits.push(d);
            }
            let mut out = ByteWriter::with_capacity(raw_len);
            let mut cursor = vec![0usize; planes];
            for _ in 0..values {
                let mut v = 0u64;
                let mut sh = 0u32;
                for k in 0..planes {
                    let sym = digits[k][cursor[k]];
                    cursor[k] += 1;
                    v |= ((sym & 0x7F) as u64) << sh;
                    sh += 7;
                    if (sym & 0x80) == 0 {
                        break;
                    }
                }
                out.write_var_u64(v);
            }
            if out.buf.len() != raw_len {
                return Err(format!("Plane aux decode mismatch: {} vs {raw_len}", out.buf.len()));
            }
            Ok(out.buf)
        }
        _ => Err(format!("Unsupported compressed aux mode {mode}")),
    }
}

#[derive(Clone, Debug)]
pub struct BlockSlice64 {
    pub start_byte: u64,
    pub end_byte: u64,
    pub first_record_idx: u64,
    pub record_count: u64,
}

pub struct FilePartition {
    pub slices: Vec<BlockSlice64>,
    pub size: u64,
    pub crc32: u32,
}

/// Single-pass 64-bit streaming FASTA partitioner (works on any >2 GiB file in ~1 MiB RAM).
pub fn partition_file(path: &Path, target_block_bytes: u64, split_tail_block: bool) -> io::Result<FilePartition> {
    let f = File::open(path)?;
    let size = f.metadata()?.len();
    if size == 0 {
        return Ok(FilePartition {
            slices: vec![BlockSlice64 {
                start_byte: 0,
                end_byte: 0,
                first_record_idx: 0,
                record_count: 1,
            }],
            size: 0,
            crc32: 0,
        });
    }

    let mut reader = BufReader::with_capacity(1 << 20, f);
    let mut slices = Vec::new();
    let mut crc = 0u32;

    let mut io_buf = vec![0u8; 1 << 20];
    let mut io_pos = 0usize;
    let mut io_lim = 0usize;

    let mut read_next = |crc_acc: &mut u32| -> io::Result<Option<u8>> {
        if io_pos >= io_lim {
            io_lim = reader.read(&mut io_buf)?;
            if io_lim == 0 {
                return Ok(None);
            }
            *crc_acc = ffi::crc32_ieee(*crc_acc, &io_buf[..io_lim]);
            io_pos = 0;
        }
        let b = io_buf[io_pos];
        io_pos += 1;
        Ok(Some(b))
    };

    let Some(first_b) = read_next(&mut crc)? else {
        return Ok(FilePartition { slices: vec![], size: 0, crc32: 0 });
    };
    let starts_with_hdr = first_b == b'>';
    let mut block_start = 0u64;
    let mut global_rec: i64 = if starts_with_hdr { -1 } else { 0 };
    let mut block_first_rec = 0u64;
    let mut block_rec_cnt = if starts_with_hdr { 0u64 } else { 1u64 };
    let mut offset = 0u64;
    let mut cur_b = first_b;

    loop {
        let line_start = offset;
        let line_first_b = cur_b;
        let mut line_end: i64 = -1;
        let mut last_content_off: i64 = -1;
        let mut last_content_b = 0u8;

        loop {
            if cur_b == b'\n' {
                line_end = offset as i64;
                break;
            }
            last_content_off = offset as i64;
            last_content_b = cur_b;
            match read_next(&mut crc)? {
                None => break,
                Some(nb) => {
                    offset += 1;
                    cur_b = nb;
                }
            }
        }
        let eof = line_end < 0;
        let content_end = if eof {
            size
        } else if last_content_off == line_end - 1 && last_content_b == b'\r' {
            (line_end - 1) as u64
        } else {
            line_end as u64
        };
        let is_hdr = content_end > line_start && line_first_b == b'>';

        if is_hdr {
            if line_start > block_start && (line_start - block_start) >= target_block_bytes && block_rec_cnt > 0 {
                slices.push(BlockSlice64 {
                    start_byte: block_start,
                    end_byte: line_start,
                    first_record_idx: block_first_rec,
                    record_count: block_rec_cnt,
                });
                block_start = line_start;
                block_first_rec = (global_rec + 1) as u64;
                block_rec_cnt = 0;
            }
            global_rec += 1;
            if block_rec_cnt == 0 {
                block_first_rec = global_rec as u64;
            }
            block_rec_cnt += 1;
        } else {
            let line_chars = content_end - line_start;
            if line_start > block_start && (line_start - block_start) >= target_block_bytes * 2 {
                slices.push(BlockSlice64 {
                    start_byte: block_start,
                    end_byte: line_start,
                    first_record_idx: block_first_rec,
                    record_count: block_rec_cnt.max(1),
                });
                block_start = line_start;
                block_first_rec = global_rec.max(0) as u64;
                block_rec_cnt = 1;
            } else if line_chars > target_block_bytes * 2 {
                if line_start > block_start {
                    slices.push(BlockSlice64 {
                        start_byte: block_start,
                        end_byte: line_start,
                        first_record_idx: block_first_rec,
                        record_count: block_rec_cnt.max(1),
                    });
                    block_start = line_start;
                    block_first_rec = global_rec.max(0) as u64;
                    block_rec_cnt = 1;
                }
                let next_p = if eof { size } else { (line_end + 1) as u64 };
                let cut_thr = if split_tail_block { target_block_bytes } else { target_block_bytes * 2 };
                while next_p - block_start > cut_thr {
                    let cut = block_start + target_block_bytes;
                    slices.push(BlockSlice64 {
                        start_byte: block_start,
                        end_byte: cut,
                        first_record_idx: block_first_rec,
                        record_count: 1,
                    });
                    block_start = cut;
                }
                if eof { break; }
                offset = (line_end + 1) as u64;
                match read_next(&mut crc)? {
                    None => break,
                    Some(nb) => { cur_b = nb; continue; }
                }
            }
        }

        if eof { break; }
        offset = (line_end + 1) as u64;
        match read_next(&mut crc)? {
            None => break,
            Some(nb) => cur_b = nb,
        }
    }

    if block_start < size {
        slices.push(BlockSlice64 {
            start_byte: block_start,
            end_byte: size,
            first_record_idx: global_rec.min(block_first_rec as i64).max(0) as u64,
            record_count: block_rec_cnt.max(1),
        });
    }

    Ok(FilePartition { slices, size, crc32: crc })
}

pub struct FactorizedChunk {
    pub rna_mode: bool,
    pub num_headers: usize,
    pub header_seq_line_deltas: Vec<u64>,
    pub header_bytes: Vec<Vec<u8>>,
    pub num_seq_lines: usize,
    pub seq_line_lengths: Vec<u64>,
    pub line_eol_types: Vec<u8>,
    pub total_seq_chars: usize,
    pub lower_upper_gaps: Vec<u64>,
    pub lower_lengths: Vec<u64>,
    pub n_acgt_gaps: Vec<u64>,
    pub n_run_lengths: Vec<u64>,
    pub exc_gaps: Vec<u64>,
    pub exc_symbols: Vec<u8>,
    pub dna_bases: Vec<u8>,
}

impl FactorizedChunk {
    pub fn record_count(&self) -> usize {
        let leading_lines = if self.num_headers > 0 { self.header_seq_line_deltas[0] as usize } else { self.num_seq_lines };
        let has_leading = self.num_headers == 0 || leading_lines > 0;
        self.num_headers + usize::from(has_leading)
    }
}

pub fn factorize(raw: &[u8]) -> FactorizedChunk {
    if raw.is_empty() {
        return FactorizedChunk {
            rna_mode: false,
            num_headers: 0,
            header_seq_line_deltas: Vec::new(),
            header_bytes: Vec::new(),
            num_seq_lines: 0,
            seq_line_lengths: Vec::new(),
            line_eol_types: Vec::new(),
            total_seq_chars: 0,
            lower_upper_gaps: Vec::new(),
            lower_lengths: Vec::new(),
            n_acgt_gaps: Vec::new(),
            n_run_lengths: Vec::new(),
            exc_gaps: Vec::new(),
            exc_symbols: Vec::new(),
            dna_bases: Vec::new(),
        };
    }

    let mut seq_starts = Vec::with_capacity(raw.len() / 60 + 4);
    let mut seq_lens = Vec::with_capacity(raw.len() / 60 + 4);
    let mut eol_types = Vec::with_capacity(raw.len() / 60 + 4);
    let mut hdr_deltas = Vec::new();
    let mut hdr_bytes = Vec::new();
    let mut seq_since_hdr = 0u64;
    let mut has_t = false;
    let mut has_u = false;
    let mut total_seq_chars = 0usize;

    let mut p = 0usize;
    while p < raw.len() {
        let line_start = p;
        let nl = raw[p..].iter().position(|&b| b == b'\n').map_or(raw.len(), |i| p + i);
        let (eol, content_end) = if nl < raw.len() {
            if nl > line_start && raw[nl - 1] == b'\r' {
                p = nl + 1;
                (2u8, nl - 1)
            } else {
                p = nl + 1;
                (1u8, nl)
            }
        } else {
            p = raw.len();
            (0u8, raw.len())
        };
        eol_types.push(eol);

        if content_end > line_start && raw[line_start] == b'>' {
            hdr_deltas.push(seq_since_hdr);
            hdr_bytes.push(raw[line_start + 1..content_end].to_vec());
            seq_since_hdr = 0;
        } else {
            let slen = content_end - line_start;
            seq_starts.push(line_start);
            seq_lens.push(slen as u64);
            seq_since_hdr += 1;
            total_seq_chars += slen;
            for &b in &raw[line_start..content_end] {
                if b == b'T' || b == b't' { has_t = true; }
                else if b == b'U' || b == b'u' { has_u = true; }
            }
        }
    }

    let rna_mode = has_u && !has_t;
    let fourth = if rna_mode { b'U' } else { b'T' };

    let mut dna = Vec::with_capacity(total_seq_chars);
    let mut low_gaps = Vec::new();
    let mut low_lens = Vec::new();
    let mut in_lower = false;
    let mut cur_case_run = 0u64;
    let mut pending_up_gap = 0u64;

    let mut n_gaps = Vec::new();
    let mut n_lens = Vec::new();
    let mut cur_acgt_since_n = 0u64;
    let mut cur_n_run = 0u64;

    let mut exc_gaps = Vec::new();
    let mut exc_syms = Vec::new();
    let mut cur_acgtn_since_exc = 0u64;

    for (&s, &sl) in seq_starts.iter().zip(seq_lens.iter()) {
        let e = s + sl as usize;
        for &orig in &raw[s..e] {
            let mut b = orig;
            if b >= b'a' && b <= b'z' {
                if !in_lower {
                    pending_up_gap = cur_case_run;
                    cur_case_run = 1;
                    in_lower = true;
                } else {
                    cur_case_run += 1;
                }
                b -= 32;
            } else {
                if in_lower {
                    low_gaps.push(pending_up_gap);
                    low_lens.push(cur_case_run);
                    cur_case_run = 1;
                    in_lower = false;
                } else {
                    cur_case_run += 1;
                }
            }

            let code = if b == b'A' { 0 } else if b == b'C' { 1 } else if b == b'G' { 2 } else if b == fourth { 3 } else { -1 };
            if code >= 0 {
                if cur_n_run > 0 {
                    n_gaps.push(cur_acgt_since_n);
                    n_lens.push(cur_n_run);
                    cur_n_run = 0;
                    cur_acgt_since_n = 0;
                }
                dna.push(code as u8);
                cur_acgt_since_n += 1;
                cur_acgtn_since_exc += 1;
            } else if b == b'N' {
                cur_n_run += 1;
                cur_acgtn_since_exc += 1;
            } else {
                exc_gaps.push(cur_acgtn_since_exc);
                exc_syms.push(b);
                cur_acgtn_since_exc = 0;
            }
        }
    }
    if in_lower {
        low_gaps.push(pending_up_gap);
        low_lens.push(cur_case_run);
    }
    if cur_n_run > 0 {
        n_gaps.push(cur_acgt_since_n);
        n_lens.push(cur_n_run);
    }

    FactorizedChunk {
        rna_mode,
        num_headers: hdr_bytes.len(),
        header_seq_line_deltas: hdr_deltas,
        header_bytes: hdr_bytes,
        num_seq_lines: seq_lens.len(),
        seq_line_lengths: seq_lens,
        line_eol_types: eol_types,
        total_seq_chars,
        lower_upper_gaps: low_gaps,
        lower_lengths: low_lens,
        n_acgt_gaps: n_gaps,
        n_run_lengths: n_lens,
        exc_gaps,
        exc_symbols: exc_syms,
        dna_bases: dna,
    }
}

fn write_rle_u64(out: &mut ByteWriter, arr: &[u64]) {
    let mut i = 0;
    while i < arr.len() {
        let v = arr[i];
        let mut run = 1usize;
        while i + run < arr.len() && arr[i + run] == v { run += 1; }
        out.write_var_u64(run as u64);
        out.write_var_u64(v);
        i += run;
    }
}

fn read_rle_u64(inp: &mut ByteReader<'_>, expected: usize) -> Result<Vec<u64>, String> {
    let mut out = Vec::with_capacity(expected);
    while inp.has_remaining() && out.len() < expected {
        let run = inp.read_var_u64()? as usize;
        let val = inp.read_var_u64()?;
        out.resize(out.len() + run, val);
    }
    if out.len() != expected {
        return Err(format!("RLE int mismatch: expected {expected}, got {}", out.len()));
    }
    Ok(out)
}

fn write_rle_u8(out: &mut ByteWriter, arr: &[u8]) {
    let mut i = 0;
    while i < arr.len() {
        let v = arr[i];
        let mut run = 1usize;
        while i + run < arr.len() && arr[i + run] == v { run += 1; }
        out.write_var_u64(run as u64);
        out.write_var_u64(v as u64);
        i += run;
    }
}

fn read_rle_u8(inp: &mut ByteReader<'_>, expected: usize) -> Result<Vec<u8>, String> {
    let mut out = Vec::with_capacity(expected);
    while inp.has_remaining() && out.len() < expected {
        let run = inp.read_var_u64()? as usize;
        let val = inp.read_var_u64()? as u8;
        out.resize(out.len() + run, val);
    }
    if out.len() != expected {
        return Err(format!("RLE byte mismatch: expected {expected}, got {}", out.len()));
    }
    Ok(out)
}

pub fn write_aux_metadata(out: &mut ByteWriter, c: &FactorizedChunk) {
    let mut flags = 0u8;
    if c.rna_mode { flags |= 1; }
    if c.num_headers > 0 { flags |= 2; }
    if !c.lower_upper_gaps.is_empty() { flags |= 4; }
    if !c.n_acgt_gaps.is_empty() { flags |= 8; }
    if !c.exc_gaps.is_empty() { flags |= 16; }

    out.write_u8(flags);
    out.write_var_u64(c.dna_bases.len() as u64);
    out.write_var_u64(c.total_seq_chars as u64);
    out.write_var_u64(c.num_seq_lines as u64);

    if c.num_headers > 0 {
        out.write_var_u64(c.num_headers as u64);
        let mut hbuf = ByteWriter::with_capacity(c.num_headers * 24);
        for &d in &c.header_seq_line_deltas { hbuf.write_var_u64(d); }
        let mut prev: &[u8] = &[];
        for cur in &c.header_bytes {
            let max_s = prev.len().min(cur.len());
            let mut sh = 0usize;
            while sh < max_s && prev[sh] == cur[sh] { sh += 1; }
            hbuf.write_var_u64(sh as u64);
            hbuf.write_var_u64((cur.len() - sh) as u64);
            hbuf.write_bytes(&cur[sh..]);
            prev = cur;
        }
        write_compressed_aux(out, &hbuf.buf);
    }

    let mut lbuf = ByteWriter::with_capacity(32);
    write_rle_u64(&mut lbuf, &c.seq_line_lengths);
    write_compressed_aux(out, &lbuf.buf);

    let mut ebuf = ByteWriter::with_capacity(16);
    write_rle_u8(&mut ebuf, &c.line_eol_types);
    write_compressed_aux(out, &ebuf.buf);

    if !c.lower_upper_gaps.is_empty() {
        out.write_var_u64(c.lower_upper_gaps.len() as u64);
        let min_low = *c.lower_lengths.iter().min().unwrap();
        out.write_var_u64(min_low);
        let mut gcol = ByteWriter::with_capacity(c.lower_upper_gaps.len() * 2);
        let mut lcol = ByteWriter::with_capacity(c.lower_lengths.len() * 2);
        for (&g, &l) in c.lower_upper_gaps.iter().zip(c.lower_lengths.iter()) {
            gcol.write_var_u64(g);
            lcol.write_var_u64(l - min_low);
        }
        write_compressed_aux(out, &gcol.buf);
        write_compressed_aux(out, &lcol.buf);
    }

    if !c.n_acgt_gaps.is_empty() {
        out.write_var_u64(c.n_acgt_gaps.len() as u64);
        let mut nbuf = ByteWriter::with_capacity(c.n_acgt_gaps.len() * 4);
        for (&g, &l) in c.n_acgt_gaps.iter().zip(c.n_run_lengths.iter()) {
            nbuf.write_var_u64(g);
            nbuf.write_var_u64(l);
        }
        write_compressed_aux(out, &nbuf.buf);
    }

    if !c.exc_gaps.is_empty() {
        out.write_var_u64(c.exc_gaps.len() as u64);
        let mut gbuf = ByteWriter::with_capacity(c.exc_gaps.len() * 2);
        for &g in &c.exc_gaps { gbuf.write_var_u64(g); }
        write_compressed_aux(out, &gbuf.buf);
        write_compressed_aux(out, &c.exc_symbols);
    }
}

pub fn read_aux_metadata(inp: &mut ByteReader<'_>, expected_dna_len_out: &mut usize) -> Result<FactorizedChunk, String> {
    let flags = inp.read_u8()?;
    let rna_mode = (flags & 1) != 0;
    let has_hdr = (flags & 2) != 0;
    let has_low = (flags & 4) != 0;
    let has_n = (flags & 8) != 0;
    let has_exc = (flags & 16) != 0;

    let dna_len = inp.read_var_u64()? as usize;
    *expected_dna_len_out = dna_len;
    let total_seq_chars = inp.read_var_u64()? as usize;
    let num_seq_lines = inp.read_var_u64()? as usize;

    let mut header_deltas = Vec::new();
    let mut header_bytes: Vec<Vec<u8>> = Vec::new();
    let num_headers = if has_hdr {
        let nh = inp.read_var_u64()? as usize;
        let hdr_raw = read_compressed_aux(inp)?;
        let mut hr = ByteReader::new(&hdr_raw);
        for _ in 0..nh { header_deltas.push(hr.read_var_u64()?); }
        let mut prev = Vec::new();
        for _ in 0..nh {
            let sh = hr.read_var_u64()? as usize;
            let sfx = hr.read_var_u64()? as usize;
            let mut cur = Vec::with_capacity(sh + sfx);
            cur.extend_from_slice(&prev[..sh]);
            cur.extend_from_slice(hr.read_bytes(sfx)?);
            prev = cur.clone();
            header_bytes.push(cur);
        }
        nh
    } else {
        0
    };

    let ll_raw = read_compressed_aux(inp)?;
    let seq_line_lengths = read_rle_u64(&mut ByteReader::new(&ll_raw), num_seq_lines)?;

    let eol_raw = read_compressed_aux(inp)?;
    let line_eol_types = read_rle_u8(&mut ByteReader::new(&eol_raw), num_headers + num_seq_lines)?;

    let (mut lower_upper_gaps, mut lower_lengths) = (Vec::new(), Vec::new());
    if has_low {
        let nl = inp.read_var_u64()? as usize;
        let min_low = inp.read_var_u64()?;
        let graw = read_compressed_aux(inp)?;
        let lraw = read_compressed_aux(inp)?;
        let (mut gr, mut lr) = (ByteReader::new(&graw), ByteReader::new(&lraw));
        for _ in 0..nl {
            lower_upper_gaps.push(gr.read_var_u64()?);
            lower_lengths.push(lr.read_var_u64()? + min_low);
        }
    }

    let (mut n_acgt_gaps, mut n_run_lengths) = (Vec::new(), Vec::new());
    if has_n {
        let nn = inp.read_var_u64()? as usize;
        let nraw = read_compressed_aux(inp)?;
        let mut nr = ByteReader::new(&nraw);
        for _ in 0..nn {
            n_acgt_gaps.push(nr.read_var_u64()?);
            n_run_lengths.push(nr.read_var_u64()?);
        }
    }

    let (mut exc_gaps, mut exc_symbols) = (Vec::new(), Vec::new());
    if has_exc {
        let ne = inp.read_var_u64()? as usize;
        let graw = read_compressed_aux(inp)?;
        let mut gr = ByteReader::new(&graw);
        for _ in 0..ne { exc_gaps.push(gr.read_var_u64()?); }
        exc_symbols = read_compressed_aux(inp)?;
    }

    Ok(FactorizedChunk {
        rna_mode,
        num_headers,
        header_seq_line_deltas: header_deltas,
        header_bytes,
        num_seq_lines,
        seq_line_lengths,
        line_eol_types,
        total_seq_chars,
        lower_upper_gaps,
        lower_lengths,
        n_acgt_gaps,
        n_run_lengths,
        exc_gaps,
        exc_symbols,
        dna_bases: Vec::new(),
    })
}

pub fn reconstruct_to(c: &FactorizedChunk, dna_bases: &[u8], out: &mut dyn Write, crc: &mut u32) -> io::Result<u64> {
    if c.num_headers == 0 && c.num_seq_lines == 0 {
        return Ok(0);
    }
    let alpha = if c.rna_mode { [b'A', b'C', b'G', b'U'] } else { [b'A', b'C', b'G', b'T'] };
    let acgtn_len = c.total_seq_chars - c.exc_gaps.len();

    let mut acgtn = vec![0u8; acgtn_len];
    if c.n_acgt_gaps.is_empty() {
        for (i, &b) in dna_bases.iter().enumerate() {
            acgtn[i] = alpha[(b & 3) as usize];
        }
    } else {
        let mut sp = 0usize;
        let mut dp = 0usize;
        for (&gap, &nlen) in c.n_acgt_gaps.iter().zip(c.n_run_lengths.iter()) {
            for _ in 0..gap {
                acgtn[dp] = alpha[(dna_bases[sp] & 3) as usize];
                dp += 1;
                sp += 1;
            }
            for _ in 0..nlen {
                acgtn[dp] = b'N';
                dp += 1;
            }
        }
        while sp < dna_bases.len() {
            acgtn[dp] = alpha[(dna_bases[sp] & 3) as usize];
            dp += 1;
            sp += 1;
        }
    }

    let mut seq_chars = if c.exc_gaps.is_empty() {
        acgtn
    } else {
        let mut sc = vec![0u8; c.total_seq_chars];
        let mut sp = 0usize;
        let mut dp = 0usize;
        for (&gap, &sym) in c.exc_gaps.iter().zip(c.exc_symbols.iter()) {
            let g = gap as usize;
            if g > 0 {
                sc[dp..dp + g].copy_from_slice(&acgtn[sp..sp + g]);
                sp += g;
                dp += g;
            }
            sc[dp] = sym;
            dp += 1;
        }
        let rem = acgtn.len() - sp;
        if rem > 0 {
            sc[dp..dp + rem].copy_from_slice(&acgtn[sp..]);
        }
        sc
    };

    if !c.lower_upper_gaps.is_empty() {
        let mut pos = 0usize;
        for (&gap, &len) in c.lower_upper_gaps.iter().zip(c.lower_lengths.iter()) {
            pos += gap as usize;
            for k in 0..len as usize {
                seq_chars[pos + k] = seq_chars[pos + k].wrapping_add(32);
            }
            pos += len as usize;
        }
    }

    let mut io_buf = Vec::with_capacity(65536);
    let mut written = 0u64;
    fn flush_buf(out: &mut dyn Write, buf: &mut Vec<u8>, crc_val: &mut u32, written: &mut u64) -> io::Result<()> {
        if !buf.is_empty() {
            out.write_all(buf)?;
            *crc_val = ffi::crc32_ieee(*crc_val, buf);
            *written += buf.len() as u64;
            buf.clear();
        }
        Ok(())
    }

    let mut h = 0usize;
    let mut sl = 0usize;
    let mut seq_since_hdr = 0u64;
    let mut line_idx = 0usize;
    let mut seq_char_pos = 0usize;

    while h < c.num_headers || sl < c.num_seq_lines {
        while h < c.num_headers && c.header_seq_line_deltas[h] == seq_since_hdr {
            let hb = &c.header_bytes[h];
            if io_buf.len() + 3 + hb.len() > 65536 {
                flush_buf(out, &mut io_buf, crc, &mut written)?;
            }
            io_buf.push(b'>');
            if hb.len() + 2 <= 65536 - io_buf.len() {
                io_buf.extend_from_slice(hb);
            } else {
                flush_buf(out, &mut io_buf, crc, &mut written)?;
                out.write_all(hb)?;
                *crc = ffi::crc32_ieee(*crc, hb);
                written += hb.len() as u64;
            }
            let eol = c.line_eol_types[line_idx];
            line_idx += 1;
            if eol == 2 { io_buf.extend_from_slice(b"\r\n"); }
            else if eol == 1 { io_buf.push(b'\n'); }
            h += 1;
            seq_since_hdr = 0;
        }
        if sl < c.num_seq_lines {
            let l = c.seq_line_lengths[sl] as usize;
            if io_buf.len() + l + 2 > 65536 {
                flush_buf(out, &mut io_buf, crc, &mut written)?;
            }
            if l > 0 {
                if l + 2 <= 65536 - io_buf.len() {
                    io_buf.extend_from_slice(&seq_chars[seq_char_pos..seq_char_pos + l]);
                } else {
                    flush_buf(out, &mut io_buf, crc, &mut written)?;
                    out.write_all(&seq_chars[seq_char_pos..seq_char_pos + l])?;
                    *crc = ffi::crc32_ieee(*crc, &seq_chars[seq_char_pos..seq_char_pos + l]);
                    written += l as u64;
                }
                seq_char_pos += l;
            }
            let eol = c.line_eol_types[line_idx];
            line_idx += 1;
            if eol == 2 { io_buf.extend_from_slice(b"\r\n"); }
            else if eol == 1 { io_buf.push(b'\n'); }
            sl += 1;
            seq_since_hdr += 1;
        }
    }
    flush_buf(out, &mut io_buf, crc, &mut written)?;
    Ok(written)
}
