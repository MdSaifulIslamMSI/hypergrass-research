use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use crate::fasta::{ByteReader, ByteWriter};
use crate::ffi::{self, SyncmerIndex};
use crate::hgs2::{
    self, CompressOptions, FLAG3_FASTQ, FLAG3_HAS_CDS, MAGIC3, MODE_REF,
    MODE_SINGLE, StreamingSha256, compress_panel_names, decompress_panel_names,
    load_panel_for_compress, verify_and_load_panel_for_decode,
};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum OmicInputKind {
    Auto,
    Fasta,
    Fastq,
}

#[derive(Clone)]
pub struct Hgs3CompressOptions {
    pub format: OmicInputKind,
    pub base: CompressOptions,
    pub emit_cds_digest: bool,
}

impl Default for Hgs3CompressOptions {
    fn default() -> Self {
        Self {
            format: OmicInputKind::Auto,
            base: CompressOptions::default(),
            emit_cds_digest: false,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct CdsSummary {
    pub orf_count: u64,
    pub longest_orf_aa: u32,
    pub total_orf_aa: u64,
    pub gc_permille: u16,
    pub proteome_sha256: [u8; 32],
}

pub struct Hgs3Report {
    pub is_fastq: bool,
    pub mode: u8,
    pub raw_bytes: u64,
    pub archive_bytes: u64,
    pub total_records: u64,
    pub dna_bases: u64,
    pub qual_bytes: u64,
    pub meta_comp_bytes: usize,
    pub mosaic_comp_bytes: usize,
    pub novelty_comp_bytes: usize,
    pub qual_comp_bytes: usize,
    pub cds: Option<CdsSummary>,
    pub seconds: f64,
}

impl Hgs3Report {
    pub fn bpb(&self) -> f64 {
        if self.raw_bytes == 0 {
            0.0
        } else {
            (self.archive_bytes as f64 * 8.0) / (self.raw_bytes as f64)
        }
    }

    pub fn ratio(&self) -> f64 {
        if self.archive_bytes == 0 {
            0.0
        } else {
            (self.raw_bytes as f64) / (self.archive_bytes as f64)
        }
    }
}

pub(crate) fn build_mismatch_mask(len: usize, variants: &[ffi::HgsVariantRecord]) -> Vec<u8> {
    let mut mask = vec![0u8; len];
    for v in variants {
        let start = v.tgt_pos as usize;
        if start >= len {
            continue;
        }
        let span = match v.kind {
            0 => 1usize,
            1 | 3 | 4 => (v.alt_len as usize).max(1),
            _ => 0usize,
        };
        let end = start.saturating_add(span).min(len);
        for b in &mut mask[start..end] {
            *b = 1;
        }
    }
    mask
}

const CODON_TABLE: [u8; 64] = *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

pub fn scan_cds_summary(two_bit: &[u8]) -> Option<CdsSummary> {
    if two_bit.is_empty() {
        return None;
    }
    let mut gc = 0u64;
    for &b in two_bit {
        if b == 1 || b == 2 {
            gc += 1;
        }
    }
    let gc_permille = ((gc * 1000) / (two_bit.len() as u64)).min(1000) as u16;
    let n = two_bit.len();

    let mut orf_count = 0u64;
    let mut longest_orf_aa = 0u32;
    let mut total_orf_aa = 0u64;
    let mut sha = StreamingSha256::new();
    let mut pep_buf: Vec<u8> = Vec::with_capacity(256);

    let mut scan_strand = |rev: bool| {
        for frame in 0..3usize {
            let mut pos = frame;
            let mut in_orf = false;
            pep_buf.clear();
            while pos + 3 <= n {
                let (b0, b1, b2) = if !rev {
                    (two_bit[pos] & 3, two_bit[pos + 1] & 3, two_bit[pos + 2] & 3)
                } else {
                    let i = n - 1 - pos;
                    (
                        (two_bit[i] ^ 3) & 3,
                        (two_bit[i - 1] ^ 3) & 3,
                        (two_bit[i - 2] ^ 3) & 3,
                    )
                };
                let idx = ((b0 as usize) << 4) | ((b1 as usize) << 2) | (b2 as usize);
                let aa = CODON_TABLE[idx];
                if !in_orf {
                    if aa == b'M' {
                        in_orf = true;
                        pep_buf.clear();
                        pep_buf.push(b'M');
                    }
                } else if aa == b'*' {
                    if pep_buf.len() >= 30 {
                        let plen = pep_buf.len() as u32;
                        orf_count += 1;
                        total_orf_aa += plen as u64;
                        if plen > longest_orf_aa {
                            longest_orf_aa = plen;
                        }
                        sha.update(&pep_buf);
                        sha.update(b"*");
                    }
                    in_orf = false;
                } else {
                    pep_buf.push(aa);
                }
                pos += 3;
            }
        }
    };

    scan_strand(false);
    scan_strand(true);

    Some(CdsSummary {
        orf_count,
        longest_orf_aa,
        total_orf_aa,
        gc_permille,
        proteome_sha256: sha.finish(),
    })
}

fn serialize_cds_summary(s: &CdsSummary) -> Result<Vec<u8>, String> {
    let mut w = ByteWriter::with_capacity(64);
    w.write_var_u64(s.orf_count);
    w.write_var_u64(s.longest_orf_aa as u64);
    w.write_var_u64(s.total_orf_aa);
    w.write_var_u64(s.gc_permille as u64);
    w.write_bytes(&s.proteome_sha256);
    ffi::zlib_compress(&w.buf).map_err(|e| e.to_string())
}

pub(crate) fn deserialize_cds_summary(comp: &[u8]) -> Result<CdsSummary, String> {
    let raw = ffi::zlib_decompress(comp).map_err(|e| e.to_string())?;
    let mut r = ByteReader::new(&raw);
    let orf_count = r.read_var_u64()?;
    let longest_orf_aa = r.read_var_u64()? as u32;
    let total_orf_aa = r.read_var_u64()?;
    let gc_permille = r.read_var_u64()? as u16;
    let mut proteome_sha256 = [0u8; 32];
    proteome_sha256.copy_from_slice(r.read_bytes(32)?);
    Ok(CdsSummary {
        orf_count,
        longest_orf_aa,
        total_orf_aa,
        gc_permille,
        proteome_sha256,
    })
}

fn scan_line(raw: &[u8], pos: &mut usize) -> (usize, usize, u8) {
    let start = *pos;
    let n = raw.len();
    while *pos < n {
        if raw[*pos] == b'\n' {
            let end = *pos;
            *pos += 1;
            return (start, end, 0);
        }
        if raw[*pos] == b'\r' && *pos + 1 < n && raw[*pos + 1] == b'\n' {
            let end = *pos;
            *pos += 2;
            return (start, end, 1);
        }
        *pos += 1;
    }
    (start, n, 2)
}

fn eol_bytes(code: u8) -> &'static [u8] {
    match code {
        0 => b"\n",
        1 => b"\r\n",
        _ => b"",
    }
}

pub fn detect_input_kind(path: &Path, explicit: OmicInputKind) -> Result<OmicInputKind, String> {
    if explicit != OmicInputKind::Auto {
        return Ok(explicit);
    }
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    if name.ends_with(".fq") || name.ends_with(".fastq") {
        return Ok(OmicInputKind::Fastq);
    }
    let mut f = File::open(path).map_err(|e| e.to_string())?;
    let mut probe = [0u8; 64];
    use std::io::Read;
    let n = f.read(&mut probe).map_err(|e| e.to_string())?;
    for &b in &probe[..n] {
        if b.is_ascii_whitespace() {
            continue;
        }
        return Ok(if b == b'@' {
            OmicInputKind::Fastq
        } else {
            OmicInputKind::Fasta
        });
    }
    Ok(OmicInputKind::Fasta)
}

pub(crate) struct ParsedHgs3Fastq<'a> {
    pub mode: u8,
    pub raw_bytes: u64,
    pub expected_crc: u32,
    pub total_records: u64,
    pub panel_names: Vec<String>,
    pub panel_sha: [u8; 32],
    pub total_dna_bases: u64,
    pub total_qual_bytes: u64,
    pub read_len_hint: u32,
    pub meta_comp: &'a [u8],
    pub mosaic_bytes: &'a [u8],
    pub nov_len: u64,
    pub ans_nov_comp: &'a [u8],
    pub qual_comp: &'a [u8],
    pub cds_comp: &'a [u8],
}

pub(crate) fn parse_hgs3_fastq<'a>(
    bytes: &'a [u8],
    decode_panel: bool,
) -> Result<ParsedHgs3Fastq<'a>, String> {
    let mut r = ByteReader::new(bytes);
    if r.read_bytes(4)? != MAGIC3 {
        return Err("Not an HGS3 archive".into());
    }
    let flags = r.read_u8()?;
    let mode = flags & 3;
    if (flags & FLAG3_FASTQ) == 0 {
        return Err("HGS3 archive does not carry FASTQ payload".into());
    }
    let has_cds = (flags & FLAG3_HAS_CDS) != 0;
    let raw_bytes = r.read_u64_le()?;
    let expected_crc = r.read_u32_le()?;
    let total_records = r.read_var_u64()?;
    let num_units = r.read_var_u64()?;

    let mut panel_names = Vec::new();
    let mut panel_sha = [0u8; 32];
    if mode == MODE_REF {
        let pcnt = r.read_var_u64()? as usize;
        let cn_len = r.read_var_u64()? as usize;
        let cn = r.read_bytes(cn_len)?;
        if decode_panel {
            panel_names = decompress_panel_names(cn, pcnt)?;
        }
        panel_sha.copy_from_slice(r.read_bytes(32)?);
    }

    if raw_bytes == 0 || num_units == 0 || total_records == 0 {
        return Ok(ParsedHgs3Fastq {
            mode,
            raw_bytes,
            expected_crc,
            total_records,
            panel_names,
            panel_sha,
            total_dna_bases: 0,
            total_qual_bytes: 0,
            read_len_hint: 1,
            meta_comp: &[],
            mosaic_bytes: &[],
            nov_len: 0,
            ans_nov_comp: &[],
            qual_comp: &[],
            cds_comp: &[],
        });
    }

    let total_dna_bases = r.read_var_u64()?;
    let total_qual_bytes = r.read_var_u64()?;
    let read_len_hint = (r.read_var_u64()? as u32).max(1);

    let meta_len = r.read_var_u64()? as usize;
    let meta_comp = r.read_bytes(meta_len)?;
    let mos_len = r.read_var_u64()? as usize;
    let mosaic_bytes = r.read_bytes(mos_len)?;
    let nov_len = r.read_var_u64()?;
    let ans_len = r.read_var_u64()? as usize;
    let ans_nov_comp = r.read_bytes(ans_len)?;
    let qual_len = r.read_var_u64()? as usize;
    let qual_comp = r.read_bytes(qual_len)?;
    let cds_comp = if has_cds && r.has_remaining() {
        let clen = r.read_var_u64()? as usize;
        r.read_bytes(clen)?
    } else {
        &[]
    };

    Ok(ParsedHgs3Fastq {
        mode,
        raw_bytes,
        expected_crc,
        total_records,
        panel_names,
        panel_sha,
        total_dna_bases,
        total_qual_bytes,
        read_len_hint,
        meta_comp,
        mosaic_bytes,
        nov_len,
        ans_nov_comp,
        qual_comp,
        cds_comp,
    })
}

pub(crate) fn extract_fastq_read_bounds(
    meta_comp: &[u8],
    n_recs: usize,
) -> Result<(usize, Vec<u64>), String> {
    if n_recs == 0 || meta_comp.is_empty() {
        return Ok((0, Vec::new()));
    }
    let meta_raw = ffi::zlib_decompress(meta_comp).map_err(|e| e.to_string())?;
    let mut mr = ByteReader::new(&meta_raw);
    let uflags = mr.read_u8()?;
    let eol_mode = uflags & 3;
    let fixed_seq_len = (uflags & 16) != 0;
    let _final_eol = mr.read_u8()?;
    let cp_len = mr.read_var_u64()? as usize;
    let _common_prefix = mr.read_bytes(cp_len)?;
    if eol_mode == 2 {
        let _eols = mr.read_bytes(n_recs * 4)?;
    }
    if fixed_seq_len {
        let fl = mr.read_var_u64()? as usize;
        Ok((fl, Vec::new()))
    } else {
        let mut ends = Vec::with_capacity(n_recs);
        let mut acc = 0u64;
        for _ in 0..n_recs {
            acc += mr.read_var_u64()?;
            ends.push(acc);
        }
        Ok((0, ends))
    }
}

fn lcp_len(a: &[u8], b: &[u8]) -> usize {
    let lim = a.len().min(b.len());
    let mut i = 0;
    while i < lim && a[i] == b[i] {
        i += 1;
    }
    i
}

fn compress_fastq_to_hgs3(
    input: &Path,
    output: &Path,
    opts: &Hgs3CompressOptions,
) -> Result<Hgs3Report, String> {
    let t0 = Instant::now();
    let raw = fs::read(input).map_err(|e| e.to_string())?;
    let raw_bytes = raw.len() as u64;
    let raw_crc = ffi::crc32_ieee(0, &raw);
    let threads = opts.base.threads.max(1);

    let mode = if opts.base.reference_panel.is_some() {
        MODE_REF
    } else {
        MODE_SINGLE
    };

    if raw.is_empty() {
        let mut out = ByteWriter::with_capacity(64);
        out.write_bytes(MAGIC3);
        out.write_u8((mode & 3) | FLAG3_FASTQ);
        out.write_u64_le(0);
        out.write_u32_le(0);
        out.write_var_u64(0);
        out.write_var_u64(0);
        if mode == MODE_REF {
            let lp = load_panel_for_compress(&opts.base)?;
            out.write_var_u64(lp.names.len() as u64);
            let cn = compress_panel_names(&lp.names)?;
            out.write_var_u64(cn.len() as u64);
            out.write_bytes(&cn);
            out.write_bytes(&lp.manifest_sha256);
        }
        fs::write(output, &out.buf).map_err(|e| e.to_string())?;
        return Ok(Hgs3Report {
            is_fastq: true,
            mode,
            raw_bytes: 0,
            archive_bytes: out.buf.len() as u64,
            total_records: 0,
            dna_bases: 0,
            qual_bytes: 0,
            meta_comp_bytes: 0,
            mosaic_comp_bytes: 0,
            novelty_comp_bytes: 0,
            qual_comp_bytes: 0,
            cds: None,
            seconds: t0.elapsed().as_secs_f64(),
        });
    }

    struct RecSpan {
        id_s: usize,
        id_e: usize,
        plus_s: usize,
        plus_e: usize,
        seq_s: usize,
        seq_e: usize,
        qual_s: usize,
        qual_e: usize,
        eols: [u8; 4],
    }

    let mut recs: Vec<RecSpan> = Vec::with_capacity((raw.len() / 120).max(16));
    let mut pos = 0usize;
    let n = raw.len();
    let mut total_dna_bases = 0usize;
    let mut total_qual_bytes = 0usize;

    while pos < n {
        if raw[pos] != b'@' {
            return Err(format!(
                "Invalid FASTQ record start byte 0x{:02X} at offset {pos}",
                raw[pos]
            ));
        }
        pos += 1;
        let (id_s, id_e, e1) = scan_line(&raw, &mut pos);
        if pos >= n && e1 == 2 {
            return Err("Truncated FASTQ record after @id line".into());
        }
        let (seq_s, seq_e, e2) = scan_line(&raw, &mut pos);
        if pos >= n || raw[pos] != b'+' {
            return Err(format!("Expected '+' FASTQ separator line at offset {pos}"));
        }
        pos += 1;
        let (plus_s, plus_e, e3) = scan_line(&raw, &mut pos);
        let (qual_s, qual_e, e4) = scan_line(&raw, &mut pos);

        total_dna_bases += seq_e - seq_s;
        total_qual_bytes += qual_e - qual_s;
        recs.push(RecSpan {
            id_s,
            id_e,
            plus_s,
            plus_e,
            seq_s,
            seq_e,
            qual_s,
            qual_e,
            eols: [e1, e2, e3, e4],
        });
    }

    let n_recs = recs.len();
    let first_eol = recs[0].eols[0];
    let uniform_interior_eol = first_eol < 2
        && recs
            .iter()
            .enumerate()
            .all(|(i, r)| {
                r.eols[0] == first_eol
                    && r.eols[1] == first_eol
                    && r.eols[2] == first_eol
                    && (i + 1 == n_recs || r.eols[3] == first_eol)
            });
    let eol_mode: u8 = if uniform_interior_eol { first_eol } else { 2 };
    let final_eol: u8 = recs[n_recs - 1].eols[3];

    let all_plus_empty = recs.iter().all(|r| r.plus_s == r.plus_e);
    let all_plus_eq_id = !all_plus_empty
        && recs
            .iter()
            .all(|r| raw[r.plus_s..r.plus_e] == raw[r.id_s..r.id_e]);
    let plus_mode: u8 = if all_plus_empty {
        0
    } else if all_plus_eq_id {
        1
    } else {
        2
    };

    let first_seq_len = recs[0].seq_e - recs[0].seq_s;
    let fixed_seq_len = recs.iter().all(|r| (r.seq_e - r.seq_s) == first_seq_len);
    let fixed_qual_len = recs
        .iter()
        .all(|r| (r.qual_e - r.qual_s) == (r.seq_e - r.seq_s));

    let mut sorted_lens: Vec<usize> = recs.iter().map(|r| r.seq_e - r.seq_s).collect();
    sorted_lens.sort_unstable();
    let read_len_hint = (sorted_lens[sorted_lens.len() / 2] as u32).max(1);

    let mut common_prefix = &raw[recs[0].id_s..recs[0].id_e];
    for r in &recs[1..] {
        let cur = &raw[r.id_s..r.id_e];
        let l = lcp_len(common_prefix, cur);
        common_prefix = &common_prefix[..l];
        if common_prefix.is_empty() {
            break;
        }
    }
    let cp_len = common_prefix.len();

    // Extract contiguous DNA 2-bit array, ambiguity/case exceptions, and quality stream
    let mut dna_2bit = Vec::with_capacity(total_dna_bases);
    let mut quals_all = Vec::with_capacity(total_qual_bytes);
    let mut exc_pos: Vec<usize> = Vec::new();
    let mut exc_byte: Vec<u8> = Vec::new();
    let acgt_up = *b"ACGT";

    for r in &recs {
        for &ch in &raw[r.seq_s..r.seq_e] {
            let gpos = dna_2bit.len();
            let tb = match ch {
                b'A' | b'a' => 0u8,
                b'C' | b'c' => 1u8,
                b'G' | b'g' => 2u8,
                b'T' | b't' | b'U' | b'u' => 3u8,
                _ => 0u8,
            };
            dna_2bit.push(tb);
            if ch != acgt_up[tb as usize] {
                exc_pos.push(gpos);
                exc_byte.push(ch);
            }
        }
        quals_all.extend_from_slice(&raw[r.qual_s..r.qual_e]);
    }

    // Serialize FASTQ metadata stream cleanly and zlib-compress once
    let mut mw = ByteWriter::with_capacity(n_recs * 12 + 256);
    let uflags = (eol_mode & 3)
        | ((plus_mode & 3) << 2)
        | (if fixed_seq_len { 16 } else { 0 })
        | (if fixed_qual_len { 32 } else { 0 });
    mw.write_u8(uflags);
    mw.write_u8(final_eol);
    mw.write_var_u64(cp_len as u64);
    mw.write_bytes(common_prefix);

    if eol_mode == 2 {
        for r in &recs {
            mw.write_bytes(&r.eols);
        }
    }
    if fixed_seq_len {
        mw.write_var_u64(first_seq_len as u64);
    } else {
        for r in &recs {
            mw.write_var_u64((r.seq_e - r.seq_s) as u64);
        }
    }
    if !fixed_qual_len {
        for r in &recs {
            mw.write_var_u64((r.qual_e - r.qual_s) as u64);
        }
    }

    let mut prev_suffix: &[u8] = &[];
    for r in &recs {
        let suf = &raw[r.id_s + cp_len..r.id_e];
        let l = lcp_len(prev_suffix, suf);
        mw.write_var_u64(l as u64);
        mw.write_var_u64((suf.len() - l) as u64);
        mw.write_bytes(&suf[l..]);
        prev_suffix = suf;
    }

    if plus_mode == 2 {
        for r in &recs {
            let plus = &raw[r.plus_s..r.plus_e];
            if plus.is_empty() {
                mw.write_u8(0);
            } else if plus == &raw[r.id_s..r.id_e] {
                mw.write_u8(1);
            } else {
                mw.write_u8(2);
                mw.write_var_u64(plus.len() as u64);
                mw.write_bytes(plus);
            }
        }
    }

    if exc_pos.is_empty() {
        mw.write_u8(0);
    } else if exc_pos.len() * 3 <= total_dna_bases {
        mw.write_u8(1);
        mw.write_var_u64(exc_pos.len() as u64);
        let mut last = 0usize;
        for i in 0..exc_pos.len() {
            mw.write_var_u64((exc_pos[i] - last) as u64);
            mw.write_u8(exc_byte[i]);
            last = exc_pos[i];
        }
    } else {
        mw.write_u8(2);
        let mut dense = vec![0u8; total_dna_bases];
        for i in 0..exc_pos.len() {
            let p = exc_pos[i];
            dense[p] = exc_byte[i] ^ acgt_up[dna_2bit[p] as usize];
        }
        mw.write_bytes(&dense);
    }

    let meta_comp = ffi::zlib_compress(&mw.buf).map_err(|e| e.to_string())?;
    drop(mw);

    // Reference-conditioned DNA Mosaic + Parallel 4-stripe rANS Novelty
    let loaded_panel = if mode == MODE_REF {
        Some(load_panel_for_compress(&opts.base)?)
    } else {
        None
    };
    let panel_refs: Vec<&[u8]> = loaded_panel
        .as_ref()
        .map_or_else(Vec::new, |lp| lp.seqs.iter().map(|s| s.as_slice()).collect());

    let encoded_mosaic = if let Some(ref lp) = loaded_panel {
        let total_panel: u64 = panel_refs.iter().map(|s| s.len() as u64).sum();
        let sub_shift = if total_panel > 300_000_000 { 2 } else { 0 };
        let build_sec = total_panel < 250_000_000;
        let idx = SyncmerIndex::build(&panel_refs, sub_shift, build_sec, threads);
        ffi::encode_mosaic(
            &dna_2bit,
            Some(&idx),
            lp.seqs.len(),
            -1,
            &[],
            opts.base.mosaic_flags,
            threads,
        )?
    } else {
        ffi::encode_mosaic(&dna_2bit, None, 0, -1, &[], opts.base.mosaic_flags, threads)?
    };

    let nov_len = encoded_mosaic.novelty_bases.len() as u64;
    let ans_nov_comp = ffi::ans_novelty_encode(&encoded_mosaic.novelty_bases, threads)?;

    // Pangenome-conditioned quality coding using deterministic Mosaic mismatches & homopolymers
    let inspect_rep = ffi::mosaic_inspect(
        &encoded_mosaic.mosaic_bytes,
        &encoded_mosaic.novelty_bases,
        &panel_refs,
        panel_refs.len(),
        true,
        false,
    )?;
    let mut mismatch_mask = build_mismatch_mask(dna_2bit.len(), &inspect_rep.variants);
    let mut homo_bins = ffi::compute_homo_bins(&dna_2bit);
    mismatch_mask.resize(quals_all.len(), 0);
    homo_bins.resize(quals_all.len(), 0);
    let qual_comp = ffi::ans4_encode_qual_mt(&quals_all, &mismatch_mask, &homo_bins, read_len_hint, threads)?;

    let cds_summary = if opts.emit_cds_digest {
        scan_cds_summary(&dna_2bit)
    } else {
        None
    };

    let mut out = ByteWriter::with_capacity(128 + meta_comp.len() + encoded_mosaic.mosaic_bytes.len());
    let mut hflags = (mode & 3) | FLAG3_FASTQ;
    if cds_summary.is_some() {
        hflags |= FLAG3_HAS_CDS;
    }
    out.write_bytes(MAGIC3);
    out.write_u8(hflags);
    out.write_u64_le(raw_bytes);
    out.write_u32_le(raw_crc);
    out.write_var_u64(n_recs as u64);
    out.write_var_u64(1);

    if let Some(ref lp) = loaded_panel {
        out.write_var_u64(lp.names.len() as u64);
        let cn = compress_panel_names(&lp.names)?;
        out.write_var_u64(cn.len() as u64);
        out.write_bytes(&cn);
        out.write_bytes(&lp.manifest_sha256);
    }

    out.write_var_u64(total_dna_bases as u64);
    out.write_var_u64(total_qual_bytes as u64);
    out.write_var_u64(read_len_hint as u64);

    out.write_var_u64(meta_comp.len() as u64);
    out.write_bytes(&meta_comp);
    out.write_var_u64(encoded_mosaic.mosaic_bytes.len() as u64);
    out.write_bytes(&encoded_mosaic.mosaic_bytes);
    out.write_var_u64(nov_len);
    out.write_var_u64(ans_nov_comp.len() as u64);
    out.write_bytes(&ans_nov_comp);
    out.write_var_u64(qual_comp.len() as u64);
    out.write_bytes(&qual_comp);
    if let Some(ref cds) = cds_summary {
        let cb = serialize_cds_summary(cds)?;
        out.write_var_u64(cb.len() as u64);
        out.write_bytes(&cb);
    }

    fs::write(output, &out.buf).map_err(|e| e.to_string())?;
    Ok(Hgs3Report {
        is_fastq: true,
        mode,
        raw_bytes,
        archive_bytes: out.buf.len() as u64,
        total_records: n_recs as u64,
        dna_bases: total_dna_bases as u64,
        qual_bytes: total_qual_bytes as u64,
        meta_comp_bytes: meta_comp.len(),
        mosaic_comp_bytes: encoded_mosaic.mosaic_bytes.len(),
        novelty_comp_bytes: ans_nov_comp.len(),
        qual_comp_bytes: qual_comp.len(),
        cds: cds_summary,
        seconds: t0.elapsed().as_secs_f64(),
    })
}

fn decompress_fastq_hgs3(
    bytes: &[u8],
    output: &Path,
    panel: Option<&Path>,
    threads: i32,
) -> Result<f64, String> {
    let t0 = Instant::now();
    let p = parse_hgs3_fastq(bytes, true)?;

    let mut panel_seqs = Vec::new();
    if p.mode == MODE_REF && !p.panel_names.is_empty() {
        panel_seqs = verify_and_load_panel_for_decode(panel, &p.panel_names, &p.panel_sha)?;
    }

    if p.raw_bytes == 0 || p.total_records == 0 {
        if p.expected_crc != 0 {
            return Err("CRC mismatch on empty HGS3 archive".into());
        }
        File::create(output).map_err(|e| e.to_string())?;
        return Ok(t0.elapsed().as_secs_f64());
    }

    let nov_dna = ffi::ans_novelty_decode(p.ans_nov_comp, p.nov_len, threads)?;
    let panel_refs: Vec<&[u8]> = panel_seqs.iter().map(|v| v.as_slice()).collect();
    let dna_2bit = ffi::decode_mosaic(p.mosaic_bytes, &nov_dna, &panel_refs)?;
    if dna_2bit.len() != p.total_dna_bases as usize {
        return Err(format!(
            "HGS3 FASTQ DNA length mismatch: {} vs {}",
            dna_2bit.len(),
            p.total_dna_bases
        ));
    }

    let inspect_rep = ffi::mosaic_inspect(
        p.mosaic_bytes,
        &nov_dna,
        &panel_refs,
        panel_refs.len(),
        true,
        false,
    )?;
    let qual_len = p.total_qual_bytes as usize;
    let mut mismatch_mask = build_mismatch_mask(dna_2bit.len(), &inspect_rep.variants);
    let mut homo_bins = ffi::compute_homo_bins(&dna_2bit);
    mismatch_mask.resize(qual_len, 0);
    homo_bins.resize(qual_len, 0);
    let quals = ffi::ans4_decode_qual_mt(
        p.qual_comp,
        &mismatch_mask,
        &homo_bins,
        qual_len,
        p.read_len_hint,
        threads,
    )?;

    let meta_raw = ffi::zlib_decompress(p.meta_comp).map_err(|e| e.to_string())?;
    let mut mr = ByteReader::new(&meta_raw);
    let uflags = mr.read_u8()?;
    let eol_mode = uflags & 3;
    let plus_mode = (uflags >> 2) & 3;
    let fixed_seq_len = (uflags & 16) != 0;
    let fixed_qual_len = (uflags & 32) != 0;
    let final_eol = mr.read_u8()?;
    let cp_len = mr.read_var_u64()? as usize;
    let common_prefix = mr.read_bytes(cp_len)?;

    let n_recs = p.total_records as usize;
    let mut per_rec_eols = Vec::new();
    if eol_mode == 2 {
        per_rec_eols.resize(n_recs, [0u8; 4]);
        for slot in &mut per_rec_eols {
            slot.copy_from_slice(mr.read_bytes(4)?);
        }
    }

    let mut seq_lens = Vec::with_capacity(n_recs);
    if fixed_seq_len {
        let l = mr.read_var_u64()? as usize;
        seq_lens.resize(n_recs, l);
    } else {
        for _ in 0..n_recs {
            seq_lens.push(mr.read_var_u64()? as usize);
        }
    }

    let mut qual_lens = Vec::with_capacity(n_recs);
    if fixed_qual_len {
        qual_lens.extend_from_slice(&seq_lens);
    } else {
        for _ in 0..n_recs {
            qual_lens.push(mr.read_var_u64()? as usize);
        }
    }

    let mut ids: Vec<Vec<u8>> = Vec::with_capacity(n_recs);
    let mut prev_suffix: Vec<u8> = Vec::new();
    for _ in 0..n_recs {
        let l = mr.read_var_u64()? as usize;
        let rem = mr.read_var_u64()? as usize;
        if l > prev_suffix.len() {
            return Err("Corrupt FASTQ ID LCP length".into());
        }
        let mut suf = Vec::with_capacity(l + rem);
        suf.extend_from_slice(&prev_suffix[..l]);
        suf.extend_from_slice(mr.read_bytes(rem)?);
        let mut full = Vec::with_capacity(common_prefix.len() + suf.len());
        full.extend_from_slice(common_prefix);
        full.extend_from_slice(&suf);
        prev_suffix = suf;
        ids.push(full);
    }

    let mut plus_lines: Vec<Vec<u8>> = Vec::with_capacity(if plus_mode == 2 { n_recs } else { 0 });
    if plus_mode == 2 {
        for id in &ids {
            let tag = mr.read_u8()?;
            match tag {
                0 => plus_lines.push(Vec::new()),
                1 => plus_lines.push(id.clone()),
                2 => {
                    let pl = mr.read_var_u64()? as usize;
                    plus_lines.push(mr.read_bytes(pl)?.to_vec());
                }
                _ => return Err("Invalid FASTQ plus line tag".into()),
            }
        }
    }

    let total_dna = dna_2bit.len();
    let acgt_up = *b"ACGT";
    let mut seq_ascii = Vec::with_capacity(total_dna);
    for &b in &dna_2bit {
        seq_ascii.push(acgt_up[(b & 3) as usize]);
    }
    let mask_mode = mr.read_u8()?;
    match mask_mode {
        0 => {}
        1 => {
            let cnt = mr.read_var_u64()? as usize;
            let mut cur = 0usize;
            for _ in 0..cnt {
                cur += mr.read_var_u64()? as usize;
                let ch = mr.read_u8()?;
                if cur >= seq_ascii.len() {
                    return Err("FASTQ ambiguity index out of range".into());
                }
                seq_ascii[cur] = ch;
            }
        }
        2 => {
            let dense = mr.read_bytes(total_dna)?;
            for (dst, &d) in seq_ascii.iter_mut().zip(dense.iter()) {
                if d != 0 {
                    *dst ^= d;
                }
            }
        }
        _ => return Err("Unknown FASTQ ambiguity mask mode".into()),
    }

    let mut out_f = BufWriter::with_capacity(65536, File::create(output).map_err(|e| e.to_string())?);
    let mut crc = 0u32;
    let mut written = 0u64;
    let mut write_slice = |s: &[u8]| -> Result<(), String> {
        if !s.is_empty() {
            out_f.write_all(s).map_err(|e| e.to_string())?;
            crc = ffi::crc32_ieee(crc, s);
            written += s.len() as u64;
        }
        Ok(())
    };

    let mut d_off = 0usize;
    let mut q_off = 0usize;
    for i in 0..n_recs {
        let eols = if eol_mode == 2 {
            per_rec_eols[i]
        } else if i + 1 == n_recs {
            [eol_mode, eol_mode, eol_mode, final_eol]
        } else {
            [eol_mode; 4]
        };

        write_slice(b"@")?;
        write_slice(&ids[i])?;
        write_slice(eol_bytes(eols[0]))?;

        let sl = seq_lens[i];
        write_slice(&seq_ascii[d_off..d_off + sl])?;
        write_slice(eol_bytes(eols[1]))?;
        d_off += sl;

        write_slice(b"+")?;
        if plus_mode == 1 {
            write_slice(&ids[i])?;
        } else if plus_mode == 2 {
            write_slice(&plus_lines[i])?;
        }
        write_slice(eol_bytes(eols[2]))?;

        let ql = qual_lens[i];
        write_slice(&quals[q_off..q_off + ql])?;
        write_slice(eol_bytes(eols[3]))?;
        q_off += ql;
    }
    out_f.flush().map_err(|e| e.to_string())?;

    if written != p.raw_bytes {
        return Err(format!(
            "HGS3 FASTQ byte length mismatch: {written} vs {}",
            p.raw_bytes
        ));
    }
    if crc != p.expected_crc {
        return Err(format!(
            "HGS3 FASTQ CRC-32 mismatch: 0x{crc:08X} vs 0x{:08X}",
            p.expected_crc
        ));
    }
    Ok(t0.elapsed().as_secs_f64())
}

pub fn compress_omics(
    input: &Path,
    output: &Path,
    opts: &Hgs3CompressOptions,
) -> Result<Hgs3Report, String> {
    let kind = detect_input_kind(input, opts.format)?;
    if kind == OmicInputKind::Fastq {
        return compress_fastq_to_hgs3(input, output, opts);
    }

    let cds_bytes_opt = if opts.emit_cds_digest {
        let raw = fs::read(input).map_err(|e| e.to_string())?;
        let mut two = Vec::with_capacity(raw.len());
        let mut in_hdr = false;
        let mut at_sol = true;
        for &b in &raw {
            if in_hdr {
                if b == b'\n' {
                    in_hdr = false;
                    at_sol = true;
                }
                continue;
            }
            if b == b'>' && at_sol {
                in_hdr = true;
                at_sol = false;
                continue;
            }
            at_sol = b == b'\n' || b == b'\r';
            match b {
                b'A' | b'a' => two.push(0),
                b'C' | b'c' => two.push(1),
                b'G' | b'g' => two.push(2),
                b'T' | b't' | b'U' | b'u' => two.push(3),
                _ => {}
            }
        }
        scan_cds_summary(&two)
    } else {
        None
    };

    let cds_comp = if let Some(ref c) = cds_bytes_opt {
        Some(serialize_cds_summary(c)?)
    } else {
        None
    };

    let rep = hgs2::compress_with_container(
        input,
        output,
        &opts.base,
        MAGIC3,
        true,
        cds_comp.as_deref(),
    )?;

    Ok(Hgs3Report {
        is_fastq: false,
        mode: rep.mode,
        raw_bytes: rep.raw_bytes,
        archive_bytes: rep.archive_bytes,
        total_records: 0,
        dna_bases: rep.target_dna_bases,
        qual_bytes: 0,
        meta_comp_bytes: 0,
        mosaic_comp_bytes: 0,
        novelty_comp_bytes: 0,
        qual_comp_bytes: 0,
        cds: cds_bytes_opt,
        seconds: rep.seconds,
    })
}

pub fn decompress_omics(
    archive: &Path,
    output: &Path,
    panel: Option<&Path>,
    threads: i32,
    bounded: bool,
) -> Result<f64, String> {
    let bytes = fs::read(archive).map_err(|e| e.to_string())?;
    if bytes.len() >= 4 && &bytes[..4] == hgs2::MAGIC {
        return hgs2::decompress_file(archive, output, panel, threads, bounded);
    }
    if bytes.len() < 5 || &bytes[..4] != MAGIC3 {
        return Err("Archive is neither HGS2 nor HGS3".into());
    }
    if (bytes[4] & FLAG3_FASTQ) != 0 {
        return decompress_fastq_hgs3(&bytes, output, panel, threads);
    }
    hgs2::decompress_with_container(archive, output, panel, threads, bounded, MAGIC3)
}
