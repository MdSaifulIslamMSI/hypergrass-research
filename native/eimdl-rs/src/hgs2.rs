use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::fasta::{self, ByteReader, ByteWriter, FactorizedChunk};
use crate::ffi::{self, MosaicFlags, SyncmerIndex};

pub const MAGIC: &[u8; 4] = b"HGS2";
pub const MODE_SINGLE: u8 = 0;
pub const MODE_REF: u8 = 1;
pub const MODE_COLLECTION: u8 = 2;
pub const FLAG_RAW_FLAT_ACGT: u8 = 4;
pub const FLAG_RAW_RNA_U: u8 = 8;

pub const TARGET_UNIT_BYTES: u64 = 4 * 1024 * 1024;

#[derive(Clone)]
pub struct CompressOptions {
    pub mode: i32, // -1 auto, 0 single, 1 ref, 2 collection
    pub single_genome: bool,
    pub bounded: bool,
    pub threads: i32,
    pub reference_panel: Option<PathBuf>,
    pub exclude_panel_names: HashSet<String>,
    pub explicit_panel_names: Option<Vec<String>>,
    pub panel_limit: usize,
    pub mosaic_flags: MosaicFlags,
}

impl Default for CompressOptions {
    fn default() -> Self {
        Self {
            mode: -1,
            single_genome: false,
            bounded: false,
            threads: std::thread::available_parallelism().map_or(4, |n| n.get() as i32),
            reference_panel: None,
            exclude_panel_names: HashSet::new(),
            explicit_panel_names: None,
            panel_limit: 0,
            mosaic_flags: MosaicFlags::default(),
        }
    }
}

pub struct CompressionReport {
    pub mode: u8,
    pub raw_bytes: u64,
    pub target_dna_bases: u64,
    pub archive_bytes: u64,
    pub panel_genomes: usize,
    pub panel_bases: u64,
    pub panel_match_bases: u64,
    pub self_match_bases: u64,
    pub tandem_match_bases: u64,
    pub novelty_bases: u64,
    pub panel_tokens: u64,
    pub self_tokens: u64,
    pub tandem_tokens: u64,
    pub literal_runs: u64,
    pub indel_splices: u64,
    pub substitutions: u64,
    pub seconds: f64,
}

impl CompressionReport {
    pub fn bpb(&self) -> f64 {
        if self.raw_bytes == 0 { 0.0 } else { (self.archive_bytes as f64 * 8.0) / (self.raw_bytes as f64) }
    }
    pub fn ratio(&self) -> f64 {
        if self.archive_bytes == 0 { 0.0 } else { (self.raw_bytes as f64) / (self.archive_bytes as f64) }
    }
    pub fn mode_name(&self) -> &'static str {
        match self.mode {
            MODE_SINGLE => "MODE_S",
            MODE_REF => "MODE_R",
            MODE_COLLECTION => "MODE_C",
            _ => "UNKNOWN",
        }
    }
}

pub struct ArchiveInfo {
    pub mode: u8,
    pub raw_flat_acgt: bool,
    pub raw_rna_u: bool,
    pub raw_bytes: u64,
    pub crc32: u32,
    pub total_records: u64,
    pub num_units: usize,
    pub total_novelty_bases: u64,
}

fn load_unit_slice(file: &mut File, scratch: &mut Vec<u8>, start: u64, end: u64) -> Result<usize, String> {
    let len = (end - start) as usize;
    if scratch.len() < len {
        scratch.resize(len, 0);
    }
    file.seek(SeekFrom::Start(start)).map_err(|e| e.to_string())?;
    file.read_exact(&mut scratch[..len]).map_err(|e| e.to_string())?;
    Ok(len)
}

fn stream_extract_panel_file(path: &Path) -> Result<(u64, [u8; 32], Vec<u8>), String> {
    let mut f = File::open(path).map_err(|e| e.to_string())?;
    let file_len = f.metadata().map_err(|e| e.to_string())?.len();
    let mut reader = BufReader::with_capacity(65536, &mut f);

    let mut buf = [0u8; 65536];
    let mut two_bit = Vec::with_capacity(file_len.min(64 * 1024 * 1024) as usize);
    let mut in_hdr = false;
    let mut at_line_start = true;
    let mut sha = StreamingSha256::new();
    let mut total = 0u64;
    loop {
        let n = reader.read(&mut buf).map_err(|e| e.to_string())?;
        if n == 0 { break; }
        sha.update(&buf[..n]);
        total += n as u64;
        for &b in &buf[..n] {
            if in_hdr {
                if b == b'\n' {
                    in_hdr = false;
                    at_line_start = true;
                }
                continue;
            }
            if b == b'>' && at_line_start {
                in_hdr = true;
                at_line_start = false;
                continue;
            }
            at_line_start = b == b'\n' || b == b'\r';
            match b {
                b'A' | b'a' => two_bit.push(0),
                b'C' | b'c' => two_bit.push(1),
                b'G' | b'g' => two_bit.push(2),
                b'T' | b't' | b'U' | b'u' => two_bit.push(3),
                _ => {}
            }
        }
    }
    let digest = sha.finish();
    Ok((total, digest, two_bit))
}

struct StreamingSha256 {
    state: [u32; 8],
    bitlen: u64,
    buf: [u8; 64],
    buf_len: usize,
}

impl StreamingSha256 {
    fn new() -> Self {
        Self {
            state: [
                0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
            ],
            bitlen: 0,
            buf: [0u8; 64],
            buf_len: 0,
        }
    }

    fn transform(&mut self, block: &[u8; 64]) {
        const K: [u32; 64] = [
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
            0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
            0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
            0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
            0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
            0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
            0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
            0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
        ];
        let mut w = [0u32; 64];
        for i in 0..16 {
            w[i] = u32::from_be_bytes([block[i * 4], block[i * 4 + 1], block[i * 4 + 2], block[i * 4 + 3]]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16].wrapping_add(s0).wrapping_add(w[i - 7]).wrapping_add(s1);
        }
        let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h] = self.state;
        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ ((!e) & g);
            let t1 = h.wrapping_add(s1).wrapping_add(ch).wrapping_add(K[i]).wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let t2 = s0.wrapping_add(maj);
            h = g; g = f; f = e; e = d.wrapping_add(t1);
            d = c; c = b; b = a; a = t1.wrapping_add(t2);
        }
        self.state[0] = self.state[0].wrapping_add(a);
        self.state[1] = self.state[1].wrapping_add(b);
        self.state[2] = self.state[2].wrapping_add(c);
        self.state[3] = self.state[3].wrapping_add(d);
        self.state[4] = self.state[4].wrapping_add(e);
        self.state[5] = self.state[5].wrapping_add(f);
        self.state[6] = self.state[6].wrapping_add(g);
        self.state[7] = self.state[7].wrapping_add(h);
    }

    fn update(&mut self, data: &[u8]) {
        self.bitlen += (data.len() as u64) * 8;
        let mut i = 0usize;
        if self.buf_len > 0 {
            let take = (64 - self.buf_len).min(data.len());
            self.buf[self.buf_len..self.buf_len + take].copy_from_slice(&data[..take]);
            self.buf_len += take;
            i += take;
            if self.buf_len == 64 {
                let blk = self.buf;
                self.transform(&blk);
                self.buf_len = 0;
            }
        }
        while i + 64 <= data.len() {
            let blk: [u8; 64] = data[i..i + 64].try_into().unwrap();
            self.transform(&blk);
            i += 64;
        }
        if i < data.len() {
            self.buf_len = data.len() - i;
            self.buf[..self.buf_len].copy_from_slice(&data[i..]);
        }
    }

    fn finish(mut self) -> [u8; 32] {
        self.buf[self.buf_len] = 0x80;
        self.buf_len += 1;
        if self.buf_len > 56 {
            self.buf[self.buf_len..].fill(0);
            let blk = self.buf;
            self.transform(&blk);
            self.buf_len = 0;
        }
        self.buf[self.buf_len..56].fill(0);
        self.buf[56..64].copy_from_slice(&self.bitlen.to_be_bytes());
        let blk = self.buf;
        self.transform(&blk);
        let mut out = [0u8; 32];
        for (i, w) in self.state.iter().enumerate() {
            out[i * 4..(i + 1) * 4].copy_from_slice(&w.to_be_bytes());
        }
        out
    }
}

fn compute_manifest_sha256(names: &[String], byte_lens: &[u64], hashes: &[[u8; 32]]) -> [u8; 32] {
    let mut sha = StreamingSha256::new();
    sha.update(&(names.len() as u32).to_be_bytes());
    for i in 0..names.len() {
        let nb = names[i].as_bytes();
        sha.update(&(nb.len() as u32).to_be_bytes());
        sha.update(nb);
        sha.update(&byte_lens[i].to_be_bytes());
        sha.update(&hashes[i]);
    }
    sha.finish()
}

fn compress_panel_names(names: &[String]) -> Result<Vec<u8>, String> {
    let mut raw = ByteWriter::with_capacity(names.len() * 8 + 16);
    let mut prev: &[u8] = &[];
    for s in names {
        let cur = s.as_bytes();
        let max_s = prev.len().min(cur.len());
        let mut lcp = 0;
        while lcp < max_s && prev[lcp] == cur[lcp] { lcp += 1; }
        raw.write_var_u64(lcp as u64);
        raw.write_var_u64((cur.len() - lcp) as u64);
        raw.write_bytes(&cur[lcp..]);
        prev = cur;
    }
    ffi::zlib_compress(&raw.buf).map_err(|e| e.to_string())
}

fn decompress_panel_names(comp: &[u8], count: usize) -> Result<Vec<String>, String> {
    let raw = ffi::zlib_decompress(comp).map_err(|e| e.to_string())?;
    let mut inp = ByteReader::new(&raw);
    let mut names = Vec::with_capacity(count);
    let mut prev = Vec::new();
    for _ in 0..count {
        let lcp = inp.read_var_u64()? as usize;
        let rem = inp.read_var_u64()? as usize;
        if lcp > prev.len() {
            return Err("Malformed panel name prefix length".into());
        }
        let mut cur = Vec::with_capacity(lcp + rem);
        cur.extend_from_slice(&prev[..lcp]);
        cur.extend_from_slice(inp.read_bytes(rem)?);
        let s = String::from_utf8(cur.clone()).map_err(|e| e.to_string())?;
        prev = cur;
        names.push(s);
    }
    if inp.has_remaining() {
        return Err("Trailing bytes in panel names stream".into());
    }
    Ok(names)
}

struct LoadedPanel {
    names: Vec<String>,
    seqs: Vec<Vec<u8>>,
    manifest_sha256: [u8; 32],
    total_bases: u64,
}

fn load_panel_for_compress(opts: &CompressOptions) -> Result<LoadedPanel, String> {
    let panel_path = opts.reference_panel.as_ref().ok_or("Missing reference panel")?;
    let mut files = Vec::new();
    if panel_path.is_file() {
        files.push(panel_path.clone());
    } else if panel_path.is_dir() {
        if let Some(ref exp) = opts.explicit_panel_names {
            for name in exp {
                if opts.exclude_panel_names.contains(name) { continue; }
                let p = panel_path.join(name);
                if !p.is_file() {
                    return Err(format!("Explicit panel file missing: {}", p.display()));
                }
                files.push(p);
            }
        } else {
            for entry in fs::read_dir(panel_path).map_err(|e| e.to_string())? {
                let p = entry.map_err(|e| e.to_string())?.path();
                if !p.is_file() { continue; }
                let fname = p.file_name().and_then(|n| n.to_str()).unwrap_or("").to_string();
                if opts.exclude_panel_names.contains(&fname) { continue; }
                let low = fname.to_lowercase();
                if low.ends_with(".seq") || low.ends_with(".fa") || low.ends_with(".fna") || low.ends_with(".fasta") {
                    files.push(p);
                }
            }
            files.sort_by(|a, b| a.file_name().cmp(&b.file_name()));
        }
    } else {
        return Err(format!("Reference panel path not found: {}", panel_path.display()));
    }
    if opts.panel_limit > 0 && files.len() > opts.panel_limit {
        files.truncate(opts.panel_limit);
    }

    let mut names = Vec::with_capacity(files.len());
    let mut seqs = Vec::with_capacity(files.len());
    let mut hashes = Vec::with_capacity(files.len());
    let mut lens = Vec::with_capacity(files.len());
    let mut total_bases = 0u64;
    for p in files {
        names.push(p.file_name().and_then(|n| n.to_str()).unwrap_or("").to_string());
        let (blen, digest, dna) = stream_extract_panel_file(&p)?;
        total_bases += dna.len() as u64;
        lens.push(blen);
        hashes.push(digest);
        seqs.push(dna);
    }
    let manifest_sha256 = compute_manifest_sha256(&names, &lens, &hashes);
    Ok(LoadedPanel { names, seqs, manifest_sha256, total_bases })
}

fn verify_and_load_panel_for_decode(panel_path: Option<&Path>, expected_names: &[String], expected_sha: &[u8; 32]) -> Result<Vec<Vec<u8>>, String> {
    let pp = panel_path.ok_or("Mode R archive requires --panel / --ref path to decode")?;
    let mut files = Vec::with_capacity(expected_names.len());
    if pp.is_file() {
        let got_name = pp.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if expected_names.len() != 1 || got_name != expected_names[0] {
            return Err("Reference file name mismatch".into());
        }
        files.push(pp.to_path_buf());
    } else if pp.is_dir() {
        for name in expected_names {
            let p = pp.join(name);
            if !p.is_file() {
                return Err(format!("Missing required reference panel file: {}", p.display()));
            }
            files.push(p);
        }
    } else {
        return Err("Reference panel path not found".into());
    }

    let mut seqs = Vec::with_capacity(files.len());
    let mut hashes = Vec::with_capacity(files.len());
    let mut lens = Vec::with_capacity(files.len());
    for p in files {
        let (blen, digest, dna) = stream_extract_panel_file(&p)?;
        lens.push(blen);
        hashes.push(digest);
        seqs.push(dna);
    }
    let actual_sha = compute_manifest_sha256(expected_names, &lens, &hashes);
    if &actual_sha != expected_sha {
        return Err("Reference panel SHA-256 manifest mismatch".into());
    }
    Ok(seqs)
}

fn check_flat_upper_file(path: &Path) -> Result<Option<(u64, u32, bool)>, String> {
    let mut f = File::open(path).map_err(|e| e.to_string())?;
    let size = f.metadata().map_err(|e| e.to_string())?.len();
    if size == 0 {
        return Ok(None);
    }
    let mut buf = vec![0u8; 65536];
    let mut crc = 0u32;
    let mut has_t = false;
    let mut has_u = false;
    loop {
        let n = f.read(&mut buf).map_err(|e| e.to_string())?;
        if n == 0 { break; }
        crc = ffi::crc32_ieee(crc, &buf[..n]);
        for &b in &buf[..n] {
            match b {
                b'A' | b'C' | b'G' => {}
                b'T' => {
                    has_t = true;
                    if has_u { return Ok(None); }
                }
                b'U' => {
                    has_u = true;
                    if has_t { return Ok(None); }
                }
                _ => return Ok(None),
            }
        }
    }
    Ok(Some((size, crc, has_u && !has_t)))
}

fn build_prim_order(
    sketches: &[Vec<u64>],
    out_rank_to_unit: &mut [usize],
    out_mst_parent: &mut [i32],
    out_sim_order: &mut [Vec<i32>],
) {
    let n = sketches.len();
    let mut sim = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        sim[i][i] = 1.0;
        for j in (i + 1)..n {
            let jac = ffi::sketch_jaccard(&sketches[i], &sketches[j]);
            sim[i][j] = jac;
            sim[j][i] = jac;
        }
    }
    let mut root = 0usize;
    let mut best_sum = -1.0f64;
    for (i, row) in sim.iter().enumerate().take(n) {
        let s: f64 = row.iter().enumerate().filter(|&(j, _)| j != i).map(|(_, &v)| v).sum();
        if s > best_sum {
            best_sum = s;
            root = i;
        }
    }
    let mut visited = vec![false; n];
    let mut max_sim = vec![-1.0f64; n];
    let mut best_tr = vec![0i32; n];
    visited[root] = true;
    out_rank_to_unit[0] = root;
    out_mst_parent[0] = -1;
    for v in 0..n {
        if !visited[v] {
            max_sim[v] = sim[root][v];
        }
    }
    for r in 1..n {
        let mut pick = 0usize;
        let mut bsim = -1.0f64;
        for v in 0..n {
            if !visited[v] && max_sim[v] > bsim {
                bsim = max_sim[v];
                pick = v;
            }
        }
        visited[pick] = true;
        out_rank_to_unit[r] = pick;
        out_mst_parent[r] = best_tr[pick];

        let mut anc: Vec<usize> = (0..r).collect();
        anc.sort_by(|&a, &b| {
            sim[pick][out_rank_to_unit[b]]
                .partial_cmp(&sim[pick][out_rank_to_unit[a]])
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.cmp(&b))
        });
        out_sim_order[r] = anc.into_iter().map(|a| a as i32).collect();

        for v in 0..n {
            if !visited[v] && sim[pick][v] > max_sim[v] {
                max_sim[v] = sim[pick][v];
                best_tr[v] = r as i32;
            }
        }
    }
}

pub fn compress_file(input: &Path, output: &Path, opts: &CompressOptions) -> Result<CompressionReport, String> {
    let t0 = Instant::now();
    let file_size = fs::metadata(input).map_err(|e| e.to_string())?.len();
    let threads = opts.threads.max(1);

    let mut mode = if opts.reference_panel.is_some() {
        MODE_REF
    } else if opts.mode >= 0 {
        opts.mode as u8
    } else if opts.single_genome || file_size < 6_500_000 {
        MODE_SINGLE
    } else {
        MODE_COLLECTION
    };

    if file_size == 0 {
        let mut out = ByteWriter::with_capacity(64);
        out.write_bytes(MAGIC);
        out.write_u8(mode & 3);
        out.write_u64_le(0);
        out.write_u32_le(0);
        out.write_var_u64(0);
        out.write_var_u64(0);
        let mut pg = 0;
        let mut pb = 0;
        if mode == MODE_REF {
            let lp = load_panel_for_compress(opts)?;
            pg = lp.names.len();
            pb = lp.total_bases;
            out.write_var_u64(pg as u64);
            let cn = compress_panel_names(&lp.names)?;
            out.write_var_u64(cn.len() as u64);
            out.write_bytes(&cn);
            out.write_bytes(&lp.manifest_sha256);
        }
        fs::write(output, &out.buf).map_err(|e| e.to_string())?;
        return Ok(CompressionReport {
            mode,
            raw_bytes: 0,
            target_dna_bases: 0,
            archive_bytes: out.buf.len() as u64,
            panel_genomes: pg,
            panel_bases: pb,
            panel_match_bases: 0,
            self_match_bases: 0,
            tandem_match_bases: 0,
            novelty_bases: 0,
            panel_tokens: 0,
            self_tokens: 0,
            tandem_tokens: 0,
            literal_runs: 0,
            indel_splices: 0,
            substitutions: 0,
            seconds: t0.elapsed().as_secs_f64(),
        });
    }

    // Check fast-path flat ACGT/U streaming
    let flat_check = check_flat_upper_file(input)?;
    let raw_flat_acgt = flat_check.is_some();
    let (raw_crc, raw_rna_u) = if let Some((_, crc, rna)) = flat_check {
        (crc, rna)
    } else {
        (0, false)
    };

    let mut unit_ranges: Vec<(u64, u64)> = Vec::new();
    let final_crc;
    if raw_flat_acgt {
        final_crc = raw_crc;
        let num_slices = if mode == MODE_COLLECTION && file_size > (TARGET_UNIT_BYTES * 3) / 2 {
            ((file_size + TARGET_UNIT_BYTES / 2) / TARGET_UNIT_BYTES).max(1)
        } else {
            1
        };
        let mut cur = 0u64;
        for i in 0..num_slices {
            let rem = file_size - cur;
            let len = if i + 1 == num_slices { rem } else { (rem / (num_slices - i)).min(TARGET_UNIT_BYTES) };
            unit_ranges.push((cur, cur + len));
            cur += len;
        }
    } else {
        let part = fasta::partition_file(input, TARGET_UNIT_BYTES, true).map_err(|e| e.to_string())?;
        final_crc = part.crc32;
        if mode != MODE_COLLECTION && file_size <= 16 * 1024 * 1024 && !opts.bounded {
            unit_ranges.push((0, file_size));
        } else {
            let mut idx = 0usize;
            while idx < part.slices.len() {
                let st = part.slices[idx].start_byte;
                let mut en = part.slices[idx].end_byte;
                idx += 1;
                while idx < part.slices.len() && (en - st) + (part.slices[idx].end_byte - part.slices[idx].start_byte) <= TARGET_UNIT_BYTES + 524_288 {
                    en = part.slices[idx].end_byte;
                    idx += 1;
                }
                unit_ranges.push((st, en));
            }
        }
        if mode == MODE_COLLECTION && unit_ranges.len() <= 1 && opts.mode < 0 {
            mode = MODE_SINGLE;
        }
    }

    let num_units = unit_ranges.len();
    let mut out_file = BufWriter::with_capacity(65536, File::create(output).map_err(|e| e.to_string())?);

    // Determine whether we stream unit DNA from disk in Stage 1 (--bounded or huge files)
    let stream_from_disk = opts.bounded || file_size > 800 * 1024 * 1024;
    let hub_cap: usize = if mode == MODE_COLLECTION && num_units > 1 && stream_from_disk {
        num_units.min(if file_size > 1_000_000_000 { 16 } else { 32 })
    } else {
        num_units
    };

    // Pass 1: Scan units one-by-one using single scratch buffer to serialize aux metadata & build sketches
    let mut fin = File::open(input).map_err(|e| e.to_string())?;
    let mut raw_scratch = Vec::new();
    let mut total_records = 0u64;
    let mut first_rec_idx = vec![0u64; num_units];
    let mut rec_counts = vec![0u64; num_units];
    let mut flat_unit_lens = vec![0u64; num_units];
    let mut aux_blobs: Vec<Vec<u8>> = if raw_flat_acgt { Vec::new() } else { Vec::with_capacity(num_units) };
    let mut sketches: Vec<Vec<u64>> = if mode == MODE_COLLECTION && num_units > 1 { Vec::with_capacity(num_units) } else { Vec::new() };
    let mut kept_dna: Vec<Option<Vec<u8>>> = vec![None; num_units];

    for (u, &(st, en)) in unit_ranges.iter().enumerate() {
        let len = load_unit_slice(&mut fin, &mut raw_scratch, st, en)?;
        let slice = &raw_scratch[..len];
        let dna = if raw_flat_acgt {
            let mut two = vec![0u8; len];
            ffi::ascii_to_2bit_checked(slice, &mut two, if raw_rna_u { b'U' } else { b'T' });
            flat_unit_lens[u] = len as u64;
            if u == 0 {
                total_records = 1;
                rec_counts[u] = 1;
            }
            two
        } else {
            let mut c = fasta::factorize(slice);
            let rc = c.record_count() as u64;
            first_rec_idx[u] = total_records;
            rec_counts[u] = rc;
            total_records += rc;
            let mut abuf = ByteWriter::with_capacity(256);
            fasta::write_aux_metadata(&mut abuf, &c);
            aux_blobs.push(abuf.buf);
            std::mem::take(&mut c.dna_bases)
        };

        if mode == MODE_COLLECTION && num_units > 1 {
            sketches.push(ffi::compute_bottom_sketch(&dna, 2048));
        }
        if !stream_from_disk {
            kept_dna[u] = Some(dna);
        }
    }

    // Write HGS2 container header directly
    let mut hdr = ByteWriter::with_capacity(64 + num_units * 16);
    let mut hflags = mode & 3;
    if raw_flat_acgt { hflags |= FLAG_RAW_FLAT_ACGT; }
    if raw_rna_u { hflags |= FLAG_RAW_RNA_U; }
    hdr.write_bytes(MAGIC);
    hdr.write_u8(hflags);
    hdr.write_u64_le(file_size);
    hdr.write_u32_le(final_crc);
    hdr.write_var_u64(total_records);
    hdr.write_var_u64(num_units as u64);

    let loaded_panel = if mode == MODE_REF {
        let lp = load_panel_for_compress(opts)?;
        hdr.write_var_u64(lp.names.len() as u64);
        let cn = compress_panel_names(&lp.names)?;
        hdr.write_var_u64(cn.len() as u64);
        hdr.write_bytes(&cn);
        hdr.write_bytes(&lp.manifest_sha256);
        Some(lp)
    } else {
        None
    };

    let mut rank_to_unit = vec![0usize; num_units];
    let mut mst_parent = vec![-1i32; num_units];
    let mut sim_orders = vec![Vec::new(); num_units];

    if mode == MODE_COLLECTION && num_units > 1 {
        build_prim_order(&sketches, &mut rank_to_unit, &mut mst_parent, &mut sim_orders);
        for &ru in &rank_to_unit { hdr.write_var_u64(ru as u64); }
        for r in 0..num_units { hdr.write_var_u64(r.min(hub_cap) as u64); }
        drop(sketches);
    } else {
        for (u, slot) in rank_to_unit.iter_mut().enumerate().take(num_units) { *slot = u; }
    }

    if !raw_flat_acgt {
        for u in 0..num_units {
            hdr.write_var_u64(first_rec_idx[u]);
            hdr.write_var_u64(rec_counts[u]);
            hdr.write_var_u64(aux_blobs[u].len() as u64);
            hdr.write_bytes(&aux_blobs[u]);
        }
        drop(aux_blobs);
    } else {
        for &len in flat_unit_lens.iter().take(num_units) {
            hdr.write_var_u64(len);
        }
    }
    out_file.write_all(&hdr.buf).map_err(|e| e.to_string())?;
    let mut total_written_bytes = hdr.buf.len() as u64;

    // Helper to get unit DNA (from kept_dna OR streamed directly from disk if stream_from_disk!)
    let mut fetch_unit_dna = |u: usize| -> Result<Vec<u8>, String> {
        if let Some(d) = kept_dna[u].take() {
            return Ok(d);
        }
        let (st, en) = unit_ranges[u];
        let len = load_unit_slice(&mut fin, &mut raw_scratch, st, en)?;
        let slice = &raw_scratch[..len];
        if raw_flat_acgt {
            let mut two = vec![0u8; len];
            ffi::ascii_to_2bit_checked(slice, &mut two, if raw_rna_u { b'U' } else { b'T' });
            Ok(two)
        } else {
            Ok(fasta::factorize(slice).dna_bases)
        }
    };

    let mut report = CompressionReport {
        mode,
        raw_bytes: file_size,
        target_dna_bases: 0,
        archive_bytes: 0,
        panel_genomes: loaded_panel.as_ref().map_or(0, |p| p.names.len()),
        panel_bases: loaded_panel.as_ref().map_or(0, |p| p.total_bases),
        panel_match_bases: 0,
        self_match_bases: 0,
        tandem_match_bases: 0,
        novelty_bases: 0,
        panel_tokens: 0,
        self_tokens: 0,
        tandem_tokens: 0,
        literal_runs: 0,
        indel_splices: 0,
        substitutions: 0,
        seconds: 0.0,
    };

    // Stage 1: Mosaic Encoding (streaming per rank directly into output + compact novelty collector)
    let mut novelty_super = Vec::new();

    let mut emit_mosaic = |em: ffi::MosaicEncoded, target_len: u64, wout: &mut BufWriter<File>| -> Result<(), String> {
        report.target_dna_bases += target_len;
        report.panel_match_bases += em.panel_match_bases;
        report.self_match_bases += em.self_match_bases;
        report.tandem_match_bases += em.tandem_match_bases;
        report.novelty_bases += em.novelty_bases.len() as u64;
        report.panel_tokens += em.panel_tokens;
        report.self_tokens += em.self_tokens;
        report.tandem_tokens += em.tandem_tokens;
        report.literal_runs += em.literal_runs;
        report.indel_splices += em.indel_splices;
        report.substitutions += em.substitutions;

        let mut lw = ByteWriter::with_capacity(8);
        lw.write_var_u64(em.mosaic_bytes.len() as u64);
        wout.write_all(&lw.buf).map_err(|e| e.to_string())?;
        wout.write_all(&em.mosaic_bytes).map_err(|e| e.to_string())?;
        total_written_bytes += (lw.buf.len() + em.mosaic_bytes.len()) as u64;
        novelty_super.extend_from_slice(&em.novelty_bases);
        Ok(())
    };

    if mode == MODE_SINGLE {
        for u in 0..num_units {
            let dna = fetch_unit_dna(u)?;
            let tl = dna.len() as u64;
            let em = ffi::encode_mosaic(&dna, None, 0, -1, &[], opts.mosaic_flags, threads)?;
            emit_mosaic(em, tl, &mut out_file)?;
        }
    } else if mode == MODE_REF {
        let lp = loaded_panel.as_ref().unwrap();
        let refs: Vec<&[u8]> = lp.seqs.iter().map(|v| v.as_slice()).collect();
        let idx = SyncmerIndex::build(&refs, 0, !opts.bounded, threads);
        let panel_sk: Vec<Vec<u64>> = lp.seqs.iter().map(|s| ffi::compute_bottom_sketch(s, 2048)).collect();
        for u in 0..num_units {
            let dna = fetch_unit_dna(u)?;
            let tl = dna.len() as u64;
            let tsk = ffi::compute_bottom_sketch(&dna, 2048);
            let mut sims: Vec<(f64, usize)> = panel_sk.iter().enumerate()
                .map(|(d, ps)| (ffi::sketch_jaccard(&tsk, ps), d))
                .collect();
            sims.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal).then_with(|| a.1.cmp(&b.1)));
            let pref = if !sims.is_empty() && sims[0].0 > 0.05 { sims[0].1 as i64 } else { -1 };
            let ord: Vec<i32> = sims.into_iter().map(|(_, d)| d as i32).collect();
            let em = ffi::encode_mosaic(&dna, Some(&idx), lp.seqs.len(), pref, &ord, opts.mosaic_flags, threads)?;
            emit_mosaic(em, tl, &mut out_file)?;
        }
    } else {
        // MODE_COLLECTION (bounded hub bank or full bank)
        let mut hub_dna: Vec<Vec<u8>> = Vec::with_capacity(hub_cap);
        for &u in rank_to_unit.iter().take(hub_cap) {
            hub_dna.push(fetch_unit_dna(u)?);
        }
        let hub_refs: Vec<&[u8]> = hub_dna.iter().map(|v| v.as_slice()).collect();
        let hub_bytes: u64 = hub_dna.iter().map(|v| (v.len() as u64 + 3) >> 2).sum();
        let sub_shift = if opts.bounded && hub_bytes > 64 * 1024 * 1024 { 1 } else { 0 };
        let build_sec = !opts.bounded && !stream_from_disk && num_units <= 16;
        let idx = SyncmerIndex::build(&hub_refs, sub_shift, build_sec, threads);

        for r in 0..num_units {
            let dna = if r < hub_cap {
                std::mem::take(&mut hub_dna[r])
            } else {
                fetch_unit_dna(rank_to_unit[r])?
            };
            let tl = dna.len() as u64;
            let em = if r == 0 {
                ffi::encode_mosaic(&dna, None, 0, -1, &[], opts.mosaic_flags, threads)?
            } else {
                let anc = r.min(hub_cap);
                let pref = if mst_parent[r] >= 0 && (mst_parent[r] as usize) < anc {
                    mst_parent[r] as i64
                } else {
                    0
                };
                ffi::encode_mosaic(&dna, Some(&idx), anc, pref, &sim_orders[r], opts.mosaic_flags, threads)?
            };
            emit_mosaic(em, tl, &mut out_file)?;
            sim_orders[r].clear();
        }
    }

    // Stage 2: Unified Novelty Superstream
    let total_nov = novelty_super.len() as u64;
    let nov_comp = ffi::novelty_encode(&novelty_super, threads, opts.bounded)?;
    drop(novelty_super);

    let mut nw = ByteWriter::with_capacity(16);
    nw.write_var_u64(total_nov);
    nw.write_var_u64(nov_comp.len() as u64);
    out_file.write_all(&nw.buf).map_err(|e| e.to_string())?;
    out_file.write_all(&nov_comp).map_err(|e| e.to_string())?;
    out_file.flush().map_err(|e| e.to_string())?;

    total_written_bytes += (nw.buf.len() + nov_comp.len()) as u64;
    report.archive_bytes = total_written_bytes;
    report.seconds = t0.elapsed().as_secs_f64();
    Ok(report)
}

pub fn decompress_file(archive: &Path, output: &Path, panel: Option<&Path>, threads: i32, bounded: bool) -> Result<f64, String> {
    let t0 = Instant::now();
    let bytes = fs::read(archive).map_err(|e| e.to_string())?;
    let mut inp = ByteReader::new(&bytes);
    let magic = inp.read_bytes(4)?;
    if magic != MAGIC {
        return Err("Not an HGS2 archive".into());
    }
    let flags = inp.read_u8()?;
    let mode = flags & 3;
    let raw_flat_acgt = (flags & FLAG_RAW_FLAT_ACGT) != 0;
    let raw_rna_u = (flags & FLAG_RAW_RNA_U) != 0;
    let raw_len = inp.read_u64_le()?;
    let expected_crc = inp.read_u32_le()?;
    let _total_recs = inp.read_var_u64()?;
    let num_units = inp.read_var_u64()? as usize;

    let mut panel_seqs = Vec::new();
    if mode == MODE_REF {
        let pcnt = inp.read_var_u64()? as usize;
        let cn_len = inp.read_var_u64()? as usize;
        let cn = inp.read_bytes(cn_len)?;
        let pnames = decompress_panel_names(cn, pcnt)?;
        let mut psha = [0u8; 32];
        psha.copy_from_slice(inp.read_bytes(32)?);
        if pcnt > 0 {
            panel_seqs = verify_and_load_panel_for_decode(panel, &pnames, &psha)?;
        }
    }

    if raw_len == 0 || num_units == 0 {
        if expected_crc != 0 { return Err("CRC mismatch on empty archive".into()); }
        File::create(output).map_err(|e| e.to_string())?;
        return Ok(t0.elapsed().as_secs_f64());
    }

    let mut rank_to_unit = vec![0usize; num_units];
    let mut max_anc_rank = vec![0usize; num_units];
    if mode == MODE_COLLECTION && num_units > 1 {
        for slot in rank_to_unit.iter_mut().take(num_units) { *slot = inp.read_var_u64()? as usize; }
        for slot in max_anc_rank.iter_mut().take(num_units) { *slot = inp.read_var_u64()? as usize; }
    } else {
        for (u, slot) in rank_to_unit.iter_mut().enumerate().take(num_units) { *slot = u; }
    }

    let mut unit_meta_slices: Vec<&[u8]> = Vec::new();
    if !raw_flat_acgt {
        for _ in 0..num_units {
            let _frec = inp.read_var_u64()?;
            let _rcnt = inp.read_var_u64()?;
            let mlen = inp.read_var_u64()? as usize;
            unit_meta_slices.push(inp.read_bytes(mlen)?);
        }
    } else {
        for _ in 0..num_units {
            let _dlen = inp.read_var_u64()?;
        }
    }

    let mut mosaic_slices = Vec::with_capacity(num_units);
    for _ in 0..num_units {
        let mlen = inp.read_var_u64()? as usize;
        mosaic_slices.push(inp.read_bytes(mlen)?);
    }

    let total_nov = inp.read_var_u64()?;
    let nov_comp_len = inp.read_var_u64()? as usize;
    let nov_comp = inp.read_bytes(nov_comp_len)?;
    let nov_super = ffi::novelty_decode(nov_comp, total_nov, threads, bounded)?;

    let mut unit_to_rank = vec![0usize; num_units];
    let mut decoded_ranks: Vec<Vec<u8>> = Vec::with_capacity(num_units);
    let mut nov_off = 0usize;

    if mode == MODE_COLLECTION && num_units > 1 {
        for r in 0..num_units {
            unit_to_rank[rank_to_unit[r]] = r;
            let (_tlen, nlen) = ffi::mosaic_peek_lengths(mosaic_slices[r])?;
            let nslice = &nov_super[nov_off..nov_off + (nlen as usize)];
            nov_off += nlen as usize;
            let anc = max_anc_rank[r];
            let donor_refs: Vec<&[u8]> = decoded_ranks[..anc].iter().map(|v| v.as_slice()).collect();
            let dna = ffi::decode_mosaic(mosaic_slices[r], nslice, &donor_refs)?;
            decoded_ranks.push(dna);
        }
    } else {
        let panel_refs: Vec<&[u8]> = panel_seqs.iter().map(|v| v.as_slice()).collect();
        for u in 0..num_units {
            unit_to_rank[u] = u;
            let (_tlen, nlen) = ffi::mosaic_peek_lengths(mosaic_slices[u])?;
            let nslice = &nov_super[nov_off..nov_off + (nlen as usize)];
            nov_off += nlen as usize;
            let dna = ffi::decode_mosaic(mosaic_slices[u], nslice, &panel_refs)?;
            decoded_ranks.push(dna);
        }
    }
    drop(nov_super);
    drop(panel_seqs);

    let mut out_f = BufWriter::with_capacity(65536, File::create(output).map_err(|e| e.to_string())?);
    let mut crc = 0u32;
    let mut written = 0u64;

    if raw_flat_acgt {
        let map = if raw_rna_u { *b"ACGU" } else { *b"ACGT" };
        let mut ascii_buf = Vec::new();
        for &r in unit_to_rank.iter().take(num_units) {
            let dna = &decoded_ranks[r];
            ascii_buf.resize(dna.len(), 0);
            for (i, &b) in dna.iter().enumerate() { ascii_buf[i] = map[(b & 3) as usize]; }
            out_f.write_all(&ascii_buf).map_err(|e| e.to_string())?;
            crc = ffi::crc32_ieee(crc, &ascii_buf);
            written += ascii_buf.len() as u64;
        }
    } else {
        for u in 0..num_units {
            let r = unit_to_rank[u];
            let mut expected_dna = 0usize;
            let chunk: FactorizedChunk = fasta::read_aux_metadata(&mut ByteReader::new(unit_meta_slices[u]), &mut expected_dna)?;
            if decoded_ranks[r].len() != expected_dna {
                return Err(format!("Unit {u} DNA length mismatch: {} vs {expected_dna}", decoded_ranks[r].len()));
            }
            written += fasta::reconstruct_to(&chunk, &decoded_ranks[r], &mut out_f, &mut crc).map_err(|e| e.to_string())?;
        }
    }
    out_f.flush().map_err(|e| e.to_string())?;

    if written != raw_len {
        return Err(format!("Decompressed size mismatch: {written} vs {raw_len}"));
    }
    if crc != expected_crc {
        return Err(format!("CRC-32 mismatch: 0x{crc:08X} vs 0x{expected_crc:08X}"));
    }
    Ok(t0.elapsed().as_secs_f64())
}

pub fn inspect_file(archive: &Path) -> Result<ArchiveInfo, String> {
    let bytes = fs::read(archive).map_err(|e| e.to_string())?;
    let mut inp = ByteReader::new(&bytes);
    if inp.read_bytes(4)? != MAGIC { return Err("Not an HGS2 archive".into()); }
    let flags = inp.read_u8()?;
    let mode = flags & 3;
    let raw_flat_acgt = (flags & FLAG_RAW_FLAT_ACGT) != 0;
    let raw_rna_u = (flags & FLAG_RAW_RNA_U) != 0;
    let raw_bytes = inp.read_u64_le()?;
    let crc32 = inp.read_u32_le()?;
    let total_records = inp.read_var_u64()?;
    let num_units = inp.read_var_u64()? as usize;
    if mode == MODE_REF {
        let _pcnt = inp.read_var_u64()?;
        let cn_len = inp.read_var_u64()? as usize;
        let _ = inp.read_bytes(cn_len)?;
        let _ = inp.read_bytes(32)?;
    }
    if num_units > 0 {
        if mode == MODE_COLLECTION && num_units > 1 {
            for _ in 0..num_units { let _ = inp.read_var_u64()?; }
            for _ in 0..num_units { let _ = inp.read_var_u64()?; }
        }
        if !raw_flat_acgt {
            for _ in 0..num_units {
                let _ = inp.read_var_u64()?;
                let _ = inp.read_var_u64()?;
                let ml = inp.read_var_u64()? as usize;
                let _ = inp.read_bytes(ml)?;
            }
        } else {
            for _ in 0..num_units { let _ = inp.read_var_u64()?; }
        }
        for _ in 0..num_units {
            let ml = inp.read_var_u64()? as usize;
            let _ = inp.read_bytes(ml)?;
        }
        let total_novelty_bases = inp.read_var_u64()?;
        return Ok(ArchiveInfo {
            mode, raw_flat_acgt, raw_rna_u, raw_bytes, crc32, total_records, num_units, total_novelty_bases,
        });
    }
    Ok(ArchiveInfo {
        mode, raw_flat_acgt, raw_rna_u, raw_bytes, crc32, total_records, num_units, total_novelty_bases: 0,
    })
}
