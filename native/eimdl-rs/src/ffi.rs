#![allow(dead_code)]

use std::ptr;
use std::slice;

#[repr(C)]
pub struct HgsOwnedBytes {
    pub data: *mut u8,
    pub len: usize,
}

#[repr(C)]
pub struct HgsMosaicResult {
    pub mosaic_data: *mut u8,
    pub mosaic_len: usize,
    pub novelty_data: *mut u8,
    pub novelty_len: usize,
    pub panel_match_bases: u64,
    pub self_match_bases: u64,
    pub tandem_match_bases: u64,
    pub panel_tokens: u64,
    pub self_tokens: u64,
    pub tandem_tokens: u64,
    pub literal_runs: u64,
    pub indel_splices: u64,
    pub substitutions: u64,
    pub status: i32,
}

#[repr(C)]
#[derive(Clone, Copy, Default, Debug)]
pub struct HgsMosaicStats {
    pub target_len: u64,
    pub novelty_len: u64,
    pub panel_match_bases: u64,
    pub self_match_bases: u64,
    pub tandem_match_bases: u64,
    pub literal_bases: u64,
    pub panel_tokens: u64,
    pub self_tokens: u64,
    pub tandem_tokens: u64,
    pub literal_runs: u64,
    pub snps: u64,
    pub transitions: u64,
    pub transversions: u64,
    pub wobble_mod3_hits: u64,
    pub splices: u64,
    pub inserted_bases: u64,
    pub deleted_bases: u64,
    pub recombination_switches: u64,
}

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct HgsVariantRecord {
    pub kind: u8,
    pub strand: u8,
    pub is_transition: u8,
    pub prev_mod3: u8,
    pub donor_id: u32,
    pub ref_pos: u64,
    pub tgt_pos: u64,
    pub ref_len: u32,
    pub alt_len: u32,
    pub ref_seq_off: u64,
    pub alt_seq_off: u64,
}

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct HgsSpanEvent {
    pub kind: u8,
    pub strand: u8,
    pub has_edits: u8,
    pub reserved: u8,
    pub donor_id: u32,
    pub ref_pos: u64,
    pub tgt_pos: u64,
    pub tgt_len: u64,
    pub ref_len: u64,
}

#[repr(C)]
pub struct HgsMosaicInspectResult {
    pub stats: HgsMosaicStats,
    pub variants: *mut HgsVariantRecord,
    pub num_variants: usize,
    pub spans: *mut HgsSpanEvent,
    pub num_spans: usize,
    pub var_seq_pool: *mut u8,
    pub var_seq_pool_len: usize,
    pub status: i32,
}

#[repr(C)]
pub struct HgsSyncmerIndex64 {
    _private: [u8; 0],
}

unsafe extern "C" {
    // C23 SIMD & Crypto Kernels (native/c_kernels/hgs_simd.c)
    fn hgs_cpu_has_avx2() -> i32;
    fn hgs_ascii_to_2bit_checked(src: *const u8, dst: *mut u8, len: usize, fourth_char: u8) -> i32;
    fn hgs_unpack_2bit_ascii(words: *const u64, offset: u64, len: usize, fourth_char: u8, out_ascii: *mut u8);
    fn hgs_crc32_ieee(prev_crc: u32, data: *const u8, len: usize) -> u32;
    fn hgs_sha256(data: *const u8, len: usize, out_digest: *mut u8);

    // C23 Deflate / Inflate Kernels (native/c_kernels/hgs_deflate.c)
    fn hgs_deflate_free(ptr: *mut u8);
    fn hgs_deflate_raw(src: *const u8, src_len: usize, out_buf: *mut *mut u8, out_len: *mut usize) -> i32;
    fn hgs_inflate_raw(src: *const u8, src_len: usize, hint: usize, out_buf: *mut *mut u8, out_len: *mut usize) -> i32;
    fn hgs_zlib_compress(src: *const u8, src_len: usize, out_buf: *mut *mut u8, out_len: *mut usize) -> i32;
    fn hgs_zlib_decompress(src: *const u8, src_len: usize, out_buf: *mut *mut u8, out_len: *mut usize) -> i32;

    // C++23 Syncmer Graph, Mosaic & Novelty Engines (native/cpp_engine/hgs_cpp_engine.cpp)
    fn hgs_cpp_free_buffer(buf: HgsOwnedBytes);
    fn hgs_cpp_compute_bottom_sketch(two_bit: *const u8, len: usize, sketch_size: usize, out_sk: *mut *mut u64, out_len: *mut usize) -> i32;
    fn hgs_cpp_free_u64(ptr: *mut u64);
    fn hgs_cpp_sketch_jaccard(a: *const u64, a_len: usize, b: *const u64, b_len: usize) -> f64;

    fn hgs_cpp_index_build(
        donor_ptrs: *const *const u8,
        donor_lens: *const u64,
        num_donors: usize,
        subsample_shift: i32,
        build_secondary: i32,
        threads: i32,
    ) -> *mut HgsSyncmerIndex64;
    fn hgs_cpp_index_free(idx: *mut HgsSyncmerIndex64);

    fn hgs_cpp_free_mosaic_result(res: *mut HgsMosaicResult);
    fn hgs_cpp_encode_mosaic(
        target: *const u8,
        tgt_len: u64,
        panel_idx: *const HgsSyncmerIndex64,
        max_donor_excl: usize,
        preferred_twin: i64,
        sim_order: *const i32,
        sim_order_len: usize,
        allow_self: i32,
        allow_tandem: i32,
        allow_approx: i32,
        allow_splice: i32,
        allow_kimura: i32,
        allow_wobble: i32,
        threads: i32,
    ) -> HgsMosaicResult;

    fn hgs_cpp_decode_mosaic(
        mosaic_bytes: *const u8,
        mosaic_len: usize,
        novelty_bytes: *const u8,
        novelty_len: usize,
        donor_ptrs: *const *const u8,
        donor_lens: *const u64,
        max_donor_excl: usize,
        out_two_bit: *mut *mut u8,
        out_len: *mut u64,
    ) -> i32;

    fn hgs_cpp_mosaic_peek_lengths(mosaic: *const u8, mosaic_len: usize, out_tgt: *mut u64, out_nov: *mut u64) -> i32;
    fn hgs_cpp_novelty_encode(bases: *const u8, len: u64, threads: i32, bounded: i32, out: *mut HgsOwnedBytes) -> i32;
    fn hgs_cpp_novelty_decode(payload: *const u8, payload_len: usize, expected_len: u64, threads: i32, bounded: i32, out: *mut HgsOwnedBytes) -> i32;

    // C23 4-Stripe SIMD rANS & Conditioned Phred Quality Kernels (native/c_kernels/hgs_ans.c)
    fn hgs_ans_free(ptr: *mut u8);
    fn hgs_ans4_encode_dna(two_bit: *const u8, len: usize, out_buf: *mut *mut u8, out_len: *mut usize) -> i32;
    fn hgs_ans4_decode_dna(payload: *const u8, payload_len: usize, expected_len: usize, out_buf: *mut *mut u8, out_len: *mut usize) -> i32;
    fn hgs_ans_compute_homo_bins(two_bit: *const u8, len: usize, out_bins: *mut u8);
    fn hgs_ans4_encode_qual(
        quals_raw: *const u8,
        mismatch_mask: *const u8,
        homo_bins: *const u8,
        len: usize,
        read_len_hint: u32,
        out_buf: *mut *mut u8,
        out_len: *mut usize,
    ) -> i32;
    fn hgs_ans4_decode_qual(
        payload: *const u8,
        payload_len: usize,
        mismatch_mask: *const u8,
        homo_bins: *const u8,
        expected_len: usize,
        read_len_hint: u32,
        out_buf: *mut *mut u8,
        out_len: *mut usize,
    ) -> i32;

    // C++23 Compressed-Domain Mosaic Inspector & Parallel rANS Novelty Wrapper
    fn hgs_cpp_free_inspect_result(res: *mut HgsMosaicInspectResult);
    fn hgs_cpp_mosaic_inspect(
        mosaic_bytes: *const u8,
        mosaic_len: usize,
        novelty_bytes: *const u8,
        novelty_len: usize,
        donor_2bit_ptrs: *const *const u8,
        donor_lens: *const u64,
        max_donor_exclusive: usize,
        collect_variants: i32,
        collect_spans: i32,
    ) -> HgsMosaicInspectResult;
    fn hgs_cpp_ans_novelty_encode(two_bit_bases: *const u8, len: u64, threads: i32, out_payload: *mut HgsOwnedBytes) -> i32;
    fn hgs_cpp_ans_novelty_decode(payload: *const u8, payload_len: usize, expected_len: u64, threads: i32, out_two_bit: *mut HgsOwnedBytes) -> i32;
}

pub fn has_avx2() -> bool {
    unsafe { hgs_cpu_has_avx2() != 0 }
}

pub fn crc32_ieee(prev: u32, data: &[u8]) -> u32 {
    unsafe { hgs_crc32_ieee(prev, data.as_ptr(), data.len()) }
}

pub fn sha256(data: &[u8]) -> [u8; 32] {
    let mut out = [0u8; 32];
    unsafe { hgs_sha256(data.as_ptr(), data.len(), out.as_mut_ptr()) };
    out
}

pub fn ascii_to_2bit_checked(src: &[u8], dst: &mut [u8], fourth_char: u8) -> bool {
    assert_eq!(src.len(), dst.len());
    unsafe { hgs_ascii_to_2bit_checked(src.as_ptr(), dst.as_mut_ptr(), src.len(), fourth_char) != 0 }
}

pub fn deflate_raw(src: &[u8]) -> Result<Vec<u8>, &'static str> {
    let mut ptr = ptr::null_mut();
    let mut len = 0usize;
    let rc = unsafe { hgs_deflate_raw(src.as_ptr(), src.len(), &mut ptr, &mut len) };
    if rc != 0 {
        return Err("hgs_deflate_raw failed");
    }
    let vec = if len == 0 || ptr.is_null() { Vec::new() } else { unsafe { slice::from_raw_parts(ptr, len) }.to_vec() };
    unsafe { hgs_deflate_free(ptr) };
    Ok(vec)
}

pub fn inflate_raw(src: &[u8], hint: usize) -> Result<Vec<u8>, &'static str> {
    let mut ptr = ptr::null_mut();
    let mut len = 0usize;
    let rc = unsafe { hgs_inflate_raw(src.as_ptr(), src.len(), hint, &mut ptr, &mut len) };
    if rc != 0 {
        return Err("hgs_inflate_raw failed");
    }
    let vec = if len == 0 || ptr.is_null() { Vec::new() } else { unsafe { slice::from_raw_parts(ptr, len) }.to_vec() };
    unsafe { hgs_deflate_free(ptr) };
    Ok(vec)
}

pub fn zlib_compress(src: &[u8]) -> Result<Vec<u8>, &'static str> {
    let mut ptr = ptr::null_mut();
    let mut len = 0usize;
    let rc = unsafe { hgs_zlib_compress(src.as_ptr(), src.len(), &mut ptr, &mut len) };
    if rc != 0 {
        return Err("hgs_zlib_compress failed");
    }
    let vec = if len == 0 || ptr.is_null() { Vec::new() } else { unsafe { slice::from_raw_parts(ptr, len) }.to_vec() };
    unsafe { hgs_deflate_free(ptr) };
    Ok(vec)
}

pub fn zlib_decompress(src: &[u8]) -> Result<Vec<u8>, &'static str> {
    let mut ptr = ptr::null_mut();
    let mut len = 0usize;
    let rc = unsafe { hgs_zlib_decompress(src.as_ptr(), src.len(), &mut ptr, &mut len) };
    if rc != 0 {
        return Err("hgs_zlib_decompress failed");
    }
    let vec = if len == 0 || ptr.is_null() { Vec::new() } else { unsafe { slice::from_raw_parts(ptr, len) }.to_vec() };
    unsafe { hgs_deflate_free(ptr) };
    Ok(vec)
}

pub fn compute_bottom_sketch(two_bit: &[u8], sketch_size: usize) -> Vec<u64> {
    let mut ptr = ptr::null_mut();
    let mut len = 0usize;
    let rc = unsafe { hgs_cpp_compute_bottom_sketch(two_bit.as_ptr(), two_bit.len(), sketch_size, &mut ptr, &mut len) };
    if rc != 0 || ptr.is_null() {
        return Vec::new();
    }
    let vec = unsafe { slice::from_raw_parts(ptr, len) }.to_vec();
    unsafe { hgs_cpp_free_u64(ptr) };
    vec
}

pub fn sketch_jaccard(a: &[u64], b: &[u64]) -> f64 {
    unsafe { hgs_cpp_sketch_jaccard(a.as_ptr(), a.len(), b.as_ptr(), b.len()) }
}

pub struct SyncmerIndex {
    ptr: *mut HgsSyncmerIndex64,
}
unsafe impl Send for SyncmerIndex {}
unsafe impl Sync for SyncmerIndex {}

impl Drop for SyncmerIndex {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            unsafe { hgs_cpp_index_free(self.ptr) };
            self.ptr = ptr::null_mut();
        }
    }
}

impl SyncmerIndex {
    pub fn build(donors: &[&[u8]], subsample_shift: i32, build_secondary: bool, threads: i32) -> Self {
        let ptrs: Vec<*const u8> = donors.iter().map(|s| s.as_ptr()).collect();
        let lens: Vec<u64> = donors.iter().map(|s| s.len() as u64).collect();
        let ptr = unsafe {
            hgs_cpp_index_build(
                ptrs.as_ptr(),
                lens.as_ptr(),
                donors.len(),
                subsample_shift,
                if build_secondary { 1 } else { 0 },
                threads,
            )
        };
        assert!(!ptr.is_null(), "hgs_cpp_index_build returned null");
        Self { ptr }
    }
}

pub struct MosaicEncoded {
    pub mosaic_bytes: Vec<u8>,
    pub novelty_bases: Vec<u8>,
    pub panel_match_bases: u64,
    pub self_match_bases: u64,
    pub tandem_match_bases: u64,
    pub panel_tokens: u64,
    pub self_tokens: u64,
    pub tandem_tokens: u64,
    pub literal_runs: u64,
    pub indel_splices: u64,
    pub substitutions: u64,
}

#[derive(Clone, Copy)]
pub struct MosaicFlags {
    pub allow_self: bool,
    pub allow_tandem: bool,
    pub allow_approx: bool,
    pub allow_splice: bool,
    pub allow_kimura: bool,
    pub allow_wobble: bool,
}

impl Default for MosaicFlags {
    fn default() -> Self {
        Self {
            allow_self: true,
            allow_tandem: true,
            allow_approx: true,
            allow_splice: true,
            allow_kimura: true,
            allow_wobble: true,
        }
    }
}

pub fn encode_mosaic(
    target: &[u8],
    panel: Option<&SyncmerIndex>,
    max_donor_excl: usize,
    preferred_twin: i64,
    sim_order: &[i32],
    flags: MosaicFlags,
    threads: i32,
) -> Result<MosaicEncoded, &'static str> {
    let idx_ptr = panel.map_or(ptr::null(), |p| p.ptr as *const HgsSyncmerIndex64);
    let mut raw = unsafe {
        hgs_cpp_encode_mosaic(
            target.as_ptr(),
            target.len() as u64,
            idx_ptr,
            max_donor_excl,
            preferred_twin,
            if sim_order.is_empty() { ptr::null() } else { sim_order.as_ptr() },
            sim_order.len(),
            if flags.allow_self { 1 } else { 0 },
            if flags.allow_tandem { 1 } else { 0 },
            if flags.allow_approx { 1 } else { 0 },
            if flags.allow_splice { 1 } else { 0 },
            if flags.allow_kimura { 1 } else { 0 },
            if flags.allow_wobble { 1 } else { 0 },
            threads,
        )
    };
    if raw.status != 0 {
        unsafe { hgs_cpp_free_mosaic_result(&mut raw) };
        return Err("hgs_cpp_encode_mosaic failed");
    }
    let mosaic_bytes = if raw.mosaic_len == 0 {
        Vec::new()
    } else {
        unsafe { slice::from_raw_parts(raw.mosaic_data, raw.mosaic_len) }.to_vec()
    };
    let novelty_bases = if raw.novelty_len == 0 {
        Vec::new()
    } else {
        unsafe { slice::from_raw_parts(raw.novelty_data, raw.novelty_len) }.to_vec()
    };
    let out = MosaicEncoded {
        mosaic_bytes,
        novelty_bases,
        panel_match_bases: raw.panel_match_bases,
        self_match_bases: raw.self_match_bases,
        tandem_match_bases: raw.tandem_match_bases,
        panel_tokens: raw.panel_tokens,
        self_tokens: raw.self_tokens,
        tandem_tokens: raw.tandem_tokens,
        literal_runs: raw.literal_runs,
        indel_splices: raw.indel_splices,
        substitutions: raw.substitutions,
    };
    unsafe { hgs_cpp_free_mosaic_result(&mut raw) };
    Ok(out)
}

pub fn decode_mosaic(
    mosaic: &[u8],
    novelty: &[u8],
    donors: &[&[u8]],
) -> Result<Vec<u8>, &'static str> {
    let ptrs: Vec<*const u8> = donors.iter().map(|s| s.as_ptr()).collect();
    let lens: Vec<u64> = donors.iter().map(|s| s.len() as u64).collect();
    let mut out_ptr = ptr::null_mut();
    let mut out_len = 0u64;
    let rc = unsafe {
        hgs_cpp_decode_mosaic(
            mosaic.as_ptr(),
            mosaic.len(),
            novelty.as_ptr(),
            novelty.len(),
            if ptrs.is_empty() { ptr::null() } else { ptrs.as_ptr() },
            if lens.is_empty() { ptr::null() } else { lens.as_ptr() },
            donors.len(),
            &mut out_ptr,
            &mut out_len,
        )
    };
    if rc != 0 {
        return Err("hgs_cpp_decode_mosaic failed");
    }
    let vec = if out_len == 0 || out_ptr.is_null() {
        Vec::new()
    } else {
        unsafe { slice::from_raw_parts(out_ptr, out_len as usize) }.to_vec()
    };
    unsafe { hgs_cpp_free_buffer(HgsOwnedBytes { data: out_ptr, len: out_len as usize }) };
    Ok(vec)
}

pub fn mosaic_peek_lengths(mosaic: &[u8]) -> Result<(u64, u64), &'static str> {
    let mut tgt = 0u64;
    let mut nov = 0u64;
    let rc = unsafe { hgs_cpp_mosaic_peek_lengths(mosaic.as_ptr(), mosaic.len(), &mut tgt, &mut nov) };
    if rc != 0 {
        return Err("hgs_cpp_mosaic_peek_lengths failed");
    }
    Ok((tgt, nov))
}

pub fn novelty_encode(bases: &[u8], threads: i32, bounded: bool) -> Result<Vec<u8>, &'static str> {
    let mut owned = HgsOwnedBytes { data: ptr::null_mut(), len: 0 };
    let rc = unsafe {
        hgs_cpp_novelty_encode(
            bases.as_ptr(),
            bases.len() as u64,
            threads,
            if bounded { 1 } else { 0 },
            &mut owned,
        )
    };
    if rc != 0 {
        return Err("hgs_cpp_novelty_encode failed");
    }
    let vec = if owned.len == 0 || owned.data.is_null() {
        Vec::new()
    } else {
        unsafe { slice::from_raw_parts(owned.data, owned.len) }.to_vec()
    };
    unsafe { hgs_cpp_free_buffer(owned) };
    Ok(vec)
}

pub fn novelty_decode(payload: &[u8], expected_len: u64, threads: i32, bounded: bool) -> Result<Vec<u8>, &'static str> {
    let mut owned = HgsOwnedBytes { data: ptr::null_mut(), len: 0 };
    let rc = unsafe {
        hgs_cpp_novelty_decode(
            payload.as_ptr(),
            payload.len(),
            expected_len,
            threads,
            if bounded { 1 } else { 0 },
            &mut owned,
        )
    };
    if rc != 0 {
        return Err("hgs_cpp_novelty_decode failed");
    }
    Ok(take_owned_cpp_bytes(owned))
}

fn take_owned_cpp_bytes(owned: HgsOwnedBytes) -> Vec<u8> {
    let vec = if owned.len == 0 || owned.data.is_null() {
        Vec::new()
    } else {
        unsafe { slice::from_raw_parts(owned.data, owned.len) }.to_vec()
    };
    unsafe { hgs_cpp_free_buffer(owned) };
    vec
}

fn take_ans_bytes(ptr: *mut u8, len: usize) -> Vec<u8> {
    let vec = if len == 0 || ptr.is_null() {
        Vec::new()
    } else {
        unsafe { slice::from_raw_parts(ptr, len) }.to_vec()
    };
    unsafe { hgs_ans_free(ptr) };
    vec
}

pub struct MosaicInspectReport {
    pub stats: HgsMosaicStats,
    pub variants: Vec<HgsVariantRecord>,
    pub spans: Vec<HgsSpanEvent>,
    pub var_seq_pool: Vec<u8>,
}

pub fn mosaic_inspect(
    mosaic: &[u8],
    novelty: &[u8],
    donors: &[&[u8]],
    max_donor_exclusive: usize,
    collect_variants: bool,
    collect_spans: bool,
) -> Result<MosaicInspectReport, &'static str> {
    let limit = max_donor_exclusive.min(donors.len());
    let ptrs: Vec<*const u8> = donors[..limit].iter().map(|s| s.as_ptr()).collect();
    let lens: Vec<u64> = donors[..limit].iter().map(|s| s.len() as u64).collect();
    let mut raw = unsafe {
        hgs_cpp_mosaic_inspect(
            mosaic.as_ptr(),
            mosaic.len(),
            novelty.as_ptr(),
            novelty.len(),
            ptrs.as_ptr(),
            lens.as_ptr(),
            limit,
            if collect_variants { 1 } else { 0 },
            if collect_spans { 1 } else { 0 },
        )
    };
    if raw.status != 0 {
        unsafe { hgs_cpp_free_inspect_result(&mut raw) };
        return Err("hgs_cpp_mosaic_inspect failed");
    }
    let variants = if raw.num_variants == 0 || raw.variants.is_null() {
        Vec::new()
    } else {
        unsafe { slice::from_raw_parts(raw.variants, raw.num_variants) }.to_vec()
    };
    let spans = if raw.num_spans == 0 || raw.spans.is_null() {
        Vec::new()
    } else {
        unsafe { slice::from_raw_parts(raw.spans, raw.num_spans) }.to_vec()
    };
    let var_seq_pool = if raw.var_seq_pool_len == 0 || raw.var_seq_pool.is_null() {
        Vec::new()
    } else {
        unsafe { slice::from_raw_parts(raw.var_seq_pool, raw.var_seq_pool_len) }.to_vec()
    };
    let stats = raw.stats;
    unsafe { hgs_cpp_free_inspect_result(&mut raw) };
    Ok(MosaicInspectReport {
        stats,
        variants,
        spans,
        var_seq_pool,
    })
}

pub fn ans_novelty_encode(bases: &[u8], threads: i32) -> Result<Vec<u8>, &'static str> {
    let mut owned = HgsOwnedBytes { data: ptr::null_mut(), len: 0 };
    let rc = unsafe { hgs_cpp_ans_novelty_encode(bases.as_ptr(), bases.len() as u64, threads, &mut owned) };
    if rc != 0 {
        return Err("hgs_cpp_ans_novelty_encode failed");
    }
    Ok(take_owned_cpp_bytes(owned))
}

pub fn ans_novelty_decode(payload: &[u8], expected_len: u64, threads: i32) -> Result<Vec<u8>, &'static str> {
    let mut owned = HgsOwnedBytes { data: ptr::null_mut(), len: 0 };
    let rc = unsafe { hgs_cpp_ans_novelty_decode(payload.as_ptr(), payload.len(), expected_len, threads, &mut owned) };
    if rc != 0 {
        return Err("hgs_cpp_ans_novelty_decode failed");
    }
    Ok(take_owned_cpp_bytes(owned))
}

pub fn ans4_encode_dna(two_bit: &[u8]) -> Result<Vec<u8>, &'static str> {
    let mut ptr = ptr::null_mut();
    let mut len = 0usize;
    let rc = unsafe { hgs_ans4_encode_dna(two_bit.as_ptr(), two_bit.len(), &mut ptr, &mut len) };
    if rc != 0 {
        return Err("hgs_ans4_encode_dna failed");
    }
    Ok(take_ans_bytes(ptr, len))
}

pub fn ans4_decode_dna(payload: &[u8], expected_len: usize) -> Result<Vec<u8>, &'static str> {
    let mut ptr = ptr::null_mut();
    let mut len = 0usize;
    let rc = unsafe { hgs_ans4_decode_dna(payload.as_ptr(), payload.len(), expected_len, &mut ptr, &mut len) };
    if rc != 0 {
        return Err("hgs_ans4_decode_dna failed");
    }
    Ok(take_ans_bytes(ptr, len))
}

pub fn compute_homo_bins(two_bit: &[u8]) -> Vec<u8> {
    let mut bins = vec![0u8; two_bit.len()];
    if !two_bit.is_empty() {
        unsafe { hgs_ans_compute_homo_bins(two_bit.as_ptr(), two_bit.len(), bins.as_mut_ptr()) };
    }
    bins
}

pub fn ans4_encode_qual(
    quals_raw: &[u8],
    mismatch_mask: &[u8],
    homo_bins: &[u8],
    read_len_hint: u32,
) -> Result<Vec<u8>, &'static str> {
    if quals_raw.len() != mismatch_mask.len() || quals_raw.len() != homo_bins.len() {
        return Err("quality conditioning array length mismatch");
    }
    let mut ptr = ptr::null_mut();
    let mut len = 0usize;
    let rc = unsafe {
        hgs_ans4_encode_qual(
            quals_raw.as_ptr(),
            mismatch_mask.as_ptr(),
            homo_bins.as_ptr(),
            quals_raw.len(),
            read_len_hint,
            &mut ptr,
            &mut len,
        )
    };
    if rc != 0 {
        return Err("hgs_ans4_encode_qual failed");
    }
    Ok(take_ans_bytes(ptr, len))
}

pub fn ans4_decode_qual(
    payload: &[u8],
    mismatch_mask: &[u8],
    homo_bins: &[u8],
    expected_len: usize,
    read_len_hint: u32,
) -> Result<Vec<u8>, &'static str> {
    if mismatch_mask.len() != expected_len || homo_bins.len() != expected_len {
        return Err("quality decode mask length mismatch");
    }
    let mut ptr = ptr::null_mut();
    let mut len = 0usize;
    let rc = unsafe {
        hgs_ans4_decode_qual(
            payload.as_ptr(),
            payload.len(),
            mismatch_mask.as_ptr(),
            homo_bins.as_ptr(),
            expected_len,
            read_len_hint,
            &mut ptr,
            &mut len,
        )
    };
    if rc != 0 {
        return Err("hgs_ans4_decode_qual failed");
    }
    Ok(take_ans_bytes(ptr, len))
}
