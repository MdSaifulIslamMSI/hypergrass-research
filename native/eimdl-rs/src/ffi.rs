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
    fn hgs_ans_qual_bwt_block_cap() -> usize;
    fn hgs_ans_qual_bwt_peek_block_len(block_payload: *const u8, avail: usize) -> usize;
    fn hgs_ans_qual_bwt_encode_block(
        blk: *const u8,
        bn: usize,
        out_buf: *mut *mut u8,
        out_len: *mut usize,
    ) -> i32;
    fn hgs_ans_qual_bwt_decode_block(
        block_payload: *const u8,
        block_avail: usize,
        expected_bn: usize,
        out_consumed: *mut usize,
        out_blk: *mut u8,
    ) -> i32;
    fn hgs_ans4_encode_qual_ctx_lane(
        quals_raw: *const u8,
        mismatch_mask: *const u8,
        homo_bins: *const u8,
        total_len: usize,
        read_len_hint: u32,
        lane: usize,
        out_state: *mut u32,
        out_fwd_buf: *mut *mut u8,
        out_fwd_len: *mut usize,
    ) -> i32;
    fn hgs_ans4_decode_qual_ctx_lane(
        lane_payload: *const u8,
        lane_payload_len: usize,
        init_state: u32,
        mismatch_mask: *const u8,
        homo_bins: *const u8,
        total_len: usize,
        read_len_hint: u32,
        lane: usize,
        out_lane_buf: *mut u8,
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

enum QualEncOut {
    Bwt(Vec<u8>),
    Ctx(u32, Vec<u8>),
}

pub fn ans4_encode_qual_mt(
    quals_raw: &[u8],
    mismatch_mask: &[u8],
    homo_bins: &[u8],
    read_len_hint: u32,
    threads: i32,
) -> Result<Vec<u8>, &'static str> {
    if quals_raw.len() != mismatch_mask.len() || quals_raw.len() != homo_bins.len() {
        return Err("quality conditioning array length mismatch");
    }
    if threads <= 1 || quals_raw.len() < 4096 {
        return ans4_encode_qual(quals_raw, mismatch_mask, homo_bins, read_len_hint);
    }

    let n1: usize = mismatch_mask.iter().map(|&m| (m & 1) as usize).sum();
    let n0 = quals_raw.len() - n1;
    let split = n0 >= 16 && n1 >= 16;

    let mut s0 = Vec::new();
    let mut s1 = Vec::new();
    if split {
        s0.reserve_exact(n0);
        s1.reserve_exact(n1);
        for (&q, &m) in quals_raw.iter().zip(mismatch_mask.iter()) {
            if (m & 1) != 0 {
                s1.push(q);
            } else {
                s0.push(q);
            }
        }
    }

    let blk_cap = unsafe { hgs_ans_qual_bwt_block_cap() }.max(1);
    let mut bwt_chunks: Vec<&[u8]> = Vec::new();
    if split {
        bwt_chunks.extend(s0.chunks(blk_cap));
        bwt_chunks.extend(s1.chunks(blk_cap));
    } else {
        bwt_chunks.extend(quals_raw.chunks(blk_cap));
    }

    let nbwt = bwt_chunks.len();
    let total_jobs = nbwt + 4;

    let mut bwt_done: Vec<Option<Vec<u8>>> = vec![None; nbwt];
    let mut ctx_done: Vec<Option<(u32, Vec<u8>)>> = vec![None; 4];

    let n_workers = (threads as usize).clamp(1, total_jobs);
    let next_idx = std::sync::atomic::AtomicUsize::new(0);
    std::thread::scope(|sc| {
        let mut handles = Vec::with_capacity(n_workers);
        for _ in 0..n_workers {
            handles.push(sc.spawn(|| -> Result<Vec<(usize, QualEncOut)>, &'static str> {
                let mut local = Vec::new();
                loop {
                    let j = next_idx.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if j >= total_jobs {
                        break;
                    }
                    if j < nbwt {
                        let ch = bwt_chunks[j];
                        let mut ptr = ptr::null_mut();
                        let mut len = 0usize;
                        let rc = unsafe {
                            hgs_ans_qual_bwt_encode_block(ch.as_ptr(), ch.len(), &mut ptr, &mut len)
                        };
                        if rc != 0 {
                            return Err("hgs_ans_qual_bwt_encode_block failed");
                        }
                        local.push((j, QualEncOut::Bwt(take_ans_bytes(ptr, len))));
                    } else {
                        let lane = j - nbwt;
                        let mut st = 0u32;
                        let mut ptr = ptr::null_mut();
                        let mut len = 0usize;
                        let rc = unsafe {
                            hgs_ans4_encode_qual_ctx_lane(
                                quals_raw.as_ptr(),
                                mismatch_mask.as_ptr(),
                                homo_bins.as_ptr(),
                                quals_raw.len(),
                                read_len_hint,
                                lane,
                                &mut st,
                                &mut ptr,
                                &mut len,
                            )
                        };
                        if rc != 0 {
                            return Err("hgs_ans4_encode_qual_ctx_lane failed");
                        }
                        local.push((j, QualEncOut::Ctx(st, take_ans_bytes(ptr, len))));
                    }
                }
                Ok(local)
            }));
        }
        for h in handles {
            let chunk_res = h.join().map_err(|_| "ans4_encode_qual_mt worker panicked")??;
            for (idx, item) in chunk_res {
                match item {
                    QualEncOut::Bwt(v) => bwt_done[idx] = Some(v),
                    QualEncOut::Ctx(st, v) => ctx_done[idx - nbwt] = Some((st, v)),
                }
            }
        }
        Ok::<(), &'static str>(())
    })?;

    let mut bwt_len = 2usize;
    for b in &bwt_done {
        bwt_len = bwt_len.saturating_add(b.as_ref().ok_or("missing BWT block")?.len());
    }
    let raw_len = 1usize + quals_raw.len();
    let mut ctx_len = 33usize;
    for c in &ctx_done {
        ctx_len = ctx_len.saturating_add(c.as_ref().ok_or("missing CTX lane")?.1.len());
    }

    if bwt_len < ctx_len && bwt_len < raw_len {
        let mut out = Vec::with_capacity(bwt_len);
        out.push(4u8);
        out.push(u8::from(split));
        for b in bwt_done.into_iter().flatten() {
            out.extend_from_slice(&b);
        }
        return Ok(out);
    }
    if ctx_len < raw_len {
        let mut out = Vec::with_capacity(ctx_len);
        out.push(3u8);
        for c in &ctx_done {
            let st = c.as_ref().ok_or("missing CTX state")?.0;
            out.extend_from_slice(&st.to_le_bytes());
        }
        for c in &ctx_done {
            let fl = c.as_ref().ok_or("missing CTX len")?.1.len() as u32;
            out.extend_from_slice(&fl.to_le_bytes());
        }
        for (_, buf) in ctx_done.into_iter().flatten() {
            out.extend_from_slice(&buf);
        }
        return Ok(out);
    }
    let mut out = Vec::with_capacity(raw_len);
    out.push(2u8);
    out.extend_from_slice(quals_raw);
    Ok(out)
}

pub fn ans4_decode_qual_mt(
    payload: &[u8],
    mismatch_mask: &[u8],
    homo_bins: &[u8],
    expected_len: usize,
    read_len_hint: u32,
    threads: i32,
) -> Result<Vec<u8>, &'static str> {
    if mismatch_mask.len() != expected_len || homo_bins.len() != expected_len {
        return Err("quality decode mask length mismatch");
    }
    if threads <= 1 || expected_len < 4096 || payload.is_empty() {
        return ans4_decode_qual(payload, mismatch_mask, homo_bins, expected_len, read_len_hint);
    }

    match payload[0] {
        2 => {
            if payload.len() < 1 + expected_len {
                return Err("hgs_ans4_decode_qual raw truncated");
            }
            Ok(payload[1..1 + expected_len].to_vec())
        }
        4 => {
            if payload.len() < 2 {
                return Err("hgs_ans4_decode_qual BWT truncated");
            }
            let split = payload[1] != 0;
            let blk_cap = unsafe { hgs_ans_qual_bwt_block_cap() }.max(1);
            let mut pos = 2usize;

            let decode_chunks = |jobs: Vec<(&[u8], &mut [u8])>| -> Result<(), &'static str> {
                if jobs.is_empty() {
                    return Ok(());
                }
                let nw = (threads as usize).clamp(1, jobs.len());
                let iter_lock = std::sync::Mutex::new(jobs.into_iter());
                std::thread::scope(|sc| {
                    let mut handles = Vec::with_capacity(nw);
                    for _ in 0..nw {
                        handles.push(sc.spawn(|| -> Result<(), &'static str> {
                            loop {
                                let next = { iter_lock.lock().map_err(|_| "lock poisoned")?.next() };
                                let Some((in_sl, out_sl)) = next else { break };
                                let mut cons = 0usize;
                                let rc = unsafe {
                                    hgs_ans_qual_bwt_decode_block(
                                        in_sl.as_ptr(),
                                        in_sl.len(),
                                        out_sl.len(),
                                        &mut cons,
                                        out_sl.as_mut_ptr(),
                                    )
                                };
                                if rc != 0 || cons != in_sl.len() {
                                    return Err("hgs_ans_qual_bwt_decode_block failed");
                                }
                            }
                            Ok(())
                        }));
                    }
                    for h in handles {
                        h.join().map_err(|_| "ans4_decode_qual_mt worker panicked")??;
                    }
                    Ok(())
                })
            };

            if !split {
                let mut out = vec![0u8; expected_len];
                let mut jobs = Vec::new();
                for ch in out.chunks_mut(blk_cap) {
                    let blen = unsafe { hgs_ans_qual_bwt_peek_block_len(payload[pos..].as_ptr(), payload.len() - pos) };
                    if blen == 0 || pos + blen > payload.len() {
                        return Err("invalid BWT quality block length");
                    }
                    jobs.push((&payload[pos..pos + blen], ch));
                    pos += blen;
                }
                if pos != payload.len() {
                    return Err("trailing bytes in BWT quality payload");
                }
                decode_chunks(jobs)?;
                Ok(out)
            } else {
                let n1: usize = mismatch_mask.iter().map(|&m| (m & 1) as usize).sum();
                let n0 = expected_len - n1;
                let mut s0 = vec![0u8; n0];
                let mut s1 = vec![0u8; n1];
                let mut jobs = Vec::new();
                for ch in s0.chunks_mut(blk_cap) {
                    let blen = unsafe { hgs_ans_qual_bwt_peek_block_len(payload[pos..].as_ptr(), payload.len() - pos) };
                    if blen == 0 || pos + blen > payload.len() {
                        return Err("invalid BWT quality s0 block length");
                    }
                    jobs.push((&payload[pos..pos + blen], ch));
                    pos += blen;
                }
                for ch in s1.chunks_mut(blk_cap) {
                    let blen = unsafe { hgs_ans_qual_bwt_peek_block_len(payload[pos..].as_ptr(), payload.len() - pos) };
                    if blen == 0 || pos + blen > payload.len() {
                        return Err("invalid BWT quality s1 block length");
                    }
                    jobs.push((&payload[pos..pos + blen], ch));
                    pos += blen;
                }
                if pos != payload.len() {
                    return Err("trailing bytes in split BWT quality payload");
                }
                decode_chunks(jobs)?;
                let mut out = Vec::with_capacity(expected_len);
                let mut i0 = 0usize;
                let mut i1 = 0usize;
                for &m in mismatch_mask {
                    if (m & 1) != 0 {
                        out.push(s1[i1]);
                        i1 += 1;
                    } else {
                        out.push(s0[i0]);
                        i0 += 1;
                    }
                }
                Ok(out)
            }
        }
        3 => {
            if payload.len() < 33 {
                return Err("hgs_ans4_decode_qual ctx truncated");
            }
            let mut states = [0u32; 4];
            let mut blens = [0usize; 4];
            for l in 0..4 {
                states[l] = u32::from_le_bytes(payload[1 + l * 4..5 + l * 4].try_into().unwrap());
                blens[l] = u32::from_le_bytes(payload[17 + l * 4..21 + l * 4].try_into().unwrap()) as usize;
            }
            let mut pos = 33usize;
            let mut in_lanes: [&[u8]; 4] = [&[]; 4];
            for l in 0..4 {
                if pos + blens[l] > payload.len() {
                    return Err("hgs_ans4_decode_qual lane out of bounds");
                }
                in_lanes[l] = &payload[pos..pos + blens[l]];
                pos += blens[l];
            }
            let base_q = expected_len >> 2;
            let rem = expected_len & 3;
            let ll0 = base_q + usize::from(0 < rem);
            let ll1 = base_q + usize::from(1 < rem);
            let ll2 = base_q + usize::from(2 < rem);

            let mut out = vec![0u8; expected_len];
            let (o0, r1) = out.split_at_mut(ll0);
            let (o1, r2) = r1.split_at_mut(ll1);
            let (o2, o3) = r2.split_at_mut(ll2);
            let out_lanes: [(usize, &mut [u8]); 4] = [(0, o0), (1, o1), (2, o2), (3, o3)];

            std::thread::scope(|sc| {
                let mut handles = Vec::with_capacity(4);
                for (lane, out_buf) in out_lanes {
                    let inp = in_lanes[lane];
                    let st = states[lane];
                    handles.push(sc.spawn(move || -> Result<(), &'static str> {
                        let rc = unsafe {
                            hgs_ans4_decode_qual_ctx_lane(
                                inp.as_ptr(),
                                inp.len(),
                                st,
                                mismatch_mask.as_ptr(),
                                homo_bins.as_ptr(),
                                expected_len,
                                read_len_hint,
                                lane,
                                out_buf.as_mut_ptr(),
                            )
                        };
                        if rc != 0 {
                            return Err("hgs_ans4_decode_qual_ctx_lane failed");
                        }
                        Ok(())
                    }));
                }
                for h in handles {
                    h.join().map_err(|_| "ans4 ctx decode panicked")??;
                }
                Ok::<(), &'static str>(())
            })?;
            Ok(out)
        }
        _ => ans4_decode_qual(payload, mismatch_mask, homo_bins, expected_len, read_len_hint),
    }
}
