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
    let vec = if owned.len == 0 || owned.data.is_null() {
        Vec::new()
    } else {
        unsafe { slice::from_raw_parts(owned.data, owned.len) }.to_vec()
    };
    unsafe { hgs_cpp_free_buffer(owned) };
    Ok(vec)
}
