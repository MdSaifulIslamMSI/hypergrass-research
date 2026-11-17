#include "hgs_simd.h"
#include <string.h>

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 202311L && defined(__has_include)
#if __has_include(<stdbit.h>)
#include <stdbit.h>
#define HGS_HAVE_C23_STDBIT 1
#endif
#endif

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 202311L
#define HGS_CONSTEXPR static constexpr
static_assert(sizeof(uint64_t) == 8);
#else
#define HGS_CONSTEXPR static const
#endif

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <cpuid.h>
#endif
#endif

static inline uint32_t hgs_ctz64(uint64_t x) {
#if defined(HGS_HAVE_C23_STDBIT)
    return (uint32_t)stdc_trailing_zeros(x);
#elif defined(__GNUC__) || defined(__clang__)
    return (uint32_t)__builtin_ctzll(x);
#elif defined(_MSC_VER)
    unsigned long idx;
    _BitScanForward64(&idx, x);
    return (uint32_t)idx;
#else
    uint32_t n = 0;
    if ((x & 0xFFFFFFFFULL) == 0) { n += 32; x >>= 32; }
    if ((x & 0xFFFFULL) == 0) { n += 16; x >>= 16; }
    if ((x & 0xFFULL) == 0) { n += 8; x >>= 8; }
    if ((x & 0xFULL) == 0) { n += 4; x >>= 4; }
    if ((x & 0x3ULL) == 0) { n += 2; x >>= 2; }
    if ((x & 0x1ULL) == 0) { n += 1; }
    return n;
#endif
}

int hgs_cpu_has_avx2(void) {
#if defined(__x86_64__) || defined(_M_X64)
    static int cached = -1;
    if (cached >= 0) return cached;
#if defined(_MSC_VER)
    int info[4] = {0};
    __cpuid(info, 0);
    if (info[0] < 7) { cached = 0; return 0; }
    __cpuidex(info, 7, 0);
    int osxsave = 0;
    int info1[4] = {0};
    __cpuid(info1, 1);
    osxsave = (info1[2] & (1 << 27)) && (info1[2] & (1 << 28));
    cached = (osxsave && ((info[1] & (1 << 5)) != 0)) ? 1 : 0;
#else
    unsigned int eax = 0, ebx = 0, ecx = 0, edx = 0;
    if (__get_cpuid_max(0, NULL) < 7) { cached = 0; return 0; }
    __cpuid(1, eax, ebx, ecx, edx);
    if ((ecx & (1U << 27)) == 0 || (ecx & (1U << 28)) == 0) { cached = 0; return 0; }
    __cpuid_count(7, 0, eax, ebx, ecx, edx);
    cached = ((ebx & (1U << 5)) != 0) ? 1 : 0;
#endif
    return cached;
#else
    return 0;
#endif
}

HGS_CONSTEXPR uint8_t ASCII_MAP_2BIT[256] = {
    ['A'] = 0, ['a'] = 0,
    ['C'] = 1, ['c'] = 1,
    ['G'] = 2, ['g'] = 2,
    ['T'] = 3, ['t'] = 3,
    ['U'] = 3, ['u'] = 3
};

#if (defined(__x86_64__) || defined(_M_X64)) && (defined(__GNUC__) || defined(__clang__))
__attribute__((target("avx2")))
static int hgs_ascii_to_2bit_avx2(const uint8_t* src, uint8_t* dst, size_t len, uint8_t fourth_char) {
    const __m256i v_A = _mm256_set1_epi8('A');
    const __m256i v_C = _mm256_set1_epi8('C');
    const __m256i v_G = _mm256_set1_epi8('G');
    const __m256i v_F = _mm256_set1_epi8((char)fourth_char);
    const __m256i c_1 = _mm256_set1_epi8(1);
    const __m256i c_2 = _mm256_set1_epi8(2);
    const __m256i c_3 = _mm256_set1_epi8(3);

    size_t i = 0;
    int all_strict_upper = 1;
    for (; i + 32 <= len; i += 32) {
        __m256i x = _mm256_loadu_si256((const __m256i*)(const void*)(src + i));
        __m256i mA = _mm256_cmpeq_epi8(x, v_A);
        __m256i mC = _mm256_cmpeq_epi8(x, v_C);
        __m256i mG = _mm256_cmpeq_epi8(x, v_G);
        __m256i mF = _mm256_cmpeq_epi8(x, v_F);
        __m256i valid = _mm256_or_si256(_mm256_or_si256(mA, mC), _mm256_or_si256(mG, mF));
        if ((uint32_t)_mm256_movemask_epi8(valid) != 0xFFFFFFFFU) {
            all_strict_upper = 0;
            for (size_t k = 0; k < 32; k++) {
                dst[i + k] = ASCII_MAP_2BIT[src[i + k]];
            }
        } else {
            __m256i code = _mm256_or_si256(
                _mm256_and_si256(mC, c_1),
                _mm256_or_si256(_mm256_and_si256(mG, c_2), _mm256_and_si256(mF, c_3)));
            _mm256_storeu_si256((__m256i*)(void*)(dst + i), code);
        }
    }
    for (; i < len; i++) {
        uint8_t ch = src[i];
        if (ch == 'A') dst[i] = 0;
        else if (ch == 'C') dst[i] = 1;
        else if (ch == 'G') dst[i] = 2;
        else if (ch == fourth_char) dst[i] = 3;
        else {
            all_strict_upper = 0;
            dst[i] = ASCII_MAP_2BIT[ch];
        }
    }
    return all_strict_upper;
}
#endif

int hgs_ascii_to_2bit_checked(const uint8_t* src, uint8_t* dst, size_t len, uint8_t fourth_char) {
#if (defined(__x86_64__) || defined(_M_X64)) && (defined(__GNUC__) || defined(__clang__))
    if (len >= 64 && hgs_cpu_has_avx2()) {
        return hgs_ascii_to_2bit_avx2(src, dst, len, fourth_char);
    }
#endif
    int all_strict_upper = 1;
    for (size_t i = 0; i < len; i++) {
        uint8_t ch = src[i];
        if (ch == 'A') dst[i] = 0;
        else if (ch == 'C') dst[i] = 1;
        else if (ch == 'G') dst[i] = 2;
        else if (ch == fourth_char) dst[i] = 3;
        else {
            all_strict_upper = 0;
            dst[i] = ASCII_MAP_2BIT[ch];
        }
    }
    return all_strict_upper;
}

void hgs_pack_2bit_words(const uint8_t* bases, size_t len, uint64_t* out_words) {
    size_t full_words = len >> 5;
    size_t idx = 0;
    for (size_t w = 0; w < full_words; w++) {
        uint64_t acc0 = 0, acc1 = 0;
        for (int j = 0; j < 16; j++) {
            acc0 |= ((uint64_t)(bases[idx + (size_t)j] & 3U)) << (j << 1);
            acc1 |= ((uint64_t)(bases[idx + 16 + (size_t)j] & 3U)) << ((j + 16) << 1);
        }
        out_words[w] = acc0 | acc1;
        idx += 32;
    }
    if (idx < len) {
        uint64_t acc = 0;
        int shift = 0;
        while (idx < len) {
            acc |= ((uint64_t)(bases[idx++] & 3U)) << shift;
            shift += 2;
        }
        out_words[full_words++] = acc;
    }
    out_words[full_words] = 0ULL;
}

void hgs_unpack_2bit_ascii(const uint64_t* words, uint64_t offset, size_t len, uint8_t fourth_char, uint8_t* out_ascii) {
    uint8_t lut[4] = {'A', 'C', 'G', fourth_char};
    size_t k = 0;
    while (k + 32 <= len) {
        uint64_t w = hgs_word32_forward(words, offset + k);
        for (int j = 0; j < 32; j++) {
            out_ascii[k + (size_t)j] = lut[w & 3ULL];
            w >>= 2;
        }
        k += 32;
    }
    while (k < len) {
        uint64_t pos = offset + k;
        uint64_t w = words[pos >> 5];
        uint32_t sh = (uint32_t)((pos & 31ULL) << 1);
        out_ascii[k] = lut[(w >> sh) & 3ULL];
        k++;
    }
}

void hgs_unpack_2bit_bytes(const uint64_t* words, uint64_t offset, size_t len, uint8_t* out_bytes) {
    size_t k = 0;
    while (k + 32 <= len) {
        uint64_t w = hgs_word32_forward(words, offset + k);
        for (int j = 0; j < 32; j++) {
            out_bytes[k + (size_t)j] = (uint8_t)(w & 3ULL);
            w >>= 2;
        }
        k += 32;
    }
    while (k < len) {
        uint64_t pos = offset + k;
        uint64_t w = words[pos >> 5];
        uint32_t sh = (uint32_t)((pos & 31ULL) << 1);
        out_bytes[k] = (uint8_t)((w >> sh) & 3ULL);
        k++;
    }
}

static inline uint8_t base_at(const uint64_t* words, uint64_t p) {
    return (uint8_t)((words[p >> 5] >> ((p & 31ULL) << 1)) & 3ULL);
}

uint64_t hgs_count_exact_match_2bit(
    const uint64_t* ref_words,
    uint64_t ref_offset,
    uint64_t ref_len,
    uint8_t strand,
    uint64_t ref_strand_pos,
    const uint64_t* tgt_words,
    uint64_t tgt_offset,
    uint64_t tgt_len,
    uint64_t tgt_pos,
    uint64_t max_len) {
    if (ref_strand_pos >= ref_len || tgt_pos >= tgt_len) return 0;
    uint64_t r_rem = ref_len - ref_strand_pos;
    uint64_t t_rem = tgt_len - tgt_pos;
    uint64_t limit = max_len;
    if (r_rem < limit) limit = r_rem;
    if (t_rem < limit) limit = t_rem;

    uint64_t k = 0;
    if (strand == 0) {
        uint64_t r_base = ref_offset + ref_strand_pos;
        uint64_t t_base = tgt_offset + tgt_pos;
        while (k + 32 <= limit) {
            uint64_t diff = hgs_word32_forward(ref_words, r_base + k) ^ hgs_word32_forward(tgt_words, t_base + k);
            if (diff != 0ULL) {
                return k + (hgs_ctz64(diff) >> 1);
            }
            k += 32;
        }
        while (k < limit) {
            if (base_at(ref_words, r_base + k) != base_at(tgt_words, t_base + k)) break;
            k++;
        }
        return k;
    } else {
        uint64_t fwd_top = ref_offset + ref_len - 1ULL - ref_strand_pos;
        uint64_t t_base = tgt_offset + tgt_pos;
        while (k + 32 <= limit) {
            uint64_t r_word = hgs_revcomp_word32(hgs_word32_forward(ref_words, fwd_top - k - 31ULL));
            uint64_t diff = r_word ^ hgs_word32_forward(tgt_words, t_base + k);
            if (diff != 0ULL) {
                return k + (hgs_ctz64(diff) >> 1);
            }
            k += 32;
        }
        while (k < limit) {
            uint8_t rb = (uint8_t)(base_at(ref_words, fwd_top - k) ^ 3U);
            uint8_t tb = base_at(tgt_words, t_base + k);
            if (rb != tb) break;
            k++;
        }
        return k;
    }
}

void hgs_scan_open_syncmers(
    const uint64_t* words,
    uint64_t start_offset,
    uint64_t seg_len,
    int k,
    int s,
    int t,
    HgsSyncmerCb cb,
    void* ctx) {
    if ((uint64_t)k > seg_len || k <= 0 || s <= 0 || k > 32) return;
    int w = k - s + 1;
    if (w <= 0 || w > 32 || t < 0 || t >= w) return;

    uint64_t k_mask = (k == 32) ? ~0ULL : ((1ULL << ((unsigned)k << 1)) - 1ULL);
    uint64_t s_mask = (1ULL << ((unsigned)s << 1)) - 1ULL;
    unsigned rc_k_shift = (unsigned)((k - 1) << 1);
    unsigned rc_s_shift = (unsigned)((s - 1) << 1);

    uint32_t win_hashes[32];
    uint64_t fwd_k = 0, rc_k = 0, fwd_s = 0, rc_s = 0;

    for (uint64_t rel_pos = 0; rel_pos < seg_len; rel_pos++) {
        uint64_t b = base_at(words, start_offset + rel_pos);
        uint64_t cb_val = b ^ 3ULL;

        fwd_k = ((fwd_k << 2) | b) & k_mask;
        rc_k = (rc_k >> 2) | (cb_val << rc_k_shift);
        fwd_s = ((fwd_s << 2) | b) & s_mask;
        rc_s = (rc_s >> 2) | (cb_val << rc_s_shift);

        if (rel_pos + 1 >= (uint64_t)s) {
            uint64_t s_idx = rel_pos + 1ULL - (uint64_t)s;
            uint64_t canon_s = (fwd_s <= rc_s) ? fwd_s : rc_s;
            win_hashes[s_idx % (uint64_t)w] = hgs_mix32(canon_s);
        }
        if (rel_pos + 1 >= (uint64_t)k) {
            uint64_t k_start = rel_pos + 1ULL - (uint64_t)k;
            uint32_t center = win_hashes[(k_start + (uint64_t)t) % (uint64_t)w];
            int ok = 1;
            for (int j = 0; j < w; j++) {
                if (j == t) continue;
                if (win_hashes[(k_start + (uint64_t)j) % (uint64_t)w] <= center) {
                    ok = 0;
                    break;
                }
            }
            if (ok) {
                uint8_t strand = (fwd_k <= rc_k) ? 0U : 1U;
                uint64_t canon_k = (strand == 0U) ? fwd_k : rc_k;
                cb(k_start, canon_k, strand, ctx);
            }
        }
    }
}

/* Slice-by-8 IEEE 802.3 CRC-32 */
static uint32_t crc_tables[8][256];
static int crc_init_done = 0;

static void init_crc32_tables(void) {
    if (crc_init_done) return;
    for (uint32_t i = 0; i < 256; i++) {
        uint32_t c = i;
        for (int j = 0; j < 8; j++) {
            c = (c & 1U) ? ((c >> 1) ^ 0xEDB88320U) : (c >> 1);
        }
        crc_tables[0][i] = c;
    }
    for (int t = 1; t < 8; t++) {
        for (uint32_t i = 0; i < 256; i++) {
            uint32_t c = crc_tables[t - 1][i];
            crc_tables[t][i] = (c >> 8) ^ crc_tables[0][c & 0xFFU];
        }
    }
    crc_init_done = 1;
}

uint32_t hgs_crc32_ieee(uint32_t prev_crc, const uint8_t* data, size_t len) {
    init_crc32_tables();
    uint32_t crc = ~prev_crc;
    size_t i = 0;
    while (i + 8 <= len) {
        uint32_t lo = (uint32_t)data[i] | ((uint32_t)data[i+1] << 8) | ((uint32_t)data[i+2] << 16) | ((uint32_t)data[i+3] << 24);
        uint32_t hi = (uint32_t)data[i+4] | ((uint32_t)data[i+5] << 8) | ((uint32_t)data[i+6] << 16) | ((uint32_t)data[i+7] << 24);
        lo ^= crc;
        crc = crc_tables[7][lo & 0xFFU] ^ crc_tables[6][(lo >> 8) & 0xFFU] ^
              crc_tables[5][(lo >> 16) & 0xFFU] ^ crc_tables[4][(lo >> 24) & 0xFFU] ^
              crc_tables[3][hi & 0xFFU] ^ crc_tables[2][(hi >> 8) & 0xFFU] ^
              crc_tables[1][(hi >> 16) & 0xFFU] ^ crc_tables[0][(hi >> 24) & 0xFFU];
        i += 8;
    }
    while (i < len) {
        crc = (crc >> 8) ^ crc_tables[0][(crc ^ data[i++]) & 0xFFU];
    }
    return ~crc;
}

/* FIPS 180-4 SHA-256 */
HGS_CONSTEXPR uint32_t SHA256_K[64] = {
    0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U, 0x3956c25bU, 0x59f111f1U, 0x923f82a4U, 0xab1c5ed5U,
    0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U, 0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U,
    0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU, 0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
    0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U, 0xc6e00bf3U, 0xd5a79147U, 0x06ca6351U, 0x14292967U,
    0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U, 0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U,
    0xa2bfe8a1U, 0xa81a664bU, 0xc24b8b70U, 0xc76c51a3U, 0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
    0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U, 0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU, 0x682e6ff3U,
    0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U, 0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U
};

static inline uint32_t rotr32(uint32_t x, uint32_t n) { return (x >> n) | (x << (32U - n)); }

static void sha256_transform(HgsSha256Ctx* ctx, const uint8_t block[64]) {
    uint32_t w[64];
    for (int i = 0; i < 16; i++) {
        w[i] = ((uint32_t)block[i*4] << 24) | ((uint32_t)block[i*4+1] << 16) |
               ((uint32_t)block[i*4+2] << 8) | (uint32_t)block[i*4+3];
    }
    for (int i = 16; i < 64; i++) {
        uint32_t s0 = rotr32(w[i-15], 7) ^ rotr32(w[i-15], 18) ^ (w[i-15] >> 3);
        uint32_t s1 = rotr32(w[i-2], 17) ^ rotr32(w[i-2], 19) ^ (w[i-2] >> 10);
        w[i] = w[i-16] + s0 + w[i-7] + s1;
    }
    uint32_t a = ctx->state[0], b = ctx->state[1], c = ctx->state[2], d = ctx->state[3];
    uint32_t e = ctx->state[4], f = ctx->state[5], g = ctx->state[6], h = ctx->state[7];
    for (int i = 0; i < 64; i++) {
        uint32_t S1 = rotr32(e, 6) ^ rotr32(e, 11) ^ rotr32(e, 25);
        uint32_t ch = (e & f) ^ ((~e) & g);
        uint32_t temp1 = h + S1 + ch + SHA256_K[i] + w[i];
        uint32_t S0 = rotr32(a, 2) ^ rotr32(a, 13) ^ rotr32(a, 22);
        uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
        uint32_t temp2 = S0 + maj;
        h = g; g = f; f = e; e = d + temp1;
        d = c; c = b; b = a; a = temp1 + temp2;
    }
    ctx->state[0] += a; ctx->state[1] += b; ctx->state[2] += c; ctx->state[3] += d;
    ctx->state[4] += e; ctx->state[5] += f; ctx->state[6] += g; ctx->state[7] += h;
}

void hgs_sha256_init(HgsSha256Ctx* ctx) {
    ctx->state[0] = 0x6a09e667U; ctx->state[1] = 0xbb67ae85U;
    ctx->state[2] = 0x3c6ef372U; ctx->state[3] = 0xa54ff53aU;
    ctx->state[4] = 0x510e527fU; ctx->state[5] = 0x9b05688cU;
    ctx->state[6] = 0x1f83d9abU; ctx->state[7] = 0x5be0cd19U;
    ctx->bitlen = 0;
    ctx->buf_len = 0;
}

void hgs_sha256_update(HgsSha256Ctx* ctx, const uint8_t* data, size_t len) {
    ctx->bitlen += (uint64_t)len * 8ULL;
    size_t i = 0;
    if (ctx->buf_len > 0) {
        size_t fill = 64 - ctx->buf_len;
        if (fill > len) fill = len;
        memcpy(ctx->buf + ctx->buf_len, data, fill);
        ctx->buf_len += fill;
        i += fill;
        if (ctx->buf_len == 64) {
            sha256_transform(ctx, ctx->buf);
            ctx->buf_len = 0;
        }
    }
    while (i + 64 <= len) {
        sha256_transform(ctx, data + i);
        i += 64;
    }
    if (i < len) {
        ctx->buf_len = len - i;
        memcpy(ctx->buf, data + i, ctx->buf_len);
    }
}

void hgs_sha256_final(HgsSha256Ctx* ctx, uint8_t out_digest[32]) {
    ctx->buf[ctx->buf_len++] = 0x80U;
    if (ctx->buf_len > 56) {
        memset(ctx->buf + ctx->buf_len, 0, 64 - ctx->buf_len);
        sha256_transform(ctx, ctx->buf);
        ctx->buf_len = 0;
    }
    memset(ctx->buf + ctx->buf_len, 0, 56 - ctx->buf_len);
    for (int i = 0; i < 8; i++) {
        ctx->buf[56 + i] = (uint8_t)(ctx->bitlen >> ((7 - i) * 8));
    }
    sha256_transform(ctx, ctx->buf);
    for (int i = 0; i < 8; i++) {
        out_digest[i*4]   = (uint8_t)(ctx->state[i] >> 24);
        out_digest[i*4+1] = (uint8_t)(ctx->state[i] >> 16);
        out_digest[i*4+2] = (uint8_t)(ctx->state[i] >> 8);
        out_digest[i*4+3] = (uint8_t)(ctx->state[i]);
    }
}

void hgs_sha256(const uint8_t* data, size_t len, uint8_t out_digest[32]) {
    HgsSha256Ctx ctx;
    hgs_sha256_init(&ctx);
    hgs_sha256_update(&ctx, data, len);
    hgs_sha256_final(&ctx, out_digest);
}
