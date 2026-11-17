#ifndef HGS_SIMD_H
#define HGS_SIMD_H

#ifdef __cplusplus
#include <cstddef>
#include <cstdint>
#define HGS_NODISCARD [[nodiscard]]
extern "C" {
#else
#include <stddef.h>
#include <stdint.h>
#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 202311L
#define HGS_NODISCARD [[nodiscard]]
#else
#define HGS_NODISCARD
#endif
#endif

/**
 * Callback invoked for each strand-symmetric canonical open syncmer (k, s, t).
 * pos: 0-based offset within scanned interval
 * canon_kmer: unsigned 2k-bit canonical k-mer min(fwd, rc)
 * strand: 0 if fwd <= rc, 1 if rc < fwd
 */
typedef void (*HgsSyncmerCb)(uint64_t pos, uint64_t canon_kmer, uint8_t strand, void* ctx);

/** Returns 1 if host CPU supports AVX2 (runtime CPUID check), 0 otherwise. */
HGS_NODISCARD int hgs_cpu_has_avx2(void);

/**
 * Vectorized (AVX2 when supported, SWAR fallback) ASCII DNA/RNA -> 2-bit byte array (0..3).
 * Returns 1 if ALL bytes in src[0..len) were valid uppercase ACGT (fourth_char='T')
 * or uppercase ACGU (fourth_char='U'), else 0 (still writes lowercase/U mapped to 0..3
 * and 0 for non-ACGTU bytes).
 */
HGS_NODISCARD int hgs_ascii_to_2bit_checked(const uint8_t* src, uint8_t* dst, size_t len, uint8_t fourth_char);

/**
 * Packs 2-bit bytes (0..3 per byte) into LSB-first 64-bit words (32 bases per word),
 * writing ceil(len/32) words plus one trailing 0 sentinel word cleanly if out_words has room.
 */
void hgs_pack_2bit_words(const uint8_t* bases, size_t len, uint64_t* out_words);

/** Unpacks LSB-first 64-bit 2-bit words [offset .. offset+len) into uppercase ASCII. */
void hgs_unpack_2bit_ascii(const uint64_t* words, uint64_t offset, size_t len, uint8_t fourth_char, uint8_t* out_ascii);

/** Unpacks LSB-first 64-bit 2-bit words [offset .. offset+len) into 2-bit bytes (0..3). */
void hgs_unpack_2bit_bytes(const uint64_t* words, uint64_t offset, size_t len, uint8_t* out_bytes);

/** Extracts 32 consecutive bases starting at bit-packed base index pos (LSB-first). */
static inline uint64_t hgs_word32_forward(const uint64_t* words, uint64_t base_pos) {
    uint64_t w_idx = base_pos >> 5;
    uint32_t shift = (uint32_t)((base_pos & 31ULL) << 1);
    if (shift == 0) {
        return words[w_idx];
    }
    return (words[w_idx] >> shift) | (words[w_idx + 1] << (64U - shift));
}

/** Exact 64-bit SWAR Watson-Crick reverse complement involution of 32 packed bases. */
static inline uint64_t hgs_revcomp_word32(uint64_t w) {
    w = ~w;
    w = ((w & 0x3333333333333333ULL) << 2) | ((w >> 2) & 0x3333333333333333ULL);
    w = ((w & 0x0F0F0F0F0F0F0F0FULL) << 4) | ((w >> 4) & 0x0F0F0F0F0F0F0F0FULL);
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_bswap64(w);
#else
    return ((w & 0x00000000000000FFULL) << 56) | ((w & 0x000000000000FF00ULL) << 40) |
           ((w & 0x0000000000FF0000ULL) << 24) | ((w & 0x00000000FF000000ULL) << 8) |
           ((w & 0x000000FF00000000ULL) >> 8) | ((w & 0x0000FF0000000000ULL) >> 24) |
           ((w & 0x00FF000000000000ULL) >> 40) | ((w & 0xFF00000000000000ULL) >> 56);
#endif
}

/** Murmur/SplitMix64 finalizers matching Packed2BitSequence.mix32 / mix64 bit-for-bit. */
static inline uint32_t hgs_mix32(uint64_t z) {
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z = z ^ (z >> 31);
    return (uint32_t)z;
}

static inline uint64_t hgs_mix64(uint64_t z) {
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

/** Counts consecutive exact matching 2-bit bases (strand 0 = fwd, 1 = rev-comp). */
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
    uint64_t max_len);

/** Linear O(n) strand-symmetric canonical open-syncmer scanner over 2-bit packed words. */
void hgs_scan_open_syncmers(
    const uint64_t* words,
    uint64_t start_offset,
    uint64_t seg_len,
    int k,
    int s,
    int t,
    HgsSyncmerCb cb,
    void* ctx);

/** Slice-by-8 IEEE 802.3 CRC-32 (pass prev_crc = 0 on initial call). */
HGS_NODISCARD uint32_t hgs_crc32_ieee(uint32_t prev_crc, const uint8_t* data, size_t len);

/** Self-contained FIPS 180-4 SHA-256 digest (32 bytes written to out_digest). */
void hgs_sha256(const uint8_t* data, size_t len, uint8_t out_digest[32]);

/** Streaming SHA-256 context. */
typedef struct {
    uint32_t state[8];
    uint64_t bitlen;
    uint8_t buf[64];
    size_t buf_len;
} HgsSha256Ctx;

void hgs_sha256_init(HgsSha256Ctx* ctx);
void hgs_sha256_update(HgsSha256Ctx* ctx, const uint8_t* data, size_t len);
void hgs_sha256_final(HgsSha256Ctx* ctx, uint8_t out_digest[32]);

#ifdef __cplusplus
}
#endif

#endif /* HGS_SIMD_H */
