#ifndef HGS_ANS_H
#define HGS_ANS_H

#ifdef __cplusplus
#include <cstddef>
#include <cstdint>
#define HGS_ANS_NODISCARD [[nodiscard]]
extern "C" {
#else
#include <stddef.h>
#include <stdint.h>
#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 202311L
#define HGS_ANS_NODISCARD [[nodiscard]]
#else
#define HGS_ANS_NODISCARD
#endif
#endif

/** Frees buffers allocated by hgs_ans4_* functions. */
void hgs_ans_free(uint8_t* ptr);

/**
 * 4-stripe SIMD-friendly interleaved rANS encoder for 2-bit DNA/novelty bases (0..3),
 * with per-block Order-0/2/4/6 adaptive probability mixing and strict Law IV raw 2-bit ceiling.
 * Returns 0 on success, -1 on error. Caller frees *out_buf with hgs_ans_free().
 */
HGS_ANS_NODISCARD int hgs_ans4_encode_dna(
    const uint8_t* two_bit_bases,
    size_t len,
    uint8_t** out_buf,
    size_t* out_len);

/**
 * 4-stripe interleaved rANS decoder matching hgs_ans4_encode_dna().
 * Writes exactly expected_len 2-bit bytes (0..3) into caller-allocated out_two_bit[0..expected_len)
 * or allocates *out_two_bit when out_len != NULL.
 */
HGS_ANS_NODISCARD int hgs_ans4_decode_dna(
    const uint8_t* payload,
    size_t payload_len,
    size_t expected_len,
    uint8_t** out_buf,
    size_t* out_len);

/**
 * Computes deterministic homopolymer run-length bin (0..3) at each base position
 * over 2-bit bases (0..3). Shared identically by encoder & decoder (0-byte overhead).
 */
void hgs_ans_compute_homo_bins(const uint8_t* two_bit_bases, size_t len, uint8_t* out_bins);

/**
 * Pangenome-conditioned 4-stripe interleaved rANS Phred quality encoder.
 * Conditions base quality scores on:
 *   - mismatch_mask[i] (0 = matched donor haplotype, 1 = SNP/insertion/novelty)
 *   - homo_bins[i]     (0..3 homopolymer context from decoded DNA)
 *   - read_len_hint    (for sequencing cycle quartile context; 0 if variable)
 * Guarantees lossless recovery for ALL byte values 0..255 with raw bypass ceiling.
 */
HGS_ANS_NODISCARD int hgs_ans4_encode_qual(
    const uint8_t* quals_raw,
    const uint8_t* mismatch_mask,
    const uint8_t* homo_bins,
    size_t len,
    uint32_t read_len_hint,
    uint8_t** out_buf,
    size_t* out_len);

/** Decodes Pangenome-conditioned 4-stripe rANS quality stream back to exact original ASCII bytes. */
HGS_ANS_NODISCARD int hgs_ans4_decode_qual(
    const uint8_t* payload,
    size_t payload_len,
    const uint8_t* mismatch_mask,
    const uint8_t* homo_bins,
    size_t expected_len,
    uint32_t read_len_hint,
    uint8_t** out_buf,
    size_t* out_len);

#ifdef __cplusplus
}
#endif

#endif /* HGS_ANS_H */
