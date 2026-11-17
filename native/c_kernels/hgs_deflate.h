#ifndef HGS_DEFLATE_H
#define HGS_DEFLATE_H

#include <stddef.h>
#include <stdint.h>

#if defined(__cplusplus) || (defined(__STDC_VERSION__) && __STDC_VERSION__ >= 202311L)
#define HGS_DEF_NODISCARD [[nodiscard]]
#else
#define HGS_DEF_NODISCARD
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** Frees a heap buffer allocated by hgs_deflate_* / hgs_inflate_* functions. */
void hgs_deflate_free(uint8_t* ptr);

/**
 * Raw RFC 1951 DEFLATE (Deflater nowrap=true compatible).
 * Caller must free *out_buf with hgs_deflate_free(). Returns 0 on success, -1 on error.
 */
HGS_DEF_NODISCARD int hgs_deflate_raw(const uint8_t* src, size_t src_len, uint8_t** out_buf, size_t* out_len);

/**
 * Raw RFC 1951 INFLATE supporting BTYPE 0 (stored), 1 (fixed), and 2 (dynamic Huffman).
 * Caller must free *out_buf with hgs_deflate_free(). Returns 0 on success, -1 on error.
 */
HGS_DEF_NODISCARD int hgs_inflate_raw(const uint8_t* src, size_t src_len, size_t expected_len_hint, uint8_t** out_buf, size_t* out_len);

/** RFC 1950 zlib-wrapped DEFLATE (for HGS2 panelNames). */
HGS_DEF_NODISCARD int hgs_zlib_compress(const uint8_t* src, size_t src_len, uint8_t** out_buf, size_t* out_len);

/** RFC 1950 zlib-wrapped INFLATE (verifies Adler-32 trailer). */
HGS_DEF_NODISCARD int hgs_zlib_decompress(const uint8_t* src, size_t src_len, uint8_t** out_buf, size_t* out_len);

#ifdef __cplusplus
}
#endif

#endif /* HGS_DEFLATE_H */
