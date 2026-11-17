#include "hgs_deflate.h"
#include <stdlib.h>
#include <string.h>

void hgs_deflate_free(uint8_t* ptr) {
    free(ptr);
}

typedef struct {
    const uint8_t* in;
    size_t in_len;
    size_t pos;
    uint32_t bitbuf;
    int bitcnt;
    uint8_t* out;
    size_t out_len;
    size_t out_cap;
} InfState;

static int inf_ensure(InfState* s, size_t extra) {
    if (s->out_len + extra <= s->out_cap) return 0;
    size_t ncap = s->out_cap ? (s->out_cap * 2) : 4096;
    while (ncap < s->out_len + extra) {
        if (ncap > (SIZE_MAX / 2)) return -1;
        ncap *= 2;
    }
    uint8_t* p = (uint8_t*)realloc(s->out, ncap);
    if (!p) return -1;
    s->out = p;
    s->out_cap = ncap;
    return 0;
}

static int inf_bits(InfState* s, int need, uint32_t* val) {
    while (s->bitcnt < need) {
        if (s->pos >= s->in_len) return -1;
        s->bitbuf |= ((uint32_t)s->in[s->pos++]) << s->bitcnt;
        s->bitcnt += 8;
    }
    *val = s->bitbuf & ((1U << need) - 1U);
    s->bitbuf >>= need;
    s->bitcnt -= need;
    return 0;
}

#define MAXBITS 15
#define MAXLCODES 288
#define MAXDCODES 30
#define MAXCODES (MAXLCODES + MAXDCODES)

typedef struct {
    int16_t count[MAXBITS + 1];
    int16_t symbol[MAXLCODES];
} HuffTree;

static int huff_build(HuffTree* h, const int16_t* lengths, int n) {
    int offs[MAXBITS + 1];
    for (int i = 0; i <= MAXBITS; i++) h->count[i] = 0;
    for (int i = 0; i < n; i++) {
        if (lengths[i] < 0 || lengths[i] > MAXBITS) return -1;
        h->count[lengths[i]]++;
    }
    if (h->count[0] == n) return 0;

    int left = 1;
    for (int len = 1; len <= MAXBITS; len++) {
        left <<= 1;
        left -= h->count[len];
        if (left < 0) return -1;
    }
    offs[1] = 0;
    for (int len = 1; len < MAXBITS; len++) {
        offs[len + 1] = offs[len] + h->count[len];
    }
    for (int sym = 0; sym < n; sym++) {
        if (lengths[sym] != 0) {
            h->symbol[offs[lengths[sym]]++] = (int16_t)sym;
        }
    }
    return 0;
}

static int huff_decode(InfState* s, const HuffTree* h, int* out_sym) {
    int code = 0, first = 0, index = 0;
    for (int len = 1; len <= MAXBITS; len++) {
        uint32_t bit = 0;
        if (inf_bits(s, 1, &bit) != 0) return -1;
        code |= (int)bit;
        int count = h->count[len];
        if (code - count < first) {
            *out_sym = h->symbol[index + (code - first)];
            return 0;
        }
        index += count;
        first += count;
        first <<= 1;
        code <<= 1;
    }
    return -1;
}

static int inf_codes(InfState* s, const HuffTree* lencode, const HuffTree* distcode) {
    static const int16_t lens[29] = {
        3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31,
        35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258
    };
    static const int16_t lext[29] = {
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,
        3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0
    };
    static const int16_t dists[30] = {
        1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193,
        257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577
    };
    static const int16_t dext[30] = {
        0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6,
        7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13
    };

    for (;;) {
        int sym = 0;
        if (huff_decode(s, lencode, &sym) != 0) return -1;
        if (sym < 256) {
            if (inf_ensure(s, 1) != 0) return -1;
            s->out[s->out_len++] = (uint8_t)sym;
        } else if (sym == 256) {
            return 0;
        } else {
            sym -= 257;
            if (sym >= 29) return -1;
            uint32_t ex = 0;
            if (lext[sym] && inf_bits(s, lext[sym], &ex) != 0) return -1;
            size_t len = (size_t)lens[sym] + ex;
            int dsym = 0;
            if (huff_decode(s, distcode, &dsym) != 0 || dsym >= 30) return -1;
            uint32_t dex = 0;
            if (dext[dsym] && inf_bits(s, dext[dsym], &dex) != 0) return -1;
            size_t dist = (size_t)dists[dsym] + dex;
            if (dist == 0 || dist > s->out_len) return -1;
            if (inf_ensure(s, len) != 0) return -1;
            for (size_t i = 0; i < len; i++) {
                s->out[s->out_len] = s->out[s->out_len - dist];
                s->out_len++;
            }
        }
    }
}

int hgs_inflate_raw(const uint8_t* src, size_t src_len, size_t expected_len_hint, uint8_t** out_buf, size_t* out_len) {
    *out_buf = NULL;
    *out_len = 0;
    InfState s;
    memset(&s, 0, sizeof(s));
    s.in = src;
    s.in_len = src_len;
    s.out_cap = expected_len_hint ? expected_len_hint : (src_len * 2 + 256);
    s.out = (uint8_t*)malloc(s.out_cap);
    if (!s.out) return -1;

    uint32_t bfinal = 0;
    do {
        uint32_t btype = 0;
        if (inf_bits(&s, 1, &bfinal) != 0) goto fail;
        if (inf_bits(&s, 2, &btype) != 0) goto fail;
        if (btype == 0) {
            s.bitbuf = 0;
            s.bitcnt = 0;
            if (s.pos + 4 > s.in_len) goto fail;
            uint16_t len = (uint16_t)s.in[s.pos] | ((uint16_t)s.in[s.pos + 1] << 8);
            uint16_t nlen = (uint16_t)s.in[s.pos + 2] | ((uint16_t)s.in[s.pos + 3] << 8);
            s.pos += 4;
            if ((uint16_t)~len != nlen || s.pos + len > s.in_len) goto fail;
            if (inf_ensure(&s, len) != 0) goto fail;
            memcpy(s.out + s.out_len, s.in + s.pos, len);
            s.out_len += len;
            s.pos += len;
        } else if (btype == 1) {
            int16_t lengths[MAXCODES];
            int idx = 0;
            for (; idx < 144; idx++) lengths[idx] = 8;
            for (; idx < 256; idx++) lengths[idx] = 9;
            for (; idx < 280; idx++) lengths[idx] = 7;
            for (; idx < 288; idx++) lengths[idx] = 8;
            HuffTree lencode, distcode;
            if (huff_build(&lencode, lengths, 288) != 0) goto fail;
            for (idx = 0; idx < 30; idx++) lengths[idx] = 5;
            if (huff_build(&distcode, lengths, 30) != 0) goto fail;
            if (inf_codes(&s, &lencode, &distcode) != 0) goto fail;
        } else if (btype == 2) {
            static const int16_t order[19] = {
                16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
            };
            uint32_t hlit_r = 0, hdist_r = 0, hclen_r = 0;
            if (inf_bits(&s, 5, &hlit_r) != 0 || inf_bits(&s, 5, &hdist_r) != 0 || inf_bits(&s, 4, &hclen_r) != 0) goto fail;
            int nlen = (int)hlit_r + 257;
            int ndist = (int)hdist_r + 1;
            int ncode = (int)hclen_r + 4;
            if (nlen > 286 || ndist > MAXDCODES) goto fail;

            int16_t cl_lens[19];
            for (int i = 0; i < 19; i++) cl_lens[i] = 0;
            for (int i = 0; i < ncode; i++) {
                uint32_t v = 0;
                if (inf_bits(&s, 3, &v) != 0) goto fail;
                cl_lens[order[i]] = (int16_t)v;
            }
            HuffTree cl_tree;
            if (huff_build(&cl_tree, cl_lens, 19) != 0) goto fail;

            int16_t lengths[MAXCODES];
            int idx = 0;
            while (idx < nlen + ndist) {
                int sym = 0;
                if (huff_decode(&s, &cl_tree, &sym) != 0) goto fail;
                if (sym < 16) {
                    lengths[idx++] = (int16_t)sym;
                } else {
                    int16_t rep_val = 0;
                    int rep_cnt = 0;
                    uint32_t ex = 0;
                    if (sym == 16) {
                        if (idx == 0 || inf_bits(&s, 2, &ex) != 0) goto fail;
                        rep_val = lengths[idx - 1];
                        rep_cnt = 3 + (int)ex;
                    } else if (sym == 17) {
                        if (inf_bits(&s, 3, &ex) != 0) goto fail;
                        rep_cnt = 3 + (int)ex;
                    } else {
                        if (inf_bits(&s, 7, &ex) != 0) goto fail;
                        rep_cnt = 11 + (int)ex;
                    }
                    if (idx + rep_cnt > nlen + ndist) goto fail;
                    while (rep_cnt-- > 0) lengths[idx++] = rep_val;
                }
            }
            if (lengths[256] == 0) goto fail;
            HuffTree lencode, distcode;
            if (huff_build(&lencode, lengths, nlen) != 0) goto fail;
            if (huff_build(&distcode, lengths + nlen, ndist) != 0) goto fail;
            if (inf_codes(&s, &lencode, &distcode) != 0) goto fail;
        } else {
            goto fail;
        }
    } while (!bfinal);

    *out_buf = s.out;
    *out_len = s.out_len;
    return 0;

fail:
    free(s.out);
    return -1;
}

/* Raw DEFLATE encoder (LZ77 + Fixed Huffman BTYPE=1, automatic stored fallback) */
typedef struct {
    uint8_t* buf;
    size_t len;
    size_t cap;
    uint32_t bitbuf;
    int bitcnt;
} DefOut;

static int def_grow(DefOut* o, size_t need) {
    if (o->len + need <= o->cap) return 0;
    size_t ncap = o->cap ? o->cap * 2 : 1024;
    while (ncap < o->len + need) ncap *= 2;
    uint8_t* p = (uint8_t*)realloc(o->buf, ncap);
    if (!p) return -1;
    o->buf = p;
    o->cap = ncap;
    return 0;
}

static int def_put_lsb(DefOut* o, uint32_t val, int nbits) {
    o->bitbuf |= (val << o->bitcnt);
    o->bitcnt += nbits;
    while (o->bitcnt >= 8) {
        if (def_grow(o, 1) != 0) return -1;
        o->buf[o->len++] = (uint8_t)(o->bitbuf & 0xFFU);
        o->bitbuf >>= 8;
        o->bitcnt -= 8;
    }
    return 0;
}

static int def_put_msb(DefOut* o, uint32_t code, int nbits) {
    uint32_t rev = 0;
    for (int i = 0; i < nbits; i++) {
        rev = (rev << 1) | ((code >> i) & 1U);
    }
    return def_put_lsb(o, rev, nbits);
}

static int emit_fixed_lit(DefOut* o, int sym) {
    if (sym <= 143) return def_put_msb(o, (uint32_t)(0x30 + sym), 8);
    if (sym <= 255) return def_put_msb(o, (uint32_t)(0x190 + (sym - 144)), 9);
    if (sym <= 279) return def_put_msb(o, (uint32_t)(sym - 256), 7);
    return def_put_msb(o, (uint32_t)(0xC0 + (sym - 280)), 8);
}

static int emit_fixed_match(DefOut* o, size_t len, size_t dist) {
    static const int lbase[29] = {
        3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258
    };
    static const int lext[29] = {
        0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,0
    };
    static const int dbase[30] = {
        1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,
        1025,1537,2049,3073,4097,6145,8193,12289,16385,24577
    };
    static const int dext[30] = {
        0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13
    };
    int li = 0;
    while (li < 28 && (int)len >= lbase[li + 1]) li++;
    if (emit_fixed_lit(o, 257 + li) != 0) return -1;
    if (lext[li] && def_put_lsb(o, (uint32_t)((int)len - lbase[li]), lext[li]) != 0) return -1;
    int di = 0;
    while (di < 29 && (int)dist >= dbase[di + 1]) di++;
    if (def_put_msb(o, (uint32_t)di, 5) != 0) return -1;
    if (dext[di] && def_put_lsb(o, (uint32_t)((int)dist - dbase[di]), dext[di]) != 0) return -1;
    return 0;
}

static int deflate_stored_fallback(const uint8_t* src, size_t src_len, uint8_t** out_buf, size_t* out_len) {
    size_t blocks = (src_len == 0) ? 1 : ((src_len + 65534) / 65535);
    uint8_t* out = (uint8_t*)malloc(src_len + blocks * 5 + 8);
    if (!out) return -1;
    size_t op = 0, sp = 0;
    do {
        size_t chunk = src_len - sp;
        if (chunk > 65535) chunk = 65535;
        int bfinal = (sp + chunk == src_len) ? 1 : 0;
        out[op++] = (uint8_t)bfinal;
        uint16_t l16 = (uint16_t)chunk, nl16 = (uint16_t)~l16;
        out[op++] = (uint8_t)(l16 & 0xFF); out[op++] = (uint8_t)(l16 >> 8);
        out[op++] = (uint8_t)(nl16 & 0xFF); out[op++] = (uint8_t)(nl16 >> 8);
        if (chunk > 0) { memcpy(out + op, src + sp, chunk); op += chunk; sp += chunk; }
    } while (sp < src_len);
    *out_buf = out;
    *out_len = op;
    return 0;
}

int hgs_deflate_raw(const uint8_t* src, size_t src_len, uint8_t** out_buf, size_t* out_len) {
    if (src_len <= 12) return deflate_stored_fallback(src, src_len, out_buf, out_len);

    DefOut o;
    memset(&o, 0, sizeof(o));
    o.cap = src_len + 256;
    o.buf = (uint8_t*)malloc(o.cap);
    if (!o.buf) return -1;

    /* Single final fixed-Huffman block: BFINAL=1, BTYPE=01 (LSB-first bits: 1, then 1,0) */
    def_put_lsb(&o, 3U, 3);

    enum { HSIZE = 32768 };
    int32_t* head = (int32_t*)malloc(HSIZE * sizeof(int32_t));
    int32_t* prev = (int32_t*)malloc(32768 * sizeof(int32_t));
    if (!head || !prev) { free(head); free(prev); free(o.buf); return -1; }
    for (int i = 0; i < HSIZE; i++) head[i] = -1;
    for (int i = 0; i < 32768; i++) prev[i] = -1;

    size_t pos = 0;
    while (pos < src_len) {
        size_t best_len = 0, best_dist = 0;
        if (pos + 3 <= src_len) {
            uint32_t h = (((uint32_t)src[pos] << 10) ^ ((uint32_t)src[pos+1] << 5) ^ src[pos+2]) & (HSIZE - 1);
            int32_t cur = head[h];
            prev[pos & 32767] = cur;
            head[h] = (int32_t)pos;
            int chain = 0;
            while (cur >= 0 && chain < 32) {
                size_t dist = pos - (size_t)cur;
                if (dist == 0 || dist > 32768) break;
                size_t maxl = src_len - pos;
                if (maxl > 258) maxl = 258;
                if (maxl >= 3 && src[cur + best_len] == src[pos + best_len]) {
                    size_t ml = 0;
                    while (ml < maxl && src[(size_t)cur + ml] == src[pos + ml]) ml++;
                    if (ml >= 3 && ml > best_len) {
                        best_len = ml;
                        best_dist = dist;
                        if (ml == 258) break;
                    }
                }
                cur = prev[cur & 32767];
                chain++;
            }
        }
        if (best_len >= 3) {
            if (emit_fixed_match(&o, best_len, best_dist) != 0) goto err;
            for (size_t k = 1; k < best_len && (pos + k + 3 <= src_len); k++) {
                size_t np = pos + k;
                uint32_t nh = (((uint32_t)src[np] << 10) ^ ((uint32_t)src[np+1] << 5) ^ src[np+2]) & (HSIZE - 1);
                prev[np & 32767] = head[nh];
                head[nh] = (int32_t)np;
            }
            pos += best_len;
        } else {
            if (emit_fixed_lit(&o, src[pos++]) != 0) goto err;
        }
    }
    free(head);
    free(prev);

    if (emit_fixed_lit(&o, 256) != 0) { free(o.buf); return -1; }
    if (o.bitcnt > 0) {
        if (def_grow(&o, 1) != 0) { free(o.buf); return -1; }
        o.buf[o.len++] = (uint8_t)(o.bitbuf & 0xFFU);
    }
    if (o.len >= src_len + ((src_len + 65534) / 65535) * 5) {
        free(o.buf);
        return deflate_stored_fallback(src, src_len, out_buf, out_len);
    }
    *out_buf = o.buf;
    *out_len = o.len;
    return 0;

err:
    free(head);
    free(prev);
    free(o.buf);
    return -1;
}

static uint32_t adler32_calc(const uint8_t* data, size_t len) {
    uint32_t s1 = 1, s2 = 0;
    for (size_t i = 0; i < len; i++) {
        s1 = (s1 + data[i]) % 65521U;
        s2 = (s2 + s1) % 65521U;
    }
    return (s2 << 16) | s1;
}

int hgs_zlib_compress(const uint8_t* src, size_t src_len, uint8_t** out_buf, size_t* out_len) {
    uint8_t* raw = NULL;
    size_t rlen = 0;
    if (hgs_deflate_raw(src, src_len, &raw, &rlen) != 0) return -1;
    uint8_t* out = (uint8_t*)malloc(rlen + 6);
    if (!out) { free(raw); return -1; }
    out[0] = 0x78; out[1] = 0x9C;
    memcpy(out + 2, raw, rlen);
    free(raw);
    uint32_t ad = adler32_calc(src, src_len);
    out[2 + rlen]     = (uint8_t)(ad >> 24);
    out[2 + rlen + 1] = (uint8_t)(ad >> 16);
    out[2 + rlen + 2] = (uint8_t)(ad >> 8);
    out[2 + rlen + 3] = (uint8_t)(ad);
    *out_buf = out;
    *out_len = rlen + 6;
    return 0;
}

int hgs_zlib_decompress(const uint8_t* src, size_t src_len, uint8_t** out_buf, size_t* out_len) {
    if (src_len < 6) return -1;
    uint8_t cmf = src[0], flg = src[1];
    if ((cmf & 0x0F) != 8 || (((uint16_t)cmf << 8 | flg) % 31) != 0) return -1;
    size_t off = (flg & 0x20) ? 6 : 2;
    if (src_len < off + 4) return -1;
    uint8_t* raw = NULL;
    size_t rlen = 0;
    if (hgs_inflate_raw(src + off, src_len - off - 4, 0, &raw, &rlen) != 0) return -1;
    uint32_t want = ((uint32_t)src[src_len-4] << 24) | ((uint32_t)src[src_len-3] << 16) |
                    ((uint32_t)src[src_len-2] << 8) | (uint32_t)src[src_len-1];
    if (adler32_calc(raw, rlen) != want) { free(raw); return -1; }
    *out_buf = raw;
    *out_len = rlen;
    return 0;
}
