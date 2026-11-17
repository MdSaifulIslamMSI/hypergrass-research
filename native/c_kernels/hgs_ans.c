#include "hgs_ans.h"
#include <stdlib.h>
#include <string.h>

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 202311L
static_assert(sizeof(uint32_t) == 4);
#endif

#define ANS_PROB_BITS 12U
#define ANS_SCALE (1U << ANS_PROB_BITS) /* 4096 */
#define ANS_L (1U << 23)                /* 8,388,608 (overflow-safe Byte-rANS) */
#define ANS_XMAX_MUL ((ANS_L >> ANS_PROB_BITS) << 8) /* 524,288 */

#define ANS_MODE_RAW_2BIT 0U
#define ANS_MODE_RANS4_DNA 1U
#define ANS_MODE_RAW_QUAL 2U
#define ANS_MODE_RANS4_QUAL 3U

#define DNA_BLK_CAP (1U << 16) /* 65536 bases per multi-lane block */
#define QUAL_BLK_CAP (1U << 16)

void hgs_ans_free(uint8_t* ptr) {
    free(ptr);
}

typedef struct {
    uint16_t cum;
    uint16_t freq;
} AnsSymRec;

typedef struct {
    uint8_t* data;
    size_t len;
    size_t cap;
} AnsVec;

static int vec_reserve(AnsVec* v, size_t need) {
    if (need <= v->cap) return 0;
    size_t ncap = v->cap ? v->cap : 1024U;
    while (ncap < need) {
        if (ncap > (SIZE_MAX >> 1)) return -1;
        ncap <<= 1;
    }
    uint8_t* p = (uint8_t*)realloc(v->data, ncap);
    if (!p) return -1;
    v->data = p;
    v->cap = ncap;
    return 0;
}

static int vec_push(AnsVec* v, uint8_t b) {
    if (v->len + 1U > v->cap && vec_reserve(v, v->len + 1U) != 0) return -1;
    v->data[v->len++] = b;
    return 0;
}

static int vec_append(AnsVec* v, const uint8_t* src, size_t n) {
    if (n == 0) return 0;
    if (v->len + n > v->cap && vec_reserve(v, v->len + n) != 0) return -1;
    memcpy(v->data + v->len, src, n);
    v->len += n;
    return 0;
}

static int vec_push_u32_le(AnsVec* v, uint32_t x) {
    uint8_t b[4] = {
        (uint8_t)(x & 0xFFU),
        (uint8_t)((x >> 8) & 0xFFU),
        (uint8_t)((x >> 16) & 0xFFU),
        (uint8_t)((x >> 24) & 0xFFU)
    };
    return vec_append(v, b, 4);
}

static inline uint32_t read_u32_le(const uint8_t* p) {
    return (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}

static inline void ans_normalize_cdf(const uint32_t* weights, size_t n, uint16_t* out_cum, uint16_t* out_freq) {
    uint64_t tot = 0;
    size_t best = 0;
    uint32_t best_w = 0;
    for (size_t i = 0; i < n; i++) {
        uint32_t w = weights[i] + 1U;
        tot += w;
        if (w > best_w) {
            best_w = w;
            best = i;
        }
    }
    uint32_t budget = ANS_SCALE - (uint32_t)n;
    uint32_t sum = 0;
    for (size_t i = 0; i < n; i++) {
        uint32_t f = 1U + (uint32_t)(((uint64_t)(weights[i] + 1U) * budget) / tot);
        out_freq[i] = (uint16_t)f;
        sum += f;
    }
    out_freq[best] = (uint16_t)(out_freq[best] + (ANS_SCALE - sum));
    uint16_t c = 0;
    for (size_t i = 0; i < n; i++) {
        out_cum[i] = c;
        c = (uint16_t)(c + out_freq[i]);
    }
}

static inline int ans_encode_one(uint32_t* state, AnsVec* rev_bytes, uint16_t cum, uint16_t freq) {
    uint32_t x = *state;
    uint32_t xmax = ANS_XMAX_MUL * (uint32_t)freq;
    while (x >= xmax) {
        if (vec_push(rev_bytes, (uint8_t)(x & 0xFFU)) != 0) return -1;
        x >>= 8;
    }
    *state = ((x / (uint32_t)freq) << ANS_PROB_BITS) + (x % (uint32_t)freq) + (uint32_t)cum;
    return 0;
}

static inline int ans_decode_find_and_step(
    uint32_t* state,
    const uint8_t** p_cur,
    const uint8_t* p_end,
    const uint16_t* cums,
    const uint16_t* freqs,
    size_t n_syms,
    size_t* out_sym) {
    uint32_t x = *state;
    uint32_t slot = x & (ANS_SCALE - 1U);
    size_t sym = 0;
    while (sym + 1U < n_syms && slot >= (uint32_t)cums[sym + 1U]) {
        sym++;
    }
    *out_sym = sym;
    uint32_t f = (uint32_t)freqs[sym];
    uint32_t c = (uint32_t)cums[sym];
    x = f * (x >> ANS_PROB_BITS) + (slot - c);
    while (x < ANS_L) {
        if (*p_cur >= p_end) return -1;
        x = (x << 8) | (uint32_t)(*(*p_cur)++);
    }
    *state = x;
    return 0;
}

/* ---------------- Adaptive 4-Stripe DNA Predictor (Order 0/2/4/6) ---------------- */

typedef struct {
    uint16_t o2[16][4];
    uint16_t o4[256][4];
    uint16_t o6[4096][4];
    uint32_t ctx2;
    uint32_t ctx4;
    uint32_t ctx6;
} DnaLaneModel;

static void dna_model_init(DnaLaneModel* m) {
    memset(m, 0, sizeof(*m));
}

static inline void dna_model_cdf(const DnaLaneModel* m, uint16_t cum[4], uint16_t freq[4]) {
    const uint16_t* p2 = m->o2[m->ctx2];
    const uint16_t* p4 = m->o4[m->ctx4];
    const uint16_t* p6 = m->o6[m->ctx6];
    uint32_t w[4];
    for (int s = 0; s < 4; s++) {
        w[s] = (uint32_t)p2[s] + ((uint32_t)p4[s] << 2) + ((uint32_t)p6[s] << 4);
    }
    ans_normalize_cdf(w, 4, cum, freq);
}

static inline void dna_model_update(DnaLaneModel* m, uint8_t b) {
    uint16_t* p2 = m->o2[m->ctx2];
    uint16_t* p4 = m->o4[m->ctx4];
    uint16_t* p6 = m->o6[m->ctx6];
    if (++p2[b] > 2048U) { for (int s = 0; s < 4; s++) p2[s] >>= 1; }
    if (++p4[b] > 2048U) { for (int s = 0; s < 4; s++) p4[s] >>= 1; }
    if (++p6[b] > 2048U) { for (int s = 0; s < 4; s++) p6[s] >>= 1; }
    m->ctx2 = ((m->ctx2 << 2) | b) & 15U;
    m->ctx4 = ((m->ctx4 << 2) | b) & 255U;
    m->ctx6 = ((m->ctx6 << 2) | b) & 4095U;
}

static int pack_raw_2bit(const uint8_t* bases, size_t len, AnsVec* out) {
    if (vec_push(out, (uint8_t)ANS_MODE_RAW_2BIT) != 0) return -1;
    size_t nbytes = (len + 3U) >> 2;
    if (vec_reserve(out, out->len + nbytes) != 0) return -1;
    size_t i = 0;
    while (i + 4U <= len) {
        uint8_t v = (uint8_t)((bases[i] & 3U) | ((bases[i + 1U] & 3U) << 2) |
                              ((bases[i + 2U] & 3U) << 4) | ((bases[i + 3U] & 3U) << 6));
        out->data[out->len++] = v;
        i += 4U;
    }
    if (i < len) {
        uint8_t v = 0;
        int sh = 0;
        while (i < len) {
            v |= (uint8_t)((bases[i++] & 3U) << sh);
            sh += 2;
        }
        out->data[out->len++] = v;
    }
    return 0;
}

int hgs_ans4_encode_dna(const uint8_t* two_bit_bases, size_t len, uint8_t** out_buf, size_t* out_len) {
    *out_buf = NULL;
    *out_len = 0;
    AnsVec out = {NULL, 0, 0};
    if (len == 0) {
        if (vec_push(&out, (uint8_t)ANS_MODE_RAW_2BIT) != 0) return -1;
        *out_buf = out.data;
        *out_len = out.len;
        return 0;
    }
    if (len < 64U) {
        if (pack_raw_2bit(two_bit_bases, len, &out) != 0) { free(out.data); return -1; }
        *out_buf = out.data;
        *out_len = out.len;
        return 0;
    }

    DnaLaneModel* models = (DnaLaneModel*)calloc(4, sizeof(DnaLaneModel));
    AnsSymRec* recs[4] = {NULL, NULL, NULL, NULL};
    AnsVec rev[4] = {{NULL,0,0},{NULL,0,0},{NULL,0,0},{NULL,0,0}};
    if (!models) return -1;

    size_t base_q = len >> 2;
    size_t rem = len & 3U;
    size_t lane_len[4];
    size_t lane_off[4];
    size_t cur = 0;
    for (size_t l = 0; l < 4U; l++) {
        lane_len[l] = base_q + (l < rem ? 1U : 0U);
        lane_off[l] = cur;
        cur += lane_len[l];
        if (lane_len[l] > 0) {
            recs[l] = (AnsSymRec*)malloc(lane_len[l] * sizeof(AnsSymRec));
            if (!recs[l]) {
                for (size_t k = 0; k < 4U; k++) free(recs[k]);
                free(models);
                return -1;
            }
        }
        dna_model_init(&models[l]);
    }

    /* Pass 1: forward causal context & CDF snapshot per stripe */
    for (size_t l = 0; l < 4U; l++) {
        const uint8_t* src = two_bit_bases + lane_off[l];
        size_t n = lane_len[l];
        uint16_t cum[4], freq[4];
        for (size_t i = 0; i < n; i++) {
            uint8_t b = src[i] & 3U;
            dna_model_cdf(&models[l], cum, freq);
            recs[l][i].cum = cum[b];
            recs[l][i].freq = freq[b];
            dna_model_update(&models[l], b);
        }
    }

    /* Pass 2: backward interleaved 4-lane rANS encode */
    uint32_t state[4] = {ANS_L, ANS_L, ANS_L, ANS_L};
    int err = 0;
    for (size_t l = 0; l < 4U && !err; l++) {
        size_t n = lane_len[l];
        for (size_t i = n; i > 0; i--) {
            if (ans_encode_one(&state[l], &rev[l], recs[l][i - 1U].cum, recs[l][i - 1U].freq) != 0) {
                err = -1;
                break;
            }
        }
    }

    for (size_t l = 0; l < 4U; l++) free(recs[l]);
    free(models);

    size_t raw_size = 1U + ((len + 3U) >> 2);
    size_t rans_size = 1U + 32U;
    for (size_t l = 0; l < 4U; l++) rans_size += rev[l].len;

    if (err != 0 || rans_size >= raw_size) {
        for (size_t l = 0; l < 4U; l++) free(rev[l].data);
        if (err != 0 || pack_raw_2bit(two_bit_bases, len, &out) != 0) {
            free(out.data);
            return -1;
        }
        *out_buf = out.data;
        *out_len = out.len;
        return 0;
    }

    if (vec_push(&out, (uint8_t)ANS_MODE_RANS4_DNA) != 0) goto fail;
    for (size_t l = 0; l < 4U; l++) {
        if (vec_push_u32_le(&out, state[l]) != 0) goto fail;
    }
    for (size_t l = 0; l < 4U; l++) {
        if (vec_push_u32_le(&out, (uint32_t)rev[l].len) != 0) goto fail;
    }
    for (size_t l = 0; l < 4U; l++) {
        if (vec_reserve(&out, out.len + rev[l].len) != 0) goto fail;
        for (size_t i = rev[l].len; i > 0; i--) {
            out.data[out.len++] = rev[l].data[i - 1U];
        }
        free(rev[l].data);
        rev[l].data = NULL;
    }
    *out_buf = out.data;
    *out_len = out.len;
    return 0;

fail:
    for (size_t l = 0; l < 4U; l++) free(rev[l].data);
    free(out.data);
    return -1;
}

int hgs_ans4_decode_dna(const uint8_t* payload, size_t payload_len, size_t expected_len, uint8_t** out_buf, size_t* out_len) {
    *out_buf = NULL;
    *out_len = 0;
    if (payload_len < 1U) return -1;
    uint8_t mode = payload[0];
    if (expected_len == 0) {
        return (mode == ANS_MODE_RAW_2BIT) ? 0 : -1;
    }
    uint8_t* out = (uint8_t*)malloc(expected_len);
    if (!out) return -1;

    if (mode == ANS_MODE_RAW_2BIT) {
        size_t need = 1U + ((expected_len + 3U) >> 2);
        if (payload_len < need) { free(out); return -1; }
        const uint8_t* src = payload + 1U;
        for (size_t i = 0; i < expected_len; i++) {
            out[i] = (uint8_t)((src[i >> 2] >> ((i & 3U) << 1)) & 3U);
        }
        *out_buf = out;
        *out_len = expected_len;
        return 0;
    }
    if (mode != ANS_MODE_RANS4_DNA || payload_len < 33U) {
        free(out);
        return -1;
    }

    uint32_t state[4];
    uint32_t blen[4];
    const uint8_t* p_cur[4];
    const uint8_t* p_end[4];
    for (size_t l = 0; l < 4U; l++) {
        state[l] = read_u32_le(payload + 1U + l * 4U);
        blen[l] = read_u32_le(payload + 17U + l * 4U);
    }
    size_t pos = 33U;
    for (size_t l = 0; l < 4U; l++) {
        if (pos + (size_t)blen[l] > payload_len) { free(out); return -1; }
        p_cur[l] = payload + pos;
        p_end[l] = payload + pos + (size_t)blen[l];
        pos += (size_t)blen[l];
    }

    DnaLaneModel* models = (DnaLaneModel*)calloc(4, sizeof(DnaLaneModel));
    if (!models) { free(out); return -1; }
    size_t base_q = expected_len >> 2;
    size_t rem = expected_len & 3U;
    size_t lane_len[4];
    size_t lane_off[4];
    size_t cur = 0;
    for (size_t l = 0; l < 4U; l++) {
        lane_len[l] = base_q + (l < rem ? 1U : 0U);
        lane_off[l] = cur;
        cur += lane_len[l];
        dna_model_init(&models[l]);
    }

    /* Decode 4 independent stripes synchronously across 4-way lanes */
    size_t max_q = lane_len[0];
    uint16_t cum[4][4];
    uint16_t freq[4][4];
    for (size_t k = 0; k < max_q; k++) {
        for (size_t l = 0; l < 4U; l++) {
            if (k < lane_len[l]) dna_model_cdf(&models[l], cum[l], freq[l]);
        }
        for (size_t l = 0; l < 4U; l++) {
            if (k < lane_len[l]) {
                size_t sym = 0;
                if (ans_decode_find_and_step(&state[l], &p_cur[l], p_end[l], cum[l], freq[l], 4U, &sym) != 0) {
                    free(models);
                    free(out);
                    return -1;
                }
                out[lane_off[l] + k] = (uint8_t)sym;
                dna_model_update(&models[l], (uint8_t)sym);
            }
        }
    }

    for (size_t l = 0; l < 4U; l++) {
        if (state[l] != ANS_L || p_cur[l] != p_end[l]) {
            free(models);
            free(out);
            return -1;
        }
    }
    free(models);
    *out_buf = out;
    *out_len = expected_len;
    return 0;
}

/* ---------------- Pangenome-Conditioned Phred Quality 4-Stripe Coder ---------------- */

void hgs_ans_compute_homo_bins(const uint8_t* two_bit_bases, size_t len, uint8_t* out_bins) {
    if (len == 0) return;
    uint32_t run = 1U;
    out_bins[0] = 0U;
    uint8_t prev = two_bit_bases[0] & 3U;
    for (size_t i = 1U; i < len; i++) {
        uint8_t b = two_bit_bases[i] & 3U;
        if (b == prev) {
            run++;
        } else {
            run = 1U;
            prev = b;
        }
        out_bins[i] = (run == 1U) ? 0U : (run <= 3U ? 1U : (run <= 6U ? 2U : 3U));
    }
}

/* 2-level hierarchical context table over (mismatch, homo_bin, prev_qbin, cycle_bin):
 *   same_cnt[128][2]: repeat-match probability
 *   coarse_cnt[128][9]: Phred-33 bins 0..7 ([33..96]) + bin 8 (raw byte escape)
 *   fine_cnt[64][8]:   within-bin offset 0..7 conditioned on (coarse, mismatch, prev_qbin>>1)
 *   nibble_uniform:    16-sym CDF for escaped non-[33..96] bytes (100% lossless) */
typedef struct {
    uint16_t same[128][2];
    uint16_t coarse[128][9];
    uint16_t fine[128][8];
} QualModel;

static inline uint32_t qual_ctx(uint8_t prev_q, uint8_t mismatch, uint8_t homo_bin, size_t pos, uint32_t read_len_hint) {
    uint32_t pq_bin = (prev_q >= 33U && prev_q <= 96U) ? (uint32_t)((prev_q - 33U) >> 4) : 3U; /* 2 bits */
    uint32_t cyc_bin = (read_len_hint > 0U) ? (uint32_t)(((pos % read_len_hint) * 2U) / read_len_hint) : 0U; /* 1 bit */
    return ((uint32_t)(mismatch & 1U) << 5) | ((uint32_t)(homo_bin & 3U) << 3) | (pq_bin << 1) | (cyc_bin & 1U);
}

static inline uint32_t qual_fine_ctx(uint32_t c_bin, uint8_t mismatch, uint8_t prev_q) {
    uint32_t pq = (prev_q >= 33U && prev_q <= 96U) ? (uint32_t)((prev_q - 33U) >> 5) : 0U;
    return ((c_bin & 7U) << 2) | ((uint32_t)(mismatch & 1U) << 1) | (pq & 1U);
}

typedef struct {
    AnsSymRec* data;
    size_t len;
    size_t cap;
} QualRecVec;

static int qrec_push(QualRecVec* v, uint16_t cum, uint16_t freq) {
    if (v->len + 1U > v->cap) {
        size_t ncap = v->cap ? (v->cap << 1) : 4096U;
        AnsSymRec* p = (AnsSymRec*)realloc(v->data, ncap * sizeof(AnsSymRec));
        if (!p) return -1;
        v->data = p;
        v->cap = ncap;
    }
    v->data[v->len].cum = cum;
    v->data[v->len].freq = freq;
    v->len++;
    return 0;
}

int hgs_ans4_encode_qual(
    const uint8_t* quals_raw,
    const uint8_t* mismatch_mask,
    const uint8_t* homo_bins,
    size_t len,
    uint32_t read_len_hint,
    uint8_t** out_buf,
    size_t* out_len) {
    *out_buf = NULL;
    *out_len = 0;
    AnsVec out = {NULL, 0, 0};
    if (len == 0) {
        if (vec_push(&out, (uint8_t)ANS_MODE_RAW_QUAL) != 0) return -1;
        *out_buf = out.data;
        *out_len = out.len;
        return 0;
    }

    QualModel* models = (QualModel*)calloc(4, sizeof(QualModel));
    QualRecVec recs[4] = {{NULL,0,0},{NULL,0,0},{NULL,0,0},{NULL,0,0}};
    AnsVec rev[4] = {{NULL,0,0},{NULL,0,0},{NULL,0,0},{NULL,0,0}};
    if (!models) return -1;

    uint32_t uniform16_w[16] = {0};
    uint16_t u16_cum[16], u16_freq[16];
    ans_normalize_cdf(uniform16_w, 16, u16_cum, u16_freq);

    size_t base_q = len >> 2;
    size_t rem = len & 3U;
    size_t lane_len[4], lane_off[4];
    size_t cur = 0;
    for (size_t l = 0; l < 4U; l++) {
        lane_len[l] = base_q + (l < rem ? 1U : 0U);
        lane_off[l] = cur;
        cur += lane_len[l];
    }

    /* Pass 1: forward hierarchical causal probability snapshot */
    for (size_t l = 0; l < 4U; l++) {
        uint8_t prev_q = 66U; /* typical Phred+33 start */
        QualModel* m = &models[l];
        for (size_t k = 0; k < lane_len[l]; k++) {
            size_t gpos = lane_off[l] + k;
            uint8_t q = quals_raw[gpos];
            uint8_t mm = mismatch_mask ? (mismatch_mask[gpos] & 1U) : 0U;
            uint8_t hb = homo_bins ? (homo_bins[gpos] & 3U) : 0U;
            uint32_t ctx = qual_ctx(prev_q, mm, hb, gpos, read_len_hint);

            uint32_t sw[2] = {m->same[ctx][0], m->same[ctx][1]};
            uint16_t scum[2], sfreq[2];
            ans_normalize_cdf(sw, 2U, scum, sfreq);
            int is_same = (q == prev_q) ? 1 : 0;
            if (qrec_push(&recs[l], scum[is_same], sfreq[is_same]) != 0) goto fail_oom;
            if (++m->same[ctx][is_same] > 4096U) {
                m->same[ctx][0] >>= 1;
                m->same[ctx][1] >>= 1;
            }

            if (!is_same) {
                uint32_t cw[9];
                for (int c = 0; c < 9; c++) cw[c] = m->coarse[ctx][c];
                uint16_t ccum[9], cfreq[9];
                ans_normalize_cdf(cw, 9U, ccum, cfreq);
                if (q >= 33U && q <= 96U) {
                    uint32_t rel = (uint32_t)(q - 33U);
                    uint32_t cbin = rel >> 3;
                    uint32_t fbin = rel & 7U;
                    if (qrec_push(&recs[l], ccum[cbin], cfreq[cbin]) != 0) goto fail_oom;
                    if (++m->coarse[ctx][cbin] > 4096U) {
                        for (int c = 0; c < 9; c++) m->coarse[ctx][c] >>= 1;
                    }
                    uint32_t fctx = qual_fine_ctx(cbin, mm, prev_q);
                    uint32_t fw[8];
                    for (int f = 0; f < 8; f++) fw[f] = m->fine[fctx][f];
                    uint16_t fcum[8], ffreq[8];
                    ans_normalize_cdf(fw, 8U, fcum, ffreq);
                    if (qrec_push(&recs[l], fcum[fbin], ffreq[fbin]) != 0) goto fail_oom;
                    if (++m->fine[fctx][fbin] > 4096U) {
                        for (int f = 0; f < 8; f++) m->fine[fctx][f] >>= 1;
                    }
                } else {
                    if (qrec_push(&recs[l], ccum[8], cfreq[8]) != 0) goto fail_oom;
                    if (++m->coarse[ctx][8] > 4096U) {
                        for (int c = 0; c < 9; c++) m->coarse[ctx][c] >>= 1;
                    }
                    uint32_t hi = (uint32_t)((q >> 4) & 0xFU);
                    uint32_t lo = (uint32_t)(q & 0xFU);
                    if (qrec_push(&recs[l], u16_cum[hi], u16_freq[hi]) != 0) goto fail_oom;
                    if (qrec_push(&recs[l], u16_cum[lo], u16_freq[lo]) != 0) goto fail_oom;
                }
            }
            prev_q = q;
        }
    }

    uint32_t state[4] = {ANS_L, ANS_L, ANS_L, ANS_L};
    int err = 0;
    for (size_t l = 0; l < 4U && !err; l++) {
        for (size_t i = recs[l].len; i > 0; i--) {
            if (ans_encode_one(&state[l], &rev[l], recs[l].data[i - 1U].cum, recs[l].data[i - 1U].freq) != 0) {
                err = -1;
                break;
            }
        }
    }
    for (size_t l = 0; l < 4U; l++) free(recs[l].data);
    free(models);

    size_t rans_size = 33U;
    for (size_t l = 0; l < 4U; l++) rans_size += rev[l].len;

    if (err != 0 || rans_size >= len + 1U) {
        for (size_t l = 0; l < 4U; l++) free(rev[l].data);
        if (err != 0 || vec_push(&out, (uint8_t)ANS_MODE_RAW_QUAL) != 0 || vec_append(&out, quals_raw, len) != 0) {
            free(out.data);
            return -1;
        }
        *out_buf = out.data;
        *out_len = out.len;
        return 0;
    }

    if (vec_push(&out, (uint8_t)ANS_MODE_RANS4_QUAL) != 0) goto fail;
    for (size_t l = 0; l < 4U; l++) { if (vec_push_u32_le(&out, state[l]) != 0) goto fail; }
    for (size_t l = 0; l < 4U; l++) { if (vec_push_u32_le(&out, (uint32_t)rev[l].len) != 0) goto fail; }
    for (size_t l = 0; l < 4U; l++) {
        if (vec_reserve(&out, out.len + rev[l].len) != 0) goto fail;
        for (size_t i = rev[l].len; i > 0; i--) out.data[out.len++] = rev[l].data[i - 1U];
        free(rev[l].data);
        rev[l].data = NULL;
    }
    *out_buf = out.data;
    *out_len = out.len;
    return 0;

fail_oom:
    for (size_t l = 0; l < 4U; l++) { free(recs[l].data); free(rev[l].data); }
    free(models);
    free(out.data);
    return -1;

fail:
    for (size_t l = 0; l < 4U; l++) free(rev[l].data);
    free(out.data);
    return -1;
}

int hgs_ans4_decode_qual(
    const uint8_t* payload,
    size_t payload_len,
    const uint8_t* mismatch_mask,
    const uint8_t* homo_bins,
    size_t expected_len,
    uint32_t read_len_hint,
    uint8_t** out_buf,
    size_t* out_len) {
    *out_buf = NULL;
    *out_len = 0;
    if (payload_len < 1U) return -1;
    uint8_t mode = payload[0];
    if (expected_len == 0) {
        return (mode == ANS_MODE_RAW_QUAL) ? 0 : -1;
    }
    uint8_t* out = (uint8_t*)malloc(expected_len);
    if (!out) return -1;

    if (mode == ANS_MODE_RAW_QUAL) {
        if (payload_len < 1U + expected_len) { free(out); return -1; }
        memcpy(out, payload + 1U, expected_len);
        *out_buf = out;
        *out_len = expected_len;
        return 0;
    }
    if (mode != ANS_MODE_RANS4_QUAL || payload_len < 33U) { free(out); return -1; }

    uint32_t state[4], blen[4];
    const uint8_t* p_cur[4];
    const uint8_t* p_end[4];
    for (size_t l = 0; l < 4U; l++) {
        state[l] = read_u32_le(payload + 1U + l * 4U);
        blen[l] = read_u32_le(payload + 17U + l * 4U);
    }
    size_t pos = 33U;
    for (size_t l = 0; l < 4U; l++) {
        if (pos + (size_t)blen[l] > payload_len) { free(out); return -1; }
        p_cur[l] = payload + pos;
        p_end[l] = payload + pos + (size_t)blen[l];
        pos += (size_t)blen[l];
    }

    QualModel* models = (QualModel*)calloc(4, sizeof(QualModel));
    if (!models) { free(out); return -1; }
    uint32_t uniform16_w[16] = {0};
    uint16_t u16_cum[16], u16_freq[16];
    ans_normalize_cdf(uniform16_w, 16U, u16_cum, u16_freq);

    size_t base_q = expected_len >> 2;
    size_t rem = expected_len & 3U;
    size_t lane_len[4], lane_off[4];
    size_t cur = 0;
    for (size_t l = 0; l < 4U; l++) {
        lane_len[l] = base_q + (l < rem ? 1U : 0U);
        lane_off[l] = cur;
        cur += lane_len[l];
    }

    for (size_t l = 0; l < 4U; l++) {
        uint8_t prev_q = 66U;
        QualModel* m = &models[l];
        for (size_t k = 0; k < lane_len[l]; k++) {
            size_t gpos = lane_off[l] + k;
            uint8_t mm = mismatch_mask ? (mismatch_mask[gpos] & 1U) : 0U;
            uint8_t hb = homo_bins ? (homo_bins[gpos] & 3U) : 0U;
            uint32_t ctx = qual_ctx(prev_q, mm, hb, gpos, read_len_hint);

            uint32_t sw[2] = {m->same[ctx][0], m->same[ctx][1]};
            uint16_t scum[2], sfreq[2];
            ans_normalize_cdf(sw, 2U, scum, sfreq);
            size_t is_same = 0;
            if (ans_decode_find_and_step(&state[l], &p_cur[l], p_end[l], scum, sfreq, 2U, &is_same) != 0) goto fail;
            if (++m->same[ctx][is_same] > 4096U) {
                m->same[ctx][0] >>= 1;
                m->same[ctx][1] >>= 1;
            }

            uint8_t q = prev_q;
            if (!is_same) {
                uint32_t cw[9];
                for (int c = 0; c < 9; c++) cw[c] = m->coarse[ctx][c];
                uint16_t ccum[9], cfreq[9];
                ans_normalize_cdf(cw, 9U, ccum, cfreq);
                size_t cbin = 0;
                if (ans_decode_find_and_step(&state[l], &p_cur[l], p_end[l], ccum, cfreq, 9U, &cbin) != 0) goto fail;
                if (++m->coarse[ctx][cbin] > 4096U) {
                    for (int c = 0; c < 9; c++) m->coarse[ctx][c] >>= 1;
                }
                if (cbin < 8U) {
                    uint32_t fctx = qual_fine_ctx((uint32_t)cbin, mm, prev_q);
                    uint32_t fw[8];
                    for (int f = 0; f < 8; f++) fw[f] = m->fine[fctx][f];
                    uint16_t fcum[8], ffreq[8];
                    ans_normalize_cdf(fw, 8U, fcum, ffreq);
                    size_t fbin = 0;
                    if (ans_decode_find_and_step(&state[l], &p_cur[l], p_end[l], fcum, ffreq, 8U, &fbin) != 0) goto fail;
                    if (++m->fine[fctx][fbin] > 4096U) {
                        for (int f = 0; f < 8; f++) m->fine[fctx][f] >>= 1;
                    }
                    q = (uint8_t)(33U + (cbin << 3) + fbin);
                } else {
                    size_t hi = 0, lo = 0;
                    if (ans_decode_find_and_step(&state[l], &p_cur[l], p_end[l], u16_cum, u16_freq, 16U, &hi) != 0) goto fail;
                    if (ans_decode_find_and_step(&state[l], &p_cur[l], p_end[l], u16_cum, u16_freq, 16U, &lo) != 0) goto fail;
                    q = (uint8_t)((hi << 4) | lo);
                }
            }
            out[gpos] = q;
            prev_q = q;
        }
    }

    for (size_t l = 0; l < 4U; l++) {
        if (state[l] != ANS_L || p_cur[l] != p_end[l]) goto fail;
    }
    free(models);
    *out_buf = out;
    *out_len = expected_len;
    return 0;

fail:
    free(models);
    free(out);
    return -1;
}
