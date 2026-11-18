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
#define ANS_MODE_RANS4_QUAL_BWT 4U

#define DNA_BLK_CAP (1U << 16) /* 65536 bases per multi-lane block */
#define QUAL_BLK_CAP (1U << 16)
#define QUAL_BWT_BLK 2000000U  /* 2 MB BWT block cap (~10 MB peak SA-IS workspace) */

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

static inline void ans_normalize_cdf_with_tot(
    const uint32_t* weights,
    size_t n,
    uint32_t tot,
    size_t best,
    uint16_t* out_cum,
    uint16_t* out_freq) {
    uint32_t budget = ANS_SCALE - (uint32_t)n;
    uint32_t sum = 0;
    if (tot <= 8192U) {
        /* Exact over 1 <= tot <= 8192 for (w + 1) * budget < 2^23 (verified zero boundary error). */
        uint64_t inv = ((1ULL << 42) + (uint64_t)tot - 1ULL) / (uint64_t)tot;
        for (size_t i = 0; i < n; i++) {
            uint64_t num = (uint64_t)(weights[i] + 1U) * (uint64_t)budget;
            uint32_t f = 1U + (uint32_t)((num * inv) >> 42);
            out_freq[i] = (uint16_t)f;
            sum += f;
        }
    } else {
        for (size_t i = 0; i < n; i++) {
            uint32_t f = 1U + (uint32_t)(((uint64_t)(weights[i] + 1U) * (uint64_t)budget) / (uint64_t)tot);
            out_freq[i] = (uint16_t)f;
            sum += f;
        }
    }
    out_freq[best] = (uint16_t)(out_freq[best] + (ANS_SCALE - sum));
    uint16_t c = 0;
    for (size_t i = 0; i < n; i++) {
        out_cum[i] = c;
        c = (uint16_t)(c + out_freq[i]);
    }
}

static inline void ans_normalize_cdf(const uint32_t* weights, size_t n, uint16_t* out_cum, uint16_t* out_freq) {
    uint32_t tot = 0;
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
    ans_normalize_cdf_with_tot(weights, n, tot, best, out_cum, out_freq);
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
    if (n_syms <= 16U) {
        while (sym + 1U < n_syms && slot >= (uint32_t)cums[sym + 1U]) {
            sym++;
        }
    } else {
        size_t lo = 0U, hi = n_syms - 1U;
        while (lo < hi) {
            size_t mid = lo + ((hi - lo + 1U) >> 1);
            if ((uint32_t)cums[mid] <= slot) {
                lo = mid;
            } else {
                hi = mid - 1U;
            }
        }
        sym = lo;
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

/* ---------------- Pangenome-Stratified BWT + MTF + RUNA/RUNB + 4x-rANS Quality Coder (Mode 4) ---------------- */

static int sais_is_lms(const uint8_t* t, int32_t i) {
    if (i <= 0) return 0;
    int t_i = (t[i >> 3] >> (i & 7)) & 1;
    int t_p = (t[(i - 1) >> 3] >> ((i - 1) & 7)) & 1;
    return t_i && !t_p;
}

static inline int sais_tget(const uint8_t* t, int32_t i) {
    return (t[i >> 3] >> (i & 7)) & 1;
}

static inline void sais_tset(uint8_t* t, int32_t i, int v) {
    if (v) t[i >> 3] |= (uint8_t)(1U << (i & 7));
    else t[i >> 3] &= (uint8_t)~(1U << (i & 7));
}

static void sais_get_buckets(const int32_t* s, int32_t* bkt, int32_t n, int32_t K, int end) {
    for (int32_t i = 0; i < K; i++) bkt[i] = 0;
    for (int32_t i = 0; i < n; i++) bkt[s[i]]++;
    int32_t sum = 0;
    for (int32_t i = 0; i < K; i++) {
        sum += bkt[i];
        bkt[i] = end ? sum : (sum - bkt[i]);
    }
}

static void sais_induce_sal(const uint8_t* t, int32_t* sa, const int32_t* s, int32_t* bkt, int32_t n, int32_t K, int end) {
    sais_get_buckets(s, bkt, n, K, end);
    for (int32_t i = 0; i < n; i++) {
        int32_t j = sa[i] - 1;
        if (j >= 0 && !sais_tget(t, j)) sa[bkt[s[j]]++] = j;
    }
}

static void sais_induce_sas(const uint8_t* t, int32_t* sa, const int32_t* s, int32_t* bkt, int32_t n, int32_t K, int end) {
    sais_get_buckets(s, bkt, n, K, end);
    for (int32_t i = n - 1; i >= 0; i--) {
        int32_t j = sa[i] - 1;
        if (j >= 0 && sais_tget(t, j)) sa[--bkt[s[j]]] = j;
    }
}

static int sais_compute(const int32_t* s, int32_t* sa, int32_t n, int32_t K) {
    size_t t_bytes = (size_t)((n + 7) >> 3);
    uint8_t* t = (uint8_t*)calloc(t_bytes, 1U);
    int32_t* bkt = (int32_t*)malloc((size_t)K * sizeof(int32_t));
    if (!t || !bkt) { free(t); free(bkt); return -1; }

    sais_tset(t, n - 1, 1);
    if (n >= 2) sais_tset(t, n - 2, 0);
    for (int32_t i = n - 3; i >= 0; i--) {
        int v = (s[i] < s[i + 1]) || (s[i] == s[i + 1] && sais_tget(t, i + 1));
        sais_tset(t, i, v);
    }

    sais_get_buckets(s, bkt, n, K, 1);
    for (int32_t i = 0; i < n; i++) sa[i] = -1;
    for (int32_t i = 1; i < n; i++) {
        if (sais_is_lms(t, i)) sa[--bkt[s[i]]] = i;
    }
    sais_induce_sal(t, sa, s, bkt, n, K, 0);
    sais_induce_sas(t, sa, s, bkt, n, K, 1);

    int32_t n1 = 0;
    for (int32_t i = 0; i < n; i++) {
        if (sais_is_lms(t, sa[i])) sa[n1++] = sa[i];
    }
    for (int32_t i = n1; i < n; i++) sa[i] = -1;
    int32_t name = 0;
    int32_t prev = -1;
    for (int32_t i = 0; i < n1; i++) {
        int32_t pos = sa[i];
        int diff = 0;
        for (int32_t d = 0; d < n; d++) {
            if (prev == -1 || s[pos + d] != s[prev + d] || sais_tget(t, pos + d) != sais_tget(t, prev + d)) {
                diff = 1;
                break;
            } else if (d > 0 && (sais_is_lms(t, pos + d) || sais_is_lms(t, prev + d))) {
                break;
            }
        }
        if (diff) { name++; prev = pos; }
        int32_t slot = pos >> 1;
        sa[n1 + slot] = name - 1;
    }
    for (int32_t i = n - 1, j = n - 1; i >= n1; i--) {
        if (sa[i] >= 0) sa[j--] = sa[i];
    }

    int32_t* sa1 = sa;
    int32_t* s1 = sa + n - n1;
    if (name < n1) {
        if (sais_compute(s1, sa1, n1, name) != 0) { free(t); free(bkt); return -1; }
    } else {
        for (int32_t i = 0; i < n1; i++) sa1[s1[i]] = i;
    }

    sais_get_buckets(s, bkt, n, K, 1);
    for (int32_t i = 1, j = 0; i < n; i++) {
        if (sais_is_lms(t, i)) s1[j++] = i;
    }
    for (int32_t i = 0; i < n1; i++) sa1[i] = s1[sa1[i]];
    for (int32_t i = n1; i < n; i++) sa[i] = -1;
    for (int32_t i = n1 - 1; i >= 0; i--) {
        int32_t j = sa[i];
        sa[i] = -1;
        sa[--bkt[s[j]]] = j;
    }
    sais_induce_sal(t, sa, s, bkt, n, K, 0);
    sais_induce_sas(t, sa, s, bkt, n, K, 1);

    free(t);
    free(bkt);
    return 0;
}

typedef struct {
    uint16_t* data;
    size_t len;
    size_t cap;
} SymU16Vec;

static int sym16_push(SymU16Vec* v, uint16_t x) {
    if (v->len + 1U > v->cap) {
        size_t ncap = v->cap ? (v->cap << 1) : 4096U;
        uint16_t* p = (uint16_t*)realloc(v->data, ncap * sizeof(uint16_t));
        if (!p) return -1;
        v->data = p;
        v->cap = ncap;
    }
    v->data[v->len++] = x;
    return 0;
}

static int flush_zero_run(SymU16Vec* out, size_t run) {
    while (run > 0U) {
        run--;
        if (sym16_push(out, (uint16_t)(run & 1U)) != 0) return -1;
        run >>= 1;
    }
    return 0;
}

static inline void bwt_update_cdf_state(
    uint32_t* w,
    size_t nsym,
    size_t sym,
    uint32_t* row_sum,
    size_t* best_idx,
    uint32_t* best_w) {
    w[sym] += 8U;
    *row_sum += 8U;
    if (*row_sum > 2048U) {
        uint32_t ns = 0U;
        size_t mb = 0U;
        uint32_t mbw = 0U;
        for (size_t s0 = 0; s0 < nsym; s0++) {
            uint32_t hv = w[s0] >> 1;
            w[s0] = hv;
            ns += hv;
            uint32_t bw = hv + 1U;
            if (bw > mbw) {
                mbw = bw;
                mb = s0;
            }
        }
        *row_sum = ns;
        *best_idx = mb;
        *best_w = mbw;
    } else {
        uint32_t nw = w[sym] + 1U;
        if (nw > *best_w || (nw == *best_w && sym < *best_idx)) {
            *best_w = nw;
            *best_idx = sym;
        }
    }
}

static int bwt_encode_block(const uint8_t* blk, size_t bn, AnsVec* dst) {
    int32_t n = (int32_t)bn;
    int32_t* s = (int32_t*)malloc((size_t)(n + 1) * sizeof(int32_t));
    int32_t* sa = (int32_t*)malloc((size_t)(n + 1) * sizeof(int32_t));
    uint8_t* bwt = (uint8_t*)malloc(bn);
    if (!s || !sa || !bwt) { free(s); free(sa); free(bwt); return -1; }

    uint8_t bitmap[32] = {0};
    for (size_t i = 0; i < bn; i++) {
        uint8_t v = blk[i];
        s[i] = (int32_t)v + 1;
        bitmap[v >> 3] |= (uint8_t)(1U << (v & 7));
    }
    s[n] = 0;

    if (sais_compute(s, sa, n + 1, 257) != 0) { free(s); free(sa); free(bwt); return -1; }
    free(s);

    uint32_t eof_idx = 0U;
    size_t bi = 0U;
    for (int32_t i = 0; i <= n; i++) {
        int32_t p = sa[i];
        if (p == 0) eof_idx = (uint32_t)i;
        else bwt[bi++] = blk[p - 1];
    }
    free(sa);

    uint8_t mtf[256];
    size_t sym_cnt = 0U;
    for (size_t c = 0; c < 256U; c++) {
        if (bitmap[c >> 3] & (1U << (c & 7))) mtf[sym_cnt++] = (uint8_t)c;
    }

    SymU16Vec syms = {NULL, 0, 0};
    size_t zrun = 0U;
    for (size_t i = 0; i < bn; i++) {
        uint8_t b = bwt[i];
        if (b == mtf[0]) {
            zrun++;
        } else {
            if (flush_zero_run(&syms, zrun) != 0) { free(syms.data); free(bwt); return -1; }
            zrun = 0U;
            size_t idx = 1U;
            while (idx < sym_cnt && mtf[idx] != b) idx++;
            for (size_t j = idx; j > 0U; j--) mtf[j] = mtf[j - 1U];
            mtf[0] = b;
            if (sym16_push(&syms, (uint16_t)(idx + 1U)) != 0) { free(syms.data); free(bwt); return -1; }
        }
    }
    free(bwt);
    if (flush_zero_run(&syms, zrun) != 0) { free(syms.data); return -1; }

    size_t nsym = sym_cnt + 1U; /* 0=RUNA, 1=RUNB, 2..sym_cnt = MTF rank 1..sym_cnt-1 */
    size_t rlen = syms.len;
    size_t base_q = rlen >> 2;
    size_t rem = rlen & 3U;
    size_t lane_len[4], lane_off[4], cur = 0U;
    for (size_t l = 0; l < 4U; l++) {
        lane_len[l] = base_q + (l < rem ? 1U : 0U);
        lane_off[l] = cur;
        cur += lane_len[l];
    }

    uint32_t (*weights)[3][257] = (uint32_t(*)[3][257])calloc(4U, sizeof(uint32_t[3][257]));
    QualRecVec recs[4] = {{NULL,0,0},{NULL,0,0},{NULL,0,0},{NULL,0,0}};
    AnsVec rev[4] = {{NULL,0,0},{NULL,0,0},{NULL,0,0},{NULL,0,0}};
    if (!weights) { free(syms.data); return -1; }

    uint16_t cums[257], freqs[257];
    uint32_t row_sum[4][3];
    size_t best_idx[4][3];
    uint32_t best_w[4][3];
    memset(row_sum, 0, sizeof(row_sum));
    for (size_t l = 0; l < 4U; l++) {
        for (size_t c = 0; c < 3U; c++) {
            best_idx[l][c] = 0U;
            best_w[l][c] = 1U;
        }
    }
    for (size_t l = 0; l < 4U; l++) {
        size_t ctx = 2U;
        for (size_t k = 0; k < lane_len[l]; k++) {
            uint16_t sym = syms.data[lane_off[l] + k];
            ans_normalize_cdf_with_tot(
                weights[l][ctx], nsym, row_sum[l][ctx] + (uint32_t)nsym, best_idx[l][ctx], cums, freqs);
            if (qrec_push(&recs[l], cums[sym], freqs[sym]) != 0) {
                free(weights); free(syms.data);
                for (size_t x = 0; x < 4U; x++) free(recs[x].data);
                return -1;
            }
            bwt_update_cdf_state(
                weights[l][ctx], nsym, (size_t)sym, &row_sum[l][ctx], &best_idx[l][ctx], &best_w[l][ctx]);
            ctx = (sym == 0U) ? 0U : (sym == 1U ? 1U : 2U);
        }
    }
    free(weights);
    free(syms.data);

    uint32_t state[4] = {ANS_L, ANS_L, ANS_L, ANS_L};
    int err = 0;
    for (size_t l = 0; l < 4U && !err; l++) {
        for (size_t i = recs[l].len; i > 0U; i--) {
            if (ans_encode_one(&state[l], &rev[l], recs[l].data[i - 1U].cum, recs[l].data[i - 1U].freq) != 0) {
                err = -1;
                break;
            }
        }
        free(recs[l].data);
        recs[l].data = NULL;
    }
    for (size_t l = 0; l < 4U; l++) free(recs[l].data);

    if (err != 0 ||
        vec_push_u32_le(dst, eof_idx) != 0 ||
        vec_append(dst, bitmap, 32U) != 0 ||
        vec_push_u32_le(dst, (uint32_t)rlen) != 0) {
        for (size_t l = 0; l < 4U; l++) free(rev[l].data);
        return -1;
    }
    for (size_t l = 0; l < 4U; l++) { if (vec_push_u32_le(dst, state[l]) != 0) { for (size_t x = 0; x < 4U; x++) free(rev[x].data); return -1; } }
    for (size_t l = 0; l < 4U; l++) { if (vec_push_u32_le(dst, (uint32_t)rev[l].len) != 0) { for (size_t x = 0; x < 4U; x++) free(rev[x].data); return -1; } }
    for (size_t l = 0; l < 4U; l++) {
        if (vec_reserve(dst, dst->len + rev[l].len) != 0) { for (size_t x = 0; x < 4U; x++) free(rev[x].data); return -1; }
        for (size_t i = rev[l].len; i > 0U; i--) dst->data[dst->len++] = rev[l].data[i - 1U];
        free(rev[l].data);
        rev[l].data = NULL;
    }
    return 0;
}

static int bwt_decode_block(const uint8_t* p_base, size_t p_avail, size_t* consumed, uint8_t* out_blk, size_t bn) {
    if (p_avail < 72U) return -1;
    uint32_t eof_idx = read_u32_le(p_base);
    if (eof_idx == 0U || (size_t)eof_idx > bn) return -1;
    const uint8_t* bitmap = p_base + 4U;
    uint32_t rlen = read_u32_le(p_base + 36U);
    uint32_t state[4], blen[4];
    for (size_t l = 0; l < 4U; l++) {
        state[l] = read_u32_le(p_base + 40U + l * 4U);
        blen[l] = read_u32_le(p_base + 56U + l * 4U);
    }
    size_t pos = 72U;
    const uint8_t* p_cur[4];
    const uint8_t* p_end[4];
    for (size_t l = 0; l < 4U; l++) {
        if (pos + (size_t)blen[l] > p_avail) return -1;
        p_cur[l] = p_base + pos;
        p_end[l] = p_base + pos + (size_t)blen[l];
        pos += (size_t)blen[l];
    }
    *consumed = pos;

    uint8_t mtf[256];
    size_t sym_cnt = 0U;
    for (size_t c = 0; c < 256U; c++) {
        if (bitmap[c >> 3] & (1U << (c & 7))) mtf[sym_cnt++] = (uint8_t)c;
    }
    if (sym_cnt == 0U) return -1;
    size_t nsym = sym_cnt + 1U;

    uint16_t* syms = (uint16_t*)malloc(((size_t)rlen + 1U) * sizeof(uint16_t));
    uint32_t (*weights)[3][257] = (uint32_t(*)[3][257])calloc(4U, sizeof(uint32_t[3][257]));
    if (!syms || !weights) { free(syms); free(weights); return -1; }

    size_t base_q = (size_t)rlen >> 2;
    size_t rem = (size_t)rlen & 3U;
    size_t lane_len[4], lane_off[4], cur = 0U;
    for (size_t l = 0; l < 4U; l++) {
        lane_len[l] = base_q + (l < rem ? 1U : 0U);
        lane_off[l] = cur;
        cur += lane_len[l];
    }

    uint16_t cums[257], freqs[257];
    uint32_t row_sum[4][3];
    size_t best_idx[4][3];
    uint32_t best_w[4][3];
    memset(row_sum, 0, sizeof(row_sum));
    for (size_t l = 0; l < 4U; l++) {
        for (size_t c = 0; c < 3U; c++) {
            best_idx[l][c] = 0U;
            best_w[l][c] = 1U;
        }
    }
    for (size_t l = 0; l < 4U; l++) {
        size_t ctx = 2U;
        for (size_t k = 0; k < lane_len[l]; k++) {
            ans_normalize_cdf_with_tot(
                weights[l][ctx], nsym, row_sum[l][ctx] + (uint32_t)nsym, best_idx[l][ctx], cums, freqs);
            size_t sym = 0U;
            if (ans_decode_find_and_step(&state[l], &p_cur[l], p_end[l], cums, freqs, nsym, &sym) != 0) {
                free(syms); free(weights); return -1;
            }
            syms[lane_off[l] + k] = (uint16_t)sym;
            bwt_update_cdf_state(
                weights[l][ctx], nsym, sym, &row_sum[l][ctx], &best_idx[l][ctx], &best_w[l][ctx]);
            ctx = (sym == 0U) ? 0U : (sym == 1U ? 1U : 2U);
        }
        if (state[l] != ANS_L || p_cur[l] != p_end[l]) { free(syms); free(weights); return -1; }
    }
    free(weights);

    uint8_t* bwt = (uint8_t*)malloc(bn);
    uint32_t* pred = (uint32_t*)malloc((bn + 1U) * sizeof(uint32_t));
    if (!bwt || !pred) { free(syms); free(bwt); free(pred); return -1; }

    size_t bi = 0U;
    size_t si = 0U;
    uint32_t counts[256] = {0};
    while (si < (size_t)rlen) {
        if (syms[si] <= 1U) {
            size_t run = 0U;
            size_t pow2 = 1U;
            while (si < (size_t)rlen && syms[si] <= 1U) {
                run += (syms[si] == 0U ? 1U : 2U) * pow2;
                pow2 <<= 1;
                si++;
            }
            if (bi + run > bn) { free(syms); free(bwt); free(pred); return -1; }
            uint8_t ch = mtf[0];
            counts[ch] += (uint32_t)run;
            for (size_t k = 0; k < run; k++) bwt[bi++] = ch;
        } else {
            size_t idx = (size_t)syms[si++] - 1U;
            if (idx >= sym_cnt) { free(syms); free(bwt); free(pred); return -1; }
            uint8_t ch = mtf[idx];
            for (size_t j = idx; j > 0U; j--) mtf[j] = mtf[j - 1U];
            mtf[0] = ch;
            if (bi >= bn) { free(syms); free(bwt); free(pred); return -1; }
            counts[ch]++;
            bwt[bi++] = ch;
        }
    }
    free(syms);
    if (bi != bn) { free(bwt); free(pred); return -1; }

    uint32_t start[256];
    uint32_t cum_c = 1U;
    for (size_t c = 0; c < 256U; c++) {
        start[c] = cum_c;
        cum_c += counts[c];
    }
    for (size_t r = 0U; r <= bn; r++) {
        if (r == (size_t)eof_idx) continue;
        size_t lpos = (r < (size_t)eof_idx) ? r : (r - 1U);
        uint8_t ch = bwt[lpos];
        pred[r] = start[ch]++;
    }
    uint32_t cur_row = 0U;
    for (size_t step = 0U; step < bn; step++) {
        size_t lpos = (cur_row < eof_idx) ? (size_t)cur_row : (size_t)(cur_row - 1U);
        out_blk[bn - 1U - step] = bwt[lpos];
        cur_row = pred[cur_row];
    }
    free(bwt);
    free(pred);
    return 0;
}

static int try_encode_qual_bwt(const uint8_t* quals_raw, const uint8_t* mismatch_mask, size_t len, AnsVec* out_bwt) {
    if (len < 256U) return -1;
    size_t n0 = 0U, n1 = 0U;
    if (mismatch_mask) {
        for (size_t i = 0; i < len; i++) {
            if (mismatch_mask[i] & 1U) n1++;
            else n0++;
        }
    }
    int split = (mismatch_mask && n0 >= 16U && n1 >= 16U) ? 1 : 0;
    if (vec_push(out_bwt, (uint8_t)ANS_MODE_RANS4_QUAL_BWT) != 0) return -1;
    if (vec_push(out_bwt, (uint8_t)split) != 0) return -1;

    if (!split) {
        for (size_t off = 0U; off < len; off += QUAL_BWT_BLK) {
            size_t bn = (len - off < QUAL_BWT_BLK) ? (len - off) : QUAL_BWT_BLK;
            if (bwt_encode_block(quals_raw + off, bn, out_bwt) != 0) return -1;
        }
    } else {
        uint8_t* s0 = (uint8_t*)malloc(n0);
        uint8_t* s1 = (uint8_t*)malloc(n1);
        if (!s0 || !s1) { free(s0); free(s1); return -1; }
        size_t i0 = 0U, i1 = 0U;
        for (size_t i = 0; i < len; i++) {
            if (mismatch_mask[i] & 1U) s1[i1++] = quals_raw[i];
            else s0[i0++] = quals_raw[i];
        }
        for (size_t off = 0U; off < n0; off += QUAL_BWT_BLK) {
            size_t bn = (n0 - off < QUAL_BWT_BLK) ? (n0 - off) : QUAL_BWT_BLK;
            if (bwt_encode_block(s0 + off, bn, out_bwt) != 0) { free(s0); free(s1); return -1; }
        }
        free(s0);
        for (size_t off = 0U; off < n1; off += QUAL_BWT_BLK) {
            size_t bn = (n1 - off < QUAL_BWT_BLK) ? (n1 - off) : QUAL_BWT_BLK;
            if (bwt_encode_block(s1 + off, bn, out_bwt) != 0) { free(s1); return -1; }
        }
        free(s1);
    }
    return 0;
}

static int decode_qual_bwt(const uint8_t* payload, size_t payload_len, const uint8_t* mismatch_mask, size_t expected_len, uint8_t* out) {
    if (payload_len < 2U) return -1;
    uint8_t split = payload[1];
    size_t pos = 2U;
    if (!split) {
        for (size_t off = 0U; off < expected_len; off += QUAL_BWT_BLK) {
            size_t bn = (expected_len - off < QUAL_BWT_BLK) ? (expected_len - off) : QUAL_BWT_BLK;
            size_t cons = 0U;
            if (bwt_decode_block(payload + pos, payload_len - pos, &cons, out + off, bn) != 0) return -1;
            pos += cons;
        }
        return (pos == payload_len) ? 0 : -1;
    }
    if (!mismatch_mask) return -1;
    size_t n0 = 0U, n1 = 0U;
    for (size_t i = 0; i < expected_len; i++) {
        if (mismatch_mask[i] & 1U) n1++;
        else n0++;
    }
    uint8_t* s0 = (uint8_t*)malloc(n0 ? n0 : 1U);
    uint8_t* s1 = (uint8_t*)malloc(n1 ? n1 : 1U);
    if (!s0 || !s1) { free(s0); free(s1); return -1; }
    for (size_t off = 0U; off < n0; off += QUAL_BWT_BLK) {
        size_t bn = (n0 - off < QUAL_BWT_BLK) ? (n0 - off) : QUAL_BWT_BLK;
        size_t cons = 0U;
        if (bwt_decode_block(payload + pos, payload_len - pos, &cons, s0 + off, bn) != 0) { free(s0); free(s1); return -1; }
        pos += cons;
    }
    for (size_t off = 0U; off < n1; off += QUAL_BWT_BLK) {
        size_t bn = (n1 - off < QUAL_BWT_BLK) ? (n1 - off) : QUAL_BWT_BLK;
        size_t cons = 0U;
        if (bwt_decode_block(payload + pos, payload_len - pos, &cons, s1 + off, bn) != 0) { free(s0); free(s1); return -1; }
        pos += cons;
    }
    if (pos != payload_len) { free(s0); free(s1); return -1; }
    size_t i0 = 0U, i1 = 0U;
    for (size_t i = 0; i < expected_len; i++) {
        out[i] = (mismatch_mask[i] & 1U) ? s1[i1++] : s0[i0++];
    }
    free(s0);
    free(s1);
    return 0;
}

static inline void compute_4_lanes(size_t total_len, size_t lane_len[4], size_t lane_off[4]) {
    size_t base_q = total_len >> 2;
    size_t rem = total_len & 3U;
    size_t cur = 0U;
    for (size_t l = 0; l < 4U; l++) {
        lane_len[l] = base_q + (l < rem ? 1U : 0U);
        lane_off[l] = cur;
        cur += lane_len[l];
    }
}

size_t hgs_ans_qual_bwt_block_cap(void) {
    return (size_t)QUAL_BWT_BLK;
}

size_t hgs_ans_qual_bwt_peek_block_len(const uint8_t* block_payload, size_t avail) {
    if (!block_payload || avail < 72U) return 0U;
    uint64_t tot = 72ULL;
    for (size_t l = 0; l < 4U; l++) {
        tot += (uint64_t)read_u32_le(block_payload + 56U + l * 4U);
    }
    return (tot <= (uint64_t)avail) ? (size_t)tot : 0U;
}

int hgs_ans_qual_bwt_encode_block(
    const uint8_t* blk,
    size_t bn,
    uint8_t** out_buf,
    size_t* out_len) {
    if (!out_buf || !out_len) return -1;
    *out_buf = NULL;
    *out_len = 0U;
    if (!blk || bn == 0U) return -1;
    AnsVec v = {NULL, 0, 0};
    if (bwt_encode_block(blk, bn, &v) != 0) {
        free(v.data);
        return -1;
    }
    *out_buf = v.data;
    *out_len = v.len;
    return 0;
}

int hgs_ans_qual_bwt_decode_block(
    const uint8_t* block_payload,
    size_t block_avail,
    size_t expected_bn,
    size_t* out_consumed,
    uint8_t* out_blk) {
    if (!block_payload || !out_consumed || !out_blk || expected_bn == 0U) return -1;
    return bwt_decode_block(block_payload, block_avail, out_consumed, out_blk, expected_bn);
}

int hgs_ans4_encode_qual_ctx_lane(
    const uint8_t* quals_raw,
    const uint8_t* mismatch_mask,
    const uint8_t* homo_bins,
    size_t total_len,
    uint32_t read_len_hint,
    size_t lane,
    uint32_t* out_state,
    uint8_t** out_fwd_buf,
    size_t* out_fwd_len) {
    if (!out_state || !out_fwd_buf || !out_fwd_len || lane >= 4U) return -1;
    *out_state = ANS_L;
    *out_fwd_buf = NULL;
    *out_fwd_len = 0U;

    size_t lane_len[4], lane_off[4];
    compute_4_lanes(total_len, lane_len, lane_off);

    QualModel* m = (QualModel*)calloc(1U, sizeof(QualModel));
    if (!m) return -1;
    QualRecVec recs = {NULL, 0, 0};
    AnsVec rev = {NULL, 0, 0};

    uint32_t uniform16_w[16] = {0};
    uint16_t u16_cum[16], u16_freq[16];
    ans_normalize_cdf(uniform16_w, 16U, u16_cum, u16_freq);

    uint8_t prev_q = 66U;
    for (size_t k = 0; k < lane_len[lane]; k++) {
        size_t gpos = lane_off[lane] + k;
        uint8_t q = quals_raw[gpos];
        uint8_t mm = mismatch_mask ? (mismatch_mask[gpos] & 1U) : 0U;
        uint8_t hb = homo_bins ? (homo_bins[gpos] & 3U) : 0U;
        uint32_t ctx = qual_ctx(prev_q, mm, hb, gpos, read_len_hint);

        uint32_t sw[2] = {m->same[ctx][0], m->same[ctx][1]};
        uint16_t scum[2], sfreq[2];
        ans_normalize_cdf(sw, 2U, scum, sfreq);
        int is_same = (q == prev_q) ? 1 : 0;
        if (qrec_push(&recs, scum[is_same], sfreq[is_same]) != 0) goto fail;
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
                if (qrec_push(&recs, ccum[cbin], cfreq[cbin]) != 0) goto fail;
                if (++m->coarse[ctx][cbin] > 4096U) {
                    for (int c = 0; c < 9; c++) m->coarse[ctx][c] >>= 1;
                }
                uint32_t fctx = qual_fine_ctx(cbin, mm, prev_q);
                uint32_t fw[8];
                for (int f = 0; f < 8; f++) fw[f] = m->fine[fctx][f];
                uint16_t fcum[8], ffreq[8];
                ans_normalize_cdf(fw, 8U, fcum, ffreq);
                if (qrec_push(&recs, fcum[fbin], ffreq[fbin]) != 0) goto fail;
                if (++m->fine[fctx][fbin] > 4096U) {
                    for (int f = 0; f < 8; f++) m->fine[fctx][f] >>= 1;
                }
            } else {
                if (qrec_push(&recs, ccum[8], cfreq[8]) != 0) goto fail;
                if (++m->coarse[ctx][8] > 4096U) {
                    for (int c = 0; c < 9; c++) m->coarse[ctx][c] >>= 1;
                }
                uint32_t hi = (uint32_t)((q >> 4) & 0xFU);
                uint32_t lo = (uint32_t)(q & 0xFU);
                if (qrec_push(&recs, u16_cum[hi], u16_freq[hi]) != 0) goto fail;
                if (qrec_push(&recs, u16_cum[lo], u16_freq[lo]) != 0) goto fail;
            }
        }
        prev_q = q;
    }

    uint32_t st = ANS_L;
    for (size_t i = recs.len; i > 0; i--) {
        if (ans_encode_one(&st, &rev, recs.data[i - 1U].cum, recs.data[i - 1U].freq) != 0) goto fail;
    }
    free(recs.data);
    free(m);

    uint8_t* fwd = (uint8_t*)malloc(rev.len ? rev.len : 1U);
    if (!fwd) {
        free(rev.data);
        return -1;
    }
    for (size_t i = 0; i < rev.len; i++) {
        fwd[i] = rev.data[rev.len - 1U - i];
    }
    free(rev.data);

    *out_state = st;
    *out_fwd_buf = fwd;
    *out_fwd_len = rev.len;
    return 0;

fail:
    free(recs.data);
    free(rev.data);
    free(m);
    return -1;
}

int hgs_ans4_decode_qual_ctx_lane(
    const uint8_t* lane_payload,
    size_t lane_payload_len,
    uint32_t init_state,
    const uint8_t* mismatch_mask,
    const uint8_t* homo_bins,
    size_t total_len,
    uint32_t read_len_hint,
    size_t lane,
    uint8_t* out_lane_buf) {
    if (lane >= 4U || (!lane_payload && lane_payload_len > 0U)) return -1;

    size_t lane_len[4], lane_off[4];
    compute_4_lanes(total_len, lane_len, lane_off);
    if (lane_len[lane] > 0U && !out_lane_buf) return -1;

    QualModel* m = (QualModel*)calloc(1U, sizeof(QualModel));
    if (!m) return -1;

    uint32_t uniform16_w[16] = {0};
    uint16_t u16_cum[16], u16_freq[16];
    ans_normalize_cdf(uniform16_w, 16U, u16_cum, u16_freq);

    uint32_t st = init_state;
    const uint8_t* p_cur = lane_payload;
    const uint8_t* p_end = lane_payload + lane_payload_len;
    uint8_t prev_q = 66U;

    for (size_t k = 0; k < lane_len[lane]; k++) {
        size_t gpos = lane_off[lane] + k;
        uint8_t mm = mismatch_mask ? (mismatch_mask[gpos] & 1U) : 0U;
        uint8_t hb = homo_bins ? (homo_bins[gpos] & 3U) : 0U;
        uint32_t ctx = qual_ctx(prev_q, mm, hb, gpos, read_len_hint);

        uint32_t sw[2] = {m->same[ctx][0], m->same[ctx][1]};
        uint16_t scum[2], sfreq[2];
        ans_normalize_cdf(sw, 2U, scum, sfreq);
        size_t is_same = 0;
        if (ans_decode_find_and_step(&st, &p_cur, p_end, scum, sfreq, 2U, &is_same) != 0) goto fail;
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
            if (ans_decode_find_and_step(&st, &p_cur, p_end, ccum, cfreq, 9U, &cbin) != 0) goto fail;
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
                if (ans_decode_find_and_step(&st, &p_cur, p_end, fcum, ffreq, 8U, &fbin) != 0) goto fail;
                if (++m->fine[fctx][fbin] > 4096U) {
                    for (int f = 0; f < 8; f++) m->fine[fctx][f] >>= 1;
                }
                q = (uint8_t)(33U + (cbin << 3) + fbin);
            } else {
                size_t hi = 0, lo = 0;
                if (ans_decode_find_and_step(&st, &p_cur, p_end, u16_cum, u16_freq, 16U, &hi) != 0) goto fail;
                if (ans_decode_find_and_step(&st, &p_cur, p_end, u16_cum, u16_freq, 16U, &lo) != 0) goto fail;
                q = (uint8_t)((hi << 4) | lo);
            }
        }
        out_lane_buf[k] = q;
        prev_q = q;
    }

    if (st != ANS_L || p_cur != p_end) goto fail;
    free(m);
    return 0;

fail:
    free(m);
    return -1;
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

    uint32_t state[4] = {ANS_L, ANS_L, ANS_L, ANS_L};
    uint8_t* fwd[4] = {NULL, NULL, NULL, NULL};
    size_t flen[4] = {0, 0, 0, 0};
    int err = 0;
    for (size_t l = 0; l < 4U; l++) {
        if (hgs_ans4_encode_qual_ctx_lane(
                quals_raw, mismatch_mask, homo_bins, len, read_len_hint, l, &state[l], &fwd[l], &flen[l]) != 0) {
            err = -1;
            break;
        }
    }

    size_t rans_size = 33U;
    for (size_t l = 0; l < 4U; l++) rans_size += flen[l];

    AnsVec bwt_out = {NULL, 0, 0};
    if (len >= 256U) {
        if (try_encode_qual_bwt(quals_raw, mismatch_mask, len, &bwt_out) != 0) {
            free(bwt_out.data);
            bwt_out.data = NULL;
            bwt_out.len = 0;
        }
    }

    size_t raw_size = len + 1U;
    size_t ctx_size = (err == 0 && rans_size < raw_size) ? rans_size : SIZE_MAX;
    size_t bwt_size = (bwt_out.data != NULL && bwt_out.len < raw_size) ? bwt_out.len : SIZE_MAX;

    if (bwt_size < ctx_size && bwt_size < raw_size) {
        for (size_t l = 0; l < 4U; l++) free(fwd[l]);
        *out_buf = bwt_out.data;
        *out_len = bwt_out.len;
        return 0;
    }
    free(bwt_out.data);

    if (err != 0 || rans_size >= raw_size) {
        for (size_t l = 0; l < 4U; l++) free(fwd[l]);
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
    for (size_t l = 0; l < 4U; l++) { if (vec_push_u32_le(&out, (uint32_t)flen[l]) != 0) goto fail; }
    for (size_t l = 0; l < 4U; l++) {
        if (vec_append(&out, fwd[l], flen[l]) != 0) goto fail;
        free(fwd[l]);
        fwd[l] = NULL;
    }
    *out_buf = out.data;
    *out_len = out.len;
    return 0;

fail:
    for (size_t l = 0; l < 4U; l++) free(fwd[l]);
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
    if (mode == ANS_MODE_RANS4_QUAL_BWT) {
        if (decode_qual_bwt(payload, payload_len, mismatch_mask, expected_len, out) != 0) {
            free(out);
            return -1;
        }
        *out_buf = out;
        *out_len = expected_len;
        return 0;
    }
    if (mode != ANS_MODE_RANS4_QUAL || payload_len < 33U) { free(out); return -1; }

    uint32_t state[4], blen[4];
    const uint8_t* p_cur[4];
    for (size_t l = 0; l < 4U; l++) {
        state[l] = read_u32_le(payload + 1U + l * 4U);
        blen[l] = read_u32_le(payload + 17U + l * 4U);
    }
    size_t pos = 33U;
    for (size_t l = 0; l < 4U; l++) {
        if (pos + (size_t)blen[l] > payload_len) { free(out); return -1; }
        p_cur[l] = payload + pos;
        pos += (size_t)blen[l];
    }

    size_t lane_len[4], lane_off[4];
    compute_4_lanes(expected_len, lane_len, lane_off);

    for (size_t l = 0; l < 4U; l++) {
        if (hgs_ans4_decode_qual_ctx_lane(
                p_cur[l], (size_t)blen[l], state[l], mismatch_mask, homo_bins, expected_len, read_len_hint, l, out + lane_off[l]) != 0) {
            free(out);
            return -1;
        }
    }

    *out_buf = out;
    *out_len = expected_len;
    return 0;
}
