#ifndef HGS_CPP_ENGINE_HPP
#define HGS_CPP_ENGINE_HPP

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
#define HGS_CPP_NODISCARD [[nodiscard]]
extern "C" {
#else
#define HGS_CPP_NODISCARD
#endif

/** Owned heap buffer returned across FFI; free strictly via hgs_cpp_free_buffer(). */
typedef struct {
    uint8_t* data;
    size_t len;
} HgsOwnedBytes;

void hgs_cpp_free_buffer(HgsOwnedBytes buf);

/** Computes bottom-k canonical primary open-syncmer hash sketch over 2-bit bytes (0..3). */
HGS_CPP_NODISCARD int hgs_cpp_compute_bottom_sketch(
    const uint8_t* two_bit_bases,
    size_t len,
    size_t sketch_size,
    uint64_t** out_sketch,
    size_t* out_len);

void hgs_cpp_free_u64(uint64_t* ptr);

/** Estimates sorted bottom-sketch Jaccard similarity in [0, 1]. */
HGS_CPP_NODISCARD double hgs_cpp_sketch_jaccard(const uint64_t* a, size_t a_len, const uint64_t* b, size_t b_len);

/** Opaque 64-bit Syncmer CSR Index supporting >2.14 Gbp banks. */
typedef struct HgsSyncmerIndex64 HgsSyncmerIndex64;

/**
 * Builds a 64-bit CSR canonical open-syncmer index over num_donors donors concatenated into bank.
 * Each donor i has length donor_lens[i] and bases packed in donor_two_bit_ptrs[i] (2-bit 0..3).
 */
HGS_CPP_NODISCARD HgsSyncmerIndex64* hgs_cpp_index_build(
    const uint8_t* const* donor_two_bit_ptrs,
    const uint64_t* donor_lens,
    size_t num_donors,
    int subsample_shift,
    int build_secondary,
    int threads);

void hgs_cpp_index_free(HgsSyncmerIndex64* idx);

/** Telemetry & encoded buffers returned from Stage-1 mosaic encoding. */
typedef struct {
    uint8_t* mosaic_data;
    size_t mosaic_len;
    uint8_t* novelty_data;
    size_t novelty_len;
    uint64_t panel_match_bases;
    uint64_t self_match_bases;
    uint64_t tandem_match_bases;
    uint64_t panel_tokens;
    uint64_t self_tokens;
    uint64_t tandem_tokens;
    uint64_t literal_runs;
    uint64_t indel_splices;
    uint64_t substitutions;
    int status;
} HgsMosaicResult;

void hgs_cpp_free_mosaic_result(HgsMosaicResult* res);

/**
 * Encodes target[0..tgt_len) (2-bit 0..3) into standard .hgs2 Mosaic token stream + extracted novelty bases.
 * Optionally uses panel_idx (for donor ids < max_donor_exclusive) and donor_sim_order[0..sim_order_len).
 */
HGS_CPP_NODISCARD HgsMosaicResult hgs_cpp_encode_mosaic(
    const uint8_t* target_two_bit,
    uint64_t tgt_len,
    const HgsSyncmerIndex64* panel_idx,
    size_t max_donor_exclusive,
    int64_t preferred_twin_donor,
    const int32_t* donor_sim_order,
    size_t sim_order_len,
    int allow_self,
    int allow_tandem,
    int allow_approx,
    int allow_splice,
    int allow_kimura,
    int allow_wobble,
    int threads);

/**
 * Decodes standard .hgs2 Mosaic stream into 2-bit bases (0..3).
 * donor_2bit_ptrs[0..num_donors) points to previously decoded 2-bit donor arrays of lengths donor_lens[].
 */
HGS_CPP_NODISCARD int hgs_cpp_decode_mosaic(
    const uint8_t* mosaic_bytes,
    size_t mosaic_len,
    const uint8_t* novelty_bytes,
    size_t novelty_len,
    const uint8_t* const* donor_2bit_ptrs,
    const uint64_t* donor_lens,
    size_t max_donor_exclusive,
    uint8_t** out_two_bit,
    uint64_t* out_len);

/** Parses [target_len, novelty_len] from a .hgs2 mosaic payload header cleanly. */
HGS_CPP_NODISCARD int hgs_cpp_mosaic_peek_lengths(const uint8_t* mosaic_bytes, size_t mosaic_len, uint64_t* out_tgt_len, uint64_t* out_nov_len);

/** Stage-2 Unified Novelty Encoder (exact .hgs2 FastNoveltyCoder wire format: FastCRM + Turbo + Raw2Bit). */
HGS_CPP_NODISCARD int hgs_cpp_novelty_encode(
    const uint8_t* two_bit_bases,
    uint64_t len,
    int threads,
    int bounded_mem,
    HgsOwnedBytes* out_payload);

/** Stage-2 Unified Novelty Decoder (exact .hgs2 FastNoveltyCoder wire format). */
HGS_CPP_NODISCARD int hgs_cpp_novelty_decode(
    const uint8_t* payload,
    size_t payload_len,
    uint64_t expected_len,
    int threads,
    int bounded_mem,
    HgsOwnedBytes* out_two_bit);

/** Compressed-domain evolutionary statistics decoded directly from Mosaic tokens. */
typedef struct {
    uint64_t target_len;
    uint64_t novelty_len;
    uint64_t panel_match_bases;
    uint64_t self_match_bases;
    uint64_t tandem_match_bases;
    uint64_t literal_bases;
    uint64_t panel_tokens;
    uint64_t self_tokens;
    uint64_t tandem_tokens;
    uint64_t literal_runs;
    uint64_t snps;
    uint64_t transitions;
    uint64_t transversions;
    uint64_t wobble_mod3_hits;
    uint64_t splices;
    uint64_t inserted_bases;
    uint64_t deleted_bases;
    uint64_t recombination_switches;
} HgsMosaicStats;

/**
 * Exact variant call extracted directly from compressed Mosaic tokens without FASTA emission.
 * kind: 0 = SNV, 1 = INS, 2 = DEL, 3 = COMPLEX_INDEL, 4 = LITERAL_SPAN
 */
typedef struct {
    uint8_t kind;
    uint8_t strand;
    uint8_t is_transition;
    uint8_t prev_mod3;
    uint32_t donor_id;
    uint64_t ref_pos;
    uint64_t tgt_pos;
    uint32_t ref_len;
    uint32_t alt_len;
    uint64_t ref_seq_off;
    uint64_t alt_seq_off;
} HgsVariantRecord;

/** Token span descriptor for O(tokens) compressed-domain k-mer and ancestry queries. */
typedef struct {
    uint8_t kind;     /* 0 = LIT, 1 = PANEL, 2 = SELF, 3 = TANDEM */
    uint8_t strand;   /* 0 = +, 1 = - */
    uint8_t has_edits;
    uint8_t reserved;
    uint32_t donor_id;
    uint64_t ref_pos;
    uint64_t tgt_pos;
    uint64_t tgt_len;
    uint64_t ref_len;
} HgsSpanEvent;

typedef struct {
    HgsMosaicStats stats;
    HgsVariantRecord* variants;
    size_t num_variants;
    HgsSpanEvent* spans;
    size_t num_spans;
    uint8_t* var_seq_pool; /* 2-bit bases (0..3) referenced by HgsVariantRecord offsets */
    size_t var_seq_pool_len;
    int status;
} HgsMosaicInspectResult;

void hgs_cpp_free_inspect_result(HgsMosaicInspectResult* res);

/**
 * Inspects a compressed .hgs2/.hgs3 Mosaic stream in compressed domain, computing exact
 * evolutionary telemetry and optionally extracting exact SNV/INDEL variants & token spans.
 */
HGS_CPP_NODISCARD HgsMosaicInspectResult hgs_cpp_mosaic_inspect(
    const uint8_t* mosaic_bytes,
    size_t mosaic_len,
    const uint8_t* novelty_bytes,
    size_t novelty_len,
    const uint8_t* const* donor_2bit_ptrs,
    const uint64_t* donor_lens,
    size_t max_donor_exclusive,
    int collect_variants,
    int collect_spans);

/** Parallel SIMD 4-stripe rANS Novelty Encoder for .hgs3 containers (with Law IV bypass). */
HGS_CPP_NODISCARD int hgs_cpp_ans_novelty_encode(
    const uint8_t* two_bit_bases,
    uint64_t len,
    int threads,
    HgsOwnedBytes* out_payload);

/** Parallel SIMD 4-stripe rANS Novelty Decoder for .hgs3 containers. */
HGS_CPP_NODISCARD int hgs_cpp_ans_novelty_decode(
    const uint8_t* payload,
    size_t payload_len,
    uint64_t expected_len,
    int threads,
    HgsOwnedBytes* out_two_bit);

#ifdef __cplusplus
}
#endif

#endif /* HGS_CPP_ENGINE_HPP */
