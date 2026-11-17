#include "hgs_cpp_engine.hpp"
#include "../c_kernels/hgs_simd.h"

#include <algorithm>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <thread>
#include <utility>
#include <vector>

namespace {

// ============================================================================
// Binary VarInt / Buffer Helpers (Unsigned LEB128 up to 64 bits)
// ============================================================================

struct ByteVecOut {
    std::vector<uint8_t> buf;

    void reserve(size_t n) { buf.reserve(n); }
    void write_byte(uint8_t b) { buf.push_back(b); }
    void write_bytes(const uint8_t* src, size_t len) {
        if (len > 0) buf.insert(buf.end(), src, src + len);
    }
    void write_var_u64(uint64_t v) {
        while (v >= 0x80ULL) {
            buf.push_back(static_cast<uint8_t>((v & 0x7FULL) | 0x80ULL));
            v >>= 7;
        }
        buf.push_back(static_cast<uint8_t>(v));
    }
};

struct ByteSpanIn {
    const uint8_t* data;
    size_t len;
    size_t pos = 0;

    bool has_remaining() const { return pos < len; }
    uint8_t read_byte(bool& ok) {
        if (pos >= len) { ok = false; return 0; }
        return data[pos++];
    }
    uint64_t read_var_u64(bool& ok) {
        uint64_t res = 0;
        int shift = 0;
        while (true) {
            uint8_t b = read_byte(ok);
            if (!ok) return 0;
            res |= static_cast<uint64_t>(b & 0x7FU) << shift;
            if ((b & 0x80U) == 0) return res;
            shift += 7;
            if (shift > 63) { ok = false; return 0; }
        }
    }
};

// ============================================================================
// Packed 2-Bit Sequence Store (64-bit offset safe)
// ============================================================================

struct Packed2Bit {
    std::vector<uint64_t> words;
    uint64_t length = 0;

    static Packed2Bit from_bytes(const uint8_t* bases, uint64_t len) {
        Packed2Bit p;
        p.length = len;
        size_t nw = static_cast<size_t>((len + 31ULL) >> 5);
        p.words.assign(nw + 2, 0ULL);
        if (len > 0) hgs_pack_2bit_words(bases, static_cast<size_t>(len), p.words.data());
        return p;
    }

    uint8_t base_at(uint64_t pos) const {
        return static_cast<uint8_t>((words[pos >> 5] >> ((pos & 31ULL) << 1)) & 3ULL);
    }

    uint8_t base_at_strand(uint64_t offset, uint64_t seg_len, uint8_t strand, uint64_t rel_pos) const {
        if (strand == 0) return base_at(offset + rel_pos);
        return static_cast<uint8_t>(base_at(offset + seg_len - 1ULL - rel_pos) ^ 3U);
    }

    uint64_t word32_fwd(uint64_t pos) const {
        return hgs_word32_forward(words.data(), pos);
    }

    uint64_t word32_strand(uint64_t offset, uint64_t seg_len, uint8_t strand, uint64_t rel_pos) const {
        if (strand == 0) return word32_fwd(offset + rel_pos);
        uint64_t fwd_top = offset + seg_len - 1ULL - rel_pos;
        if (fwd_top >= offset + 31ULL) {
            return hgs_revcomp_word32(word32_fwd(fwd_top - 31ULL));
        }
        uint64_t acc = 0;
        uint64_t lim = std::min<uint64_t>(32ULL, fwd_top + 1ULL - offset);
        for (uint64_t j = 0; j < lim; j++) {
            uint8_t b = static_cast<uint8_t>(base_at(fwd_top - j) ^ 3U);
            acc |= static_cast<uint64_t>(b) << (j << 1);
        }
        return acc;
    }
};

constexpr int K_PRIMARY = 16;
constexpr int S_PRIMARY = 10;
constexpr int W_PRIMARY = 7;
constexpr int T_PRIMARY = 3;

constexpr int K_SECONDARY = 13;
constexpr int S_SECONDARY = 9;
constexpr int W_SECONDARY = 5;
constexpr int T_SECONDARY = 2;

constexpr uint64_t BLOCK_SHIFT = 14; // 16384-base blocks
constexpr uint8_t MAX_WITHIN_DONOR_BUCKET_OCC = 16;

// ============================================================================
// Exact Mosaic Codec Binary Range Coder & Probability Models (wire-compatible)
// ============================================================================

constexpr int MOS_PROB_BITS = 12;
constexpr int MOS_PROB_SCALE = 1 << MOS_PROB_BITS; // 4096
constexpr int MOS_PROB_HALF = 1 << (MOS_PROB_BITS - 1); // 2048
constexpr uint32_t MOS_TOP = 1U << 24;

constexpr int MIN_MATCH_LEN = 20;
constexpr int MIN_SPLICE_ANCHOR = 12;
constexpr int MIN_TANDEM_LEN = 18;
constexpr int MAX_TANDEM_PERIOD = 64;
constexpr int MAX_BUCKET_PROBES = 48;
constexpr int MAX_SPLICE_CHAIN = 1024;

struct BitRangeEnc {
    ByteVecOut& out;
    uint64_t low = 0;
    uint32_t range = 0xFFFFFFFFU;
    uint8_t cache = 0;
    uint64_t cache_size = 1;

    explicit BitRangeEnc(ByteVecOut& o) : out(o) {}

    void shift_low() {
        uint64_t lo = low;
        if ((lo >> 32) != 0 || static_cast<uint32_t>(lo) < 0xFF000000U) {
            uint32_t carry = static_cast<uint32_t>(lo >> 32);
            uint32_t temp = cache;
            do {
                out.write_byte(static_cast<uint8_t>(temp + carry));
                temp = 0xFF;
            } while (--cache_size != 0);
            cache = static_cast<uint8_t>((lo >> 24) & 0xFFU);
        }
        cache_size++;
        low = (lo << 8) & 0xFFFFFFFFULL;
    }

    void encode_bit(int* probs, int idx, int bit) {
        int p = probs[idx];
        uint32_t bound = (range >> MOS_PROB_BITS) * static_cast<uint32_t>(p);
        if (bit == 1) {
            range = bound;
            probs[idx] = p + ((MOS_PROB_SCALE - p + 20) >> 5);
        } else {
            low += bound;
            range -= bound;
            probs[idx] = p - ((p + 20) >> 5);
        }
        while (range < MOS_TOP) {
            range <<= 8;
            shift_low();
        }
    }

    void encode_direct_bit(int bit) {
        range >>= 1;
        if (bit == 1) low += range;
        if (range < MOS_TOP) {
            range <<= 8;
            shift_low();
        }
    }

    void finish() {
        for (int i = 0; i < 5; i++) shift_low();
    }
};

struct BitRangeDec {
    ByteSpanIn& in;
    uint32_t range = 0xFFFFFFFFU;
    uint32_t code = 0;

    explicit BitRangeDec(ByteSpanIn& i) : in(i) {
        for (int k = 0; k < 5; k++) code = (code << 8) | next_byte();
    }

    uint32_t next_byte() {
        bool ok = true;
        return in.has_remaining() ? in.read_byte(ok) : 0U;
    }

    int decode_bit(int* probs, int idx) {
        int p = probs[idx];
        uint32_t bound = (range >> MOS_PROB_BITS) * static_cast<uint32_t>(p);
        int bit;
        if (code < bound) {
            range = bound;
            probs[idx] = p + ((MOS_PROB_SCALE - p + 20) >> 5);
            bit = 1;
        } else {
            code -= bound;
            range -= bound;
            probs[idx] = p - ((p + 20) >> 5);
            bit = 0;
        }
        while (range < MOS_TOP) {
            range <<= 8;
            code = (code << 8) | next_byte();
        }
        return bit;
    }

    int decode_direct_bit() {
        range >>= 1;
        int bit;
        if (code >= range) {
            code -= range;
            bit = 1;
        } else {
            bit = 0;
        }
        if (range < MOS_TOP) {
            range <<= 8;
            code = (code << 8) | next_byte();
        }
        return bit;
    }
};

inline int clz32_pos(uint32_t v) {
    return std::countl_zero(v);
}

struct MosaicModels {
    std::vector<int> is_literal = std::vector<int>(4, MOS_PROB_HALF);
    std::vector<int> is_panel_match = std::vector<int>(2, MOS_PROB_HALF);
    std::vector<int> is_tandem = std::vector<int>(2, MOS_PROB_HALF);
    std::vector<int> self_strand = std::vector<int>(2, MOS_PROB_HALF);
    std::vector<int> panel_strand = std::vector<int>(2, MOS_PROB_HALF);

    std::vector<int> mtf_hit0 = std::vector<int>(2, MOS_PROB_HALF);
    std::vector<int> mtf_hit13 = std::vector<int>(2, MOS_PROB_HALF);
    std::vector<int> mtf_choice = std::vector<int>(4, MOS_PROB_HALF);
    std::vector<int> donor_tree = std::vector<int>(4096, MOS_PROB_HALF);

    std::vector<int> use_pos_delta = std::vector<int>(2, MOS_PROB_HALF);
    std::vector<int> use_collinear_pred = std::vector<int>(2, MOS_PROB_HALF);

    std::vector<int> has_snps = std::vector<int>(8, 600);
    std::vector<int> more_snps = std::vector<int>(16, MOS_PROB_HALF);
    std::vector<int> has_splice = std::vector<int>(8, 96);
    std::vector<int> splice_is_pure_del = std::vector<int>(2, MOS_PROB_HALF);
    std::vector<int> splice_is_pure_ins = std::vector<int>(2, 3800);

    std::vector<int> kimura_ti = std::vector<int>(256, 2700);
    std::vector<int> kimura_tv = std::vector<int>(256, MOS_PROB_HALF);
    int ablate_base4[4][4];

    std::vector<int> wobble_is_third = std::vector<int>(32, 1800);
    std::vector<int> wobble_rem01 = std::vector<int>(32, MOS_PROB_HALF);

    int lit_len_slots[1][64];
    int match_len_slots[2][64];
    int self_dist_slots[2][64];
    int pos_delta_slots[2][64];
    int raw_pos_slots[1][64];
    int gap_quot_slots[4][64];
    int splice_ins_slots[1][64];
    int splice_del_slots[1][64];
    int low_bit_probs[8][128];
    int tandem_period_tree[64];

    MosaicModels() {
        auto fill_arr = [](int* a, size_t n) { std::fill(a, a + n, MOS_PROB_HALF); };
        for (int r = 0; r < 4; r++) fill_arr(ablate_base4[r], 4);
        fill_arr(lit_len_slots[0], 64);
        for (int r = 0; r < 2; r++) {
            fill_arr(match_len_slots[r], 64);
            fill_arr(self_dist_slots[r], 64);
            fill_arr(pos_delta_slots[r], 64);
        }
        fill_arr(raw_pos_slots[0], 64);
        for (int r = 0; r < 4; r++) fill_arr(gap_quot_slots[r], 64);
        fill_arr(splice_ins_slots[0], 64);
        fill_arr(splice_del_slots[0], 64);
        for (int r = 0; r < 8; r++) fill_arr(low_bit_probs[r], 128);
        fill_arr(tandem_period_tree, 64);
    }

    void write_fixed_tree(BitRangeEnc& enc, int* tree, int sym, int bits) {
        int node = 1;
        for (int b = bits - 1; b >= 0; b--) {
            int bit = (sym >> b) & 1;
            enc.encode_bit(tree, node, bit);
            node = (node << 1) | bit;
        }
    }

    int read_fixed_tree(BitRangeDec& dec, int* tree, int bits) {
        int node = 1, val = 0;
        for (int b = 0; b < bits; b++) {
            int bit = dec.decode_bit(tree, node);
            val = (val << 1) | bit;
            node = (node << 1) | bit;
        }
        return val;
    }

    void write_bounded_int(BitRangeEnc& enc, int* tree, uint32_t val, uint32_t max_excl) {
        if (max_excl <= 1) return;
        int bits = 32 - clz32_pos(max_excl - 1U);
        if (bits > 11) {
            write_count(enc, raw_pos_slots[0], 4, val);
            return;
        }
        int node = 1;
        uint32_t prefix = 0;
        for (int b = bits - 1; b >= 0; b--) {
            uint32_t zero_upper = (prefix << (b + 1)) | (1U << b);
            if (zero_upper >= max_excl) {
                prefix <<= 1;
                node <<= 1;
                continue;
            }
            int bit = static_cast<int>((val >> b) & 1U);
            enc.encode_bit(tree, node, bit);
            prefix = (prefix << 1) | static_cast<uint32_t>(bit);
            node = (node << 1) | bit;
        }
    }

    uint32_t read_bounded_int(BitRangeDec& dec, int* tree, uint32_t max_excl) {
        if (max_excl <= 1) return 0;
        int bits = 32 - clz32_pos(max_excl - 1U);
        if (bits > 11) {
            return read_count(dec, raw_pos_slots[0], 4);
        }
        int node = 1;
        uint32_t prefix = 0;
        for (int b = bits - 1; b >= 0; b--) {
            uint32_t zero_upper = (prefix << (b + 1)) | (1U << b);
            if (zero_upper >= max_excl) {
                prefix <<= 1;
                node <<= 1;
                continue;
            }
            int bit = dec.decode_bit(tree, node);
            prefix = (prefix << 1) | static_cast<uint32_t>(bit);
            node = (node << 1) | bit;
        }
        return prefix;
    }

    void write_count(BitRangeEnc& enc, int* tree, int low_ctx, uint32_t value) {
        if (value < 4) {
            write_fixed_tree(enc, tree, static_cast<int>(value), 6);
            return;
        }
        int msb = 31 - clz32_pos(value);
        int extra_bits = msb - 1;
        int second_bit = static_cast<int>((value >> extra_bits) & 1U);
        int slot = ((msb - 1) << 1) + 2 + second_bit;
        write_fixed_tree(enc, tree, slot, 6);
        uint32_t rem = value & ((1U << extra_bits) - 1U);
        int* low_row = low_bit_probs[low_ctx];
        for (int b = extra_bits - 1; b >= 0; b--) {
            int bit = static_cast<int>((rem >> b) & 1U);
            if (b < 4 && slot < 32) {
                enc.encode_bit(low_row, (slot << 2) | b, bit);
            } else {
                enc.encode_direct_bit(bit);
            }
        }
    }

    uint32_t read_count(BitRangeDec& dec, int* tree, int low_ctx) {
        int slot = read_fixed_tree(dec, tree, 6);
        if (slot < 4) return static_cast<uint32_t>(slot);
        int extra_bits = (slot - 2) >> 1;
        uint32_t second_bit = static_cast<uint32_t>((slot - 2) & 1);
        uint32_t rem = 0;
        int* low_row = low_bit_probs[low_ctx];
        for (int b = extra_bits - 1; b >= 0; b--) {
            int bit = (b < 4 && slot < 32) ? dec.decode_bit(low_row, (slot << 2) | b) : dec.decode_direct_bit();
            rem = (rem << 1) | static_cast<uint32_t>(bit);
        }
        return ((2U | second_bit) << extra_bits) | rem;
    }
};

struct MosaicConfigC {
    bool has_panel = false;
    bool allow_approx = true;
    bool allow_splice = true;
    bool allow_kimura = true;
    bool allow_wobble = true;
    bool allow_self = true;
    bool allow_tandem = true;

    uint8_t to_header_byte() const {
        uint8_t f = 0;
        if (has_panel) f |= 1;
        if (allow_approx) f |= 2;
        if (allow_splice) f |= 4;
        if (allow_kimura) f |= 8;
        if (allow_wobble) f |= 16;
        if (allow_self) f |= 32;
        if (allow_tandem) f |= 64;
        return f;
    }

    static MosaicConfigC from_header_byte(uint8_t f) {
        MosaicConfigC c;
        c.has_panel = (f & 1) != 0;
        c.allow_approx = (f & 2) != 0;
        c.allow_splice = (f & 4) != 0;
        c.allow_kimura = (f & 8) != 0;
        c.allow_wobble = (f & 16) != 0;
        c.allow_self = (f & 32) != 0;
        c.allow_tandem = (f & 64) != 0;
        return c;
    }
};

inline int gap_regime(int ema_gap) {
    if (ema_gap < 16) return 0;
    if (ema_gap < 96) return 1;
    if (ema_gap < 512) return 2;
    return 3;
}

inline void init_mtf(std::vector<int>& mtf, size_t max_donors) {
    size_t states = std::max<size_t>(4, max_donors << 1);
    mtf.resize(4);
    for (size_t i = 0; i < 4; i++) mtf[i] = static_cast<int>(i % states);
}

inline void promote_mtf(std::vector<int>& mtf, int val) {
    if (mtf[0] == val) return;
    int pos = 3;
    for (int i = 1; i < 3; i++) {
        if (mtf[i] == val) { pos = i; break; }
    }
    for (int i = pos; i > 0; i--) mtf[i] = mtf[i - 1];
    mtf[0] = val;
}

} // namespace

// ============================================================================
// 64-bit Syncmer CSR Index (HgsSyncmerIndex64)
// ============================================================================

struct HgsSyncmerIndex64 {
    Packed2Bit bank;
    size_t num_donors = 0;
    std::vector<uint64_t> donor_offsets;
    std::vector<size_t> pos_to_donor_floor;

    int primary_subsample_shift = 0;
    uint32_t primary_subsample_mask = 0;
    int bucket_bits = 12;
    uint32_t bucket_mask = 0;
    std::vector<uint64_t> bucket_offsets;
    std::vector<uint64_t> pos_and_strand; // bit 63 = strand, bits 0..62 = global_pos
    std::vector<uint16_t> sig16;

    int sec_bucket_bits = 0;
    uint32_t sec_bucket_mask = 0;
    std::vector<uint64_t> sec_bucket_offsets;
    std::vector<uint64_t> sec_pos_and_strand;
    std::vector<uint16_t> sec_sig16;

    size_t donor_of_global_pos(uint64_t global_pos) const {
        size_t blk = static_cast<size_t>(global_pos >> BLOCK_SHIFT);
        size_t d = (blk < pos_to_donor_floor.size()) ? pos_to_donor_floor[blk] : 0;
        while (d + 1 < num_donors && global_pos >= donor_offsets[d + 1]) d++;
        return d;
    }

    bool accepts_primary(uint64_t canon_k16) const {
        return primary_subsample_mask == 0 || (hgs_mix32(canon_k16) & primary_subsample_mask) == 0;
    }

    uint32_t bucket_for_primary(uint64_t canon_k16) const {
        return (hgs_mix32(canon_k16) >> primary_subsample_shift) & bucket_mask;
    }

    uint16_t sig_for_primary(uint64_t canon_k16) const {
        uint32_t h = hgs_mix32(canon_k16);
        return static_cast<uint16_t>((canon_k16 >> 16) ^ (h >> 16));
    }

    bool has_secondary() const { return !sec_bucket_offsets.empty(); }

    int64_t sec_bucket_for(uint64_t canon_k13) const {
        uint32_t h = hgs_mix32(canon_k13);
        if ((h & 1U) != 0U) return -1;
        return static_cast<int64_t>((h >> 1) & sec_bucket_mask);
    }

    uint16_t sec_sig_for(uint64_t canon_k13) const {
        uint32_t h = hgs_mix32(canon_k13);
        return static_cast<uint16_t>((canon_k13 >> 16) ^ (h >> 16));
    }
};

namespace {

int choose_bucket_bits(uint64_t est_entries, int min_bits, int max_bits) {
    int bits = min_bits;
    while (bits < max_bits && (1ULL << bits) < (est_entries >> 1)) bits++;
    return bits;
}

void build_csr_pass(
    const Packed2Bit& bank,
    const std::vector<uint64_t>& donor_offsets,
    int k,
    int s,
    int t,
    int bucket_bits,
    int subsample_shift,
    std::vector<uint64_t>& out_bucket_offsets,
    std::vector<uint64_t>& out_pos_strand,
    std::vector<uint16_t>& out_sig16) {
    size_t n = donor_offsets.size() - 1;
    size_t num_buckets = 1ULL << bucket_bits;
    uint32_t mask = static_cast<uint32_t>(num_buckets - 1);
    uint32_t sub_mask = (1U << subsample_shift) - 1U;

    std::vector<uint64_t> counts(num_buckets, 0);
    std::vector<int64_t> last_donor(num_buckets, -1);
    std::vector<uint8_t> occ(num_buckets, 0);

    struct Pass1Ctx {
        uint64_t donor_id;
        uint32_t sub_mask;
        int subsample_shift;
        uint32_t mask;
        std::vector<uint64_t>* counts;
        std::vector<int64_t>* last_donor;
        std::vector<uint8_t>* occ;
    };

    auto cb1 = [](uint64_t /*rel*/, uint64_t canon_k, uint8_t /*strand*/, void* vctx) {
        auto* ctx = static_cast<Pass1Ctx*>(vctx);
        uint32_t h = hgs_mix32(canon_k);
        if ((h & ctx->sub_mask) != 0U) return;
        uint32_t b = (h >> ctx->subsample_shift) & ctx->mask;
        if ((*ctx->last_donor)[b] != static_cast<int64_t>(ctx->donor_id)) {
            (*ctx->last_donor)[b] = static_cast<int64_t>(ctx->donor_id);
            (*ctx->occ)[b] = 1;
            (*ctx->counts)[b]++;
        } else if ((*ctx->occ)[b] < MAX_WITHIN_DONOR_BUCKET_OCC) {
            (*ctx->occ)[b]++;
            (*ctx->counts)[b]++;
        }
    };

    Pass1Ctx c1{0, sub_mask, subsample_shift, mask, &counts, &last_donor, &occ};
    for (size_t d = 0; d < n; d++) {
        c1.donor_id = d;
        uint64_t start = donor_offsets[d];
        uint64_t len = donor_offsets[d + 1] - start;
        hgs_scan_open_syncmers(bank.words.data(), start, len, k, s, t, cb1, &c1);
    }

    out_bucket_offsets.resize(num_buckets + 1);
    uint64_t running = 0;
    for (size_t b = 0; b < num_buckets; b++) {
        out_bucket_offsets[b] = running;
        running += counts[b];
    }
    out_bucket_offsets[num_buckets] = running;

    out_pos_strand.resize(static_cast<size_t>(running));
    out_sig16.resize(static_cast<size_t>(running));
    std::vector<uint64_t> write_cur(out_bucket_offsets.begin(), out_bucket_offsets.begin() + num_buckets);
    std::fill(last_donor.begin(), last_donor.end(), -1);

    struct Pass2Ctx {
        uint64_t donor_id;
        uint64_t d_start;
        uint32_t sub_mask;
        int subsample_shift;
        uint32_t mask;
        std::vector<int64_t>* last_donor;
        std::vector<uint8_t>* occ;
        std::vector<uint64_t>* write_cur;
        std::vector<uint64_t>* pos_strand;
        std::vector<uint16_t>* sig16;
    };

    auto cb2 = [](uint64_t rel_pos, uint64_t canon_k, uint8_t strand, void* vctx) {
        auto* ctx = static_cast<Pass2Ctx*>(vctx);
        uint32_t h = hgs_mix32(canon_k);
        if ((h & ctx->sub_mask) != 0U) return;
        uint32_t b = (h >> ctx->subsample_shift) & ctx->mask;
        bool accept = false;
        if ((*ctx->last_donor)[b] != static_cast<int64_t>(ctx->donor_id)) {
            (*ctx->last_donor)[b] = static_cast<int64_t>(ctx->donor_id);
            (*ctx->occ)[b] = 1;
            accept = true;
        } else if ((*ctx->occ)[b] < MAX_WITHIN_DONOR_BUCKET_OCC) {
            (*ctx->occ)[b]++;
            accept = true;
        }
        if (accept) {
            size_t dst = static_cast<size_t>((*ctx->write_cur)[b]++);
            (*ctx->pos_strand)[dst] = (ctx->d_start + rel_pos) | (static_cast<uint64_t>(strand) << 63);
            (*ctx->sig16)[dst] = static_cast<uint16_t>((canon_k >> 16) ^ (h >> 16));
        }
    };

    Pass2Ctx c2{0, 0, sub_mask, subsample_shift, mask, &last_donor, &occ, &write_cur, &out_pos_strand, &out_sig16};
    for (size_t d = 0; d < n; d++) {
        c2.donor_id = d;
        c2.d_start = donor_offsets[d];
        uint64_t len = donor_offsets[d + 1] - c2.d_start;
        hgs_scan_open_syncmers(bank.words.data(), c2.d_start, len, k, s, t, cb2, &c2);
    }
}

// Self-repeat 4-way set-associative table (64-bit positions, native cache sizing)
struct SelfSyncmerTable64 {
    static constexpr size_t WAYS = 4;
    uint32_t set_mask = 0;
    std::vector<uint64_t> pos_strand; // UINT64_MAX = empty, bit 63 = strand, bits 0..62 = pos
    std::vector<uint16_t> sigs;
    std::vector<uint8_t> clock_hand;

    explicit SelfSyncmerTable64(int set_bits) {
        size_t num_sets = 1ULL << set_bits;
        set_mask = static_cast<uint32_t>(num_sets - 1);
        pos_strand.assign(num_sets * WAYS, std::numeric_limits<uint64_t>::max());
        sigs.assign(num_sets * WAYS, 0);
        clock_hand.assign(num_sets, 0);
    }

    void insert(uint64_t canon_k, uint64_t pos, uint8_t strand) {
        uint32_t h = hgs_mix32(canon_k);
        uint32_t set = h & set_mask;
        uint16_t sig = static_cast<uint16_t>((canon_k >> 16) ^ (h >> 16));
        size_t base = static_cast<size_t>(set) << 2;
        for (size_t w = 0; w < WAYS; w++) {
            size_t slot = base + w;
            uint64_t prev = pos_strand[slot];
            if (prev != std::numeric_limits<uint64_t>::max() && sigs[slot] == sig) {
                uint64_t ppos = prev & ((1ULL << 63) - 1ULL);
                if (pos - ppos < 64) {
                    pos_strand[slot] = pos | (static_cast<uint64_t>(strand) << 63);
                    return;
                }
            }
        }
        uint8_t hand = clock_hand[set] & 3U;
        size_t slot = base + hand;
        pos_strand[slot] = pos | (static_cast<uint64_t>(strand) << 63);
        sigs[slot] = sig;
        clock_hand[set] = static_cast<uint8_t>((hand + 1U) & 3U);
    }
};

struct Subblock {
    uint64_t len = 0;
    std::vector<uint64_t> snp_offsets;
    uint64_t d_tgt_after = 0;
    uint64_t d_ref_after = 0;
};

struct ParsedToken {
    int kind = 0; // 0 = LIT, 1 = PANEL, 2 = SELF, 3 = TANDEM
    uint64_t lit_len = 0;
    uint64_t donor_id = 0; // or tandem period
    uint8_t strand = 0;
    uint64_t start_ref_pos = 0;
    uint64_t total_tgt_span = 0;
    uint64_t total_ref_span = 0;
    std::vector<Subblock> subblocks;
};

double approx_count_bits(uint64_t value) {
    if (value <= 1) return 2.5;
    int msb = 63;
    while (msb > 0 && ((value >> msb) & 1ULL) == 0ULL) msb--;
    return 3.5 + msb;
}

double estimate_edit_bits(uint64_t len, size_t snp_count, const MosaicConfigC& cfg) {
    if (snp_count == 0) return 1.0;
    double avg_gap = static_cast<double>(len) / static_cast<double>(snp_count + 1);
    double gap_bits = 1.8 + std::max(1.5, std::log2(avg_gap + 1.0));
    double sub_bits = cfg.allow_kimura ? 1.32 : 2.00;
    return 1.5 + static_cast<double>(snp_count) * (gap_bits + sub_bits + 0.8);
}

bool extend_subblock_mdl(
    const Packed2Bit& ref_bank,
    uint64_t ref_offset,
    uint64_t ref_len,
    uint8_t strand,
    uint64_t ref_start,
    const Packed2Bit& target,
    uint64_t tgt_start,
    uint64_t tgt_limit,
    uint64_t min_len,
    bool allow_approx,
    bool clonal_gear,
    Subblock& out_sb) {
    if (ref_start >= ref_len || tgt_start >= tgt_limit) return false;
    uint64_t max_len = std::min(ref_len - ref_start, tgt_limit - tgt_start);
    if (max_len < min_len) return false;
    if (!allow_approx) {
        uint64_t ex = hgs_count_exact_match_2bit(
            ref_bank.words.data(), ref_offset, ref_len, strand, ref_start,
            target.words.data(), 0, target.length, tgt_start, max_len);
        if (ex < min_len) return false;
        out_sb.len = ex;
        out_sb.snp_offsets.clear();
        out_sb.d_tgt_after = 0;
        out_sb.d_ref_after = 0;
        return true;
    }

    uint64_t k = 0;
    std::vector<uint64_t> snps;
    int consec_mm = 0;
    int window_mm = 0;
    const double base_gain = clonal_gear ? 0.008 : 1.865;
    const double snp_pen = clonal_gear ? 11.2 : 7.8;
    const double x_drop = clonal_gear ? -16.0 : -22.0;

    double cur_score = 0.0;
    double best_score = -1e9;
    uint64_t best_len = 0;
    size_t best_snp_cnt = 0;

    while (k < max_len) {
        uint64_t ex = hgs_count_exact_match_2bit(
            ref_bank.words.data(), ref_offset, ref_len, strand, ref_start + k,
            target.words.data(), 0, target.length, tgt_start + k, max_len - k);
        if (ex > 0) {
            k += ex;
            consec_mm = 0;
            window_mm = clonal_gear
                ? std::max(0, window_mm - static_cast<int>(std::min<uint64_t>(ex, 1000000)))
                : std::max(0, window_mm - static_cast<int>(std::min<uint64_t>(ex >> 2, 1000000)));
            cur_score += static_cast<double>(ex) * base_gain;
            if (k >= min_len && cur_score > best_score) {
                best_score = cur_score;
                best_len = k;
                best_snp_cnt = snps.size();
            }
            if (k >= max_len) break;
        }

        if (clonal_gear && k < 16) break;
        consec_mm++;
        if (clonal_gear) {
            if (consec_mm >= 2 || window_mm > 0) break;
            window_mm = 512;
        } else {
            window_mm += 4;
            if (consec_mm >= 4 || window_mm >= 28 || (k >= 32 && snps.size() * 7 > k)) break;
        }
        if (cur_score - best_score < x_drop && best_len >= min_len) break;
        snps.push_back(k);
        cur_score -= snp_pen;
        k++;
    }

    if (best_len < min_len) return false;
    uint64_t req_tail = clonal_gear ? 128 : 1;
    while (best_snp_cnt > 0 && (best_len - 1 - snps[best_snp_cnt - 1]) < (req_tail - 1)) {
        best_len = snps[best_snp_cnt - 1];
        best_snp_cnt--;
    }
    while (best_snp_cnt > 0 && snps[best_snp_cnt - 1] >= best_len - 1) {
        best_len = snps[best_snp_cnt - 1];
        best_snp_cnt--;
    }
    if (best_len < min_len) return false;
    snps.resize(best_snp_cnt);
    out_sb.len = best_len;
    out_sb.snp_offsets = std::move(snps);
    out_sb.d_tgt_after = 0;
    out_sb.d_ref_after = 0;
    return true;
}

bool try_find_syntenic_splice(
    const Packed2Bit& ref_bank,
    uint64_t ref_offset,
    uint64_t ref_len,
    uint8_t strand,
    uint64_t end_ref,
    const Packed2Bit& target,
    uint64_t end_tgt,
    uint64_t tgt_limit,
    uint64_t& out_dtgt,
    uint64_t& out_dref) {
    if (tgt_limit <= end_tgt || ref_len <= end_ref) return false;
    uint64_t rem_t = tgt_limit - end_tgt;
    uint64_t rem_r = ref_len - end_ref;
    if (rem_t < MIN_SPLICE_ANCHOR || rem_r < MIN_SPLICE_ANCHOR) return false;
    const uint32_t mask12 = 0xFFFFFFU;

    uint64_t tgt_w0 = target.word32_fwd(end_tgt);
    uint32_t tgt_a12 = static_cast<uint32_t>(tgt_w0) & mask12;
    uint64_t max_del = std::min<uint64_t>(192, rem_r - MIN_SPLICE_ANCHOR);
    for (uint64_t br = 1; br <= max_del; br += 21) {
        uint64_t span = std::min<uint64_t>(21, max_del - br + 1);
        uint64_t rc = ref_bank.word32_strand(ref_offset, ref_len, strand, end_ref + br);
        for (uint64_t j = 0; j < span; j++) {
            if ((static_cast<uint32_t>(rc) & mask12) == tgt_a12) {
                out_dtgt = 0;
                out_dref = br + j;
                return true;
            }
            rc >>= 2;
        }
    }

    uint64_t ref_w0 = ref_bank.word32_strand(ref_offset, ref_len, strand, end_ref);
    uint64_t max_ins = std::min<uint64_t>(12, rem_t - MIN_SPLICE_ANCHOR);
    if (max_ins >= 1) {
        uint32_t ref_a12 = static_cast<uint32_t>(ref_w0) & mask12;
        uint64_t tw = tgt_w0 >> 2;
        for (uint64_t dt = 1; dt <= max_ins; dt++) {
            if ((static_cast<uint32_t>(tw) & mask12) == ref_a12) {
                out_dtgt = dt;
                out_dref = 0;
                return true;
            }
            tw >>= 2;
        }
    }

    uint64_t max_ct = std::min<uint64_t>(6, rem_t - MIN_SPLICE_ANCHOR);
    uint64_t max_cr = std::min<uint64_t>(12, rem_r - MIN_SPLICE_ANCHOR);
    if (max_ct >= 1 && max_cr >= 1) {
        uint64_t tw = tgt_w0 >> 2;
        for (uint64_t dt = 1; dt <= max_ct; dt++) {
            uint32_t t12 = static_cast<uint32_t>(tw) & mask12;
            uint64_t rw = ref_w0 >> 2;
            for (uint64_t dr = 1; dr <= max_cr; dr++) {
                if ((static_cast<uint32_t>(rw) & mask12) == t12) {
                    out_dtgt = dt;
                    out_dref = dr;
                    return true;
                }
                rw >>= 2;
            }
            tw >>= 2;
        }
    }
    return false;
}

bool build_chain(
    int kind,
    uint64_t donor_id,
    uint8_t strand,
    uint64_t start_ref_pos,
    const Packed2Bit& ref_bank,
    uint64_t ref_offset,
    uint64_t ref_len,
    const Packed2Bit& target,
    uint64_t start_tgt_pos,
    uint64_t tgt_limit,
    const MosaicConfigC& cfg,
    bool clonal_gear,
    double header_bits,
    ParsedToken& out_tok,
    double& out_gain) {
    Subblock sb0;
    if (!extend_subblock_mdl(
            ref_bank, ref_offset, ref_len, strand, start_ref_pos,
            target, start_tgt_pos, tgt_limit, MIN_MATCH_LEN,
            cfg.allow_approx, clonal_gear, sb0)) {
        return false;
    }

    std::vector<Subblock> chain;
    chain.push_back(std::move(sb0));
    uint64_t cur_ref = start_ref_pos + chain[0].len;
    uint64_t cur_tgt = start_tgt_pos + chain[0].len;
    double total_edit = estimate_edit_bits(chain[0].len, chain[0].snp_offsets.size(), cfg);
    double total_splice_cost = 0.0;
    uint64_t total_nov_splice = 0;

    while (cfg.allow_splice && chain.size() < MAX_SPLICE_CHAIN
           && cur_tgt + MIN_SPLICE_ANCHOR <= tgt_limit && cur_ref + MIN_SPLICE_ANCHOR <= ref_len) {
        uint64_t dt = 0, dr = 0;
        if (!try_find_syntenic_splice(ref_bank, ref_offset, ref_len, strand, cur_ref, target, cur_tgt, tgt_limit, dt, dr)) {
            break;
        }
        Subblock next_sb;
        if (!extend_subblock_mdl(
                ref_bank, ref_offset, ref_len, strand, cur_ref + dr,
                target, cur_tgt + dt, tgt_limit, MIN_SPLICE_ANCHOR,
                cfg.allow_approx, clonal_gear, next_sb)) {
            break;
        }
        double scost = 6.0 + approx_count_bits(dt) + approx_count_bits(dr) + static_cast<double>(dt) * 1.9;
        double ecost = estimate_edit_bits(next_sb.len, next_sb.snp_offsets.size(), cfg);
        double step_gain = static_cast<double>(next_sb.len) * 1.865 - (scost + ecost + 3.0);
        if (step_gain <= 0.0) break;

        chain.back().d_tgt_after = dt;
        chain.back().d_ref_after = dr;
        total_edit += ecost;
        total_splice_cost += scost + 3.0;
        total_nov_splice += dt;
        cur_ref += dr + next_sb.len;
        cur_tgt += dt + next_sb.len;
        chain.push_back(std::move(next_sb));
    }

    uint64_t span_tgt = cur_tgt - start_tgt_pos;
    uint64_t span_ref = cur_ref - start_ref_pos;
    uint64_t matched = span_tgt - total_nov_splice;
    double total_cost = header_bits + approx_count_bits(chain[0].len) + total_edit + total_splice_cost + 2.5;
    double net = static_cast<double>(matched) * 1.865 - total_cost;
    if (net <= 6.0) return false;

    out_tok.kind = kind;
    out_tok.donor_id = donor_id;
    out_tok.strand = strand;
    out_tok.start_ref_pos = start_ref_pos;
    out_tok.total_tgt_span = span_tgt;
    out_tok.total_ref_span = span_ref;
    out_tok.subblocks = std::move(chain);
    out_gain = net;
    return true;
}

double estimate_panel_hdr_bits(
    size_t donor_id,
    uint8_t strand,
    uint64_t ref_pos,
    uint64_t tgt_pos,
    const std::vector<int>& mtf,
    int64_t last_ref,
    int64_t last_tgt,
    size_t max_donors) {
    int ds = static_cast<int>((donor_id << 1) | strand);
    double donor_bits;
    if (max_donors <= 1) donor_bits = 0.0;
    else if (ds == mtf[0]) donor_bits = 0.25;
    else if (ds == mtf[1] || ds == mtf[2] || ds == mtf[3]) donor_bits = 3.0;
    else donor_bits = 5.0 + std::log2(static_cast<double>(max_donors));

    double pos_bits = approx_count_bits(ref_pos) + 1.0;
    if (last_ref >= 0) {
        int64_t p0 = last_ref;
        int64_t p1 = last_ref + (static_cast<int64_t>(tgt_pos) - last_tgt);
        int64_t bd = std::min(std::llabs(static_cast<int64_t>(ref_pos) - p0),
                              std::llabs(static_cast<int64_t>(ref_pos) - p1));
        if (bd <= 65536) pos_bits = 2.0 + approx_count_bits(static_cast<uint64_t>(bd));
    }
    return 2.0 + donor_bits + pos_bits;
}

void encode_snp_gap(
    BitRangeEnc& enc,
    MosaicModels& m,
    bool wobble,
    uint32_t gap,
    uint32_t max_gap,
    int reg,
    int prev_mod3_zero) {
    if (!wobble) {
        m.write_count(enc, m.gap_quot_slots[reg], 5, gap);
        return;
    }
    uint32_t q = gap / 3U;
    uint32_t rem = gap - 3U * q;
    m.write_count(enc, m.gap_quot_slots[reg], 5, q);
    uint32_t max_rem = std::min<uint32_t>(2U, max_gap - 3U * q);
    if (max_rem == 0U) return;
    int qbin = (q == 0) ? 0 : ((q < 8) ? 1 : 2);
    int wctx = (reg << 3) | (qbin << 1) | prev_mod3_zero;
    if (max_rem == 1U) {
        enc.encode_bit(m.wobble_rem01.data(), wctx, static_cast<int>(rem));
        return;
    }
    if (rem == 2U) {
        enc.encode_bit(m.wobble_is_third.data(), wctx, 1);
    } else {
        enc.encode_bit(m.wobble_is_third.data(), wctx, 0);
        enc.encode_bit(m.wobble_rem01.data(), wctx, static_cast<int>(rem));
    }
}

uint32_t decode_snp_gap(
    BitRangeDec& dec,
    MosaicModels& m,
    bool wobble,
    uint32_t max_gap,
    int reg,
    int prev_mod3_zero,
    bool& ok) {
    if (!wobble) return m.read_count(dec, m.gap_quot_slots[reg], 5);
    uint32_t q = m.read_count(dec, m.gap_quot_slots[reg], 5);
    if (3ULL * q > max_gap) { ok = false; return 0; }
    uint32_t max_rem = std::min<uint32_t>(2U, max_gap - 3U * q);
    if (max_rem == 0U) return 3U * q;
    int qbin = (q == 0) ? 0 : ((q < 8) ? 1 : 2);
    int wctx = (reg << 3) | (qbin << 1) | prev_mod3_zero;
    if (max_rem == 1U) {
        return 3U * q + static_cast<uint32_t>(dec.decode_bit(m.wobble_rem01.data(), wctx));
    }
    if (dec.decode_bit(m.wobble_is_third.data(), wctx) == 1) return 3U * q + 2U;
    return 3U * q + static_cast<uint32_t>(dec.decode_bit(m.wobble_rem01.data(), wctx));
}

void encode_mtf_ds(BitRangeEnc& enc, MosaicModels& m, std::vector<int>& mtf, size_t max_donors, int ds, int after_lit) {
    size_t total_states = max_donors << 1;
    if (total_states <= 1) return;
    if (ds == mtf[0]) {
        enc.encode_bit(m.mtf_hit0.data(), after_lit, 1);
        return;
    }
    enc.encode_bit(m.mtf_hit0.data(), after_lit, 0);
    if (total_states == 2) { promote_mtf(mtf, ds); return; }
    int slot = -1;
    int max_slot = static_cast<int>(std::min<size_t>(4, total_states));
    for (int s = 1; s < max_slot; s++) {
        if (mtf[s] == ds) { slot = s; break; }
    }
    if (slot >= 1) {
        if (total_states > 4) enc.encode_bit(m.mtf_hit13.data(), 0, 1);
        if (max_slot == 3) enc.encode_bit(m.mtf_choice.data(), 0, slot - 1);
        else if (max_slot >= 4) {
            if (slot == 1) enc.encode_bit(m.mtf_choice.data(), 0, 0);
            else {
                enc.encode_bit(m.mtf_choice.data(), 0, 1);
                enc.encode_bit(m.mtf_choice.data(), 1, slot - 2);
            }
        }
        promote_mtf(mtf, ds);
        return;
    }
    if (total_states > 4) enc.encode_bit(m.mtf_hit13.data(), 0, 0);
    uint32_t did = static_cast<uint32_t>(ds >> 1);
    int st = ds & 1;
    m.write_bounded_int(enc, m.donor_tree.data(), did, static_cast<uint32_t>(max_donors));
    enc.encode_bit(m.panel_strand.data(), 0, st);
    promote_mtf(mtf, ds);
}

int decode_mtf_ds(BitRangeDec& dec, MosaicModels& m, std::vector<int>& mtf, size_t max_donors, int after_lit) {
    size_t total_states = max_donors << 1;
    if (total_states <= 1) return 0;
    if (dec.decode_bit(m.mtf_hit0.data(), after_lit) == 1) return mtf[0];
    if (total_states == 2) {
        int ds = mtf[0] ^ 1;
        promote_mtf(mtf, ds);
        return ds;
    }
    int max_slot = static_cast<int>(std::min<size_t>(4, total_states));
    bool in13 = (total_states <= 4) || (dec.decode_bit(m.mtf_hit13.data(), 0) == 1);
    if (in13) {
        int slot = (max_slot == 3) ? (1 + dec.decode_bit(m.mtf_choice.data(), 0))
                   : ((dec.decode_bit(m.mtf_choice.data(), 0) == 0) ? 1 : (2 + dec.decode_bit(m.mtf_choice.data(), 1)));
        int ds = mtf[slot];
        promote_mtf(mtf, ds);
        return ds;
    }
    uint32_t did = m.read_bounded_int(dec, m.donor_tree.data(), static_cast<uint32_t>(max_donors));
    int st = dec.decode_bit(m.panel_strand.data(), 0);
    int ds = static_cast<int>((did << 1) | static_cast<uint32_t>(st));
    promote_mtf(mtf, ds);
    return ds;
}

void encode_panel_ref_pos(
    BitRangeEnc& enc,
    MosaicModels& m,
    uint64_t ref_pos,
    uint64_t tgt_pos,
    uint64_t donor_len,
    int64_t last_ref,
    int64_t last_tgt,
    int is_mtf0) {
    if (last_ref >= 0 && donor_len > 0) {
        int64_t dl = static_cast<int64_t>(donor_len);
        int64_t pred0 = std::clamp<int64_t>(last_ref, 0, dl - 1);
        int64_t pred1 = std::clamp<int64_t>(last_ref + (static_cast<int64_t>(tgt_pos) - last_tgt), 0, dl - 1);
        int64_t d0 = static_cast<int64_t>(ref_pos) - pred0;
        int64_t d1 = static_cast<int64_t>(ref_pos) - pred1;
        int64_t best_d = (std::llabs(d1) <= std::llabs(d0)) ? d1 : d0;
        int use_col = (std::llabs(d1) <= std::llabs(d0)) ? 1 : 0;
        if (std::llabs(best_d) <= 65536) {
            enc.encode_bit(m.use_pos_delta.data(), is_mtf0, 1);
            if (pred0 != pred1) enc.encode_bit(m.use_collinear_pred.data(), is_mtf0, use_col);
            uint32_t zz = static_cast<uint32_t>((best_d << 1) ^ (best_d >> 63));
            m.write_count(enc, m.pos_delta_slots[use_col], 3, zz);
            return;
        }
        enc.encode_bit(m.use_pos_delta.data(), is_mtf0, 0);
    }
    m.write_count(enc, m.raw_pos_slots[0], 4, static_cast<uint32_t>(ref_pos));
}

uint64_t decode_panel_ref_pos(
    BitRangeDec& dec,
    MosaicModels& m,
    uint64_t tgt_pos,
    uint64_t donor_len,
    int64_t last_ref,
    int64_t last_tgt,
    int is_mtf0,
    bool& ok) {
    if (donor_len == 0) { ok = false; return 0; }
    if (last_ref >= 0 && dec.decode_bit(m.use_pos_delta.data(), is_mtf0) == 1) {
        int64_t dl = static_cast<int64_t>(donor_len);
        int64_t pred0 = std::clamp<int64_t>(last_ref, 0, dl - 1);
        int64_t pred1 = std::clamp<int64_t>(last_ref + (static_cast<int64_t>(tgt_pos) - last_tgt), 0, dl - 1);
        int64_t pred = pred0;
        int use_col = 1;
        if (pred0 != pred1) {
            use_col = dec.decode_bit(m.use_collinear_pred.data(), is_mtf0);
            pred = (use_col == 1) ? pred1 : pred0;
        }
        uint32_t zz = m.read_count(dec, m.pos_delta_slots[use_col], 3);
        int32_t sdelta = static_cast<int32_t>((zz >> 1) ^ -static_cast<int32_t>(zz & 1U));
        int64_t rp = pred + static_cast<int64_t>(sdelta);
        if (rp < 0 || static_cast<uint64_t>(rp) >= donor_len) { ok = false; return 0; }
        return static_cast<uint64_t>(rp);
    }
    uint64_t rp = m.read_count(dec, m.raw_pos_slots[0], 4);
    if (rp >= donor_len) { ok = false; return 0; }
    return rp;
}

void encode_subblocks(
    BitRangeEnc& enc,
    MosaicModels& m,
    const MosaicConfigC& cfg,
    const Packed2Bit& ref_bank,
    uint64_t ref_offset,
    uint64_t ref_len,
    uint8_t strand,
    uint64_t start_ref_pos,
    const Packed2Bit& target,
    uint64_t start_tgt_pos,
    const std::vector<Subblock>& subblocks) {
    uint64_t cur_ref = start_ref_pos;
    uint64_t cur_tgt = start_tgt_pos;
    int ema_gap = 128;

    for (size_t b = 0; b < subblocks.size(); b++) {
        const Subblock& sb = subblocks[b];
        uint64_t min_l = (b == 0) ? MIN_MATCH_LEN : MIN_SPLICE_ANCHOR;
        m.write_count(enc, m.match_len_slots[b == 0 ? 0 : 1], 1, static_cast<uint32_t>(sb.len - min_l));

        if (cfg.allow_approx) {
            int len_bin = sb.len < 64 ? 0 : (sb.len < 512 ? 1 : (sb.len < 4096 ? 2 : 3));
            int has_snps = sb.snp_offsets.empty() ? 0 : 1;
            enc.encode_bit(m.has_snps.data(), ((b > 0 ? 1 : 0) << 2) | len_bin, has_snps);
            if (has_snps == 1) {
                uint64_t local_off = 0;
                int prev_mod3 = 1;
                for (size_t s = 0; s < sb.snp_offsets.size(); s++) {
                    uint64_t snp_off = sb.snp_offsets[s];
                    uint32_t gap = static_cast<uint32_t>(snp_off - local_off);
                    uint32_t max_gap = static_cast<uint32_t>((sb.len - 1) - local_off);
                    int reg = gap_regime(ema_gap);
                    if (s > 0) {
                        int rem_bin = max_gap < 16 ? 0 : (max_gap < 128 ? 1 : 2);
                        enc.encode_bit(m.more_snps.data(), (reg << 2) | rem_bin, 1);
                    }
                    encode_snp_gap(enc, m, cfg.allow_wobble, gap, max_gap, reg, prev_mod3);
                    prev_mod3 = ((gap + 1U) % 3U == 0U) ? 1 : 0;
                    ema_gap = (ema_gap * 3 + std::min<int>(2048, static_cast<int>(gap))) >> 2;

                    uint64_t mr = cur_ref + snp_off;
                    uint64_t mt = cur_tgt + snp_off;
                    uint8_t rb = ref_bank.base_at_strand(ref_offset, ref_len, strand, mr);
                    uint8_t tb = target.base_at(mt);
                    if (cfg.allow_kimura) {
                        uint8_t pr = (mr > 0) ? ref_bank.base_at_strand(ref_offset, ref_len, strand, mr - 1) : 0;
                        uint8_t nr = (mr + 1 < ref_len) ? ref_bank.base_at_strand(ref_offset, ref_len, strand, mr + 1) : 0;
                        int ctx = (rb << 6) | (pr << 4) | (nr << 2) | ((reg <= 1 ? 1 : 0) << 1) | (cfg.allow_wobble ? prev_mod3 : 0);
                        uint8_t delta = static_cast<uint8_t>(tb ^ rb);
                        if (delta == 2) {
                            enc.encode_bit(m.kimura_ti.data(), ctx, 1);
                        } else {
                            enc.encode_bit(m.kimura_ti.data(), ctx, 0);
                            enc.encode_bit(m.kimura_tv.data(), ctx, delta == 3 ? 1 : 0);
                        }
                    } else {
                        uint8_t pt = (mt > 0) ? target.base_at(mt - 1) : 0;
                        m.write_fixed_tree(enc, m.ablate_base4[pt], tb, 2);
                    }
                    local_off = snp_off + 1;
                }
                if (local_off + 1 < sb.len) {
                    uint32_t max_gap = static_cast<uint32_t>((sb.len - 1) - local_off);
                    int rem_bin = max_gap < 16 ? 0 : (max_gap < 128 ? 1 : 2);
                    int reg = gap_regime(ema_gap);
                    enc.encode_bit(m.more_snps.data(), (reg << 2) | rem_bin, 0);
                }
            }
        }

        int has_snps_flag = (cfg.allow_approx && !sb.snp_offsets.empty()) ? 1 : 0;
        cur_ref += sb.len;
        cur_tgt += sb.len;

        bool can_splice = cfg.allow_splice && (cur_tgt + MIN_SPLICE_ANCHOR <= target.length) && (cur_ref + MIN_SPLICE_ANCHOR <= ref_len);
        if (can_splice) {
            bool has_sp = (b + 1 < subblocks.size());
            int sp_ctx = ((b > 0 ? 1 : 0) << 2) | (has_snps_flag << 1) | (sb.len >= 256 ? 1 : 0);
            enc.encode_bit(m.has_splice.data(), sp_ctx, has_sp ? 1 : 0);
            if (has_sp) {
                if (sb.d_tgt_after == 0) {
                    enc.encode_bit(m.splice_is_pure_del.data(), 0, 1);
                    m.write_count(enc, m.splice_del_slots[0], 7, static_cast<uint32_t>(sb.d_ref_after - 1));
                } else {
                    enc.encode_bit(m.splice_is_pure_del.data(), 0, 0);
                    if (sb.d_ref_after == 0) {
                        enc.encode_bit(m.splice_is_pure_ins.data(), 0, 1);
                        m.write_count(enc, m.splice_ins_slots[0], 6, static_cast<uint32_t>(sb.d_tgt_after - 1));
                    } else {
                        enc.encode_bit(m.splice_is_pure_ins.data(), 0, 0);
                        m.write_count(enc, m.splice_ins_slots[0], 6, static_cast<uint32_t>(sb.d_tgt_after - 1));
                        m.write_count(enc, m.splice_del_slots[0], 7, static_cast<uint32_t>(sb.d_ref_after - 1));
                    }
                }
                cur_tgt += sb.d_tgt_after;
                cur_ref += sb.d_ref_after;
            }
        }
    }
}

void copy_donor_span_to_out(
    const Packed2Bit& bank, uint64_t offset, uint64_t ref_len, uint8_t strand, uint64_t start_ref_pos,
    uint8_t* out, uint64_t dst_pos, uint64_t len) {
    if (strand == 0) {
        hgs_unpack_2bit_bytes(bank.words.data(), offset + start_ref_pos, static_cast<size_t>(len), out + dst_pos);
    } else {
        uint64_t fwd_top = offset + ref_len - 1ULL - start_ref_pos;
        uint64_t k = 0;
        while (k + 32 <= len) {
            uint64_t w = hgs_revcomp_word32(hgs_word32_forward(bank.words.data(), fwd_top - k - 31ULL));
            for (int j = 0; j < 32; j++) {
                out[dst_pos + k + static_cast<uint64_t>(j)] = static_cast<uint8_t>(w & 3ULL);
                w >>= 2;
            }
            k += 32;
        }
        while (k < len) {
            out[dst_pos + k] = static_cast<uint8_t>(bank.base_at(fwd_top - k) ^ 3U);
            k++;
        }
    }
}

bool decode_subblocks_donor(
    BitRangeDec& dec, MosaicModels& m, const MosaicConfigC& cfg,
    const Packed2Bit* bank_opt, uint64_t ref_offset, uint64_t ref_len,
    bool is_self, uint8_t strand, uint64_t start_ref_pos,
    uint8_t* out, uint64_t start_tgt_pos, uint64_t target_len,
    const uint8_t* novelty, uint64_t& nov_cur, uint64_t novelty_len,
    uint64_t& out_end_ref, uint64_t& out_end_tgt) {
    uint64_t cur_ref = start_ref_pos;
    uint64_t cur_tgt = start_tgt_pos;
    int ema_gap = 128;
    size_t b = 0;
    bool ok = true;

    auto read_ref_base = [&](uint64_t mr) -> uint8_t {
        if (!is_self) return bank_opt->base_at_strand(ref_offset, ref_len, strand, mr);
        return (strand == 0) ? (out[mr] & 3U) : static_cast<uint8_t>((out[ref_len - 1ULL - mr] ^ 3U) & 3U);
    };

    while (true) {
        uint64_t min_l = (b == 0) ? MIN_MATCH_LEN : MIN_SPLICE_ANCHOR;
        uint64_t sub_len = m.read_count(dec, m.match_len_slots[b == 0 ? 0 : 1], 1) + min_l;
        if (cur_ref + sub_len > ref_len || cur_tgt + sub_len > target_len) return false;

        if (!is_self) {
            copy_donor_span_to_out(*bank_opt, ref_offset, ref_len, strand, cur_ref, out, cur_tgt, sub_len);
        } else if (strand == 0) {
            for (uint64_t k = 0; k < sub_len; k++) out[cur_tgt + k] = out[cur_ref + k];
        } else {
            uint64_t fwd_top = ref_len - 1ULL - cur_ref;
            for (uint64_t k = 0; k < sub_len; k++) out[cur_tgt + k] = static_cast<uint8_t>((out[fwd_top - k] ^ 3U) & 3U);
        }

        int has_snps_flag = 0;
        if (cfg.allow_approx) {
            int len_bin = sub_len < 64 ? 0 : (sub_len < 512 ? 1 : (sub_len < 4096 ? 2 : 3));
            int has_snps = dec.decode_bit(m.has_snps.data(), ((b > 0 ? 1 : 0) << 2) | len_bin);
            has_snps_flag = has_snps;
            if (has_snps == 1) {
                uint64_t local_off = 0;
                int prev_mod3 = 1;
                size_t s = 0;
                while (true) {
                    if (local_off + 1 >= sub_len) break;
                    uint32_t max_gap = static_cast<uint32_t>((sub_len - 1) - local_off);
                    int reg = gap_regime(ema_gap);
                    if (s > 0) {
                        int rem_bin = max_gap < 16 ? 0 : (max_gap < 128 ? 1 : 2);
                        if (dec.decode_bit(m.more_snps.data(), (reg << 2) | rem_bin) == 0) break;
                    }
                    uint32_t gap = decode_snp_gap(dec, m, cfg.allow_wobble, max_gap, reg, prev_mod3, ok);
                    if (!ok) return false;
                    uint64_t snp_off = local_off + gap;
                    if (snp_off >= sub_len) return false;
                    prev_mod3 = ((gap + 1U) % 3U == 0U) ? 1 : 0;
                    ema_gap = (ema_gap * 3 + std::min<int>(2048, static_cast<int>(gap))) >> 2;

                    uint64_t mr = cur_ref + snp_off;
                    uint64_t mt = cur_tgt + snp_off;
                    uint8_t rb = out[mt] & 3U;
                    uint8_t tb;
                    if (cfg.allow_kimura) {
                        uint8_t pr = (mr > 0) ? read_ref_base(mr - 1) : 0;
                        uint8_t nr = (mr + 1 < ref_len) ? read_ref_base(mr + 1) : 0;
                        int ctx = (rb << 6) | (pr << 4) | (nr << 2) | ((reg <= 1 ? 1 : 0) << 1) | (cfg.allow_wobble ? prev_mod3 : 0);
                        uint8_t delta = (dec.decode_bit(m.kimura_ti.data(), ctx) == 1)
                            ? 2
                            : (dec.decode_bit(m.kimura_tv.data(), ctx) == 1 ? 3 : 1);
                        tb = static_cast<uint8_t>(rb ^ delta);
                    } else {
                        uint8_t pt = (mt > 0) ? (out[mt - 1] & 3U) : 0;
                        tb = static_cast<uint8_t>(m.read_fixed_tree(dec, m.ablate_base4[pt], 2));
                    }
                    out[mt] = tb;
                    local_off = snp_off + 1;
                    s++;
                    if (local_off + 1 >= sub_len) break;
                }
            }
        }

        cur_ref += sub_len;
        cur_tgt += sub_len;
        bool can_splice = cfg.allow_splice && (cur_tgt + MIN_SPLICE_ANCHOR <= target_len) && (cur_ref + MIN_SPLICE_ANCHOR <= ref_len);
        int sp_ctx = ((b > 0 ? 1 : 0) << 2) | (has_snps_flag << 1) | (sub_len >= 256 ? 1 : 0);
        if (!can_splice || dec.decode_bit(m.has_splice.data(), sp_ctx) == 0) break;

        uint64_t dt = 0, dr = 0;
        if (dec.decode_bit(m.splice_is_pure_del.data(), 0) == 1) {
            dr = m.read_count(dec, m.splice_del_slots[0], 7) + 1ULL;
        } else if (dec.decode_bit(m.splice_is_pure_ins.data(), 0) == 1) {
            dt = m.read_count(dec, m.splice_ins_slots[0], 6) + 1ULL;
        } else {
            dt = m.read_count(dec, m.splice_ins_slots[0], 6) + 1ULL;
            dr = m.read_count(dec, m.splice_del_slots[0], 7) + 1ULL;
        }
        if (cur_tgt + dt > target_len || nov_cur + dt > novelty_len) return false;
        if (dt > 0) {
            std::memcpy(out + cur_tgt, novelty + nov_cur, static_cast<size_t>(dt));
            nov_cur += dt;
            cur_tgt += dt;
        }
        cur_ref += dr;
        b++;
    }

    out_end_ref = cur_ref;
    out_end_tgt = cur_tgt;
    return true;
}

// ============================================================================
// Stage-2 Novelty Coder (Exact wire format for FastCRM + Turbo + Raw 2-Bit)
// ============================================================================

constexpr int NOV_PROB_BITS = 15;
constexpr int NOV_TOT = 1 << NOV_PROB_BITS; // 32768
constexpr uint32_t NOV_TOP_VALUE = 1U << 24;
constexpr int P_BITS = 14;
constexpr int P_SCALE = 1 << P_BITS; // 16384
constexpr int EXP_MAX_X = 65536;
constexpr int PAN_SUPER_THRESHOLD = 6'000'000;
constexpr uint64_t TURBO_SEG_BASES = 4'000'000ULL;
constexpr uint8_t FLAG_TURBO_NOVELTY = 4;
constexpr uint8_t SEG_RAW_2BIT = 0;
constexpr uint8_t SEG_CRM = 1;

struct NoveltyLuts {
    std::vector<int32_t> ln_lut;
    std::vector<int32_t> exp_neg_lut;
    std::vector<uint32_t> recip;

    static const NoveltyLuts& get() {
        static const NoveltyLuts inst;
        return inst;
    }

    int32_t div_round_small(int64_t num, int denom) const {
        int64_t n = num + (denom >> 1);
        if (denom < static_cast<int>(recip.size())) {
            int64_t q = (n * static_cast<int64_t>(recip[denom])) >> 30;
            while (n - q * denom >= denom) q++;
            return static_cast<int32_t>(q);
        }
        return static_cast<int32_t>(n / denom);
    }

private:
    NoveltyLuts() : ln_lut(P_SCALE + 1), exp_neg_lut(EXP_MAX_X + 1), recip(65600) {
        const int64_t LN2_Q16 = 45426;
        for (int p = 1; p <= P_SCALE; p++) {
            int int_log2 = 0;
            int64_t y = p;
            while (y < P_SCALE) { y <<= 1; int_log2--; }
            int frac_bits = 0;
            for (int b = 15; b >= 0; b--) {
                y = (y * y) >> P_BITS;
                if (y >= (P_SCALE << 1)) {
                    frac_bits |= (1 << b);
                    y >>= 1;
                }
            }
            int64_t log2_q16 = (static_cast<int64_t>(int_log2) << 16) + frac_bits;
            ln_lut[p] = static_cast<int32_t>((log2_q16 * LN2_Q16) >> 20);
        }
        ln_lut[0] = ln_lut[1] - 2839;

        const int64_t INV_LN2_Q16 = 94548;
        for (int x = 0; x <= EXP_MAX_X; x++) {
            int64_t uq16 = (static_cast<int64_t>(x) * INV_LN2_Q16) >> 12;
            int k = static_cast<int>(uq16 >> 16);
            if (k >= 16) { exp_neg_lut[x] = 1; continue; }
            int64_t r = uq16 & 0xFFFF;
            int64_t s = (r * LN2_Q16) >> 16;
            int64_t s2 = (s * s) >> 16;
            int64_t s3 = (s2 * s) >> 16;
            int64_t s4 = (s3 * s) >> 16;
            int64_t s5 = (s4 * s) >> 16;
            int64_t exp_neg_s = 65536 - s + (s2 >> 1) - (s3 / 6) + (s4 / 24) - (s5 / 120);
            exp_neg_lut[x] = static_cast<int32_t>(std::max<int64_t>(1, exp_neg_s >> k));
        }

        for (size_t d = 1; d < recip.size(); d++) {
            recip[d] = static_cast<uint32_t>((1ULL << 30) / d);
        }
    }
};

struct NovRangeEnc {
    ByteVecOut& out;
    uint64_t low = 0;
    uint32_t range = 0xFFFFFFFFU;
    int32_t cache = -1;
    uint64_t cache_size = 0;

    explicit NovRangeEnc(ByteVecOut& o) : out(o) {}

    void shift_low() {
        uint32_t carry = static_cast<uint32_t>(low >> 32);
        if (carry != 0 || static_cast<uint32_t>(low) < 0xFF000000U) {
            if (cache >= 0) out.write_byte(static_cast<uint8_t>(cache + carry));
            uint8_t fill = static_cast<uint8_t>(0xFFU + carry);
            while (cache_size > 0) { out.write_byte(fill); cache_size--; }
            cache = static_cast<int32_t>((low >> 24) & 0xFFU);
        } else {
            cache_size++;
        }
        low = (low << 8) & 0xFFFFFFFFULL;
    }

    void encode_symbol(uint32_t cum_low, uint32_t freq) {
        uint64_t r = range >> NOV_PROB_BITS;
        low += cum_low * r;
        range = static_cast<uint32_t>(freq * r);
        while (range < NOV_TOP_VALUE) {
            shift_low();
            range <<= 8;
        }
    }

    void flush() {
        for (int i = 0; i < 5; i++) shift_low();
    }
};

struct NovRangeDec {
    ByteSpanIn& in;
    uint32_t range = 0xFFFFFFFFU;
    uint32_t code = 0;

    explicit NovRangeDec(ByteSpanIn& i) : in(i) {
        for (int k = 0; k < 4; k++) code = (code << 8) | next_byte();
    }

    uint32_t next_byte() {
        bool ok = true;
        return in.has_remaining() ? in.read_byte(ok) : 0U;
    }

    int get_threshold() const {
        uint32_t r = range >> NOV_PROB_BITS;
        int t = static_cast<int>(code / r);
        return t >= NOV_TOT ? (NOV_TOT - 1) : t;
    }

    void decode_symbol(uint32_t cum_low, uint32_t freq) {
        uint32_t r = range >> NOV_PROB_BITS;
        code -= cum_low * r;
        range = freq * r;
        while (range < NOV_TOP_VALUE) {
            code = (code << 8) | next_byte();
            range <<= 8;
        }
    }
};

void pack_raw_2bit_seg(const uint8_t* bases, size_t len, std::vector<uint8_t>& out) {
    size_t raw_bytes = (len + 3) >> 2;
    out.resize(1 + raw_bytes);
    out[0] = SEG_RAW_2BIT;
    size_t full = len >> 2;
    size_t idx = 0;
    for (size_t b = 0; b < full; b++) {
        out[1 + b] = static_cast<uint8_t>((bases[idx] & 3U) | ((bases[idx + 1] & 3U) << 2) |
                                          ((bases[idx + 2] & 3U) << 4) | ((bases[idx + 3] & 3U) << 6));
        idx += 4;
    }
    if (idx < len) {
        uint8_t v = 0;
        int sh = 0;
        while (idx < len) {
            v |= static_cast<uint8_t>((bases[idx++] & 3U) << sh);
            sh += 2;
        }
        out[1 + full] = v;
    }
}

void encode_turbo_segment(const uint8_t* bases, size_t len, std::vector<uint8_t>& out_bytes) {
    const auto& luts = NoveltyLuts::get();
    ByteVecOut crm;
    size_t raw_bytes = (len + 3) >> 2;
    crm.reserve(raw_bytes + 64);
    crm.write_byte(SEG_CRM);
    NovRangeEnc rc(crm);

    std::vector<uint32_t> tab4(256, 0), tab6(4096, 0), tab8(65536, 0);
    std::vector<int32_t> rep_tab(32768, 0);
    uint32_t ctx4 = 0, ctx6 = 0, ctx8 = 0, ctx12 = 0;
    int64_t rep_pos = -1;
    int rep_streak = 0, rep_misses = 0;

    for (size_t r = 0; r < len; r++) {
        uint32_t c4 = tab4[ctx4];
        uint32_t c6 = tab6[ctx6];
        uint32_t c8 = tab8[ctx8];
        int w0 = 2 + static_cast<int>(c4 & 0x7FU) + (static_cast<int>(c6 & 0x7FU) << 1) + (static_cast<int>(c8 & 0x7FU) << 2);
        int w1 = 2 + static_cast<int>((c4 >> 8) & 0x7FU) + (static_cast<int>((c6 >> 8) & 0x7FU) << 1) + (static_cast<int>((c8 >> 8) & 0x7FU) << 2);
        int w2 = 2 + static_cast<int>((c4 >> 16) & 0x7FU) + (static_cast<int>((c6 >> 16) & 0x7FU) << 1) + (static_cast<int>((c8 >> 16) & 0x7FU) << 2);
        int w3 = 2 + static_cast<int>((c4 >> 24) & 0x7FU) + (static_cast<int>((c6 >> 24) & 0x7FU) << 1) + (static_cast<int>((c8 >> 24) & 0x7FU) << 2);
        int rep_pred = -1;
        if (rep_pos >= 0) {
            rep_pred = bases[rep_pos] & 3U;
            int boost = std::min(2048, rep_streak << 6);
            if (rep_pred == 0) w0 += boost;
            else if (rep_pred == 1) w1 += boost;
            else if (rep_pred == 2) w2 += boost;
            else w3 += boost;
        }
        int tot = w0 + w1 + w2 + w3;
        uint64_t inv = luts.recip[tot];
        int f0 = static_cast<int>((static_cast<uint64_t>(w0) * inv) >> 15);
        int f1 = static_cast<int>((static_cast<uint64_t>(w1) * inv) >> 15);
        int f2 = static_cast<int>((static_cast<uint64_t>(w2) * inv) >> 15);
        int f3 = NOV_TOT - f0 - f1 - f2;

        int sym = bases[r] & 3U;
        uint32_t cum = 0, freq = static_cast<uint32_t>(f0);
        if (sym == 1) { cum = static_cast<uint32_t>(f0); freq = static_cast<uint32_t>(f1); }
        else if (sym == 2) { cum = static_cast<uint32_t>(f0 + f1); freq = static_cast<uint32_t>(f2); }
        else if (sym == 3) { cum = static_cast<uint32_t>(f0 + f1 + f2); freq = static_cast<uint32_t>(f3); }
        rc.encode_symbol(cum, freq);

        uint32_t inc = 1U << (sym << 3);
        uint32_t nv4 = c4 + inc; tab4[ctx4] = (nv4 & 0x80808080U) ? ((nv4 & 0xFEFEFEFEU) >> 1) : nv4;
        uint32_t nv6 = c6 + inc; tab6[ctx6] = (nv6 & 0x80808080U) ? ((nv6 & 0xFEFEFEFEU) >> 1) : nv6;
        uint32_t nv8 = c8 + inc; tab8[ctx8] = (nv8 & 0x80808080U) ? ((nv8 & 0xFEFEFEFEU) >> 1) : nv8;

        if (rep_pos >= 0) {
            if (sym == rep_pred) { rep_streak++; rep_misses = 0; rep_pos++; }
            else {
                rep_misses++;
                if (rep_misses > 2) { rep_pos = -1; rep_streak = 0; }
                else rep_pos++;
            }
        }
        if (r >= 12) {
            uint32_t h12 = ((ctx12 * 0x9E3779B1U) >> 17) & 0x7FFFU;
            if (rep_pos < 0) {
                int32_t cand = rep_tab[h12] - 1;
                if (cand >= 12 && (static_cast<int64_t>(r) - cand) > 16
                    && (bases[cand] & 3U) == sym
                    && bases[cand - 1] == bases[r - 1]
                    && bases[cand - 2] == bases[r - 2]
                    && bases[cand - 3] == bases[r - 3]) {
                    rep_pos = cand + 1;
                    rep_streak = 4;
                    rep_misses = 0;
                }
            }
            rep_tab[h12] = static_cast<int32_t>(r + 1);
        }
        ctx4 = ((ctx4 << 2) | sym) & 0xFFU;
        ctx6 = ((ctx6 << 2) | sym) & 0xFFFU;
        ctx8 = ((ctx8 << 2) | sym) & 0xFFFFU;
        ctx12 = ((ctx12 << 2) | sym) & 0xFFFFFFU;
    }
    rc.flush();
    if (crm.buf.size() <= 1 + raw_bytes) {
        out_bytes = std::move(crm.buf);
    } else {
        pack_raw_2bit_seg(bases, len, out_bytes);
    }
}

bool decode_turbo_segment(const uint8_t* payload, size_t payload_len, uint8_t* dst, size_t len) {
    if (payload_len == 0) return false;
    if (payload[0] == SEG_RAW_2BIT) {
        size_t expected = 1 + ((len + 3) >> 2);
        if (payload_len != expected) return false;
        size_t idx = 0;
        for (size_t b = 1; b < payload_len && idx < len; b++) {
            uint8_t v = payload[b];
            for (int j = 0; j < 4 && idx < len; j++) {
                dst[idx++] = static_cast<uint8_t>(v & 3U);
                v >>= 2;
            }
        }
        return true;
    }
    if (payload[0] != SEG_CRM) return false;

    const auto& luts = NoveltyLuts::get();
    ByteSpanIn in{payload + 1, payload_len - 1, 0};
    NovRangeDec rd(in);

    std::vector<uint32_t> tab4(256, 0), tab6(4096, 0), tab8(65536, 0);
    std::vector<int32_t> rep_tab(32768, 0);
    uint32_t ctx4 = 0, ctx6 = 0, ctx8 = 0, ctx12 = 0;
    int64_t rep_pos = -1;
    int rep_streak = 0, rep_misses = 0;

    for (size_t r = 0; r < len; r++) {
        uint32_t c4 = tab4[ctx4];
        uint32_t c6 = tab6[ctx6];
        uint32_t c8 = tab8[ctx8];
        int w0 = 2 + static_cast<int>(c4 & 0x7FU) + (static_cast<int>(c6 & 0x7FU) << 1) + (static_cast<int>(c8 & 0x7FU) << 2);
        int w1 = 2 + static_cast<int>((c4 >> 8) & 0x7FU) + (static_cast<int>((c6 >> 8) & 0x7FU) << 1) + (static_cast<int>((c8 >> 8) & 0x7FU) << 2);
        int w2 = 2 + static_cast<int>((c4 >> 16) & 0x7FU) + (static_cast<int>((c6 >> 16) & 0x7FU) << 1) + (static_cast<int>((c8 >> 16) & 0x7FU) << 2);
        int w3 = 2 + static_cast<int>((c4 >> 24) & 0x7FU) + (static_cast<int>((c6 >> 24) & 0x7FU) << 1) + (static_cast<int>((c8 >> 24) & 0x7FU) << 2);
        int rep_pred = -1;
        if (rep_pos >= 0) {
            rep_pred = dst[rep_pos] & 3U;
            int boost = std::min(2048, rep_streak << 6);
            if (rep_pred == 0) w0 += boost;
            else if (rep_pred == 1) w1 += boost;
            else if (rep_pred == 2) w2 += boost;
            else w3 += boost;
        }
        int tot = w0 + w1 + w2 + w3;
        uint64_t inv = luts.recip[tot];
        int f0 = static_cast<int>((static_cast<uint64_t>(w0) * inv) >> 15);
        int f1 = static_cast<int>((static_cast<uint64_t>(w1) * inv) >> 15);
        int f2 = static_cast<int>((static_cast<uint64_t>(w2) * inv) >> 15);
        int f3 = NOV_TOT - f0 - f1 - f2;

        int target = rd.get_threshold();
        int sym = 0;
        uint32_t cum = 0, freq = static_cast<uint32_t>(f0);
        if (target >= f0 + f1 + f2) { sym = 3; cum = static_cast<uint32_t>(f0 + f1 + f2); freq = static_cast<uint32_t>(f3); }
        else if (target >= f0 + f1) { sym = 2; cum = static_cast<uint32_t>(f0 + f1); freq = static_cast<uint32_t>(f2); }
        else if (target >= f0) { sym = 1; cum = static_cast<uint32_t>(f0); freq = static_cast<uint32_t>(f1); }
        rd.decode_symbol(cum, freq);
        dst[r] = static_cast<uint8_t>(sym);

        uint32_t inc = 1U << (sym << 3);
        uint32_t nv4 = c4 + inc; tab4[ctx4] = (nv4 & 0x80808080U) ? ((nv4 & 0xFEFEFEFEU) >> 1) : nv4;
        uint32_t nv6 = c6 + inc; tab6[ctx6] = (nv6 & 0x80808080U) ? ((nv6 & 0xFEFEFEFEU) >> 1) : nv6;
        uint32_t nv8 = c8 + inc; tab8[ctx8] = (nv8 & 0x80808080U) ? ((nv8 & 0xFEFEFEFEU) >> 1) : nv8;

        if (rep_pos >= 0) {
            if (sym == rep_pred) { rep_streak++; rep_misses = 0; rep_pos++; }
            else {
                rep_misses++;
                if (rep_misses > 2) { rep_pos = -1; rep_streak = 0; }
                else rep_pos++;
            }
        }
        if (r >= 12) {
            uint32_t h12 = ((ctx12 * 0x9E3779B1U) >> 17) & 0x7FFFU;
            if (rep_pos < 0) {
                int32_t cand = rep_tab[h12] - 1;
                if (cand >= 12 && (static_cast<int64_t>(r) - cand) > 16
                    && (dst[cand] & 3U) == sym
                    && dst[cand - 1] == dst[r - 1]
                    && dst[cand - 2] == dst[r - 2]
                    && dst[cand - 3] == dst[r - 3]) {
                    rep_pos = cand + 1;
                    rep_streak = 4;
                    rep_misses = 0;
                }
            }
            rep_tab[h12] = static_cast<int32_t>(r + 1);
        }
        ctx4 = ((ctx4 << 2) | sym) & 0xFFU;
        ctx6 = ((ctx6 << 2) | sym) & 0xFFFU;
        ctx8 = ((ctx8 << 2) | sym) & 0xFFFFU;
        ctx12 = ((ctx12 << 2) | sym) & 0xFFFFFFU;
    }
    return true;
}

// ============================================================================
// Full FastCrmModel (19/22-Expert Chargaff + 6F-ORF-HMM + Spaced Seeds + 6 Repeat Heads + SSE)
// ============================================================================

inline void normalize_p(int p0, int p1, int p2, int p3, int out[4]) {
    if (p0 < 1) p0 = 1;
    if (p1 < 1) p1 = 1;
    if (p2 < 1) p2 = 1;
    if (p3 < 1) p3 = 1;
    int diff = P_SCALE - (p0 + p1 + p2 + p3);
    if (diff != 0) {
        if (p0 >= p1 && p0 >= p2 && p0 >= p3) p0 += diff;
        else if (p1 >= p2 && p1 >= p3) p1 += diff;
        else if (p2 >= p3) p2 += diff;
        else p3 += diff;
    }
    out[0] = p0; out[1] = p1; out[2] = p2; out[3] = p3;
}

struct MicroRepeatHead {
    bool active = false;
    int64_t cursor = 0;
    int dir = 1;
    int misses = 0;
    int max_misses = 6;
    int streak = 0;
    int total_matched = 0;
    int xor_cnt[4] = {28, 2, 4, 2};

    void activate(int64_t c, int d, int init_match, int mm) {
        active = true; cursor = c; dir = d; misses = 0; max_misses = mm;
        streak = 0; total_matched = 0;
        xor_cnt[0] = init_match; xor_cnt[1] = 2; xor_cnt[2] = 4; xor_cnt[3] = 2;
    }

    bool cursor_accessible(int64_t seg_start, int64_t global_pos, int64_t warm_limit) const {
        if (cursor < 0 || cursor >= global_pos) return false;
        if (seg_start == 0) return true;
        return cursor < warm_limit || cursor >= seg_start;
    }

    bool predict(const uint8_t* buf, int64_t seg_start, int64_t global_pos, int64_t warm_limit, int out_p[4]) {
        if (!active || !cursor_accessible(seg_start, global_pos, warm_limit)) {
            active = false;
            return false;
        }
        int tpl = buf[cursor] & 3U;
        if (dir < 0) tpl ^= 3;
        int cm = xor_cnt[0] + std::min(96, (streak << 1) + (misses == 0 ? std::min(48, total_matched) : 0));
        int c1 = xor_cnt[1], c2 = xor_cnt[2], c3 = xor_cnt[3];
        int sum = cm + c1 + c2 + c3;
        const auto& luts = NoveltyLuts::get();
        int px[4];
        px[0] = std::max(1, luts.div_round_small(static_cast<int64_t>(cm) * P_SCALE, sum));
        px[1] = std::max(1, luts.div_round_small(static_cast<int64_t>(c1) * P_SCALE, sum));
        px[2] = std::max(1, luts.div_round_small(static_cast<int64_t>(c2) * P_SCALE, sum));
        px[3] = std::max(1, P_SCALE - (px[0] + px[1] + px[2]));
        normalize_p(px[tpl], px[1 ^ tpl], px[2 ^ tpl], px[3 ^ tpl], out_p);
        return true;
    }

    void step(const uint8_t* buf, int64_t seg_start, int64_t global_pos, int64_t warm_limit, int sym) {
        if (!active || !cursor_accessible(seg_start, global_pos, warm_limit)) {
            active = false; return;
        }
        int tpl = buf[cursor] & 3U;
        if (dir < 0) tpl ^= 3;
        int delta = (sym ^ tpl) & 3;
        xor_cnt[delta] += 8;
        if (xor_cnt[0] + xor_cnt[1] + xor_cnt[2] + xor_cnt[3] >= 256) {
            xor_cnt[0] = std::max(2, xor_cnt[0] >> 1);
            xor_cnt[1] = std::max(1, xor_cnt[1] >> 1);
            xor_cnt[2] = std::max(1, xor_cnt[2] >> 1);
            xor_cnt[3] = std::max(1, xor_cnt[3] >> 1);
        }
        if (delta == 0) { misses = 0; streak++; total_matched++; }
        else {
            misses++; streak = 0;
            if (misses > max_misses) { active = false; return; }
        }
        cursor += dir;
        if (!cursor_accessible(seg_start, global_pos + 1, warm_limit)) active = false;
    }
};

struct FastCrmSnapshot {
    std::vector<uint16_t> ord0, ord1, ord2, ord3, ord4, ord5, ord6, ord7;
    std::vector<uint32_t> ord8, ord9;
    std::vector<uint16_t> orf_plus2, orf_minus2, orf_plus3, orf_minus3, orf_plus5, orf_minus5;
    std::vector<uint16_t> wobble8, gapped8;
    std::vector<uint64_t> sparse12;
    std::vector<uint64_t> rep_fine_fwd, rep_fine_rev, rep_long_fwd, rep_long_rev;
    std::vector<uint16_t> sse_counts;
    int32_t weights[8][24];
    std::vector<int64_t> bayes_scores;
    int64_t meta_bayes = 0;
    int64_t meta_neural = -2048;
};

class FastCrmModel {
public:
    bool pan_super = false;
    bool use_ord8 = false;
    int64_t warm_prefix_limit = 0;
    int idx_first_head = 13;
    int num_experts = 19;
    int64_t rep_fine_mask = 0;
    int64_t rep_long_mask = 0;
    int64_t sparse_mask = 0;

    std::vector<uint16_t> ord0, ord1, ord2, ord3, ord4, ord5, ord6, ord7;
    std::vector<uint32_t> ord8, ord9;
    std::vector<uint16_t> orf_plus2, orf_minus2, orf_plus3, orf_minus3, orf_plus5, orf_minus5;
    std::vector<uint16_t> wobble8, gapped8;
    std::vector<uint64_t> sparse12;
    std::vector<uint64_t> rep_fine_fwd, rep_fine_rev, rep_long_fwd, rep_long_rev;
    std::vector<uint16_t> sse_counts;

    MicroRepeatHead heads[6];
    int32_t exp_p[24][4] {};
    int32_t exp_ln_p[24][4] {};
    int32_t bayes_p[4] {};
    int32_t neural_p[4] {};
    int32_t weights[8][24] {};
    std::vector<int64_t> bayes_scores;
    int64_t meta_bayes = 0;
    int64_t meta_neural = -2048;

    int64_t frame_log[6] {};
    int32_t frame_w[6] {};
    int frame_wsum = 6;
    int best_frame = 0;
    int p_plus2[3][4] {}, p_minus2[3][4] {}, p_plus3[3][4] {}, p_minus3[3][4] {}, p_plus5[3][4] {}, p_minus5[3][4] {};
    int last_sse_state = -1;
    uint64_t fwd_hist = 0, rc_hist = 0;
    int active_regime = 0;
    int64_t cached_sparse_idx = 0;
    uint32_t cached_sparse_tag = 0;

    FastCrmModel(const FastCrmSnapshot* snap, bool bounded_mem, bool pan, int64_t warm_lim, uint64_t total_len, int total_segs)
        : pan_super(pan), use_ord8(pan), warm_prefix_limit(warm_lim) {
        idx_first_head = pan_super ? 16 : (use_ord8 ? 14 : 13);
        num_experts = idx_first_head + 6;
        size_t sp_size = pan_super ? (1ULL << 18) : 0ULL;
        sparse_mask = static_cast<int64_t>(sp_size) - 1;
        size_t ring_buckets = (pan_super && total_segs >= 4)
            ? (1ULL << 16)
            : (bounded_mem ? (1ULL << 14) : (total_len >= 500000ULL ? (1ULL << 19) : (1ULL << 16)));
        rep_fine_mask = static_cast<int64_t>(ring_buckets - 1);
        rep_long_mask = static_cast<int64_t>(ring_buckets - 1);

        if (snap == nullptr) {
            ord0.assign(4, 0); ord1.assign(16, 0); ord2.assign(64, 0); ord3.assign(256, 0);
            ord4.assign(1024, 0); ord5.assign(4096, 0); ord6.assign(16384, 0); ord7.assign(65536, 0);
            if (use_ord8) ord8.assign(65536, 0);
            if (pan_super) { ord9.assign(262144, 0); sparse12.assign(sp_size, 0); }
            orf_plus2.assign(16 * 12, 0); orf_minus2.assign(16 * 12, 0);
            orf_plus3.assign(64 * 12, 0); orf_minus3.assign(64 * 12, 0);
            orf_plus5.assign(1024 * 12, 0); orf_minus5.assign(1024 * 12, 0);
            seed_universal_stop_priors();
            wobble8.assign(65536 * 4, 0);
            gapped8.assign(65536 * 4, 0);
            rep_fine_fwd.assign(ring_buckets << 1, 0); rep_fine_rev.assign(ring_buckets, 0);
            rep_long_fwd.assign(ring_buckets << 1, 0); rep_long_rev.assign(ring_buckets, 0);
            sse_counts.assign(1536 * 4, 0);

            bayes_scores.resize(num_experts);
            for (int k = 0; k <= 7; k++) bayes_scores[k] = (k <= 3) ? 0 : -6144;
            bayes_scores[8] = -2048; bayes_scores[9] = -3276; bayes_scores[10] = -4500;
            bayes_scores[11] = -4096; bayes_scores[12] = -4096;
            if (use_ord8) bayes_scores[13] = -6144;
            if (pan_super) { bayes_scores[14] = -6144; bayes_scores[15] = -6144; }
            for (int h = 0; h < 6; h++) bayes_scores[idx_first_head + h] = -4096;

            for (int r = 0; r < 8; r++) {
                bool rep_r = (r & 4) != 0, hi_r = (r & 1) != 0;
                weights[r][0] = 150;
                for (int k = 1; k <= 4; k++) weights[r][k] = hi_r ? 350 : 750;
                for (int k = 5; k <= 7; k++) weights[r][k] = hi_r ? 650 : 250;
                weights[r][8] = 220; weights[r][9] = 320; weights[r][10] = 360;
                weights[r][11] = 300; weights[r][12] = 300;
                if (use_ord8) weights[r][13] = hi_r ? 650 : 220;
                if (pan_super) { weights[r][14] = hi_r ? 780 : 150; weights[r][15] = hi_r ? 650 : 150; }
                for (int h = 0; h < 6; h++) {
                    weights[r][idx_first_head + h] = rep_r ? ((h == 0 || h == 2) ? 1200 : 600) : 100;
                }
            }
        } else {
            ord0 = snap->ord0; ord1 = snap->ord1; ord2 = snap->ord2; ord3 = snap->ord3;
            ord4 = snap->ord4; ord5 = snap->ord5; ord6 = snap->ord6; ord7 = snap->ord7;
            ord8 = snap->ord8; ord9 = snap->ord9; sparse12 = snap->sparse12;
            orf_plus2 = snap->orf_plus2; orf_minus2 = snap->orf_minus2;
            orf_plus3 = snap->orf_plus3; orf_minus3 = snap->orf_minus3;
            orf_plus5 = snap->orf_plus5; orf_minus5 = snap->orf_minus5;
            wobble8 = snap->wobble8; gapped8 = snap->gapped8;
            rep_fine_fwd = snap->rep_fine_fwd; rep_fine_rev = snap->rep_fine_rev;
            rep_long_fwd = snap->rep_long_fwd; rep_long_rev = snap->rep_long_rev;
            sse_counts = snap->sse_counts;
            for (int r = 0; r < 8; r++) std::memcpy(weights[r], snap->weights[r], sizeof(weights[r]));
            bayes_scores = snap->bayes_scores;
            meta_bayes = snap->meta_bayes;
            meta_neural = snap->meta_neural;
        }
    }

    std::unique_ptr<FastCrmSnapshot> snapshot() const {
        auto s = std::make_unique<FastCrmSnapshot>();
        s->ord0 = ord0; s->ord1 = ord1; s->ord2 = ord2; s->ord3 = ord3;
        s->ord4 = ord4; s->ord5 = ord5; s->ord6 = ord6; s->ord7 = ord7;
        s->ord8 = ord8; s->ord9 = ord9; s->sparse12 = sparse12;
        s->orf_plus2 = orf_plus2; s->orf_minus2 = orf_minus2;
        s->orf_plus3 = orf_plus3; s->orf_minus3 = orf_minus3;
        s->orf_plus5 = orf_plus5; s->orf_minus5 = orf_minus5;
        s->wobble8 = wobble8; s->gapped8 = gapped8;
        s->rep_fine_fwd = rep_fine_fwd; s->rep_fine_rev = rep_fine_rev;
        s->rep_long_fwd = rep_long_fwd; s->rep_long_rev = rep_long_rev;
        s->sse_counts = sse_counts;
        for (int r = 0; r < 8; r++) std::memcpy(s->weights[r], weights[r], sizeof(weights[r]));
        s->bayes_scores = bayes_scores;
        s->meta_bayes = meta_bayes;
        s->meta_neural = meta_neural;
        return s;
    }

    void predict(const uint8_t* buf, int64_t seg_start, int64_t rel_pos, int out_freqs[4]) {
        const auto& luts = NoveltyLuts::get();
        int64_t global_pos = seg_start + rel_pos;
        spawn_repeat_candidates(seg_start, rel_pos, global_pos);
        static const int UNIFORM[4] = {4096, 4096, 4096, 4096};
        int high_sum = 0;

        short_predict(ord0.data(), 0, UNIFORM, 32, exp_p[0]);
        if (rel_pos >= 1) short_predict(ord1.data(), static_cast<size_t>((fwd_hist & 3ULL) << 2), exp_p[0], 24, exp_p[1]);
        else std::memcpy(exp_p[1], exp_p[0], 16);
        if (rel_pos >= 2) short_predict(ord2.data(), static_cast<size_t>((fwd_hist & 15ULL) << 2), exp_p[1], 24, exp_p[2]);
        else std::memcpy(exp_p[2], exp_p[1], 16);
        if (rel_pos >= 3) short_predict(ord3.data(), static_cast<size_t>((fwd_hist & 63ULL) << 2), exp_p[2], 16, exp_p[3]);
        else std::memcpy(exp_p[3], exp_p[2], 16);
        if (rel_pos >= 4) short_predict(ord4.data(), static_cast<size_t>((fwd_hist & 255ULL) << 2), exp_p[3], 16, exp_p[4]);
        else std::memcpy(exp_p[4], exp_p[3], 16);
        if (rel_pos >= 5) short_predict(ord5.data(), static_cast<size_t>((fwd_hist & 1023ULL) << 2), exp_p[4], 10, exp_p[5]);
        else std::memcpy(exp_p[5], exp_p[4], 16);
        if (rel_pos >= 6) short_predict(ord6.data(), static_cast<size_t>((fwd_hist & 4095ULL) << 2), exp_p[5], 8, exp_p[6]);
        else std::memcpy(exp_p[6], exp_p[5], 16);
        if (rel_pos >= 7) high_sum += short_predict(ord7.data(), static_cast<size_t>((fwd_hist & 16383ULL) << 2), exp_p[6], 4, exp_p[7]);
        else std::memcpy(exp_p[7], exp_p[6], 16);

        int64_t max_fl = frame_log[0];
        for (int s = 1; s < 6; s++) max_fl = std::max(max_fl, frame_log[s]);
        int fw_sum = 0, best_s = 0;
        for (int s = 0; s < 6; s++) {
            if (max_fl - frame_log[s] > 26624) frame_log[s] = max_fl - 26624;
            int w = std::max(1, luts.exp_neg_lut[static_cast<size_t>(max_fl - frame_log[s])] >> 3);
            frame_w[s] = w;
            fw_sum += w;
            if (w > frame_w[best_s]) best_s = s;
        }
        best_frame = best_s;
        frame_wsum = fw_sum;
        int rmod3 = static_cast<int>(rel_pos % 3);
        int phase = rmod3;
        if (frame_w[best_s] * 2 >= fw_sum) {
            int cb = rmod3 - (best_s < 3 ? best_s : best_s - 3);
            if (cb < 0) cb += 3;
            phase = (best_s < 3) ? cb : (2 - cb);
        }

        size_t c2o = static_cast<size_t>(fwd_hist & 15ULL);
        size_t c3o = static_cast<size_t>(fwd_hist & 63ULL);
        size_t c5o = static_cast<size_t>(fwd_hist & 1023ULL);
        for (int c = 0; c < 3; c++) {
            if (rel_pos >= 2) {
                short_predict(orf_plus2.data(), ((c2o * 3) + c) << 2, exp_p[1], 16, p_plus2[c]);
                short_predict(orf_minus2.data(), ((c2o * 3) + c) << 2, exp_p[1], 16, p_minus2[c]);
            } else {
                std::memcpy(p_plus2[c], exp_p[1], 16);
                std::memcpy(p_minus2[c], exp_p[1], 16);
            }
            if (rel_pos >= 3) {
                short_predict(orf_plus3.data(), ((c3o * 3) + c) << 2, p_plus2[c], 12, p_plus3[c]);
                short_predict(orf_minus3.data(), ((c3o * 3) + c) << 2, p_minus2[c], 12, p_minus3[c]);
            } else {
                std::memcpy(p_plus3[c], p_plus2[c], 16);
                std::memcpy(p_minus3[c], p_minus2[c], 16);
            }
            if (rel_pos >= 5) {
                short_predict(orf_plus5.data(), ((c5o * 3) + c) << 2, p_plus3[c], 8, p_plus5[c]);
                short_predict(orf_minus5.data(), ((c5o * 3) + c) << 2, p_minus3[c], 8, p_minus5[c]);
            } else {
                std::memcpy(p_plus5[c], p_plus3[c], 16);
                std::memcpy(p_minus5[c], p_minus3[c], 16);
            }
        }
        int64_t a2[4]{}, a3[4]{}, a5[4]{};
        for (int s = 0; s < 6; s++) {
            int64_t ws = frame_w[s];
            int c = rmod3 - (s < 3 ? s : s - 3);
            if (c < 0) c += 3;
            const int* q2 = (s < 3) ? p_plus2[c] : p_minus2[c];
            const int* q3 = (s < 3) ? p_plus3[c] : p_minus3[c];
            const int* q5 = (s < 3) ? p_plus5[c] : p_minus5[c];
            for (int k = 0; k < 4; k++) {
                a2[k] += ws * q2[k]; a3[k] += ws * q3[k]; a5[k] += ws * q5[k];
            }
        }
        normalize_p(luts.div_round_small(a2[0], fw_sum), luts.div_round_small(a2[1], fw_sum), luts.div_round_small(a2[2], fw_sum), luts.div_round_small(a2[3], fw_sum), exp_p[8]);
        normalize_p(luts.div_round_small(a3[0], fw_sum), luts.div_round_small(a3[1], fw_sum), luts.div_round_small(a3[2], fw_sum), luts.div_round_small(a3[3], fw_sum), exp_p[9]);
        normalize_p(luts.div_round_small(a5[0], fw_sum), luts.div_round_small(a5[1], fw_sum), luts.div_round_small(a5[2], fw_sum), luts.div_round_small(a5[3], fw_sum), exp_p[10]);

        if (rel_pos >= 11) short_predict(wobble8.data(), extract_wobble_key(fwd_hist) << 2, exp_p[4], 8, exp_p[11]);
        else std::memcpy(exp_p[11], exp_p[4], 16);
        if (rel_pos >= 9) short_predict(gapped8.data(), extract_gapped9_key(fwd_hist) << 2, exp_p[4], 8, exp_p[12]);
        else std::memcpy(exp_p[12], exp_p[4], 16);

        if (use_ord8) {
            if (rel_pos >= 8) {
                uint32_t p8 = ord8[static_cast<size_t>(fwd_hist & 65535ULL)];
                high_sum += swar_sum(p8);
                swar_predict(p8, exp_p[7], 3, exp_p[13]);
            } else std::memcpy(exp_p[13], exp_p[7], 16);
        }
        if (pan_super) {
            if (rel_pos >= 9) {
                uint32_t p9 = ord9[static_cast<size_t>(fwd_hist & 262143ULL)];
                high_sum += swar_sum(p9);
                swar_predict(p9, exp_p[13], 2, exp_p[14]);
            } else std::memcpy(exp_p[14], exp_p[13], 16);
            if (rel_pos >= 12) high_sum += predict_sparse12(fwd_hist & 0xFFFFFFULL, exp_p[14], exp_p[15]);
            else std::memcpy(exp_p[15], exp_p[14], 16);
        }

        bool strong_rep = false;
        for (int h = 0; h < 6; h++) {
            if (heads[h].predict(buf, seg_start, global_pos, warm_prefix_limit, exp_p[idx_first_head + h])) {
                if (heads[h].streak >= 2 || heads[h].total_matched >= 6) strong_rep = true;
            } else {
                std::memcpy(exp_p[idx_first_head + h], exp_p[4], 16);
            }
        }
        active_regime = (strong_rep ? 4 : 0) | ((phase == 2) ? 2 : 0) | ((high_sum >= 4) ? 1 : 0);

        int64_t max_b = bayes_scores[0];
        for (int e = 1; e < num_experts; e++) max_b = std::max(max_b, bayes_scores[e]);
        int64_t ba[4]{}, bw_sum = 0, logit[4]{};
        const int32_t* w_reg = weights[active_regime];
        for (int e = 0; e < num_experts; e++) {
            for (int k = 0; k < 4; k++) exp_ln_p[e][k] = luts.ln_lut[exp_p[e][k]];
            int64_t diff = max_b - bayes_scores[e];
            int64_t bw = luts.exp_neg_lut[static_cast<size_t>(diff >= EXP_MAX_X ? EXP_MAX_X : diff)];
            for (int k = 0; k < 4; k++) {
                ba[k] += bw * exp_p[e][k];
                logit[k] += static_cast<int64_t>(w_reg[e]) * exp_ln_p[e][k];
            }
            bw_sum += bw;
        }
        int64_t b_half = bw_sum >> 1, b_inv = (1LL << 36) / bw_sum;
        normalize_p(
            div_fast(ba[0], bw_sum, b_half, b_inv, 36),
            div_fast(ba[1], bw_sum, b_half, b_inv, 36),
            div_fast(ba[2], bw_sum, b_half, b_inv, 36),
            div_fast(ba[3], bw_sum, b_half, b_inv, 36),
            bayes_p);

        int64_t l0 = logit[0] >> 12, l1 = logit[1] >> 12, l2 = logit[2] >> 12, l3 = logit[3] >> 12;
        int64_t max_l = std::max({l0, l1, l2, l3});
        int e0 = luts.exp_neg_lut[static_cast<size_t>(std::min<int64_t>(EXP_MAX_X, max_l - l0))];
        int e1 = luts.exp_neg_lut[static_cast<size_t>(std::min<int64_t>(EXP_MAX_X, max_l - l1))];
        int e2 = luts.exp_neg_lut[static_cast<size_t>(std::min<int64_t>(EXP_MAX_X, max_l - l2))];
        int e3 = luts.exp_neg_lut[static_cast<size_t>(std::min<int64_t>(EXP_MAX_X, max_l - l3))];
        int esum = e0 + e1 + e2 + e3;
        int64_t ehalf = esum >> 1, einv = (1LL << 32) / esum;
        normalize_p(
            div_fast(static_cast<int64_t>(e0) * P_SCALE, esum, ehalf, einv, 32),
            div_fast(static_cast<int64_t>(e1) * P_SCALE, esum, ehalf, einv, 32),
            div_fast(static_cast<int64_t>(e2) * P_SCALE, esum, ehalf, einv, 32),
            div_fast(static_cast<int64_t>(e3) * P_SCALE, esum, ehalf, einv, 32),
            neural_p);

        int64_t max_m = std::max(meta_bayes, meta_neural);
        int64_t mb = luts.exp_neg_lut[static_cast<size_t>(std::min<int64_t>(EXP_MAX_X, max_m - meta_bayes))];
        int64_t mn = luts.exp_neg_lut[static_cast<size_t>(std::min<int64_t>(EXP_MAX_X, max_m - meta_neural))];
        int64_t msum = mb + mn, mhalf = msum >> 1, minv = (1LL << 34) / msum;
        int f0 = std::max(1, div_fast((mb * bayes_p[0] + mn * neural_p[0]) << 1, msum, mhalf, minv, 34));
        int f1 = std::max(1, div_fast((mb * bayes_p[1] + mn * neural_p[1]) << 1, msum, mhalf, minv, 34));
        int f2 = std::max(1, div_fast((mb * bayes_p[2] + mn * neural_p[2]) << 1, msum, mhalf, minv, 34));
        int f3 = std::max(1, div_fast((mb * bayes_p[3] + mn * neural_p[3]) << 1, msum, mhalf, minv, 34));
        int fix = NOV_TOT - (f0 + f1 + f2 + f3);
        if (fix != 0) {
            if (f0 >= f1 && f0 >= f2 && f0 >= f3) f0 += fix;
            else if (f1 >= f2 && f1 >= f3) f1 += fix;
            else if (f2 >= f3) f2 += fix;
            else f3 += fix;
        }

        if (rel_pos >= 2) {
            int sbest = 0, fmax = f0;
            if (f1 > fmax) { sbest = 1; fmax = f1; }
            if (f2 > fmax) { sbest = 2; fmax = f2; }
            if (f3 > fmax) { sbest = 3; fmax = f3; }
            int c2 = static_cast<int>(fwd_hist & 15ULL);
            int bin = std::min(7, std::max(0, (fmax - 8192) / 3072));
            int sstate = (((phase * 16 + c2) << 2) | sbest) * 8 + bin;
            last_sse_state = sstate;
            size_t sb = static_cast<size_t>(sstate) << 2;
            int sc0 = sse_counts[sb], sc1 = sse_counts[sb+1], sc2 = sse_counts[sb+2], sc3 = sse_counts[sb+3];
            int scsum = sc0 + sc1 + sc2 + sc3;
            if (scsum > 0) {
                int w = std::min(8, scsum >> 1);
                int sf0 = luts.div_round_small(static_cast<int64_t>(sc0) * NOV_TOT, scsum);
                int sf1 = luts.div_round_small(static_cast<int64_t>(sc1) * NOV_TOT, scsum);
                int sf2 = luts.div_round_small(static_cast<int64_t>(sc2) * NOV_TOT, scsum);
                int sf3 = luts.div_round_small(static_cast<int64_t>(sc3) * NOV_TOT, scsum);
                f0 = std::max(1, ((32 - w) * f0 + w * sf0 + 16) >> 5);
                f1 = std::max(1, ((32 - w) * f1 + w * sf1 + 16) >> 5);
                f2 = std::max(1, ((32 - w) * f2 + w * sf2 + 16) >> 5);
                f3 = std::max(1, ((32 - w) * f3 + w * sf3 + 16) >> 5);
                fix = NOV_TOT - (f0 + f1 + f2 + f3);
                if (fix != 0) {
                    if (f0 >= f1 && f0 >= f2 && f0 >= f3) f0 += fix;
                    else if (f1 >= f2 && f1 >= f3) f1 += fix;
                    else if (f2 >= f3) f2 += fix;
                    else f3 += fix;
                }
            }
        } else {
            last_sse_state = -1;
        }
        out_freqs[0] = f0; out_freqs[1] = f1; out_freqs[2] = f2; out_freqs[3] = f3;
    }

    void update(const uint8_t* buf, int64_t seg_start, int64_t rel_pos, int sym) {
        const auto& luts = NoveltyLuts::get();
        int64_t global_pos = seg_start + rel_pos;

        meta_bayes = ((meta_bayes * 255) >> 8) + ((3 * static_cast<int64_t>(luts.ln_lut[bayes_p[sym]])) >> 1);
        meta_neural = ((meta_neural * 255) >> 8) + ((3 * static_cast<int64_t>(luts.ln_lut[neural_p[sym]])) >> 1);
        for (int e = 0; e < num_experts; e++) {
            bayes_scores[e] = ((bayes_scores[e] * 255) >> 8) + ((13 * static_cast<int64_t>(exp_ln_p[e][sym])) >> 3);
        }
        int r0 = (sym == 0 ? P_SCALE : 0) - neural_p[0];
        int r1 = (sym == 1 ? P_SCALE : 0) - neural_p[1];
        int r2 = (sym == 2 ? P_SCALE : 0) - neural_p[2];
        int r3 = (sym == 3 ? P_SCALE : 0) - neural_p[3];
        int32_t* rw = weights[active_regime];
        for (int e = 0; e < num_experts; e++) {
            int64_t grad = static_cast<int64_t>(r0) * exp_ln_p[e][0] + static_cast<int64_t>(r1) * exp_ln_p[e][1] +
                           static_cast<int64_t>(r2) * exp_ln_p[e][2] + static_cast<int64_t>(r3) * exp_ln_p[e][3];
            int upd = rw[e] - (rw[e] >> 11) + static_cast<int>(grad >> 20);
            rw[e] = std::clamp(upd, 0, 8192);
        }
        if (last_sse_state >= 0) {
            size_t sb = static_cast<size_t>(last_sse_state) << 2;
            sse_counts[sb + sym] = static_cast<uint16_t>(sse_counts[sb + sym] + 2);
            int ssum = sse_counts[sb] + sse_counts[sb+1] + sse_counts[sb+2] + sse_counts[sb+3];
            if (ssum >= 128) {
                for (int k = 0; k < 4; k++) sse_counts[sb + k] >>= 1;
            }
        }
        for (int h = 0; h < 6; h++) heads[h].step(buf, seg_start, global_pos, warm_prefix_limit, sym);

        uint64_t pf = fwd_hist;
        rc_hist = (rc_hist >> 2) | (static_cast<uint64_t>((sym ^ 3) & 3) << 60);
        if (rel_pos >= 13) {
            insert_ring_fwd(rep_fine_fwd, rep_fine_mask, pf & ((1ULL << 26) - 1ULL), global_pos);
            insert_ring_rev(rep_fine_rev, rep_fine_mask, (rc_hist >> (62 - 26)) & ((1ULL << 26) - 1ULL), global_pos - 13);
        }
        if (rel_pos >= 18) {
            insert_ring_fwd(rep_long_fwd, rep_long_mask, pf & ((1ULL << 36) - 1ULL), global_pos);
            insert_ring_rev(rep_long_rev, rep_long_mask, (rc_hist >> (62 - 36)) & ((1ULL << 36) - 1ULL), global_pos - 18);
        }

        short_bump(ord0.data(), 0, sym, 2, 16384);
        if (rel_pos >= 1) short_bump(ord1.data(), static_cast<size_t>((pf & 3ULL) << 2), sym, 2, 16384);
        if (rel_pos >= 2) short_bump(ord2.data(), static_cast<size_t>((pf & 15ULL) << 2), sym, 2, 16384);
        if (rel_pos >= 3) short_bump(ord3.data(), static_cast<size_t>((pf & 63ULL) << 2), sym, 2, 8192);
        if (rel_pos >= 4) {
            short_bump(ord4.data(), static_cast<size_t>((pf & 255ULL) << 2), sym, 2, 8192);
            short_bump(ord4.data(), static_cast<size_t>(((rc_hist >> (62 - 8)) & 255ULL) << 2), static_cast<int>(((pf >> 6) & 3ULL) ^ 3ULL), 1, 8192);
        }
        if (rel_pos >= 5) {
            short_bump(ord5.data(), static_cast<size_t>((pf & 1023ULL) << 2), sym, 2, 4096);
            short_bump(ord5.data(), static_cast<size_t>(((rc_hist >> (62 - 10)) & 1023ULL) << 2), static_cast<int>(((pf >> 8) & 3ULL) ^ 3ULL), 1, 4096);
        }
        if (rel_pos >= 6) {
            short_bump(ord6.data(), static_cast<size_t>((pf & 4095ULL) << 2), sym, 2, 4096);
            short_bump(ord6.data(), static_cast<size_t>(((rc_hist >> (62 - 12)) & 4095ULL) << 2), static_cast<int>(((pf >> 10) & 3ULL) ^ 3ULL), 1, 4096);
        }
        if (rel_pos >= 7) {
            short_bump(ord7.data(), static_cast<size_t>((pf & 16383ULL) << 2), sym, 2, 2048);
            short_bump(ord7.data(), static_cast<size_t>(((rc_hist >> (62 - 14)) & 16383ULL) << 2), static_cast<int>(((pf >> 12) & 3ULL) ^ 3ULL), 1, 2048);
        }

        int rmod3 = static_cast<int>(rel_pos % 3);
        for (int s = 0; s < 6; s++) {
            int c = rmod3 - (s < 3 ? s : s - 3);
            if (c < 0) c += 3;
            int p3s = (s < 3) ? p_plus3[c][sym] : p_minus3[c][sym];
            int p5s = (s < 3) ? p_plus5[c][sym] : p_minus5[c][sym];
            int pm = (3 * p3s + 5 * p5s) >> 3;
            frame_log[s] = ((frame_log[s] * 255) >> 8) + luts.ln_lut[pm];
        }
        if (frame_w[best_frame] * 9 >= frame_wsum * 4) {
            int c = rmod3 - (best_frame < 3 ? best_frame : best_frame - 3);
            if (c < 0) c += 3;
            bool plus = best_frame < 3;
            if (rel_pos >= 2) short_bump((plus ? orf_plus2 : orf_minus2).data(), static_cast<size_t>(((((pf & 15ULL) * 3) + c)) << 2), sym, 2, 4096);
            if (rel_pos >= 3) short_bump((plus ? orf_plus3 : orf_minus3).data(), static_cast<size_t>(((((pf & 63ULL) * 3) + c)) << 2), sym, 2, 4096);
            if (rel_pos >= 5) short_bump((plus ? orf_plus5 : orf_minus5).data(), static_cast<size_t>(((((pf & 1023ULL) * 3) + c)) << 2), sym, 2, 4096);
        }
        if (rel_pos >= 11) short_bump(wobble8.data(), extract_wobble_key(pf) << 2, sym, 2, 1024);
        if (rel_pos >= 9) {
            short_bump(gapped8.data(), extract_gapped9_key(pf) << 2, sym, 2, 1024);
            short_bump(gapped8.data(), extract_gapped9_key(rc_hist >> (62 - 18)) << 2, static_cast<int>(((pf >> 16) & 3ULL) ^ 3ULL), 1, 1024);
        }
        if (use_ord8 && rel_pos >= 8) {
            size_t i8 = static_cast<size_t>(pf & 65535ULL);
            ord8[i8] = swar_bump(ord8[i8], sym, 2);
            size_t rc8 = static_cast<size_t>((rc_hist >> (62 - 16)) & 65535ULL);
            ord8[rc8] = swar_bump(ord8[rc8], static_cast<int>(((pf >> 14) & 3ULL) ^ 3ULL), 1);
        }
        if (pan_super) {
            if (rel_pos >= 9) {
                size_t i9 = static_cast<size_t>(pf & 262143ULL);
                ord9[i9] = swar_bump(ord9[i9], sym, 2);
                size_t rc9 = static_cast<size_t>((rc_hist >> (62 - 18)) & 262143ULL);
                ord9[rc9] = swar_bump(ord9[rc9], static_cast<int>(((pf >> 16) & 3ULL) ^ 3ULL), 1);
            }
            if (rel_pos >= 12) {
                update_sparse12_slot(cached_sparse_idx, cached_sparse_tag, sym, 2);
                update_sparse12_key((rc_hist >> (62 - 24)) & 0xFFFFFFULL, static_cast<int>(((pf >> 22) & 3ULL) ^ 3ULL), 1);
            }
        }
        fwd_hist = (pf << 2) | static_cast<uint64_t>(sym & 3);
    }

private:
    void seed_universal_stop_priors() {
        for (size_t ctx2 = 0; ctx2 < 16; ctx2++) {
            size_t bp2 = ((ctx2 * 3) + 2) << 2;
            if (ctx2 == 12) { orf_plus2[bp2+1] += 4; orf_plus2[bp2+3] += 4; }
            else if (ctx2 == 14) { orf_plus2[bp2+1] += 3; orf_plus2[bp2+2] += 3; orf_plus2[bp2+3] += 3; }
            if (ctx2 == 15 || ctx2 == 7 || ctx2 == 13) { orf_minus2[bp2+1] += 3; orf_minus2[bp2+2] += 3; orf_minus2[bp2+3] += 3; }
            for (size_t b0 = 0; b0 < 4; b0++) {
                size_t ctx3 = (b0 << 4) | ctx2;
                size_t bp3 = ((ctx3 * 3) + 2) << 2;
                if (ctx2 == 12) { orf_plus3[bp3+1] += 4; orf_plus3[bp3+3] += 4; }
                else if (ctx2 == 14) { orf_plus3[bp3+1] += 3; orf_plus3[bp3+2] += 3; orf_plus3[bp3+3] += 3; }
                if (ctx2 == 15 || ctx2 == 7 || ctx2 == 13) { orf_minus3[bp3+1] += 3; orf_minus3[bp3+2] += 3; orf_minus3[bp3+3] += 3; }
            }
        }
    }

    static uint64_t mix_seed(uint64_t x) {
        x ^= (x >> 23);
        x *= 0xFF51AFD7ED558CCDULL;
        x ^= (x >> 29);
        return x;
    }

    bool is_cursor_tracked(int64_t cand, int dir) const {
        for (int h = 0; h < 6; h++) {
            if (heads[h].active && heads[h].dir == dir && heads[h].cursor == cand) return true;
        }
        return false;
    }

    void promote_veteran(const MicroRepeatHead& ev) {
        if (ev.active && ev.total_matched >= 20) {
            if (!heads[5].active || ev.total_matched > heads[5].total_matched) {
                heads[5] = ev;
                heads[5].max_misses = 9;
            }
        }
    }

    void allocate_head_slot(int a, int b, int64_t cur, int dir, int initw, int maxm) {
        if (!heads[a].active) heads[a].activate(cur, dir, initw, maxm);
        else if (!heads[b].active) heads[b].activate(cur, dir, initw, maxm);
        else if (heads[b].misses >= 2 && heads[b].streak == 0) {
            promote_veteran(heads[b]);
            heads[b].activate(cur, dir, initw, maxm);
        }
    }

    void spawn_repeat_candidates(int64_t seg_start, int64_t rel_pos, int64_t global_pos) {
        auto valid = [&](int64_t c) {
            if (c < 0 || c >= global_pos) return false;
            return (seg_start == 0) ? true : (c < warm_prefix_limit || c >= seg_start);
        };
        if (rel_pos >= 18) {
            uint64_t h = mix_seed(fwd_hist & ((1ULL << 36) - 1ULL));
            int64_t slot = static_cast<int64_t>(h) & rep_long_mask;
            int64_t fslot = slot << 1;
            uint32_t tag = static_cast<uint32_t>(h >> 32); if (tag == 0) tag = 1;
            for (int k = 0; k < 2; k++) {
                uint64_t e = rep_long_fwd[static_cast<size_t>(fslot + k)];
                if (static_cast<uint32_t>(e >> 32) == tag) {
                    int64_t c = static_cast<int32_t>(e);
                    if (valid(c) && !is_cursor_tracked(c, 1)) allocate_head_slot(0, 1, c, 1, 56, 8);
                }
            }
            uint64_t er = rep_long_rev[static_cast<size_t>(slot)];
            if (static_cast<uint32_t>(er >> 32) == tag) {
                int64_t cr = static_cast<int32_t>(er);
                if (cr + 18 < global_pos && valid(cr) && !is_cursor_tracked(cr, -1)) allocate_head_slot(1, 0, cr, -1, 56, 8);
            }
        }
        if (rel_pos >= 13) {
            uint64_t h = mix_seed(fwd_hist & ((1ULL << 26) - 1ULL));
            int64_t slot = static_cast<int64_t>(h) & rep_fine_mask;
            int64_t fslot = slot << 1;
            uint32_t tag = static_cast<uint32_t>(h >> 32); if (tag == 0) tag = 1;
            for (int k = 0; k < 2; k++) {
                uint64_t e = rep_fine_fwd[static_cast<size_t>(fslot + k)];
                if (static_cast<uint32_t>(e >> 32) == tag) {
                    int64_t c = static_cast<int32_t>(e);
                    if (valid(c) && !is_cursor_tracked(c, 1)) allocate_head_slot(2, 3, c, 1, 28, 6);
                }
            }
            uint64_t er = rep_fine_rev[static_cast<size_t>(slot)];
            if (static_cast<uint32_t>(er >> 32) == tag) {
                int64_t cr = static_cast<int32_t>(er);
                if (cr + 13 < global_pos && valid(cr) && !is_cursor_tracked(cr, -1)) {
                    if (!heads[4].active || (heads[4].streak == 0 && heads[4].misses >= 2)) {
                        promote_veteran(heads[4]);
                        heads[4].activate(cr, -1, 28, 6);
                    }
                }
            }
        }
    }

    static void insert_ring_fwd(std::vector<uint64_t>& ring, int64_t mask, uint64_t kmer, int64_t next_gp) {
        uint64_t h = mix_seed(kmer);
        size_t s0 = static_cast<size_t>((static_cast<int64_t>(h) & mask) << 1);
        uint32_t tag = static_cast<uint32_t>(h >> 32); if (tag == 0) tag = 1;
        uint64_t packed = (static_cast<uint64_t>(tag) << 32) | (static_cast<uint64_t>(static_cast<uint32_t>(next_gp)));
        uint64_t e0 = ring[s0];
        if (static_cast<uint32_t>(e0 >> 32) == tag) {
            if (static_cast<int32_t>(next_gp) - static_cast<int32_t>(e0) >= 4) ring[s0 + 1] = e0;
            ring[s0] = packed;
        } else {
            ring[s0 + 1] = e0;
            ring[s0] = packed;
        }
    }

    static void insert_ring_rev(std::vector<uint64_t>& ring, int64_t mask, uint64_t rc_kmer, int64_t pred_gp) {
        uint64_t h = mix_seed(rc_kmer);
        size_t slot = static_cast<size_t>(static_cast<int64_t>(h) & mask);
        uint32_t tag = static_cast<uint32_t>(h >> 32); if (tag == 0) tag = 1;
        ring[slot] = (static_cast<uint64_t>(tag) << 32) | static_cast<uint64_t>(static_cast<uint32_t>(pred_gp));
    }

    static size_t extract_wobble_key(uint64_t h) {
        return static_cast<size_t>((h & 0xFULL) | (((h >> 6) & 0xFULL) << 4) | (((h >> 12) & 0xFULL) << 8) | (((h >> 18) & 0xFULL) << 12));
    }
    static size_t extract_gapped9_key(uint64_t h) {
        return static_cast<size_t>((h & 0xFFULL) | (((h >> 10) & 0xFFULL) << 8));
    }

    static int div_fast(int64_t num, int64_t denom, int64_t half, int64_t inv, int shift) {
        int64_t n = num + half;
        int64_t q = (n * inv) >> shift;
        if (n - q * denom >= denom) q++;
        return static_cast<int>(q);
    }

    static int div_recip(int64_t num, int denom, int half, int64_t r) {
        int64_t n = num + half;
        int64_t q = (n * r) >> 30;
        if (n - q * denom >= denom) q++;
        return static_cast<int>(q);
    }

    static int short_predict(const uint16_t* tab, size_t base, const int prior[4], int tau, int out_p[4]) {
        int c0 = tab[base], c1 = tab[base+1], c2 = tab[base+2], c3 = tab[base+3];
        int sum = c0 + c1 + c2 + c3;
        if (sum == 0) { std::memcpy(out_p, prior, 16); return 0; }
        int denom = (sum << 2) + tau, half = denom >> 1;
        int64_t r = NoveltyLuts::get().recip[denom];
        int p0 = std::max(1, div_recip((static_cast<int64_t>(c0) << 2) * P_SCALE + static_cast<int64_t>(tau) * prior[0], denom, half, r));
        int p1 = std::max(1, div_recip((static_cast<int64_t>(c1) << 2) * P_SCALE + static_cast<int64_t>(tau) * prior[1], denom, half, r));
        int p2 = std::max(1, div_recip((static_cast<int64_t>(c2) << 2) * P_SCALE + static_cast<int64_t>(tau) * prior[2], denom, half, r));
        int p3 = std::max(1, P_SCALE - (p0 + p1 + p2));
        normalize_p(p0, p1, p2, p3, out_p);
        return sum;
    }

    static void short_bump(uint16_t* tab, size_t base, int sym, int inc, int max_sum) {
        tab[base + sym] = static_cast<uint16_t>(tab[base + sym] + inc);
        int sum = tab[base] + tab[base+1] + tab[base+2] + tab[base+3];
        if (sum >= max_sum) {
            for (int k = 0; k < 4; k++) tab[base + k] >>= 1;
        }
    }

    static void swar_predict(uint32_t packed, const int prior[4], int tau, int out_p[4]) {
        if (packed == 0) { std::memcpy(out_p, prior, 16); return; }
        int c0 = packed & 0xFFU, c1 = (packed >> 8) & 0xFFU, c2 = (packed >> 16) & 0xFFU, c3 = (packed >> 24) & 0xFFU;
        int sum = c0 + c1 + c2 + c3;
        int denom = (sum << 2) + tau, half = denom >> 1;
        int64_t r = NoveltyLuts::get().recip[denom];
        int p0 = std::max(1, div_recip((static_cast<int64_t>(c0) << 2) * P_SCALE + static_cast<int64_t>(tau) * prior[0], denom, half, r));
        int p1 = std::max(1, div_recip((static_cast<int64_t>(c1) << 2) * P_SCALE + static_cast<int64_t>(tau) * prior[1], denom, half, r));
        int p2 = std::max(1, div_recip((static_cast<int64_t>(c2) << 2) * P_SCALE + static_cast<int64_t>(tau) * prior[2], denom, half, r));
        int p3 = std::max(1, P_SCALE - (p0 + p1 + p2));
        normalize_p(p0, p1, p2, p3, out_p);
    }

    static int swar_sum(uint32_t p) {
        return static_cast<int>((p & 0xFFU) + ((p >> 8) & 0xFFU) + ((p >> 16) & 0xFFU) + ((p >> 24) & 0xFFU));
    }

    static uint32_t swar_bump(uint32_t packed, int sym, int inc) {
        int sh = sym << 3;
        uint32_t c = (packed >> sh) & 0xFFU;
        if (c + static_cast<uint32_t>(inc) >= 250U) packed = (packed >> 1) & 0x7F7F7F7FU;
        return packed + (static_cast<uint32_t>(inc) << sh);
    }

    int predict_sparse12(uint64_t ctx24, const int prior[4], int out_p[4]) {
        uint64_t mixed = hgs_mix64(ctx24);
        int64_t idx0 = (static_cast<int64_t>(mixed) & sparse_mask) & ~1LL;
        uint32_t tag = static_cast<uint32_t>(mixed >> 32); if (tag == 0) tag = 1;
        cached_sparse_idx = idx0;
        cached_sparse_tag = tag;
        uint64_t e = sparse12[static_cast<size_t>(idx0)];
        if (static_cast<uint32_t>(e >> 32) != tag) {
            e = sparse12[static_cast<size_t>(idx0 + 1)];
            if (static_cast<uint32_t>(e >> 32) != tag) {
                std::memcpy(out_p, prior, 16);
                return 0;
            }
        }
        uint32_t low = static_cast<uint32_t>(e);
        swar_predict(low, prior, 2, out_p);
        return swar_sum(low);
    }

    void update_sparse12_key(uint64_t ctx24, int sym, int inc) {
        uint64_t mixed = hgs_mix64(ctx24);
        int64_t idx0 = (static_cast<int64_t>(mixed) & sparse_mask) & ~1LL;
        uint32_t tag = static_cast<uint32_t>(mixed >> 32); if (tag == 0) tag = 1;
        update_sparse12_slot(idx0, tag, sym, inc);
    }

    void update_sparse12_slot(int64_t idx0, uint32_t tag, int sym, int inc) {
        size_t i0 = static_cast<size_t>(idx0), i1 = static_cast<size_t>(idx0 + 1);
        uint64_t e0 = sparse12[i0];
        uint32_t t0 = static_cast<uint32_t>(e0 >> 32);
        if (t0 == tag || t0 == 0) {
            sparse12[i0] = (static_cast<uint64_t>(tag) << 32) | swar_bump(t0 == 0 ? 0U : static_cast<uint32_t>(e0), sym, inc);
            return;
        }
        uint64_t e1 = sparse12[i1];
        uint32_t t1 = static_cast<uint32_t>(e1 >> 32);
        if (t1 == tag || t1 == 0) {
            sparse12[i1] = (static_cast<uint64_t>(tag) << 32) | swar_bump(t1 == 0 ? 0U : static_cast<uint32_t>(e1), sym, inc);
            return;
        }
        int s0 = swar_sum(static_cast<uint32_t>(e0)), s1 = swar_sum(static_cast<uint32_t>(e1));
        if (s0 <= s1 && s0 <= 6) sparse12[i0] = (static_cast<uint64_t>(tag) << 32) | swar_bump(0U, sym, inc);
        else if (s1 < s0 && s1 <= 6) sparse12[i1] = (static_cast<uint64_t>(tag) << 32) | swar_bump(0U, sym, inc);
    }
};

void encode_crm_segment(const uint8_t* bases, size_t start, size_t len, FastCrmModel& model, std::vector<uint8_t>& out_bytes) {
    size_t raw_bytes = (len + 3) >> 2;
    ByteVecOut crm;
    crm.reserve(raw_bytes + 64);
    crm.write_byte(SEG_CRM);
    NovRangeEnc rc(crm);
    int freqs[4];

    for (size_t i = 0; i < len; i++) {
        model.predict(bases, static_cast<int64_t>(start), static_cast<int64_t>(i), freqs);
        int sym = bases[start + i] & 3U;
        uint32_t cum = 0;
        for (int s = 0; s < sym; s++) cum += static_cast<uint32_t>(freqs[s]);
        rc.encode_symbol(cum, static_cast<uint32_t>(freqs[sym]));
        model.update(bases, static_cast<int64_t>(start), static_cast<int64_t>(i), sym);
    }
    rc.flush();
    if (crm.buf.size() <= 1 + raw_bytes) {
        out_bytes = std::move(crm.buf);
    } else {
        pack_raw_2bit_seg(bases + start, len, out_bytes);
    }
}

bool decode_crm_segment(const uint8_t* payload, size_t payload_len, uint8_t* dst, size_t dst_start, size_t len, FastCrmModel& model, bool build_warm) {
    if (payload_len == 0) return false;
    if (payload[0] == SEG_RAW_2BIT) {
        size_t expected = 1 + ((len + 3) >> 2);
        if (payload_len != expected) return false;
        size_t idx = dst_start, end = dst_start + len;
        for (size_t b = 1; b < payload_len && idx < end; b++) {
            uint8_t v = payload[b];
            for (int j = 0; j < 4 && idx < end; j++) {
                dst[idx++] = static_cast<uint8_t>(v & 3U);
                v >>= 2;
            }
        }
        if (build_warm) {
            int dummy[4];
            for (size_t r = 0; r < len; r++) {
                model.predict(dst, 0, static_cast<int64_t>(r), dummy);
                model.update(dst, 0, static_cast<int64_t>(r), dst[r] & 3U);
            }
        }
        return true;
    }
    if (payload[0] != SEG_CRM) return false;
    ByteSpanIn in{payload + 1, payload_len - 1, 0};
    NovRangeDec rd(in);
    int freqs[4];
    for (size_t r = 0; r < len; r++) {
        model.predict(dst, static_cast<int64_t>(dst_start), static_cast<int64_t>(r), freqs);
        int target = rd.get_threshold();
        int cum = 0, sym = 0;
        while (sym < 3 && cum + freqs[sym] <= target) { cum += freqs[sym]; sym++; }
        rd.decode_symbol(static_cast<uint32_t>(cum), static_cast<uint32_t>(freqs[sym]));
        dst[dst_start + r] = static_cast<uint8_t>(sym);
        model.update(dst, static_cast<int64_t>(dst_start), static_cast<int64_t>(r), sym);
    }
    return true;
}

} // namespace

// ============================================================================
// Extern "C" Public FFI Implementation
// ============================================================================

extern "C" {

void hgs_cpp_free_buffer(HgsOwnedBytes buf) {
    std::free(buf.data);
}

int hgs_cpp_compute_bottom_sketch(
    const uint8_t* two_bit_bases, size_t len, size_t sketch_size, uint64_t** out_sketch, size_t* out_len) {
    *out_sketch = nullptr;
    *out_len = 0;
    if (sketch_size == 0 || len < static_cast<size_t>(K_PRIMARY)) return 0;

    Packed2Bit seq = Packed2Bit::from_bytes(two_bit_bases, len);
    std::vector<uint64_t> heap;
    heap.reserve(sketch_size);

    struct Ctx {
        std::vector<uint64_t>* heap;
        size_t cap;
    } ctx{&heap, sketch_size};

    auto cb = [](uint64_t /*pos*/, uint64_t canon_k, uint8_t /*st*/, void* v) {
        auto* c = static_cast<Ctx*>(v);
        uint64_t h = hgs_mix64(canon_k); // Unsigned uint64_t comparison matches Java (h ^ Long.MIN_VALUE) signed ordering!
        if (c->heap->size() < c->cap) {
            c->heap->push_back(h);
            if (c->heap->size() == c->cap) std::make_heap(c->heap->begin(), c->heap->end());
        } else if (h < c->heap->front()) {
            std::pop_heap(c->heap->begin(), c->heap->end());
            c->heap->back() = h;
            std::push_heap(c->heap->begin(), c->heap->end());
        }
    };
    hgs_scan_open_syncmers(seq.words.data(), 0, len, K_PRIMARY, S_PRIMARY, T_PRIMARY, cb, &ctx);
    std::sort(heap.begin(), heap.end());
    heap.erase(std::unique(heap.begin(), heap.end()), heap.end());
    if (heap.empty()) return 0;

    uint64_t* arr = static_cast<uint64_t*>(std::malloc(heap.size() * sizeof(uint64_t)));
    if (!arr) return -1;
    std::memcpy(arr, heap.data(), heap.size() * sizeof(uint64_t));
    *out_sketch = arr;
    *out_len = heap.size();
    return 0;
}

void hgs_cpp_free_u64(uint64_t* ptr) {
    std::free(ptr);
}

double hgs_cpp_sketch_jaccard(const uint64_t* a, size_t a_len, const uint64_t* b, size_t b_len) {
    if (!a_len || !b_len) return 0.0;
    size_t i = 0, j = 0, inter = 0, uni = 0;
    size_t lim = std::min(a_len, b_len);
    while (i < a_len && j < b_len && uni < lim) {
        uint64_t va = a[i], vb = b[j];
        if (va == vb) { inter++; i++; j++; }
        else if (va < vb) { i++; }
        else { j++; }
        uni++;
    }
    return uni == 0 ? 0.0 : static_cast<double>(inter) / static_cast<double>(uni);
}

HgsSyncmerIndex64* hgs_cpp_index_build(
    const uint8_t* const* donor_two_bit_ptrs,
    const uint64_t* donor_lens,
    size_t num_donors,
    int subsample_shift,
    int build_secondary,
    int /*threads*/) {
    auto idx = std::make_unique<HgsSyncmerIndex64>();
    idx->num_donors = num_donors;
    idx->donor_offsets.resize(num_donors + 1);
    uint64_t total_bases = 0;
    for (size_t d = 0; d < num_donors; d++) {
        idx->donor_offsets[d] = total_bases;
        total_bases += donor_lens[d];
    }
    idx->donor_offsets[num_donors] = total_bases;

    // Concatenate into unified 64-bit packed word bank directly
    idx->bank.length = total_bases;
    size_t total_words = static_cast<size_t>((total_bases + 31ULL) >> 5);
    idx->bank.words.assign(total_words + 2, 0ULL);

    uint64_t dst_bit = 0;
    for (size_t d = 0; d < num_donors; d++) {
        uint64_t dl = donor_lens[d];
        const uint8_t* ptr = donor_two_bit_ptrs[d];
        for (uint64_t k = 0; k < dl; k++) {
            uint64_t b = ptr[k] & 3U;
            size_t wi = static_cast<size_t>(dst_bit >> 6);
            uint32_t sh = static_cast<uint32_t>(dst_bit & 63ULL);
            idx->bank.words[wi] |= (b << sh);
            dst_bit += 2;
        }
    }

    size_t num_blocks = static_cast<size_t>((total_bases >> BLOCK_SHIFT) + 2);
    idx->pos_to_donor_floor.resize(num_blocks);
    size_t cur_d = 0;
    for (size_t b = 0; b < num_blocks; b++) {
        uint64_t bs = static_cast<uint64_t>(b) << BLOCK_SHIFT;
        while (cur_d + 1 < num_donors && idx->donor_offsets[cur_d + 1] <= bs) cur_d++;
        idx->pos_to_donor_floor[b] = cur_d;
    }

    idx->primary_subsample_shift = subsample_shift;
    idx->primary_subsample_mask = (1U << subsample_shift) - 1U;
    uint64_t est_prim = std::max<uint64_t>(1024, (total_bases / W_PRIMARY) >> subsample_shift);
    idx->bucket_bits = choose_bucket_bits(est_prim, 12, subsample_shift == 0 ? 24 : std::max(16, 20 - subsample_shift));
    idx->bucket_mask = (1U << idx->bucket_bits) - 1U;

    build_csr_pass(
        idx->bank, idx->donor_offsets, K_PRIMARY, S_PRIMARY, T_PRIMARY,
        idx->bucket_bits, subsample_shift, idx->bucket_offsets, idx->pos_and_strand, idx->sig16);

    if (build_secondary && subsample_shift == 0 && total_bases > 0 && total_bases <= 80'000'000ULL) {
        uint64_t est_sec = std::max<uint64_t>(1024, total_bases / (W_SECONDARY * 2));
        idx->sec_bucket_bits = choose_bucket_bits(est_sec, 12, 19);
        idx->sec_bucket_mask = (1U << idx->sec_bucket_bits) - 1U;
        build_csr_pass(
            idx->bank, idx->donor_offsets, K_SECONDARY, S_SECONDARY, T_SECONDARY,
            idx->sec_bucket_bits, 1, idx->sec_bucket_offsets, idx->sec_pos_and_strand, idx->sec_sig16);
    }

    return idx.release();
}

void hgs_cpp_index_free(HgsSyncmerIndex64* idx) {
    delete idx;
}

void hgs_cpp_free_mosaic_result(HgsMosaicResult* res) {
    if (!res) return;
    std::free(res->mosaic_data);
    std::free(res->novelty_data);
    std::memset(res, 0, sizeof(*res));
}

HgsMosaicResult hgs_cpp_encode_mosaic(
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
    int /*threads*/) {
    HgsMosaicResult res{};
    MosaicConfigC cfg{
        panel_idx != nullptr && max_donor_exclusive > 0,
        allow_approx != 0, allow_splice != 0, allow_kimura != 0,
        allow_wobble != 0, allow_self != 0, allow_tandem != 0};

    ByteVecOut hdr_body;
    hdr_body.write_byte(cfg.to_header_byte());
    hdr_body.write_var_u64(tgt_len);
    if (tgt_len == 0) {
        hdr_body.write_var_u64(0);
        res.mosaic_len = hdr_body.buf.size();
        res.mosaic_data = static_cast<uint8_t*>(std::malloc(res.mosaic_len));
        if (!res.mosaic_data) { res.status = -1; return res; }
        std::memcpy(res.mosaic_data, hdr_body.buf.data(), res.mosaic_len);
        res.status = 0;
        return res;
    }

    Packed2Bit target = Packed2Bit::from_bytes(target_two_bit, tgt_len);

    // Precompute primary (and optional secondary) target syncmers
    struct TgtSync {
        std::vector<uint64_t> pos;
        std::vector<uint32_t> canon;
        std::vector<uint8_t> strand;
    } ts16, ts13;

    auto extract_cb = [](uint64_t p, uint64_t ck, uint8_t st, void* v) {
        auto* ts = static_cast<TgtSync*>(v);
        ts->pos.push_back(p);
        ts->canon.push_back(static_cast<uint32_t>(ck));
        ts->strand.push_back(st);
    };
    ts16.pos.reserve(static_cast<size_t>(tgt_len / 6 + 32));
    hgs_scan_open_syncmers(target.words.data(), 0, tgt_len, K_PRIMARY, S_PRIMARY, T_PRIMARY, extract_cb, &ts16);
    if (cfg.has_panel && panel_idx->has_secondary()) {
        hgs_scan_open_syncmers(target.words.data(), 0, tgt_len, K_SECONDARY, S_SECONDARY, T_SECONDARY, extract_cb, &ts13);
    }

    std::vector<int> donor_rank;
    if (cfg.has_panel) {
        donor_rank.assign(max_donor_exclusive, static_cast<int>(max_donor_exclusive));
        int rk = 0;
        for (size_t i = 0; i < sim_order_len; i++) {
            int32_t d = donor_sim_order[i];
            if (d >= 0 && static_cast<size_t>(d) < max_donor_exclusive && donor_rank[d] == static_cast<int>(max_donor_exclusive)) {
                donor_rank[d] = rk++;
            }
        }
        if (preferred_twin_donor >= 0 && static_cast<size_t>(preferred_twin_donor) < max_donor_exclusive) {
            donor_rank[static_cast<size_t>(preferred_twin_donor)] = -1;
        }
    }

    // Native C++20 uses higher 2^20 4-way self-repeat cache (4M slots) in Mode S!
    int self_bits = (!cfg.has_panel && tgt_len >= 1'000'000ULL) ? 19 : 15;
    std::unique_ptr<SelfSyncmerTable64> self_table = cfg.allow_self ? std::make_unique<SelfSyncmerTable64>(self_bits) : nullptr;

    std::vector<int> mtf;
    init_mtf(mtf, cfg.has_panel ? max_donor_exclusive : 0);
    std::vector<int64_t> last_ref_end(cfg.has_panel ? max_donor_exclusive * 2 : 0, -1);
    std::vector<int64_t> last_tgt_end(cfg.has_panel ? max_donor_exclusive * 2 : 0, -1);
    if (cfg.has_panel && preferred_twin_donor >= 0 && static_cast<size_t>(preferred_twin_donor) < max_donor_exclusive) {
        promote_mtf(mtf, static_cast<int>(preferred_twin_donor << 1));
    }

    std::vector<ParsedToken> tokens;
    std::vector<uint8_t> nov_buf;
    nov_buf.reserve(static_cast<size_t>(std::min<uint64_t>(tgt_len, 65536)));

    size_t sync16_cur = 0, sync13_cur = 0, self_ins_cur = 0;
    uint64_t pos = 0, lit_start = 0;

    auto append_nov = [&](uint64_t s, uint64_t count) {
        for (uint64_t k = 0; k < count; k++) nov_buf.push_back(target.base_at(s + k));
    };

    auto flush_lit = [&](uint64_t up_to) {
        if (up_to > lit_start) {
            uint64_t l = up_to - lit_start;
            append_nov(lit_start, l);
            ParsedToken lt; lt.kind = 0; lt.lit_len = l; lt.total_tgt_span = l;
            if (!tokens.empty() && tokens.back().kind == 0) {
                tokens.back().lit_len += l;
                tokens.back().total_tgt_span += l;
            } else {
                tokens.push_back(std::move(lt));
                res.literal_runs++;
            }
            lit_start = up_to;
        }
    };

    auto eval_chain_at = [&](uint64_t qpos, ParsedToken& out_best, double& out_gain) -> bool {
        out_gain = -1e18;
        bool has_best = false;
        if (qpos + MIN_TANDEM_LEN > tgt_len) return false;

        // 1. Tandem Periodicity Check (P in 1..64)
        if (cfg.allow_tandem && qpos >= 1) {
            uint64_t max_p = std::min<uint64_t>(MAX_TANDEM_PERIOD, qpos);
            uint64_t best_p = 0, best_tl = 0;
            uint32_t t16 = static_cast<uint32_t>(target.word32_fwd(qpos));
            for (uint64_t p = 1; p <= max_p; p++) {
                if (static_cast<uint32_t>(target.word32_fwd(qpos - p)) != t16) continue;
                uint64_t tl = hgs_count_exact_match_2bit(
                    target.words.data(), 0, tgt_len, 0, qpos - p,
                    target.words.data(), 0, tgt_len, qpos, tgt_len - qpos);
                if (tl >= MIN_TANDEM_LEN && tl > best_tl) {
                    best_tl = tl; best_p = p;
                    if (tl >= 4096) break;
                }
            }
            uint64_t min_tan = cfg.has_panel ? MIN_TANDEM_LEN : 36ULL;
            if (best_tl >= min_tan) {
                double gain = static_cast<double>(best_tl) * 1.865 - (11.0 + approx_count_bits(best_tl - MIN_TANDEM_LEN));
                if (gain > 6.5) {
                    out_best = ParsedToken{};
                    out_best.kind = 3; out_best.donor_id = best_p;
                    out_best.start_ref_pos = qpos - best_p;
                    out_best.total_tgt_span = best_tl; out_best.total_ref_span = best_tl;
                    Subblock sb; sb.len = best_tl;
                    out_best.subblocks.push_back(std::move(sb));
                    out_gain = gain;
                    has_best = true;
                    if (best_tl >= 512) return true;
                }
            }
        }

        // Advance cursor views locally without mutating shared cursors
        size_t s16 = sync16_cur;
        while (s16 < ts16.pos.size() && ts16.pos[s16] < qpos) s16++;
        size_t s13 = sync13_cur;
        while (s13 < ts13.pos.size() && ts13.pos[s13] < qpos) s13++;

        // 2. Panel Donors (Primary k=16 syncmer probe)
        if (cfg.has_panel && s16 < ts16.pos.size() && ts16.pos[s16] <= qpos + 16) {
            uint64_t sync_off = ts16.pos[s16] - qpos;
            uint64_t sp = qpos + sync_off;
            if (sp + K_PRIMARY <= tgt_len) {
                uint64_t qcanon = ts16.canon[s16];
                if (panel_idx->accepts_primary(qcanon)) {
                    uint32_t bkt = panel_idx->bucket_for_primary(qcanon);
                    uint16_t qsig = panel_idx->sig_for_primary(qcanon);
                    uint64_t bs = panel_idx->bucket_offsets[bkt], be = panel_idx->bucket_offsets[bkt + 1];
                    struct Cand { uint64_t gpos; uint8_t cst; int rank; };
                    std::vector<Cand> cands;
                    for (uint64_t b = bs; b < be && cands.size() < MAX_BUCKET_PROBES; b++) {
                        if (panel_idx->sig16[b] != qsig) continue;
                        uint64_t pk = panel_idx->pos_and_strand[b];
                        uint64_t gp = pk & ((1ULL << 63) - 1ULL);
                        size_t did = panel_idx->donor_of_global_pos(gp);
                        if (did >= max_donor_exclusive) break;
                        int rk = donor_rank[did];
                        uint8_t cst = static_cast<uint8_t>(pk >> 63);
                        if (((did << 1) == static_cast<size_t>(mtf[0])) || (((did << 1) | 1) == static_cast<size_t>(mtf[0]))) rk = -2;
                        cands.push_back({gp, cst, rk});
                    }
                    std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) { return a.rank < b.rank; });
                    for (const auto& c : cands) {
                        size_t did = panel_idx->donor_of_global_pos(c.gpos);
                        uint64_t dstart = panel_idx->donor_offsets[did];
                        uint64_t dlen = panel_idx->donor_offsets[did + 1] - dstart;
                        uint64_t rel_fwd = c.gpos - dstart;
                        uint8_t st = ts16.strand[s16] ^ c.cst;
                        uint64_t ref_sp = (st == 0) ? rel_fwd : (dlen - K_PRIMARY - rel_fwd);
                        if (ref_sp < sync_off) continue;
                        uint64_t crp = ref_sp - sync_off;
                        if (crp + MIN_MATCH_LEN > dlen) continue;
                        if (hgs_count_exact_match_2bit(panel_idx->bank.words.data(), dstart, dlen, st, ref_sp,
                                                       target.words.data(), 0, tgt_len, sp, K_PRIMARY) < K_PRIMARY) continue;
                        if (sync_off > 0 && panel_idx->bank.base_at_strand(dstart, dlen, st, crp) != target.base_at(qpos)) continue;

                        size_t ds = (did << 1) | st;
                        bool clonal = false;
                        if (cfg.allow_approx) {
                            uint64_t ex = hgs_count_exact_match_2bit(
                                panel_idx->bank.words.data(), dstart, dlen, st, crp,
                                target.words.data(), 0, tgt_len, qpos, 160);
                            if (ex >= 128) clonal = true;
                        }
                        double hb = estimate_panel_hdr_bits(did, st, crp, qpos, mtf, last_ref_end[ds], last_tgt_end[ds], max_donor_exclusive);
                        ParsedToken cand_tok; double cg = 0;
                        if (build_chain(1, did, st, crp, panel_idx->bank, dstart, dlen, target, qpos, tgt_len, cfg, clonal, hb, cand_tok, cg)) {
                            if (!has_best || cg > out_gain) { out_best = std::move(cand_tok); out_gain = cg; has_best = true; }
                            if (out_best.total_tgt_span >= 4096) return true;
                        }
                    }
                }
            }
        }

        // 3. Secondary k=13 diverged island lookup (when panel has secondary & primary missed large match)
        if (cfg.has_panel && panel_idx->has_secondary() && (!has_best || out_best.total_tgt_span < 128)
            && s13 < ts13.pos.size() && ts13.pos[s13] <= qpos + 16) {
            uint64_t sync_off = ts13.pos[s13] - qpos;
            uint64_t sp = qpos + sync_off;
            uint64_t qcanon = ts13.canon[s13];
            int64_t sbkt = panel_idx->sec_bucket_for(qcanon);
            if (sbkt >= 0 && sp + K_SECONDARY <= tgt_len) {
                uint16_t qsig = panel_idx->sec_sig_for(qcanon);
                uint64_t bs = panel_idx->sec_bucket_offsets[sbkt], be = panel_idx->sec_bucket_offsets[sbkt + 1];
                int probes = 0;
                for (uint64_t b = bs; b < be && probes < 24; b++) {
                    if (panel_idx->sec_sig16[b] != qsig) continue;
                    uint64_t pk = panel_idx->sec_pos_and_strand[b];
                    uint64_t gp = pk & ((1ULL << 63) - 1ULL);
                    size_t did = panel_idx->donor_of_global_pos(gp);
                    if (did >= max_donor_exclusive) break;
                    probes++;
                    uint64_t dstart = panel_idx->donor_offsets[did];
                    uint64_t dlen = panel_idx->donor_offsets[did + 1] - dstart;
                    uint64_t rel_fwd = gp - dstart;
                    uint8_t st = ts13.strand[s13] ^ static_cast<uint8_t>(pk >> 63);
                    uint64_t ref_sp = (st == 0) ? rel_fwd : (dlen - K_SECONDARY - rel_fwd);
                    if (ref_sp < sync_off) continue;
                    uint64_t crp = ref_sp - sync_off;
                    if (crp + MIN_MATCH_LEN > dlen) continue;
                    if (hgs_count_exact_match_2bit(panel_idx->bank.words.data(), dstart, dlen, st, ref_sp,
                                                   target.words.data(), 0, tgt_len, sp, K_SECONDARY) < K_SECONDARY) continue;
                    size_t ds = (did << 1) | st;
                    double hb = estimate_panel_hdr_bits(did, st, crp, qpos, mtf, last_ref_end[ds], last_tgt_end[ds], max_donor_exclusive);
                    ParsedToken cand_tok; double cg = 0;
                    if (build_chain(1, did, st, crp, panel_idx->bank, dstart, dlen, target, qpos, tgt_len, cfg, false, hb, cand_tok, cg)) {
                        if (!has_best || cg > out_gain) { out_best = std::move(cand_tok); out_gain = cg; has_best = true; }
                    }
                }
            }
        }

        // 4. Bidirectional Self-Repeat Match
        if (self_table != nullptr && (!has_best || out_best.total_tgt_span < 256)
            && s16 < ts16.pos.size() && ts16.pos[s16] <= qpos + 16) {
            uint64_t sync_off = ts16.pos[s16] - qpos;
            uint64_t sp = qpos + sync_off;
            uint64_t qcanon = ts16.canon[s16];
            uint32_t h = hgs_mix32(qcanon);
            size_t base = static_cast<size_t>(h & self_table->set_mask) << 2;
            uint16_t sig = static_cast<uint16_t>((qcanon >> 16) ^ (h >> 16));
            for (size_t w = 0; w < SelfSyncmerTable64::WAYS; w++) {
                uint64_t ps = self_table->pos_strand[base + w];
                if (ps == std::numeric_limits<uint64_t>::max() || self_table->sigs[base + w] != sig) continue;
                uint64_t ref_fwd = ps & ((1ULL << 63) - 1ULL);
                if (ref_fwd + K_PRIMARY > qpos) continue;
                uint8_t st = ts16.strand[s16] ^ static_cast<uint8_t>(ps >> 63);
                uint64_t rsp = (st == 0) ? ref_fwd : (qpos - K_PRIMARY - ref_fwd);
                if (rsp < sync_off) continue;
                uint64_t crp = rsp - sync_off;
                if (crp + MIN_MATCH_LEN > qpos) continue;
                if (hgs_count_exact_match_2bit(target.words.data(), 0, qpos, st, rsp,
                                               target.words.data(), 0, tgt_len, sp, K_PRIMARY) < K_PRIMARY) continue;
                double hb = 5.0 + approx_count_bits((st == 0) ? (qpos - 1 - crp) : crp);
                ParsedToken cand_tok; double cg = 0;
                if (build_chain(2, 0, st, crp, target, 0, qpos, target, qpos, tgt_len, cfg, false, hb, cand_tok, cg)) {
                    if (!has_best || cg > out_gain) { out_best = std::move(cand_tok); out_gain = cg; has_best = true; }
                }
            }
        }
        return has_best;
    };

    while (pos < tgt_len) {
        if (self_table != nullptr) {
            while (self_ins_cur < ts16.pos.size() && ts16.pos[self_ins_cur] + K_PRIMARY <= pos) {
                self_table->insert(ts16.canon[self_ins_cur], ts16.pos[self_ins_cur], ts16.strand[self_ins_cur]);
                self_ins_cur++;
            }
        }
        while (sync16_cur < ts16.pos.size() && ts16.pos[sync16_cur] < pos) sync16_cur++;
        while (sync13_cur < ts13.pos.size() && ts13.pos[sync13_cur] < pos) sync13_cur++;

        ParsedToken best;
        double best_gain = -1e18;
        bool found = eval_chain_at(pos, best, best_gain);
        if (found && best.total_tgt_span < 256 && pos + 8 + MIN_MATCH_LEN < tgt_len) {
            for (size_t i = sync16_cur; i < ts16.pos.size() && ts16.pos[i] <= pos + 7; i++) {
                if (ts16.pos[i] > pos) {
                    ParsedToken lazy; double lg = -1e18;
                    if (eval_chain_at(ts16.pos[i], lazy, lg) && lg > best_gain + 14.0) {
                        found = false;
                    }
                    break;
                }
            }
        }

        if (!found) {
            uint64_t next_pos = pos + 1;
            if (cfg.has_panel && sync16_cur < ts16.pos.size() && ts16.pos[sync16_cur] > pos + 16) {
                uint64_t jump = ts16.pos[sync16_cur] - 8;
                if (jump > next_pos && jump < tgt_len) next_pos = jump;
            }
            pos = next_pos;
            continue;
        }

        flush_lit(pos);

        if (best.kind == 1) {
            res.panel_tokens++;
            res.panel_match_bases += best.total_tgt_span;
            size_t ds = static_cast<size_t>((best.donor_id << 1) | best.strand);
            promote_mtf(mtf, static_cast<int>(ds));
            last_ref_end[ds] = static_cast<int64_t>(best.start_ref_pos + best.total_ref_span);
            last_tgt_end[ds] = static_cast<int64_t>(pos + best.total_tgt_span);
        } else if (best.kind == 2) {
            res.self_tokens++;
            res.self_match_bases += best.total_tgt_span;
        } else if (best.kind == 3) {
            res.tandem_tokens++;
            res.tandem_match_bases += best.total_tgt_span;
        }

        uint64_t cpos = pos;
        for (size_t b = 0; b < best.subblocks.size(); b++) {
            const auto& sb = best.subblocks[b];
            res.substitutions += sb.snp_offsets.size();
            cpos += sb.len;
            if (b + 1 < best.subblocks.size()) {
                res.indel_splices++;
                if (sb.d_tgt_after > 0) {
                    append_nov(cpos, sb.d_tgt_after);
                    res.panel_match_bases -= (best.kind == 1 ? sb.d_tgt_after : 0);
                    res.self_match_bases -= (best.kind == 2 ? sb.d_tgt_after : 0);
                    cpos += sb.d_tgt_after;
                }
            }
        }
        pos += best.total_tgt_span;
        lit_start = pos;
        tokens.push_back(std::move(best));
    }
    flush_lit(tgt_len);

    // Entropy encode token stream into exact .hgs2 format
    uint64_t novelty_len = nov_buf.size();
    hdr_body.write_var_u64(novelty_len);
    BitRangeEnc enc(hdr_body);
    MosaicModels models;
    std::vector<int> emtf;
    init_mtf(emtf, cfg.has_panel ? max_donor_exclusive : 0);
    std::vector<int64_t> elr(cfg.has_panel ? max_donor_exclusive * 2 : 0, -1);
    std::vector<int64_t> elt(cfg.has_panel ? max_donor_exclusive * 2 : 0, -1);
    bool prev_lit = false, prev_pan = true;
    uint64_t cur_pos = 0;

    for (const auto& tok : tokens) {
        if (tok.kind == 0) {
            enc.encode_bit(models.is_literal.data(), prev_pan ? 1 : 2, 1);
            models.write_count(enc, models.lit_len_slots[0], 0, static_cast<uint32_t>(tok.lit_len - 1));
            cur_pos += tok.lit_len;
            prev_lit = true;
            continue;
        }
        bool after_lit = prev_lit;
        if (!prev_lit) enc.encode_bit(models.is_literal.data(), prev_pan ? 1 : 2, 0);
        prev_lit = false;
        bool is_pan = (tok.kind == 1);
        if (cfg.has_panel && (cfg.allow_self || cfg.allow_tandem)) {
            enc.encode_bit(models.is_panel_match.data(), prev_pan ? 1 : 0, is_pan ? 1 : 0);
        }
        prev_pan = is_pan;

        if (is_pan) {
            int ds = static_cast<int>((tok.donor_id << 1) | tok.strand);
            int is_m0 = (ds == emtf[0]) ? 1 : 0;
            encode_mtf_ds(enc, models, emtf, max_donor_exclusive, ds, after_lit ? 1 : 0);
            uint64_t dstart = panel_idx->donor_offsets[tok.donor_id];
            uint64_t dlen = panel_idx->donor_offsets[tok.donor_id + 1] - dstart;
            encode_panel_ref_pos(enc, models, tok.start_ref_pos, cur_pos, dlen, elr[ds], elt[ds], is_m0);
            encode_subblocks(enc, models, cfg, panel_idx->bank, dstart, dlen, tok.strand, tok.start_ref_pos, target, cur_pos, tok.subblocks);
            elr[ds] = static_cast<int64_t>(tok.start_ref_pos + tok.total_ref_span);
            elt[ds] = static_cast<int64_t>(cur_pos + tok.total_tgt_span);
            cur_pos += tok.total_tgt_span;
        } else if (tok.kind == 3) {
            if (cfg.allow_self) enc.encode_bit(models.is_tandem.data(), 0, 1);
            models.write_fixed_tree(enc, models.tandem_period_tree, static_cast<int>(tok.donor_id - 1), 6);
            models.write_count(enc, models.match_len_slots[0], 1, static_cast<uint32_t>(tok.total_tgt_span - MIN_TANDEM_LEN));
            cur_pos += tok.total_tgt_span;
        } else {
            if (cfg.allow_tandem) enc.encode_bit(models.is_tandem.data(), 0, 0);
            enc.encode_bit(models.self_strand.data(), 0, tok.strand);
            uint32_t dist = static_cast<uint32_t>((tok.strand == 0) ? (cur_pos - 1 - tok.start_ref_pos) : tok.start_ref_pos);
            models.write_count(enc, models.self_dist_slots[tok.strand], 2, dist);
            encode_subblocks(enc, models, cfg, target, 0, cur_pos, tok.strand, tok.start_ref_pos, target, cur_pos, tok.subblocks);
            cur_pos += tok.total_tgt_span;
        }
    }
    enc.finish();

    res.mosaic_len = hdr_body.buf.size();
    res.mosaic_data = static_cast<uint8_t*>(std::malloc(res.mosaic_len));
    if (!res.mosaic_data) { res.status = -1; return res; }
    std::memcpy(res.mosaic_data, hdr_body.buf.data(), res.mosaic_len);

    res.novelty_len = nov_buf.size();
    if (res.novelty_len > 0) {
        res.novelty_data = static_cast<uint8_t*>(std::malloc(res.novelty_len));
        if (!res.novelty_data) { std::free(res.mosaic_data); res.mosaic_data = nullptr; res.status = -1; return res; }
        std::memcpy(res.novelty_data, nov_buf.data(), res.novelty_len);
    }
    res.status = 0;
    return res;
}

int hgs_cpp_decode_mosaic(
    const uint8_t* mosaic_bytes,
    size_t mosaic_len,
    const uint8_t* novelty_bytes,
    size_t novelty_len,
    const uint8_t* const* donor_2bit_ptrs,
    const uint64_t* donor_lens,
    size_t max_donor_exclusive,
    uint8_t** out_two_bit,
    uint64_t* out_len) {
    *out_two_bit = nullptr;
    *out_len = 0;

    ByteSpanIn in{mosaic_bytes, mosaic_len, 0};
    bool ok = true;
    uint8_t flags = in.read_byte(ok);
    if (!ok) return -1;
    MosaicConfigC cfg = MosaicConfigC::from_header_byte(flags);
    uint64_t target_len = in.read_var_u64(ok);
    uint64_t nov_hdr_len = in.read_var_u64(ok);
    if (!ok || nov_hdr_len > novelty_len) return -1;
    *out_len = target_len;
    if (target_len == 0) return 0;
    if (cfg.has_panel && max_donor_exclusive == 0) return -1;

    uint8_t* out = static_cast<uint8_t*>(std::malloc(static_cast<size_t>(target_len)));
    if (!out) return -1;

    // Lazy packed wrapper only for referenced donors
    std::vector<Packed2Bit> donor_banks(cfg.has_panel ? max_donor_exclusive : 0);
    std::vector<uint8_t> donor_ready(cfg.has_panel ? max_donor_exclusive : 0, 0);

    BitRangeDec dec(in);
    MosaicModels models;
    std::vector<int> mtf;
    init_mtf(mtf, cfg.has_panel ? max_donor_exclusive : 0);
    std::vector<int64_t> last_ref(cfg.has_panel ? max_donor_exclusive * 2 : 0, -1);
    std::vector<int64_t> last_tgt(cfg.has_panel ? max_donor_exclusive * 2 : 0, -1);

    bool prev_lit = false, prev_pan = true;
    uint64_t pos = 0, nov_cur = 0;

    while (pos < target_len) {
        bool is_lit = (!prev_lit) && (dec.decode_bit(models.is_literal.data(), prev_pan ? 1 : 2) == 1);
        if (is_lit) {
            uint64_t lit_len = models.read_count(dec, models.lit_len_slots[0], 0) + 1ULL;
            if (pos + lit_len > target_len || nov_cur + lit_len > nov_hdr_len) { std::free(out); return -1; }
            std::memcpy(out + pos, novelty_bytes + nov_cur, static_cast<size_t>(lit_len));
            nov_cur += lit_len;
            pos += lit_len;
            prev_lit = true;
            continue;
        }
        bool after_lit = prev_lit;
        prev_lit = false;
        bool is_pan = cfg.has_panel && ((cfg.allow_self || cfg.allow_tandem)
            ? (dec.decode_bit(models.is_panel_match.data(), prev_pan ? 1 : 0) == 1) : true);
        prev_pan = is_pan;

        if (is_pan) {
            int prev_m0 = mtf[0];
            int ds = decode_mtf_ds(dec, models, mtf, max_donor_exclusive, after_lit ? 1 : 0);
            int is_m0 = (ds == prev_m0) ? 1 : 0;
            size_t did = static_cast<size_t>(ds >> 1);
            uint8_t st = static_cast<uint8_t>(ds & 1);
            if (did >= max_donor_exclusive) { std::free(out); return -1; }
            if (!donor_ready[did]) {
                donor_banks[did] = Packed2Bit::from_bytes(donor_2bit_ptrs[did], donor_lens[did]);
                donor_ready[did] = 1;
            }
            uint64_t srp = decode_panel_ref_pos(dec, models, pos, donor_lens[did], last_ref[ds], last_tgt[ds], is_m0, ok);
            if (!ok) { std::free(out); return -1; }
            uint64_t er = 0, et = 0;
            if (!decode_subblocks_donor(dec, models, cfg, &donor_banks[did], 0, donor_lens[did], false, st, srp,
                                        out, pos, target_len, novelty_bytes, nov_cur, nov_hdr_len, er, et)) {
                std::free(out); return -1;
            }
            last_ref[ds] = static_cast<int64_t>(er);
            last_tgt[ds] = static_cast<int64_t>(et);
            pos = et;
        } else {
            bool is_tan = (cfg.allow_self && cfg.allow_tandem)
                ? (dec.decode_bit(models.is_tandem.data(), 0) == 1) : cfg.allow_tandem;
            if (is_tan) {
                uint64_t period = static_cast<uint64_t>(models.read_fixed_tree(dec, models.tandem_period_tree, 6)) + 1ULL;
                uint64_t tlen = models.read_count(dec, models.match_len_slots[0], 1) + MIN_TANDEM_LEN;
                if (period > pos || pos + tlen > target_len) { std::free(out); return -1; }
                uint64_t src = pos - period;
                for (uint64_t k = 0; k < tlen; k++) out[pos + k] = out[src + k];
                pos += tlen;
            } else {
                uint8_t st = static_cast<uint8_t>(dec.decode_bit(models.self_strand.data(), 0));
                uint64_t dc = models.read_count(dec, models.self_dist_slots[st], 2);
                uint64_t srp = (st == 0) ? (pos - 1ULL - dc) : dc;
                if (srp >= pos) { std::free(out); return -1; }
                uint64_t er = 0, et = 0;
                if (!decode_subblocks_donor(dec, models, cfg, nullptr, 0, pos, true, st, srp,
                                            out, pos, target_len, novelty_bytes, nov_cur, nov_hdr_len, er, et)) {
                    std::free(out); return -1;
                }
                pos = et;
            }
        }
    }
    if (nov_cur != nov_hdr_len) { std::free(out); return -1; }
    *out_two_bit = out;
    return 0;
}

int hgs_cpp_mosaic_peek_lengths(const uint8_t* mosaic_bytes, size_t mosaic_len, uint64_t* out_tgt_len, uint64_t* out_nov_len) {
    ByteSpanIn in{mosaic_bytes, mosaic_len, 0};
    bool ok = true;
    (void)in.read_byte(ok);
    *out_tgt_len = in.read_var_u64(ok);
    *out_nov_len = in.read_var_u64(ok);
    return ok ? 0 : -1;
}

int hgs_cpp_novelty_encode(const uint8_t* two_bit_bases, uint64_t len, int threads, int bounded_mem, HgsOwnedBytes* out_payload) {
    out_payload->data = nullptr;
    out_payload->len = 0;
    if (len == 0) {
        uint8_t* b = static_cast<uint8_t*>(std::malloc(1));
        if (!b) return -1;
        b[0] = 0;
        out_payload->data = b; out_payload->len = 1;
        return 0;
    }

    int th = std::max(1, threads);
    uint64_t turbo_cutoff = bounded_mem ? 8'000'000ULL : 64'000'000ULL;
    if (len >= turbo_cutoff) {
        size_t total_segs = static_cast<size_t>((len + TURBO_SEG_BASES - 1ULL) / TURBO_SEG_BASES);
        std::vector<size_t> seg_starts(total_segs), seg_lens(total_segs);
        for (size_t i = 0; i < total_segs; i++) {
            seg_starts[i] = static_cast<size_t>(i * TURBO_SEG_BASES);
            seg_lens[i] = static_cast<size_t>(std::min<uint64_t>(TURBO_SEG_BASES, len - i * TURBO_SEG_BASES));
        }
        std::vector<std::vector<uint8_t>> encoded(total_segs);
        int workers = std::max(1, std::min<int>(th, static_cast<int>(std::min<size_t>(total_segs, bounded_mem ? 4 : 16))));
        if (workers <= 1) {
            for (size_t i = 0; i < total_segs; i++) encode_turbo_segment(two_bit_bases + seg_starts[i], seg_lens[i], encoded[i]);
        } else {
            std::vector<std::thread> pool;
            std::atomic<size_t> next_i{0};
            for (int w = 0; w < workers; w++) {
                pool.emplace_back([&]() {
                    while (true) {
                        size_t i = next_i.fetch_add(1);
                        if (i >= total_segs) break;
                        encode_turbo_segment(two_bit_bases + seg_starts[i], seg_lens[i], encoded[i]);
                    }
                });
            }
            for (auto& t : pool) t.join();
        }
        ByteVecOut out;
        out.write_var_u64(total_segs);
        out.write_byte(static_cast<uint8_t>((bounded_mem ? 1U : 0U) | FLAG_TURBO_NOVELTY));
        for (size_t i = 0; i < total_segs; i++) out.write_var_u64(seg_lens[i]);
        for (size_t i = 0; i < total_segs; i++) out.write_var_u64(encoded[i].size());
        for (size_t i = 0; i < total_segs; i++) out.write_bytes(encoded[i].data(), encoded[i].size());
        out_payload->len = out.buf.size();
        out_payload->data = static_cast<uint8_t*>(std::malloc(out_payload->len));
        if (!out_payload->data) return -1;
        std::memcpy(out_payload->data, out.buf.data(), out_payload->len);
        return 0;
    }

    bool pan_super = len > PAN_SUPER_THRESHOLD;
    uint8_t model_flags = static_cast<uint8_t>((bounded_mem ? 1U : 0U) | (pan_super ? 2U : 0U));
    size_t warm_bases = pan_super ? 131072ULL : (len >= 1'500'000ULL ? 393216ULL : 65536ULL);
    size_t min_tail = pan_super ? 524288ULL : (len >= 1'500'000ULL ? 2'000'000ULL : 131072ULL);
    size_t tail = (len > warm_bases) ? static_cast<size_t>(len - warm_bases) : 0;
    int max_workers = pan_super ? std::min(th, 8) : th;
    int worker_cnt = (max_workers <= 1 || tail < min_tail) ? 0 : static_cast<int>(std::min<size_t>(max_workers, tail / min_tail));

    if (worker_cnt <= 1) {
        FastCrmModel model(nullptr, bounded_mem != 0, pan_super, 0, len, 1);
        std::vector<uint8_t> seg;
        encode_crm_segment(two_bit_bases, 0, static_cast<size_t>(len), model, seg);
        ByteVecOut out;
        out.write_var_u64(1);
        out.write_byte(model_flags);
        out.write_var_u64(len);
        out.write_var_u64(seg.size());
        out.write_bytes(seg.data(), seg.size());
        out_payload->len = out.buf.size();
        out_payload->data = static_cast<uint8_t*>(std::malloc(out_payload->len));
        if (!out_payload->data) return -1;
        std::memcpy(out_payload->data, out.buf.data(), out_payload->len);
        return 0;
    }

    int total_segs = 1 + worker_cnt;
    FastCrmModel warm_m(nullptr, bounded_mem != 0, pan_super, static_cast<int64_t>(warm_bases), len, total_segs);
    std::vector<std::vector<uint8_t>> encoded(total_segs);
    encode_crm_segment(two_bit_bases, 0, warm_bases, warm_m, encoded[0]);
    auto snap = warm_m.snapshot();

    std::vector<size_t> seg_starts(total_segs, 0), seg_lens(total_segs, 0);
    seg_lens[0] = warm_bases;
    size_t cursor = warm_bases;
    for (int w = 0; w < worker_cnt; w++) {
        size_t sl = tail / static_cast<size_t>(worker_cnt) + (static_cast<size_t>(w) < (tail % static_cast<size_t>(worker_cnt)) ? 1 : 0);
        seg_starts[w + 1] = cursor;
        seg_lens[w + 1] = sl;
        cursor += sl;
    }
    int pool_th = bounded_mem ? std::min(worker_cnt, 6) : worker_cnt;
    std::vector<std::thread> pool;
    std::atomic<int> next_w{0};
    for (int p = 0; p < pool_th; p++) {
        pool.emplace_back([&]() {
            while (true) {
                int w = next_w.fetch_add(1);
                if (w >= worker_cnt) break;
                int idx = w + 1;
                FastCrmModel m(snap.get(), bounded_mem != 0, pan_super, static_cast<int64_t>(warm_bases), len, total_segs);
                encode_crm_segment(two_bit_bases, seg_starts[idx], seg_lens[idx], m, encoded[idx]);
            }
        });
    }
    for (auto& t : pool) t.join();

    ByteVecOut out;
    out.write_var_u64(static_cast<uint64_t>(total_segs));
    out.write_byte(model_flags);
    for (int i = 0; i < total_segs; i++) out.write_var_u64(seg_lens[i]);
    for (int i = 0; i < total_segs; i++) out.write_var_u64(encoded[i].size());
    for (int i = 0; i < total_segs; i++) out.write_bytes(encoded[i].data(), encoded[i].size());
    out_payload->len = out.buf.size();
    out_payload->data = static_cast<uint8_t*>(std::malloc(out_payload->len));
    if (!out_payload->data) return -1;
    std::memcpy(out_payload->data, out.buf.data(), out_payload->len);
    return 0;
}

int hgs_cpp_novelty_decode(const uint8_t* payload, size_t payload_len, uint64_t expected_len, int threads, int bounded_mem, HgsOwnedBytes* out_two_bit) {
    out_two_bit->data = nullptr;
    out_two_bit->len = 0;
    ByteSpanIn in{payload, payload_len, 0};
    bool ok = true;
    uint64_t total_segs_u64 = in.read_var_u64(ok);
    if (!ok) return -1;
    if (expected_len == 0) return (total_segs_u64 == 0) ? 0 : -1;
    if (total_segs_u64 == 0 || total_segs_u64 > expected_len / 1024 + 2) return -1;

    size_t total_segs = static_cast<size_t>(total_segs_u64);
    uint8_t model_flags = in.read_byte(ok);
    if (!ok || (model_flags & ~7U) != 0U) return -1;
    bool enc_bounded = (model_flags & 1U) != 0;
    bool enc_pan = (model_flags & 2U) != 0;
    bool enc_turbo = (model_flags & FLAG_TURBO_NOVELTY) != 0;

    std::vector<size_t> seg_lens(total_segs), seg_starts(total_segs), comp_lens(total_segs);
    uint64_t sum_len = 0;
    for (size_t i = 0; i < total_segs; i++) {
        uint64_t sl = in.read_var_u64(ok);
        if (!ok || sum_len + sl > expected_len) return -1;
        seg_starts[i] = static_cast<size_t>(sum_len);
        seg_lens[i] = static_cast<size_t>(sl);
        sum_len += sl;
    }
    if (sum_len != expected_len) return -1;
    for (size_t i = 0; i < total_segs; i++) {
        uint64_t cl = in.read_var_u64(ok);
        if (!ok || in.pos + cl > payload_len) return -1;
        comp_lens[i] = static_cast<size_t>(cl);
    }
    std::vector<const uint8_t*> ptrs(total_segs);
    for (size_t i = 0; i < total_segs; i++) {
        if (in.pos + comp_lens[i] > payload_len) return -1;
        ptrs[i] = payload + in.pos;
        in.pos += comp_lens[i];
    }

    uint8_t* dst = static_cast<uint8_t*>(std::malloc(static_cast<size_t>(expected_len)));
    if (!dst) return -1;
    int th = std::max(1, threads);

    if (enc_turbo) {
        int workers = std::max(1, std::min<int>(th, static_cast<int>(std::min<size_t>(total_segs, bounded_mem ? 4 : 16))));
        std::atomic<bool> fail{false};
        if (workers <= 1) {
            for (size_t i = 0; i < total_segs; i++) {
                if (!decode_turbo_segment(ptrs[i], comp_lens[i], dst + seg_starts[i], seg_lens[i])) fail = true;
            }
        } else {
            std::vector<std::thread> pool;
            std::atomic<size_t> next_i{0};
            for (int w = 0; w < workers; w++) {
                pool.emplace_back([&]() {
                    while (true) {
                        size_t i = next_i.fetch_add(1);
                        if (i >= total_segs) break;
                        if (!decode_turbo_segment(ptrs[i], comp_lens[i], dst + seg_starts[i], seg_lens[i])) fail = true;
                    }
                });
            }
            for (auto& t : pool) t.join();
        }
        if (fail) { std::free(dst); return -1; }
        out_two_bit->data = dst;
        out_two_bit->len = static_cast<size_t>(expected_len);
        return 0;
    }

    if (total_segs == 1) {
        FastCrmModel model(nullptr, enc_bounded, enc_pan, 0, expected_len, 1);
        if (!decode_crm_segment(ptrs[0], comp_lens[0], dst, 0, seg_lens[0], model, false)) { std::free(dst); return -1; }
        out_two_bit->data = dst;
        out_two_bit->len = static_cast<size_t>(expected_len);
        return 0;
    }

    size_t warm_bases = seg_lens[0];
    FastCrmModel warm_m(nullptr, enc_bounded, enc_pan, static_cast<int64_t>(warm_bases), expected_len, static_cast<int>(total_segs));
    if (!decode_crm_segment(ptrs[0], comp_lens[0], dst, 0, warm_bases, warm_m, true)) { std::free(dst); return -1; }
    auto snap = warm_m.snapshot();

    int max_w = bounded_mem ? std::min(th, 6) : (enc_pan ? std::min(th, 8) : th);
    int workers = std::min<int>(max_w, static_cast<int>(total_segs - 1));
    std::atomic<bool> fail{false};
    if (workers <= 1) {
        for (size_t i = 1; i < total_segs; i++) {
            FastCrmModel m(snap.get(), enc_bounded, enc_pan, static_cast<int64_t>(warm_bases), expected_len, static_cast<int>(total_segs));
            if (!decode_crm_segment(ptrs[i], comp_lens[i], dst, seg_starts[i], seg_lens[i], m, false)) fail = true;
        }
    } else {
        std::vector<std::thread> pool;
        std::atomic<size_t> next_i{1};
        for (int w = 0; w < workers; w++) {
            pool.emplace_back([&]() {
                while (true) {
                    size_t i = next_i.fetch_add(1);
                    if (i >= total_segs) break;
                    FastCrmModel m(snap.get(), enc_bounded, enc_pan, static_cast<int64_t>(warm_bases), expected_len, static_cast<int>(total_segs));
                    if (!decode_crm_segment(ptrs[i], comp_lens[i], dst, seg_starts[i], seg_lens[i], m, false)) fail = true;
                }
            });
        }
        for (auto& t : pool) t.join();
    }
    if (fail) { std::free(dst); return -1; }
    out_two_bit->data = dst;
    out_two_bit->len = static_cast<size_t>(expected_len);
    return 0;
}

} // extern "C"
