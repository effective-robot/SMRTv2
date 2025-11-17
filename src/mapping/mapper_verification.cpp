// mapper_verification.cpp - Alignment verification
#include "mapper.h"

namespace rnamapper {

std::optional<Alignment> Mapper::verify_one(const std::string &read,
                                            const std::vector<uint8_t> &read2b,
                                            const Cand &c,
                                            bool strand_is_rev,
                                            double mm_frac_tx,
                                            double mm_frac_jx) const {
    const int L = (int)read.size();
    const int bin = 1 << bin_shift;
    const int off0 = (c.off_q << bin_shift);

    auto try_target = [&](const uint8_t *seq2b, uint32_t len_bp, double mm_frac) -> std::optional<Alignment> {
        int max_mm = (int)(L * mm_frac + 0.5);
        int best_mm = INT_MAX, best_pos = -1;
        for (int d = 0; d < bin; ++d) {
            int pos = off0 + d;
            if (pos < 0 || pos + L > (int)len_bp) continue;
            int mm = ham.run(read2b, seq2b, (uint32_t)pos, L, max_mm);
            if (mm < best_mm) {
                best_mm = mm;
                best_pos = pos;
                if (mm == 0) break;
            }
        }
        if (best_pos < 0 || best_mm > max_mm) return std::nullopt;

        Alignment A;
        A.tid = c.tid;
        A.pos = best_pos;
        A.mm = best_mm;
        A.score = L - 2 * best_mm;
        A.is_reverse = strand_is_rev;
        return A;
    };

    if (IX.is_jx_tid(c.tid)) {
        const auto &J = IX.jx[IX.jx_index(c.tid)];
        return try_target(J.seq2b, J.len_bp, mm_frac_jx);
    } else {
        const auto &T = IX.tx[c.tid];
        return try_target(T.seq2b, T.len_bp, mm_frac_tx);
    }
}

int Mapper::choose_verify_budget(const std::vector<Cand> &C, uint16_t best_votes, uint16_t best_runner) const {
    if (C.empty()) return 0;
    int topk = verify_topk;
    if (best_votes >= (uint16_t)min_votes && (int)best_votes - (int)best_runner >= lead_margin) {
        topk = std::min(verify_topk_lo, verify_topk);
    }
    return std::min<int>(topk, (int)C.size());
}

bool Mapper::forward_is_clear(uint16_t best_votes, uint16_t best_runner) const {
    return (best_votes >= (uint16_t)min_votes) && ((int)best_votes - (int)best_runner) >= (lead_margin + 1);
}

} // namespace rnamapper
