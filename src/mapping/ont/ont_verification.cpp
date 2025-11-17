// ont_verification.cpp - ONT alignment verification using extended Hamming
#include "ont_mapper.h"
#include "../../core/dna_utils.h"
#include <algorithm>
#include <climits>

namespace rnamapper {

std::optional<Alignment> ONTMapper::verify_one_ont(
    const std::string &read,
    const std::vector<uint8_t> &read_2bit,
    const Cand &c,
    bool strand_is_rev) const
{
    const int L = (int)read.size();
    const int bin = 1 << bin_shift;
    const int off0 = (c.off_q << bin_shift);

    // ONT-specific verification parameters
    // Use 20% error tolerance (vs Illumina's 10%)
    double mm_frac = max_error_rate;  // 0.20 from header

    // Determine target sequence and length
    const uint8_t *target_seq = nullptr;
    uint32_t target_len = 0;

    if (IX.is_jx_tid(c.tid)) {
        const auto &J = IX.jx[IX.jx_index(c.tid)];
        target_seq = J.seq2b;
        target_len = J.len_bp;
    } else {
        const auto &T = IX.tx[c.tid];
        target_seq = T.seq2b;
        target_len = T.len_bp;
    }

    // Compute maximum allowed mismatches
    int max_mm = (int)(L * mm_frac);

    // Try multiple positions around the seed hit
    // ONT needs larger search window due to indels
    int best_mm = INT_MAX;
    int best_pos = -1;

    // Search window: ±bin around seed hit
    int search_start = std::max(0, off0 - bin/2);
    int search_end = std::min((int)target_len - L, off0 + bin + bin/2);

    for (int pos = search_start; pos <= search_end; pos++) {
        if (pos < 0 || pos + L > (int)target_len) {
            continue;
        }

        // Run Hamming distance
        int mm = ham.run(read_2bit, target_seq, (uint32_t)pos, L, max_mm);

        if (mm < best_mm) {
            best_mm = mm;
            best_pos = pos;

            // Early exit if very good match found
            if (mm <= (int)(L * 0.05)) {  // <5% error
                break;
            }
        }
    }

    // Check if alignment is acceptable
    if (best_pos < 0 || best_mm > max_mm) {
        return std::nullopt;
    }

    // Create alignment result
    Alignment A;
    A.tid = c.tid;
    A.pos = best_pos;
    A.mm = best_mm;
    A.score = L - 2 * best_mm;  // Simple scoring: match=1, mismatch=-1
    A.is_reverse = strand_is_rev;
    A.is_ambiguous = false;

    return A;
}

int ONTMapper::choose_verify_budget_ont(
    const std::vector<Cand> &candidates,
    uint16_t best_votes,
    uint16_t best_runner) const
{
    if (candidates.empty()) {
        return 0;
    }

    // ONT: verify more candidates due to high error rate
    // Higher error means more seed failures, so cast wider net
    int base_budget = 15;  // vs Illumina's 12

    // If clear winner, reduce budget
    if (best_votes >= (uint16_t)min_votes &&
        (int)best_votes - (int)best_runner >= lead_margin) {
        base_budget = 8;
    }

    return std::min(base_budget, (int)candidates.size());
}

} // namespace rnamapper
