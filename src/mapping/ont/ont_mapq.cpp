// ont_mapq.cpp - MAPQ calculation for ONT reads
#include "ont_mapper.h"
#include <algorithm>
#include <cmath>

namespace rnamapper {

uint8_t ONTMapper::calculate_mapq(
    const Alignment *best,
    const std::vector<Alignment> &all_hits,
    int read_len) const
{
    if (!best) {
        return 0;
    }

    // Base MAPQ from alignment quality
    int base_mapq = 0;

    if (best->mm == 0) {
        base_mapq = 60;  // Perfect match
    } else {
        double err_rate = (double)best->mm / std::max(1, read_len);
        if (err_rate < 0.01) err_rate = 0.01;  // Avoid log(0)
        base_mapq = (int)(-10.0 * std::log10(err_rate));
        base_mapq = std::min(60, std::max(0, base_mapq));
    }

    // Adjust for ONT error rate (10-20% baseline)
    // Long reads with 10% error should still get decent MAPQ
    if (read_len >= 1000) {
        base_mapq = std::min(60, base_mapq + 10);  // Bonus for long reads
    }

    // Ambiguity penalty
    if (all_hits.size() >= 2) {
        int score_diff = best->score - all_hits[1].score;

        // Check if different chromosomes
        GenomicPos gp1 = transcript_to_genomic(best->tid, best->pos, read_len);
        GenomicPos gp2 = transcript_to_genomic(all_hits[1].tid, all_hits[1].pos, read_len);

        bool diff_chr = (gp1.valid && gp2.valid && gp1.chrom != gp2.chrom);

        if (score_diff == 0 && diff_chr) {
            return 0;  // Perfect tie, different chromosomes
        } else if (score_diff <= 5 && diff_chr) {
            return std::min(base_mapq, 10);  // Near-tie, different chromosomes
        } else if (score_diff <= 10) {
            return std::min(base_mapq, 30);  // Some ambiguity
        }
    }

    return (uint8_t)base_mapq;
}

} // namespace rnamapper
