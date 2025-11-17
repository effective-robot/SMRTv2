// illumina_mapq.cpp - MAPQ calculation
#include "illumina_mapper.h"
#include <cmath>
#include <algorithm>

namespace rnamapper {

uint8_t IlluminaMapper::calculate_mapq(const Alignment *best, const std::vector<Alignment> &all_hits, int read_len) const {
    if (!best) return 0;

    // Step 1: Base MAPQ from alignment quality
    int base_mapq = 0;
    if (best->mm == 0) {
        base_mapq = 60; // Perfect match
    } else {
        double err_rate = (double)best->mm / std::max(1, read_len);
        if (err_rate < 0.01) err_rate = 0.01; // Avoid log(0)
        base_mapq = (int)(-10.0 * log10(err_rate));
        base_mapq = std::min(60, std::max(0, base_mapq));
    }

    // Step 2: Check for alternative hits (ambiguity penalty)
    if (all_hits.size() >= 2) {
        int score_diff = best->score - all_hits[1].score;

        // Get chromosomes of best and 2nd-best
        GenomicPos gp1 = transcript_to_genomic(best->tid, best->pos, read_len);
        GenomicPos gp2 = transcript_to_genomic(all_hits[1].tid, all_hits[1].pos, read_len);

        bool diff_chr = (gp1.valid && gp2.valid && gp1.chrom != gp2.chrom);

        // Apply penalties based on ambiguity
        if (score_diff == 0 && diff_chr) {
            // Perfect tie on different chromosomes → AMBIGUOUS
            return 0;
        } else if (score_diff <= 2 && diff_chr) {
            // Near-tie on different chromosomes → LOW CONFIDENCE
            return std::min(base_mapq, 3);
        } else if (score_diff <= 5 && diff_chr) {
            // Some ambiguity across chromosomes
            return std::min(base_mapq, 15);
        } else if (score_diff <= 2) {
            // Near-tie on same chromosome (isoform ambiguity)
            return std::min(base_mapq, 20);
        } else if (score_diff <= 5) {
            // Minor ambiguity
            return std::min(base_mapq, 40);
        }
    }

    // Step 3: Single hit but marked ambiguous during discovery
    if (all_hits.size() == 1 && best->is_ambiguous) {
        return std::min(base_mapq, 10);
    }

    return base_mapq;
}

} // namespace rnamapper
