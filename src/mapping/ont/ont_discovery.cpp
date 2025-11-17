// ont_discovery.cpp - ONT candidate transcript discovery
#include "ont_mapper.h"

namespace rnamapper {

// TODO Step 5: Implement candidate discovery for ONT reads
//
// This file will contain:
// - discover_ont() implementation
// - Relaxed voting thresholds (min_votes=2 vs Illumina's 3)
// - Candidate ranking with high error tolerance
// - Seed aggregation across long reads
//
// Key differences from Illumina:
// - Fewer seeds expected due to sparse sampling
// - Higher tolerance for seed mismatches
// - May need chaining for very long reads (>10kb)

} // namespace rnamapper
