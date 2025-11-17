// rescue_v16.h
#pragma once
#include <optional>
#include <string>
#include "mapper_v17_profile_detailed_forward.h" // forward-declare Mapper, IndexVX, TruthRec, Mapper::Alignment etc.

namespace rescue_v16 {

struct ReprocessResult {
    bool recovered = false;
    Mapper::Alignment aln;
    std::string tier_name;
    std::string notes;
};

// attempt_rescue: Attempts tiered rescue for a single read.
// - mapper: reference to your Mapper instance (uses public methods/fields)
// - qname: read id (for logs if you want)
// - seq: read sequence
// - truth: pointer to TruthRec (may be null); if given, success requires matching truth (same as v15 behavior)
// - orig_aln: optional pointer to original alignment (may be null)
ReprocessResult attempt_rescue(Mapper &mapper, const std::string &qname, const std::string &seq,
                               const TruthRec *truth = nullptr, const Mapper::Alignment *orig_aln = nullptr,
                               int truth_slack = 5);

} // namespace rescue_v16
