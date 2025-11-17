// mapper_v15_profile_detailed_forward.h
#pragma once
#include <string>

struct TruthRec; // from v15
struct IndexVX;
struct Mapper {
    // expose the nested types we need:
    struct Alignment {
        uint32_t tid;
        int pos;
        int mm;
        int score;
        bool is_reverse;
    };

    // public fields used by rescue
    int min_votes;
    bool skip_common;
    int sieve_seeds;
    int verify_topk;
    int verify_topk_lo;
    unsigned int df1_bytes;
    unsigned int df2_bytes;
    bool skip_rev_on_clear;

    // methods used
    void next_read_epoch();
    struct MapperDiscoverRes { /* not required here, we call discover via existing signature */ };
    // declare discover and verify_one signatures (they are defined in v15)
    struct DiscoverRes;
    DiscoverRes discover(const std::string &seq, bool rev);
    std::optional<Alignment> verify_one(const std::string &read,
                                        const std::vector<uint8_t> &read2b,
                                        const DiscoverRes::Cand &c, bool strand_is_rev,
                                        double mm_frac_tx = 0.10, double mm_frac_jx = 0.05) const;

    // helpers used
    std::vector<uint8_t> pack2(const std::string &seq) const;
};
