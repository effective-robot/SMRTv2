// mapper.cpp - Mapper implementation (part 1: initialization and seed generation)
#include "mapper.h"
#include <algorithm>
#include <sstream>

namespace rnamapper {

Mapper::HH4::HH4() {
    gen = 0;
    for (int i = 0; i < 4; i++) {
        off_q[i] = INT32_MIN;
        cnt[i] = 0;
        flg[i] = 0;
    }
}

Mapper::Mapper(const IndexVX &ix, FailureStats &st) : IX(ix), ST(st) {
    per_tid.resize(IX.n_tx + IX.n_jx);
    hh.resize(IX.n_tx + IX.n_jx);
    hh_touched.reserve(512);
    sieve_cnt.resize(IX.n_targets, 0);
    sieve_gen.resize(IX.n_targets, 0);
    tid_list.reserve(256);
}

void Mapper::next_read_epoch() {
    per_tid_gen++;
    if (per_tid_gen == 0) {
        per_tid_gen = 1;
        for (auto &t : per_tid) {
            t.gen = 0;
            t.best = t.second = 0;
            t.best_off_q = INT32_MIN;
        }
        for (auto &h : hh) {
            h.gen = 0;
            for (int i = 0; i < 4; i++) {
                h.off_q[i] = INT32_MIN;
                h.cnt[i] = 0;
                h.flg[i] = 0;
            }
        }
    }
    hh_touched.clear();
}

uint16_t Mapper::bump_hh(uint32_t tid, int32_t off_q, uint8_t flags) {
    HH4 &H = hh[tid];
    if (H.gen != per_tid_gen) {
        H.gen = per_tid_gen;
        for (int i = 0; i < 4; ++i) {
            H.off_q[i] = INT32_MIN;
            H.cnt[i] = 0;
            H.flg[i] = 0;
        }
    }
    for (int i = 0; i < 4; ++i) {
        if (H.off_q[i] == off_q) {
            H.cnt[i] += 1;
            H.flg[i] |= flags;
            return H.cnt[i];
        }
    }
    for (int i = 0; i < 4; ++i) {
        if (H.off_q[i] == INT32_MIN) {
            H.off_q[i] = off_q;
            H.cnt[i] = 1;
            H.flg[i] = flags;
            hh_touched.push_back(tid);
            return 1;
        }
    }
    int worst = 0;
    for (int i = 1; i < 4; ++i) {
        if (H.cnt[i] < H.cnt[worst])
            worst = i;
    }
    H.off_q[worst] = off_q;
    H.cnt[worst] = 1;
    H.flg[worst] = flags;
    hh_touched.push_back(tid);
    return 1;
}

void Mapper::on_vote(uint32_t tid, int32_t off_q, uint16_t cnt) {
    auto &t = per_tid[tid];
    if (t.gen != per_tid_gen) {
        t.gen = per_tid_gen;
        t.best = cnt;
        t.second = 0;
        t.best_off_q = off_q;
        return;
    }
    if (cnt >= t.best) {
        if (off_q != t.best_off_q)
            t.second = t.best;
        t.best = cnt;
        t.best_off_q = off_q;
    } else if (cnt > t.second) {
        t.second = cnt;
    }
}

std::vector<Mapper::Seed> Mapper::make_seeds(const std::string &seq) const {
    const int L = (int)seq.size(), K = (int)IX.k;
    std::vector<Seed> seeds;
    if (L < K) return seeds;
    seeds.reserve(std::max(0, (L - K + 1 + std::max(1, stride) - 1) / std::max(1, stride)));

    const uint64_t KMASK = (K == 32) ? 0xFFFFFFFFFFFFFFFFULL : ((1ULL << (2 * K)) - 1ULL);
    uint64_t roll = 0;
    int good = 0;
    std::unordered_set<uint64_t> seen;
    if (dedup_seeds) seen.reserve(256);

    auto push_seed = [&](int start, uint64_t kmer) {
        if (IX.uses_canonical) {
            // Canonical k-mer
            uint64_t rc = 0;
            uint64_t kx = kmer;
            for (int i = 0; i < (int)IX.k; i++) {
                uint64_t base = kx & 3ULL;
                rc = (rc << 2) | (3ULL - base);
                kx >>= 2;
            }
            if (rc < kmer) kmer = rc;
        }
        if (dedup_seeds && !seen.insert(kmer).second) return;
        if (!IX.bloom_may_contain(kmer)) return;
        auto [ok, idx] = IX.find_key(kmer);
        if (!ok) return;
        uint32_t byte_span = IX.postings_span_bytes(idx);
        uint32_t df_bytes = (byte_span <= 6u) ? 6u : byte_span;
        seeds.push_back(Seed{kmer, start, idx, df_bytes, IX.is_common(idx)});
    };

    for (int i = 0; i < L; i++) {
        int b = nt2b(seq[i]);
        if (b < 0) {
            roll = 0;
            good = 0;
            continue;
        }
        roll = ((roll << 2) | (uint64_t)b) & KMASK;
        if (++good >= K) {
            int start = i - K + 1;
            if ((start % std::max(1, stride)) == 0)
                push_seed(start, roll);
        }
    }

    std::sort(seeds.begin(), seeds.end(), [](const Seed &a, const Seed &b) {
        if (a.common != b.common) return a.common < b.common;
        if (a.df_bytes != b.df_bytes) return a.df_bytes < b.df_bytes;
        return a.pos < b.pos;
    });
    return seeds;
}

} // namespace rnamapper
