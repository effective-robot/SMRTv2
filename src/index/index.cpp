// index.cpp - Index loading implementation
#include "index.h"
#include <fstream>
#include <cstring>
#include <stdexcept>
#include <algorithm>
#include <cstdio>

namespace rnamapper {

uint32_t IndexVX::rd32(const uint8_t *&p, const uint8_t *end) {
    if (end - p < 4) throw std::runtime_error("rd32");
    uint32_t v;
    memcpy(&v, p, 4);
    p += 4;
    return v;
}

uint64_t IndexVX::rd64(const uint8_t *&p, const uint8_t *end) {
    if (end - p < 8) throw std::runtime_error("rd64");
    uint64_t v;
    memcpy(&v, p, 8);
    p += 8;
    return v;
}

std::string IndexVX::rd_str(const uint8_t *&p, const uint8_t *end) {
    uint32_t n = rd32(p, end);
    if (n > (uint64_t)(end - p)) throw std::runtime_error("rd_str");
    std::string s;
    s.resize(n);
    if (n) memcpy(&s[0], p, n);
    p += n;
    return s;
}

void IndexVX::load(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        fprintf(stderr, "ERROR: cannot open index %s\n", path.c_str());
        exit(1);
    }
    in.seekg(0, std::ios::end);
    size_t sz = (size_t)in.tellg();
    in.seekg(0, std::ios::beg);
    blob.resize(sz);
    in.read((char *)blob.data(), sz);
    in.close();

    const uint8_t *p = blob.data();
    const uint8_t *end = blob.data() + blob.size();

    bool is_v15 = false, is_v14 = false, is_v13 = false, is_v12 = false,
         is_v10a = false, is_v10 = false;

    if (memcmp(p, "RMIDX15\n", 8) == 0) {
        is_v15 = true;
        has_bloom = false;
    } else if (memcmp(p, "RMIDX14\n", 8) == 0) {
        is_v14 = true;
        has_bloom = false;
    } else if (memcmp(p, "RMIDX13\n", 8) == 0) {
        is_v13 = true;
        has_bloom = false;
    } else if (memcmp(p, "RMIDX12\n", 8) == 0) {
        is_v12 = true;
        has_bloom = false;
    } else if (memcmp(p, "RMIDX10A", 8) == 0) {
        is_v10a = true;
        has_bloom = true;
    } else if (memcmp(p, "RMIDXv10", 8) == 0) {
        is_v10 = true;
        has_bloom = false;
    } else {
        throw std::runtime_error("Unsupported index");
    }
    p += 8;

    k = rd32(p, end);
    flank = rd32(p, end);
    prefix_bits = rd32(p, end);

    if ((is_v13 || is_v14 || is_v15) && (size_t)(end - p) >= 4 && memcmp(p, "SQTB", 4) == 0) {
        p += 4;
        uint32_t n_contigs = rd32(p, end);
        contigs.resize(n_contigs);
        for (uint32_t i = 0; i < n_contigs; ++i) {
            std::string nm = rd_str(p, end);
            uint32_t ln = rd32(p, end);
            contigs[i] = {nm, ln};
            contig_len[nm] = ln;
        }
    }

    n_tx = rd32(p, end);
    n_jx = rd32(p, end);
    n_targets = rd32(p, end);
    if (n_tx + n_jx != n_targets) throw std::runtime_error("n_targets mismatch");

    tx.resize(n_tx);
    std::vector<uint64_t> tx_off(n_tx, 0);
    for (uint32_t i = 0; i < n_tx; i++) {
        auto &t = tx[i];
        t.id = rd_str(p, end);
        t.chrom = rd_str(p, end);
        if (end - p < 1) throw std::runtime_error("tx strand");
        t.strand = (char)*p++;
        t.len_bp = rd32(p, end);
        t.bytes = (t.len_bp + 3) / 4;
        tx_off[i] = rd64(p, end);

        if (is_v15 || is_v14 || is_v13 || is_v12) {
            uint32_t n_exons = rd32(p, end);
            t.genomic_exons.reserve(n_exons);
            for (uint32_t e = 0; e < n_exons; ++e) {
                uint32_t s0 = rd32(p, end);
                uint32_t e0 = rd32(p, end);
                t.genomic_exons.emplace_back(s0, e0);
            }
        }
    }

    jx.resize(n_jx);
    std::vector<uint64_t> jx_off(n_jx, 0);
    for (uint32_t i = 0; i < n_jx; i++) {
        auto &x = jx[i];
        x.id = rd_str(p, end);
        x.gene_id = rd_str(p, end);
        x.chrom = rd_str(p, end);
        if (end - p < 1) throw std::runtime_error("jx strand");
        x.strand = (char)*p++;
        x.donor_pos = rd32(p, end);
        x.acceptor_pos = rd32(p, end);
        x.len_bp = rd32(p, end);
        x.bytes = (x.len_bp + 3) / 4;
        jx_off[i] = rd64(p, end);
    }

    seq_blob_size = rd64(p, end);
    if (seq_blob_size > (uint64_t)(end - p)) throw std::runtime_error("seq blob");
    seq_blob = p;
    p += seq_blob_size;

    for (uint32_t i = 0; i < n_tx; i++) {
        if (tx_off[i] + tx[i].bytes > seq_blob_size) throw std::runtime_error("tx seq");
        tx[i].seq2b = seq_blob + tx_off[i];
    }
    for (uint32_t i = 0; i < n_jx; i++) {
        if (jx_off[i] + jx[i].bytes > seq_blob_size) throw std::runtime_error("jx seq");
        jx[i].seq2b = seq_blob + jx_off[i];
    }

    n_buckets = rd32(p, end);
    n_keys = rd32(p, end);

    if ((size_t)(end - p) < (size_t)(n_buckets + 1) * 4) throw std::runtime_error("bucket_cum");
    bucket_cum_storage.resize(n_buckets + 1);
    memcpy(bucket_cum_storage.data(), p, (size_t)(n_buckets + 1) * 4);
    p += (size_t)(n_buckets + 1) * 4;
    bucket_cum = bucket_cum_storage.data();

    if ((uint64_t)n_keys * 6 > (uint64_t)(end - p)) throw std::runtime_error("suffix");
    suffix_storage.resize((size_t)n_keys * 6);
    memcpy(suffix_storage.data(), p, suffix_storage.size());
    p += suffix_storage.size();
    key_suffixes = suffix_storage.data();

    if ((size_t)(end - p) < (size_t)(n_keys + 1) * 4) throw std::runtime_error("post offsets");
    post_off_storage.resize(n_keys + 1);
    memcpy(post_off_storage.data(), p, (size_t)(n_keys + 1) * 4);
    p += (size_t)(n_keys + 1) * 4;
    post_offsets = post_off_storage.data();

    uint32_t cb_size = rd32(p, end);
    if ((uint64_t)cb_size > (uint64_t)(end - p)) throw std::runtime_error("common bits");
    common_storage.resize(cb_size);
    if (cb_size) {
        memcpy(common_storage.data(), p, cb_size);
        p += cb_size;
    }
    common_bits = common_storage.data();

    uint64_t postings_size = rd64(p, end);
    if (postings_size > (uint64_t)(end - p)) throw std::runtime_error("postings");
    postings = p;
    p += postings_size;

    if (has_bloom && (size_t)(end - p) >= 20) {
        bloom_n_bits = rd64(p, end);
        bloom_k_hashes = (int)rd32(p, end);
        uint64_t vec_size = rd64(p, end);
        if (vec_size > 0 && vec_size * 8 <= (uint64_t)(end - p)) {
            bloom_bits.resize(vec_size);
            memcpy(bloom_bits.data(), p, vec_size * sizeof(uint64_t));
            p += vec_size * sizeof(uint64_t);
        } else {
            has_bloom = false;
        }
    }

    tx_index_by_id.reserve(n_tx * 1.1);
    for (uint32_t i = 0; i < n_tx; ++i)
        tx_index_by_id[tx[i].id] = i;

    is_v12_plus = (is_v12 || is_v13 || is_v14 || is_v15);
    uses_canonical = (is_v15);
    fprintf(stderr, "Loaded index: k=%u, tx=%u, jx=%u, canonical=%d\n",
            k, n_tx, n_jx, (int)uses_canonical);
}

std::pair<bool, uint32_t> IndexVX::find_key(uint64_t kmer) const {
    uint32_t pre = prefix_of(kmer);
    auto [L, R] = bucket_range(pre);
    if (L == R) return {false, L};
    uint64_t want = suffix_of_kmer(kmer);
    uint32_t lo = L, hi = R;
    while (hi - lo > 8) {
        uint32_t m = (lo + hi) >> 1;
        __builtin_prefetch(&key_suffixes[m * 6], 0, 1);
        uint64_t v = suffix_value_at(m);
        if (v < want)
            lo = m + 1;
        else
            hi = m;
    }
    for (uint32_t i = lo; i < hi; i++) {
        uint64_t v = suffix_value_at(i);
        if (v == want) return {true, i};
        if (v > want) return {false, L};
    }
    return {false, L};
}

} // namespace rnamapper
