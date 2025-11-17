// mapper_v17_profile_detailed.cpp
// Production-ready mapper v17 profile (NO truth DB) with rescue/reprocess and detailed logging.
// Compile: g++ -O3 -std=c++17 -march=native -o mapper_v17_profile_detailed mapper_v17_profile_detailed.cpp -lz
// Note: no external JSON dependency used; all logs are CSV.

#include <algorithm>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <optional>
#include <string>
#include <vector>
#include <zlib.h>
#include <stdexcept>
#include <unordered_set>
#include <unordered_map>
#include <sstream>
#include <iostream>

using namespace std;

// ======================= PROFILING STRUCT ==========================
struct FailureStats
{
    uint64_t reprocess_attempts = 0;
    uint64_t reprocess_recovered = 0;
    uint64_t reprocess_failed = 0;

    uint64_t total_reads = 0;
    uint64_t total_pairs = 0;
    uint64_t ok_mapped_reads = 0;
    uint64_t ok_pairs = 0;
    uint64_t f_len_lt_k = 0;
    uint64_t f_no_seeds = 0;
    uint64_t f_all_seeds_filtered = 0;
    uint64_t f_discover_no_shortlist = 0;
    uint64_t f_verify_failed = 0;
    uint64_t reverse_attempted = 0;
    uint64_t reverse_skipped_clear = 0;
    uint64_t f_genomic_conv = 0;
    uint64_t f_pair_both_unmapped = 0;
    uint64_t f_pair_one_unmapped = 0;
    uint64_t f_pair_incompatible = 0;
    uint64_t reads_with_shortlist = 0;
    uint64_t reads_with_verified_hit = 0;
    uint64_t reads_with_tx_hit_but_bad_genomic = 0;
};

// ======================= SMALL HELPERS =============================
static inline bool ends_with(const string &s, const string &suf)
{
    return s.size() >= suf.size() && equal(suf.rbegin(), suf.rend(), s.rbegin());
}
static inline uint32_t parse_u32(const char *s, uint32_t dflt)
{
    if (!s)
        return dflt;
    try
    {
        return (uint32_t)stoul(string(s));
    }
    catch (...)
    {
        return dflt;
    }
}
static inline int32_t parse_i32(const char *s, int32_t dflt)
{
    if (!s)
        return dflt;
    try
    {
        return (int32_t)stol(string(s));
    }
    catch (...)
    {
        return dflt;
    }
}
static inline uint64_t parse_u64(const char *s, uint64_t dflt)
{
    if (!s)
        return dflt;
    try
    {
        return (uint64_t)stoull(string(s));
    }
    catch (...)
    {
        return dflt;
    }
}

// ======================= FASTQ READER ===============================
struct LineReader
{
    bool gz;
    FILE *fp;
    gzFile gzfp;
    vector<char> buf;
    explicit LineReader(const string &path, size_t bufsz = 1 << 18)
        : gz(ends_with(path, ".gz")), fp(nullptr), gzfp(nullptr), buf(bufsz)
    {
        if (gz)
        {
            gzfp = gzopen(path.c_str(), "rb");
            if (!gzfp)
            {
                fprintf(stderr, "ERROR: cannot open %s\n", path.c_str());
                exit(1);
            }
            gzbuffer(gzfp, (unsigned int)bufsz);
        }
        else
        {
            fp = fopen(path.c_str(), "rb");
            if (!fp)
            {
                fprintf(stderr, "ERROR: cannot open %s\n", path.c_str());
                exit(1);
            }
            setvbuf(fp, nullptr, _IOFBF, bufsz);
        }
    }
    bool getline(string &out)
    {
        out.clear();
        if (gz)
        {
            for (;;)
            {
                char *r = gzgets(gzfp, buf.data(), (int)buf.size());
                if (!r)
                    return !out.empty();
                size_t len = strlen(r);
                if (len && r[len - 1] == '\n')
                {
                    out.append(r, len - 1);
                    return true;
                }
                out.append(r, len);
                if (len < buf.size() - 1)
                    return true;
            }
        }
        else
        {
            int ch;
            while ((ch = fgetc(fp)) != EOF)
            {
                if (ch == '\n')
                    return true;
                out.push_back((char)ch);
            }
            return !out.empty();
        }
    }
    ~LineReader()
    {
        if (gz && gzfp)
            gzclose(gzfp);
        if (!gz && fp)
            fclose(fp);
    }
};

struct ReadRecord
{
    string id, seq, plus, qual;
};

struct FastqReader
{
    LineReader lr;
    explicit FastqReader(const string &path) : lr(path) {}
    bool next(string &id, string &seq, string &plus, string &qual)
    {
        if (!lr.getline(id))
            return false;
        if (!lr.getline(seq))
            return false;
        if (!lr.getline(plus))
            return false;
        if (!lr.getline(qual))
            return false;
        if (!id.empty() && id[0] == '@')
            id = id.substr(1);
        auto sp = id.find_first_of(" \t");
        if (sp != string::npos)
            id = id.substr(0, sp);
        return true;
    }
    size_t next_batch(vector<ReadRecord> &batch, size_t max_reads)
    {
        batch.clear();
        batch.reserve(max_reads);
        for (size_t i = 0; i < max_reads; ++i)
        {
            ReadRecord r;
            if (!next(r.id, r.seq, r.plus, r.qual))
                break;
            batch.push_back(move(r));
        }
        return batch.size();
    }
};

// ======================= DNA HELPERS ================================
static inline string upper(const string &s)
{
    string out(s);
    for (char &c : out)
    {
        unsigned char u = toupper((unsigned char)c);
        c = (u == 'U') ? 'T' : u;
    }
    return out;
}
static inline string revcomp(const string &s)
{
    auto rc = [](char c) -> char
    {
        switch (c)
        {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        default:
            return 'N';
        }
    };
    string out(s.size(), 'N');
    for (size_t i = 0; i < s.size(); ++i)
        out[i] = rc(s[s.size() - 1 - i]);
    return out;
}
static inline int nt2b(char c)
{
    switch (c)
    {
    case 'A':
    case 'a':
        return 0;
    case 'C':
    case 'c':
        return 1;
    case 'G':
    case 'g':
        return 2;
    case 'T':
    case 't':
        return 3;
    default:
        return -1;
    }
}
static inline vector<uint8_t> pack2(const string &seq)
{
    vector<uint8_t> out((seq.size() + 3) / 4, 0);
    size_t o = 0;
    int sh = 6;
    for (char c : seq)
    {
        int b = nt2b(c);
        if (b < 0)
            b = 0;
        out[o] |= (uint8_t)(b << sh);
        sh -= 2;
        if (sh < 0)
        {
            sh = 6;
            ++o;
        }
    }
    return out;
}
static inline uint8_t get2b(const uint8_t *p, uint32_t idx)
{
    uint32_t byte = idx >> 2, sh = 6 - 2 * (idx & 3);
    return (p[byte] >> sh) & 0x3;
}

// ======================= INDEX STRUCT ===============================
struct TranscriptMeta
{
    string id, chrom;
    char strand = '+';
    uint32_t len_bp = 0;
    uint32_t bytes = 0;
    const uint8_t *seq2b = nullptr;
    vector<pair<uint32_t, uint32_t> > genomic_exons;
};
struct JunctionMeta
{
    string id, gene_id, chrom;
    char strand = '+';
    uint32_t donor_pos = 0, acceptor_pos = 0;
    uint32_t len_bp = 0;
    uint32_t bytes = 0;
    const uint8_t *seq2b = nullptr;
};

struct IndexVX
{
    uint32_t k = 31, flank = 50, prefix_bits = 20;
    uint32_t n_tx = 0, n_jx = 0, n_targets = 0;
    uint32_t n_buckets = 0, n_keys = 0;
    bool is_v12_plus = false;
    bool has_bloom = false;
    bool uses_canonical = false;

    vector<TranscriptMeta> tx;
    vector<JunctionMeta> jx;

    const uint32_t *bucket_cum = nullptr;
    const uint8_t *key_suffixes = nullptr;
    const uint32_t *post_offsets = nullptr;
    const uint8_t *common_bits = nullptr;
    const uint8_t *postings = nullptr;
    const uint8_t *seq_blob = nullptr;
    uint64_t seq_blob_size = 0;

    vector<uint8_t> blob;
    vector<uint32_t> bucket_cum_storage;
    vector<uint8_t> suffix_storage;
    vector<uint32_t> post_off_storage;
    vector<uint8_t> common_storage;
    vector<uint64_t> bloom_bits;
    uint64_t bloom_n_bits = 0;
    int bloom_k_hashes = 3;

    vector<pair<string, uint32_t> > contigs;
    unordered_map<string, uint32_t> contig_len;
    unordered_map<string, uint32_t> tx_index_by_id;

    inline uint64_t bloom_hash1(uint64_t key) const
    {
        key ^= key >> 33;
        key *= 0xff51afd7ed558ccdULL;
        key ^= key >> 33;
        key *= 0xc4ceb9fe1a85ec53ULL;
        key ^= key >> 33;
        return key;
    }
    inline uint64_t bloom_hash2(uint64_t key) const
    {
        key = ~key + (key << 21);
        key = key ^ (key >> 24);
        key = (key + (key << 3)) + (key << 8);
        key = key ^ (key >> 14);
        key = (key + (key << 2)) + (key << 4);
        key = key ^ (key >> 28);
        key = key + (key << 31);
        return key;
    }
    inline bool bloom_may_contain(uint64_t kmer) const
    {
        if (!has_bloom)
            return true;
        uint64_t h1 = bloom_hash1(kmer), h2 = bloom_hash2(kmer);
        for (int i = 0; i < bloom_k_hashes; ++i)
        {
            uint64_t h = (h1 + i * h2) % bloom_n_bits;
            if (!(bloom_bits[h >> 6] & (1ULL << (h & 63))))
                return false;
        }
        return true;
    }

    inline uint32_t postings_span_bytes(uint32_t key_idx) const
    {
        return post_offsets[key_idx + 1] - post_offsets[key_idx];
    }
    static inline uint32_t rd32(const uint8_t *&p, const uint8_t *end)
    {
        if (end - p < 4)
            throw runtime_error("rd32");
        uint32_t v;
        memcpy(&v, p, 4);
        p += 4;
        return v;
    }
    static inline uint64_t rd64(const uint8_t *&p, const uint8_t *end)
    {
        if (end - p < 8)
            throw runtime_error("rd64");
        uint64_t v;
        memcpy(&v, p, 8);
        p += 8;
        return v;
    }
    static inline string rd_str(const uint8_t *&p, const uint8_t *end)
    {
        uint32_t n = rd32(p, end);
        if (n > (uint64_t)(end - p))
            throw runtime_error("rd_str");
        string s;
        s.resize(n);
        if (n)
            memcpy(&s[0], p, n);
        p += n;
        return s;
    }

    void load(const string &path)
    {
        ifstream in(path, ios::binary);
        if (!in)
        {
            fprintf(stderr, "ERROR: cannot open index %s\n", path.c_str());
            exit(1);
        }
        in.seekg(0, ios::end);
        size_t sz = (size_t)in.tellg();
        in.seekg(0, ios::beg);
        blob.resize(sz);
        in.read((char *)blob.data(), sz);
        in.close();

        const uint8_t *p = blob.data();
        const uint8_t *end = blob.data() + blob.size();

        bool is_v15 = false, is_v14 = false, is_v13 = false, is_v12 = false, is_v10a = false, is_v10 = false;

        if (memcmp(p, "RMIDX15\n", 8) == 0)
        {
            is_v15 = true;
            has_bloom = false;
        }
        else if (memcmp(p, "RMIDX14\n", 8) == 0)
        {
            is_v14 = true;
            has_bloom = false;
        }
        else if (memcmp(p, "RMIDX13\n", 8) == 0)
        {
            is_v13 = true;
            has_bloom = false;
        }
        else if (memcmp(p, "RMIDX12\n", 8) == 0)
        {
            is_v12 = true;
            has_bloom = false;
        }
        else if (memcmp(p, "RMIDX10A", 8) == 0)
        {
            is_v10a = true;
            has_bloom = true;
        }
        else if (memcmp(p, "RMIDXv10", 8) == 0)
        {
            is_v10 = true;
            has_bloom = false;
        }
        else
            throw runtime_error("Unsupported index");
        p += 8;

        k = rd32(p, end);
        flank = rd32(p, end);
        prefix_bits = rd32(p, end);

        if ((is_v13 || is_v14 || is_v15) && (size_t)(end - p) >= 4 && memcmp(p, "SQTB", 4) == 0)
        {
            p += 4;
            uint32_t n_contigs = rd32(p, end);
            contigs.resize(n_contigs);
            for (uint32_t i = 0; i < n_contigs; ++i)
            {
                string nm = rd_str(p, end);
                uint32_t ln = rd32(p, end);
                contigs[i] = {nm, ln};
                contig_len[nm] = ln;
            }
        }

        n_tx = rd32(p, end);
        n_jx = rd32(p, end);
        n_targets = rd32(p, end);
        if (n_tx + n_jx != n_targets)
            throw runtime_error("n_targets mismatch");

        tx.resize(n_tx);
        vector<uint64_t> tx_off(n_tx, 0);
        for (uint32_t i = 0; i < n_tx; i++)
        {
            auto &t = tx[i];
            t.id = rd_str(p, end);
            t.chrom = rd_str(p, end);
            if (end - p < 1)
                throw runtime_error("tx strand");
            t.strand = (char)*p++;
            t.len_bp = rd32(p, end);
            t.bytes = (t.len_bp + 3) / 4;
            tx_off[i] = rd64(p, end);

            if (is_v15 || is_v14 || is_v13 || is_v12)
            {
                uint32_t n_exons = rd32(p, end);
                t.genomic_exons.reserve(n_exons);
                for (uint32_t e = 0; e < n_exons; ++e)
                {
                    uint32_t s0 = rd32(p, end);
                    uint32_t e0 = rd32(p, end);
                    t.genomic_exons.emplace_back(s0, e0);
                }
            }
        }

        jx.resize(n_jx);
        vector<uint64_t> jx_off(n_jx, 0);
        for (uint32_t i = 0; i < n_jx; i++)
        {
            auto &x = jx[i];
            x.id = rd_str(p, end);
            x.gene_id = rd_str(p, end);
            x.chrom = rd_str(p, end);
            if (end - p < 1)
                throw runtime_error("jx strand");
            x.strand = (char)*p++;
            x.donor_pos = rd32(p, end);
            x.acceptor_pos = rd32(p, end);
            x.len_bp = rd32(p, end);
            x.bytes = (x.len_bp + 3) / 4;
            jx_off[i] = rd64(p, end);
        }

        seq_blob_size = rd64(p, end);
        if (seq_blob_size > (uint64_t)(end - p))
            throw runtime_error("seq blob");
        seq_blob = p;
        p += seq_blob_size;

        for (uint32_t i = 0; i < n_tx; i++)
        {
            if (tx_off[i] + tx[i].bytes > seq_blob_size)
                throw runtime_error("tx seq");
            tx[i].seq2b = seq_blob + tx_off[i];
        }
        for (uint32_t i = 0; i < n_jx; i++)
        {
            if (jx_off[i] + jx[i].bytes > seq_blob_size)
                throw runtime_error("jx seq");
            jx[i].seq2b = seq_blob + jx_off[i];
        }

        n_buckets = rd32(p, end);
        n_keys = rd32(p, end);

        if ((size_t)(end - p) < (size_t)(n_buckets + 1) * 4)
            throw runtime_error("bucket_cum");
        bucket_cum_storage.resize(n_buckets + 1);
        memcpy(bucket_cum_storage.data(), p, (size_t)(n_buckets + 1) * 4);
        p += (size_t)(n_buckets + 1) * 4;
        bucket_cum = bucket_cum_storage.data();

        if ((uint64_t)n_keys * 6 > (uint64_t)(end - p))
            throw runtime_error("suffix");
        suffix_storage.resize((size_t)n_keys * 6);
        memcpy(suffix_storage.data(), p, suffix_storage.size());
        p += suffix_storage.size();
        key_suffixes = suffix_storage.data();

        if ((size_t)(end - p) < (size_t)(n_keys + 1) * 4)
            throw runtime_error("post offsets");
        post_off_storage.resize(n_keys + 1);
        memcpy(post_off_storage.data(), p, (size_t)(n_keys + 1) * 4);
        p += (size_t)(n_keys + 1) * 4;
        post_offsets = post_off_storage.data();

        uint32_t cb_size = rd32(p, end);
        if ((uint64_t)cb_size > (uint64_t)(end - p))
            throw runtime_error("common bits");
        common_storage.resize(cb_size);
        if (cb_size)
        {
            memcpy(common_storage.data(), p, cb_size);
            p += cb_size;
        }
        common_bits = common_storage.data();

        uint64_t postings_size = rd64(p, end);
        if (postings_size > (uint64_t)(end - p))
            throw runtime_error("postings");
        postings = p;
        p += postings_size;

        if (has_bloom && (size_t)(end - p) >= 20)
        {
            bloom_n_bits = rd64(p, end);
            bloom_k_hashes = (int)rd32(p, end);
            uint64_t vec_size = rd64(p, end);
            if (vec_size > 0 && vec_size * 8 <= (uint64_t)(end - p))
            {
                bloom_bits.resize(vec_size);
                memcpy(bloom_bits.data(), p, vec_size * sizeof(uint64_t));
                p += vec_size * sizeof(uint64_t);
            }
            else
            {
                has_bloom = false;
            }
        }

        tx_index_by_id.reserve(n_tx * 1.1);
        for (uint32_t i = 0; i < n_tx; ++i)
            tx_index_by_id[tx[i].id] = i;

        is_v12_plus = (is_v12 || is_v13 || is_v14 || is_v15);
        uses_canonical = (is_v15);
        fprintf(stderr, "Loaded index: k=%u, tx=%u, jx=%u, canonical=%d\n", k, n_tx, n_jx, (int)uses_canonical);
    }

    inline uint32_t prefix_of(uint64_t kmer) const
    {
        return (uint32_t)(kmer >> (62 - prefix_bits));
    }
    inline pair<uint32_t, uint32_t> bucket_range(uint32_t pre) const
    {
        return {bucket_cum[pre], bucket_cum[pre + 1]};
    }
    inline bool is_common(uint32_t key_idx) const
    {
        if (!common_bits)
            return false;
        return (common_bits[key_idx >> 3] >> (key_idx & 7)) & 1u;
    }
    inline uint64_t suffix_value_at(uint32_t idx) const
    {
        const uint8_t *a = key_suffixes + (size_t)idx * 6u;
        return (uint64_t)a[0] | ((uint64_t)a[1] << 8) | ((uint64_t)a[2] << 16) | ((uint64_t)a[3] << 24) | ((uint64_t)a[4] << 32) | ((uint64_t)a[5] << 40);
    }
    inline uint64_t suffix_of_kmer(uint64_t kmer) const
    {
        int bits = 62 - (int)prefix_bits;
        if (bits <= 0)
            return 0ULL;
        return kmer & ((bits == 64) ? ~0ULL : ((1ULL << bits) - 1ULL));
    }

    pair<bool, uint32_t> find_key(uint64_t kmer) const
    {
        uint32_t pre = prefix_of(kmer);
        auto [L, R] = bucket_range(pre);
        if (L == R)
            return {false, L};
        uint64_t want = suffix_of_kmer(kmer);
        uint32_t lo = L, hi = R;
        while (hi - lo > 8)
        {
            uint32_t m = (lo + hi) >> 1;
            __builtin_prefetch(&key_suffixes[m * 6], 0, 1);
            uint64_t v = suffix_value_at(m);
            if (v < want)
                lo = m + 1;
            else
                hi = m;
        }
        for (uint32_t i = lo; i < hi; i++)
        {
            uint64_t v = suffix_value_at(i);
            if (v == want)
                return {true, i};
            if (v > want)
                return {false, L};
        }
        return {false, L};
    }

    template <typename CB>
    inline void for_each_posting_breakable(uint32_t key_idx, CB cb, bool do_prefetch) const
    {
        const uint32_t off0 = post_offsets[key_idx];
        const uint32_t off1 = post_offsets[key_idx + 1];
        const uint8_t *p = postings + off0;
        const uint8_t *e = postings + off1;

        auto zigzag_decode = [](uint64_t n) -> int64_t
        {
            return (int64_t)((n >> 1) ^ (-(int64_t)(n & 1)));
        };
        auto rd = [&](uint64_t &x) -> bool
        {
            x = 0;
            uint32_t s = 0;
            while (p < e)
            {
                if (do_prefetch && ((uintptr_t)p & 63u) == 0)
                    __builtin_prefetch(p + 128, 0, 1);
                uint8_t b = *p++;
                x |= (uint64_t)(b & 0x7FUL) << s;
                if (!(b & 0x80u))
                    return true;
                s += 7;
            }
            return false;
        };

        if (p < e && *p == 0xFF)
        {
            p++;
            if (e - p < 6)
                return;
            uint32_t tid = ((uint32_t)p[0]) | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16);
            uint32_t pos = ((uint32_t)p[3]) | ((uint32_t)p[4] << 8) | ((uint32_t)p[5] << 16);
            cb(tid, pos);
            return;
        }

        if (p >= e)
            return;
        uint64_t G = 0;
        if (!rd(G))
            return;
        uint64_t prev_tid = 0;
        for (uint64_t g = 0; g < G && p < e; ++g)
        {
            uint64_t dt = 0, M = 0;
            if (!rd(dt))
                return;
            uint32_t tid = (uint32_t)(prev_tid + dt);
            prev_tid = tid;
            if (!rd(M))
                return;
            if (M == 0)
                continue;
            uint64_t first_pos_val = 0;
            if (!rd(first_pos_val))
                return;
            uint32_t pos = (uint32_t)first_pos_val;
            if (!cb(tid, pos))
                return;
            if (M > 1)
            {
                int64_t prev_delta = 0;
                for (uint64_t i = 1; i < M && p < e; ++i)
                {
                    uint64_t dod_enc = 0;
                    if (!rd(dod_enc))
                        return;
                    int64_t dod = zigzag_decode(dod_enc);
                    int64_t delta = prev_delta + dod;
                    pos = (uint32_t)((int64_t)pos + delta);
                    if (!cb(tid, pos))
                        return;
                    prev_delta = delta;
                }
            }
        }
    }

    inline bool is_jx_tid(uint32_t tid) const { return tid >= n_tx; }
    inline uint32_t jx_index(uint32_t tid) const { return tid - n_tx; }
};

// ======================= HAMMING ===============================
struct Hamming2B
{
    uint8_t byte_pop4[256];
    Hamming2B()
    {
        for (int x = 0; x < 256; x++)
        {
            int m = 0;
            for (int g = 0; g < 4; g++)
            {
                int sh = 6 - 2 * g;
                if (((x >> sh) & 0x3) != 0)
                    m++;
            }
            byte_pop4[x] = (uint8_t)m;
        }
    }
    int run(const vector<uint8_t> &read_packed,
            const uint8_t *target_packed,
            uint32_t target_start,
            int read_len,
            int max_mm) const
    {
        int mm = 0;
        int full = read_len >> 2, rem = read_len & 3;
        uint32_t tbyte = target_start >> 2;
        int phase = target_start & 3;
        if (phase == 0)
        {
            int full_bytes = full;
            int n64 = full_bytes >> 3;
            const uint64_t *rp = reinterpret_cast<const uint64_t *>(&read_packed[0]);
            const uint64_t *tp = reinterpret_cast<const uint64_t *>(&target_packed[tbyte]);
            for (int i = 0; i < n64; i++)
            {
                uint64_t x = rp[i] ^ tp[i];
                uint64_t y = (x | (x >> 1)) & 0x5555555555555555ULL;
                mm += __builtin_popcountll(y);
                if (mm > max_mm)
                    return mm;
            }
            int consumed = n64 << 3;
            for (int i = consumed; i < full_bytes; i++)
            {
                uint8_t x = read_packed[i] ^ target_packed[tbyte + i];
                uint8_t y = (uint8_t)((x | (x >> 1)) & 0x55);
                mm += __builtin_popcount((unsigned)y);
                if (mm > max_mm)
                    return mm;
            }
        }
        else
        {
            int sh = phase * 2;
            for (int i = 0; i < full; i++)
            {
                uint8_t t = (uint8_t)((target_packed[tbyte + i] << sh) |
                                      (target_packed[tbyte + i + 1] >> (8 - sh)));
                uint8_t x = read_packed[i] ^ t;
                mm += byte_pop4[x];
                if (mm > max_mm)
                    return mm;
            }
        }
        for (int r = 0; r < rem; r++)
        {
            uint32_t ridx = (full << 2) + r;
            uint8_t rb = (read_packed[ridx >> 2] >> (6 - 2 * (ridx & 3))) & 0x3;
            uint8_t tb = get2b(target_packed, target_start + ridx);
            if (rb != tb)
            {
                mm++;
                if (mm > max_mm)
                    return mm;
            }
        }
        return mm;
    }
};

// ======================= OUTPUT BUFFER ==========================
struct OutputBuffer
{
    vector<char> buf;
    size_t capacity;
    FILE *file = nullptr;

    explicit OutputBuffer(const string &path, size_t cap = 32 * 1024 * 1024)
    {
        capacity = cap;
        buf.reserve(capacity);
        file = fopen(path.c_str(), "wb");
        if (!file)
        {
            fprintf(stderr, "ERROR: cannot open %s\n", path.c_str());
            exit(1);
        }
        setvbuf(file, nullptr, _IOFBF, 1 << 20);
    }

    void append(const string &s)
    {
        if (buf.size() + s.size() > capacity)
            flush();
        buf.insert(buf.end(), s.begin(), s.end());
    }

    void flush()
    {
        if (!buf.empty() && file)
        {
            fwrite(buf.data(), 1, buf.size(), file);
            buf.clear();
        }
    }

    ~OutputBuffer()
    {
        flush();
        if (file)
        {
            fclose(file);
            file = nullptr;
        }
    }
};

// ======================= MAPPER ================================
struct Mapper
{
    const IndexVX &IX;
    FailureStats &ST;
    Hamming2B ham;

    int stride = 1, bin_shift = 3, min_votes = 3, lead_margin = 2, min_seeds_before_early = 6;
    bool skip_common = true;
    int verify_topk = 12, verify_topk_lo = 5;
    uint32_t df1_bytes = 128, df2_bytes = 512, post_budget = 0;
    bool prefetch = true, dedup_seeds = false;
    int sieve_seeds = 10;
    bool skip_rev_on_clear = true;
    int min_ins = 100, max_ins = 800;
    bool output_nm_md = false;

    struct HH4
    {
        int32_t off_q[4];
        uint16_t cnt[4];
        uint8_t flg[4];
        uint32_t gen;
        HH4()
        {
            gen = 0;
            for (int i = 0; i < 4; i++)
            {
                off_q[i] = INT32_MIN;
                cnt[i] = 0;
                flg[i] = 0;
            }
        }
    };
    vector<HH4> hh;
    vector<uint32_t> hh_touched;
    struct TBest
    {
        uint16_t best = 0, second = 0;
        uint32_t gen = 0;
        int32_t best_off_q = INT32_MIN;
    };
    vector<TBest> per_tid;
    uint32_t per_tid_gen = 1;
    vector<uint16_t> sieve_cnt;
    vector<uint32_t> sieve_gen;
    vector<uint32_t> tid_list;

    explicit Mapper(const IndexVX &ix, FailureStats &st) : IX(ix), ST(st)
    {
        per_tid.resize(IX.n_tx + IX.n_jx);
        hh.resize(IX.n_tx + IX.n_jx);
        hh_touched.reserve(512);
        sieve_cnt.resize(IX.n_targets, 0);
        sieve_gen.resize(IX.n_targets, 0);
        tid_list.reserve(256);
    }

    inline void next_read_epoch()
    {
        per_tid_gen++;
        if (per_tid_gen == 0)
        {
            per_tid_gen = 1;
            for (auto &t : per_tid)
            {
                t.gen = 0;
                t.best = t.second = 0;
                t.best_off_q = INT32_MIN;
            }
            for (auto &h : hh)
            {
                h.gen = 0;
                for (int i = 0; i < 4; i++)
                {
                    h.off_q[i] = INT32_MIN;
                    h.cnt[i] = 0;
                    h.flg[i] = 0;
                }
            }
        }
        hh_touched.clear();
    }

    struct Seed
    {
        uint64_t kmer;
        int pos;
        uint32_t key_idx, df_bytes;
        bool common;
    };
    inline uint16_t bump_hh(uint32_t tid, int32_t off_q, uint8_t flags)
    {
        HH4 &H = hh[tid];
        if (H.gen != per_tid_gen)
        {
            H.gen = per_tid_gen;
            for (int i = 0; i < 4; ++i)
            {
                H.off_q[i] = INT32_MIN;
                H.cnt[i] = 0;
                H.flg[i] = 0;
            }
        }
        for (int i = 0; i < 4; ++i)
        {
            if (H.off_q[i] == off_q)
            {
                H.cnt[i] += 1;
                H.flg[i] |= flags;
                return H.cnt[i];
            }
        }
        for (int i = 0; i < 4; ++i)
        {
            if (H.off_q[i] == INT32_MIN)
            {
                H.off_q[i] = off_q;
                H.cnt[i] = 1;
                H.flg[i] = flags;
                hh_touched.push_back(tid);
                return 1;
            }
        }
        int worst = 0;
        for (int i = 1; i < 4; ++i)
        {
            if (H.cnt[i] < H.cnt[worst])
                worst = i;
        }
        H.off_q[worst] = off_q;
        H.cnt[worst] = 1;
        H.flg[worst] = flags;
        hh_touched.push_back(tid);
        return 1;
    }

    vector<Seed> make_seeds(const string &seq) const
    {
        const int L = (int)seq.size(), K = (int)IX.k;
        vector<Seed> seeds;
        if (L < K)
            return seeds;
        seeds.reserve(max(0, (L - K + 1 + max(1, stride) - 1) / max(1, stride)));
        const uint64_t KMASK = (K == 32) ? 0xFFFFFFFFFFFFFFFFULL : ((1ULL << (2 * K)) - 1ULL);
        uint64_t roll = 0;
        int good = 0;
        unordered_set<uint64_t> seen;
        if (dedup_seeds)
            seen.reserve(256);

        auto push_seed = [&](int start, uint64_t kmer)
        {
            if (IX.uses_canonical)
            {
                // canonical
                uint64_t rc = 0;
                uint64_t kx = kmer;
                for (int i = 0; i < (int)IX.k; i++)
                {
                    uint64_t base = kx & 3ULL;
                    rc = (rc << 2) | (3ULL - base);
                    kx >>= 2;
                }
                if (rc < kmer)
                    kmer = rc;
            }
            if (dedup_seeds && !seen.insert(kmer).second)
                return;
            if (!IX.bloom_may_contain(kmer))
                return;
            auto [ok, idx] = IX.find_key(kmer);
            if (!ok)
                return;
            uint32_t byte_span = IX.postings_span_bytes(idx);
            uint32_t df_bytes = (byte_span <= 6u) ? 6u : byte_span;
            seeds.push_back(Seed{kmer, start, idx, df_bytes, IX.is_common(idx)});
        };

        for (int i = 0; i < L; i++)
        {
            int b = nt2b(seq[i]);
            if (b < 0)
            {
                roll = 0;
                good = 0;
                continue;
            }
            roll = ((roll << 2) | (uint64_t)b) & KMASK;
            if (++good >= K)
            {
                int start = i - K + 1;
                if ((start % max(1, stride)) == 0)
                    push_seed(start, roll);
            }
        }

        sort(seeds.begin(), seeds.end(), [](const Seed &a, const Seed &b)
             {
            if(a.common != b.common) return a.common < b.common;
            if(a.df_bytes != b.df_bytes) return a.df_bytes < b.df_bytes;
            return a.pos < b.pos; });
        return seeds;
    }

    inline void on_vote(uint32_t tid, int32_t off_q, uint16_t cnt)
    {
        auto &t = per_tid[tid];
        if (t.gen != per_tid_gen)
        {
            t.gen = per_tid_gen;
            t.best = cnt;
            t.second = 0;
            t.best_off_q = off_q;
            return;
        }
        if (cnt >= t.best)
        {
            if (off_q != t.best_off_q)
                t.second = t.best;
            t.best = cnt;
            t.best_off_q = off_q;
        }
        else if (cnt > t.second)
        {
            t.second = cnt;
        }
    }

    struct Cand
    {
        uint32_t tid;
        int32_t off_q;
        uint16_t votes;
        uint8_t flags;
    };

    struct DiscoverRes
    {
        vector<Cand> shortlist;
        uint16_t best_votes = 0, best_runner = 0;
        bool had_any_seed = false;
        bool all_seeds_filtered = false;
    };

    DiscoverRes discover(const string &seq, bool rev)
    {
        next_read_epoch();
        DiscoverRes R;
        auto seeds = make_seeds(seq);
        R.had_any_seed = !seeds.empty();

        uint32_t posts_seen = 0;
        uint16_t best_votes = 0, best_runner = 0;
        uint32_t best_tid = UINT32_MAX;
        const uint32_t KEEP_TOP_TIDS = 32;
        tid_list.clear();

        auto mark_tid = [&](uint32_t tid, uint32_t add)
        {
            if (sieve_gen[tid] != per_tid_gen)
            {
                sieve_gen[tid] = per_tid_gen;
                sieve_cnt[tid] = (uint16_t)min<uint32_t>(65535u, add);
                tid_list.push_back(tid);
            }
            else
            {
                uint32_t v = (uint32_t)sieve_cnt[tid] + add;
                sieve_cnt[tid] = (uint16_t)min<uint32_t>(65535u, v);
            }
        };

        bool saw_effective_posting = false;
        int processed = 0;
        // int s1 = (sieve_seeds > 0) ? min<int>(sieve_seeds, (int)seeds.size()) : 0;

        int s1;
        if (sieve_seeds <= 0)
            s1 = 0;
        else if ((int)seeds.size() <= 1)
            s1 = 0;
        else
            s1 = min<int>(sieve_seeds, (int)seeds.size() - 1);

        for (int i = 0; i < s1; i++)
        {
            const auto &s = seeds[i];
            if (skip_common && s.common)
                continue;
            if (s.df_bytes >= df1_bytes)
                continue;
            IX.for_each_posting_breakable(s.key_idx, [&](uint32_t tid, uint32_t)
                                          {
                posts_seen++; mark_tid(tid, 1); saw_effective_posting = true; return true; }, prefetch);
        }

        vector<uint8_t> allow;
        bool have_allow = false;
        if (s1 > 0 && !tid_list.empty())
        {
            struct Pair
            {
                uint16_t c;
                uint32_t tid;
            };
            auto cmp_min = [](const Pair &a, const Pair &b)
            { return a.c > b.c; };
            vector<Pair> heap;
            heap.reserve(KEEP_TOP_TIDS + 8);
            for (uint32_t tid : tid_list)
            {
                uint16_t c = sieve_cnt[tid];
                if ((int)heap.size() < (int)KEEP_TOP_TIDS)
                {
                    heap.push_back({c, tid});
                    push_heap(heap.begin(), heap.end(), cmp_min);
                }
                else if (c > heap.front().c)
                {
                    pop_heap(heap.begin(), heap.end(), cmp_min);
                    heap.back() = {c, tid};
                    push_heap(heap.begin(), heap.end(), cmp_min);
                }
            }
            sort(heap.begin(), heap.end(), [](const Pair &a, const Pair &b)
                 { return a.c > b.c; });
            allow.assign((IX.n_targets + 7) >> 3, 0);
            auto set_allow = [&](uint32_t tid)
            { allow[tid >> 3] |= (uint8_t)(1u << (tid & 7)); };
            for (const auto &p : heap)
                set_allow(p.tid);
            have_allow = !heap.empty();
        }

        for (size_t i = s1; i < seeds.size(); ++i)
        {
            const auto &s = seeds[i];
            if (skip_common && s.common)
                continue;
            if (s.df_bytes >= df2_bytes)
                continue;

            IX.for_each_posting_breakable(s.key_idx, [&](uint32_t tid, uint32_t pos) -> bool
                                          {
                posts_seen++;
                if(have_allow){
                    if(((allow[tid>>3]>>(tid&7))&1u)==0){
                        if(post_budget && posts_seen>=post_budget && best_tid!=UINT32_MAX) return false;
                        return true;
                    }
                }
                saw_effective_posting = true;

                int32_t off = (int32_t)pos - (int32_t)s.pos;
                int32_t off_q = off >> bin_shift;
                uint8_t flags=0;
                if(IX.is_jx_tid(tid)){
                    uint32_t jxid=IX.jx_index(tid);
                    const auto& J=IX.jx[jxid];
                    uint32_t mid = J.len_bp>>1;
                    flags |= (pos<mid)?1:2;
                }
                uint16_t c=bump_hh(tid, off_q, flags);
                on_vote(tid, off_q, c);

                const auto& tb=per_tid[tid];
                if(tb.gen==per_tid_gen && tb.best==c){
                    if(c>best_votes){ best_votes=c; best_tid=tid; best_runner=tb.second; }
                    else if(tid==best_tid){ best_runner=tb.second; }
                }

                if((processed>=min_seeds_before_early) &&
                   best_tid!=UINT32_MAX &&
                   ((int)best_votes - (int)best_runner) >= lead_margin){ return false; }
                if(post_budget && posts_seen>=post_budget && best_tid!=UINT32_MAX) return false;
                return true; }, prefetch);

            processed++;
            int need = min_votes;
            if (best_tid != UINT32_MAX &&
                best_votes >= need &&
                (int)best_votes - (int)best_runner >= lead_margin &&
                processed >= min_seeds_before_early)
            {
                break;
            }
            if (post_budget && posts_seen >= post_budget && best_tid != UINT32_MAX)
                break;
        }

        vector<Cand> tmp;
        tmp.reserve(hh_touched.size());
        for (uint32_t tid : hh_touched)
        {
            auto &tb = per_tid[tid];
            if (tb.gen != per_tid_gen)
                continue;
            const auto &H = hh[tid];
            int win = -1;
            for (int i = 0; i < 4; i++)
            {
                if (H.cnt[i] == 0)
                    continue;
                if (H.off_q[i] == tb.best_off_q && H.cnt[i] == tb.best)
                {
                    win = i;
                    break;
                }
            }
            if (win == -1)
            {
                for (int i = 0; i < 4; i++)
                {
                    if (H.cnt[i] == tb.best)
                    {
                        win = i;
                        break;
                    }
                }
            }
            if (win != -1)
            {
                tmp.push_back({tid, H.off_q[win], H.cnt[win], H.flg[win]});
            }
        }

        tmp.erase(remove_if(tmp.begin(), tmp.end(), [&](const Cand &c)
                            {
            if(!IX.is_jx_tid(c.tid)) return false;
            return ((c.flags & 3)!=3); }),
                  tmp.end());

        sort(tmp.begin(), tmp.end(), [](const Cand &a, const Cand &b)
             {
            if(a.votes!=b.votes) return a.votes>b.votes;
            if(a.tid!=b.tid) return a.tid<b.tid;
            return a.off_q<b.off_q; });

        R.shortlist = move(tmp);
        R.best_votes = best_votes;
        R.best_runner = best_runner;

        if (R.had_any_seed && !saw_effective_posting)
            R.all_seeds_filtered = true;

        return R;
    }

    struct Alignment
    {
        uint32_t tid;
        int pos;
        int mm;
        int score;
        bool is_reverse;
        bool is_ambiguous = false;
    };

    optional<Alignment> verify_one(const string &read,
                                   const vector<uint8_t> &read2b,
                                   const Cand &c,
                                   bool strand_is_rev,
                                   double mm_frac_tx = 0.10,
                                   double mm_frac_jx = 0.05) const
    {
        const int L = (int)read.size();
        const int bin = 1 << bin_shift;
        const int off0 = (c.off_q << bin_shift);

        auto try_target = [&](const uint8_t *seq2b, uint32_t len_bp, double mm_frac) -> optional<Alignment>
        {
            int max_mm = (int)(L * mm_frac + 0.5);
            int best_mm = INT_MAX, best_pos = -1;
            for (int d = 0; d < bin; ++d)
            {
                int pos = off0 + d;
                if (pos < 0 || pos + L > (int)len_bp)
                    continue;
                int mm = ham.run(read2b, seq2b, (uint32_t)pos, L, max_mm);
                if (mm < best_mm)
                {
                    best_mm = mm;
                    best_pos = pos;
                    if (mm == 0)
                        break;
                }
            }
            if (best_pos < 0 || best_mm > max_mm)
                return nullopt;

            Alignment A;
            A.tid = c.tid;
            A.pos = best_pos;
            A.mm = best_mm;
            A.score = L - 2 * best_mm;
            A.is_reverse = strand_is_rev;
            return A;
        };

        if (IX.is_jx_tid(c.tid))
        {
            const auto &J = IX.jx[IX.jx_index(c.tid)];
            return try_target(J.seq2b, J.len_bp, mm_frac_jx);
        }
        else
        {
            const auto &T = IX.tx[c.tid];
            return try_target(T.seq2b, T.len_bp, mm_frac_tx);
        }
    }

    inline int choose_verify_budget(const vector<Cand> &C, uint16_t best_votes, uint16_t best_runner) const
    {
        if (C.empty())
            return 0;
        int topk = verify_topk;
        if (best_votes >= (uint16_t)min_votes && (int)best_votes - (int)best_runner >= lead_margin)
        {
            topk = min(verify_topk_lo, verify_topk);
        }
        return min<int>(topk, (int)C.size());
    }

    inline bool forward_is_clear(uint16_t best_votes, uint16_t best_runner) const
    {
        return (best_votes >= (uint16_t)min_votes) && ((int)best_votes - (int)best_runner) >= (lead_margin + 1);
    }

    vector<Alignment> map_read_single(const string &qname, const string &seq)
    {
        ST.total_reads++;
        vector<Alignment> hits;
        const int L = (int)seq.size();
        if (L < (int)IX.k)
        {
            ST.f_len_lt_k++;
            return hits;
        }
        string fwd = upper(seq);
        string rev = revcomp(fwd);
        auto f2b = pack2(fwd);
        auto r2b = pack2(rev);

        auto F = discover(fwd, false);
        if (F.all_seeds_filtered)
            ST.f_all_seeds_filtered++;
        if (!F.shortlist.empty())
            ST.reads_with_shortlist++;

        int budgetF = choose_verify_budget(F.shortlist, F.best_votes, F.best_runner);
        for (int i = 0; i < budgetF; i++)
        {
            auto a = verify_one(fwd, f2b, F.shortlist[i], false);
            if (a)
                hits.push_back(*a);
        }

        bool do_reverse = true;
        if (skip_rev_on_clear && forward_is_clear(F.best_votes, F.best_runner) && !hits.empty())
        {
            do_reverse = false;
            ST.reverse_skipped_clear++;
        }

        if (do_reverse)
        {
            ST.reverse_attempted++;
            auto R = discover(rev, true);
            if (R.all_seeds_filtered)
                ST.f_all_seeds_filtered++;

            int budgetR = choose_verify_budget(R.shortlist, R.best_votes, R.best_runner);
            for (int i = 0; i < budgetR; i++)
            {
                auto a = verify_one(rev, r2b, R.shortlist[i], true);
                if (a)
                    hits.push_back(*a);
            }

            if (!F.shortlist.empty() || !R.shortlist.empty())
                ST.reads_with_shortlist++;
        }

        sort(hits.begin(), hits.end(), [](const Alignment &a, const Alignment &b)
             {
            if(a.score!=b.score) return a.score>b.score;
            return a.mm<b.mm; });
        hits.erase(unique(hits.begin(), hits.end(), [](const Alignment &a, const Alignment &b)
                          { return a.tid == b.tid && a.pos == b.pos && a.is_reverse == b.is_reverse; }),
                   hits.end());

        if (hits.empty())
        {
            bool f_empty = F.shortlist.empty();
            if (f_empty)
                ST.f_discover_no_shortlist++;
            else
                ST.f_verify_failed++;
        }
        else
        {
            ST.reads_with_verified_hit++;
        }

        // ===== ADD TIE DETECTION LOGIC HERE (BEFORE return hits;) =====
        if (hits.size() >= 2)
        {
            // Check if top 2 hits have very close scores
            int score_diff = hits[0].score - hits[1].score;

            // Check if they're on different chromosomes
            GenomicPos gp1 = transcript_to_genomic(hits[0].tid, hits[0].pos, (int)seq.size());
            GenomicPos gp2 = transcript_to_genomic(hits[1].tid, hits[1].pos, (int)seq.size());

            if (score_diff <= 2 && gp1.valid && gp2.valid && gp1.chrom != gp2.chrom)
            {
                hits[0].is_ambiguous = true;
            }
        }
        else if (hits.size() == 1)
        {
            // Only one hit - mark as potentially ambiguous
            hits[0].is_ambiguous = true;
        }
        // ===== END TIE DETECTION =====

        return hits;
    }

    static string strip_isoform(const string &s)
    {
        auto d = s.find('.');
        return d == string::npos ? s : s.substr(0, d);
    }

    inline bool tid_compatible_same_gene(uint32_t a, uint32_t b) const
    {
        if (a == b)
            return true;
        bool aj = IX.is_jx_tid(a), bj = IX.is_jx_tid(b);
        if (!aj && !bj)
        {
            auto g1 = strip_isoform(IX.tx[a].id);
            auto g2 = strip_isoform(IX.tx[b].id);
            return g1 == g2;
        }
        if (aj && !bj)
        {
            auto gT = strip_isoform(IX.tx[b].id);
            auto gJ = IX.jx[IX.jx_index(a)].gene_id;
            return gT == gJ;
        }
        if (!aj && bj)
        {
            auto gT = strip_isoform(IX.tx[a].id);
            auto gJ = IX.jx[IX.jx_index(b)].gene_id;
            return gT == gJ;
        }
        const auto &J1 = IX.jx[IX.jx_index(a)];
        const auto &J2 = IX.jx[IX.jx_index(b)];
        return J1.gene_id == J2.gene_id;
    }

    struct GenomicPos
    {
        string chrom;
        uint32_t pos;
        char strand;
        string cigar;
        uint32_t ref_span;
        bool valid;
    };

    GenomicPos transcript_to_genomic(uint32_t tid, int tx_pos, int read_len) const
    {
        GenomicPos gp{"", 0, '+', "", 0, false};
        if (read_len <= 0)
            return gp;

        auto build_blocks = [&](const vector<pair<uint32_t, uint32_t> > &tx_exons,
                                char strand,
                                int tx_pos0,
                                int len) -> vector<pair<uint32_t, uint32_t> >
        {
            vector<pair<uint32_t, uint32_t> > blocks;
            int tx_cursor = 0;
            int remaining = len;
            for (size_t i = 0; i < tx_exons.size() && remaining > 0; ++i)
            {
                uint32_t ex_s = tx_exons[i].first;
                uint32_t ex_e = tx_exons[i].second;
                int ex_len = (int)(ex_e - ex_s);
                if (tx_pos0 >= tx_cursor + ex_len)
                {
                    tx_cursor += ex_len;
                    continue;
                }
                int off_in_ex = max(0, tx_pos0 - tx_cursor);
                int avail = ex_len - off_in_ex;
                int take = min(remaining, avail);
                uint32_t block_start0;
                if (strand == '+')
                {
                    block_start0 = ex_s + (uint32_t)off_in_ex;
                }
                else
                {
                    block_start0 = (ex_e - (uint32_t)off_in_ex) - (uint32_t)take;
                }
                blocks.push_back({block_start0, (uint32_t)take});
                remaining -= take;
                tx_pos0 += take;
                tx_cursor += ex_len;
            }
            sort(blocks.begin(), blocks.end(),
                 [](const auto &a, const auto &b)
                 { return a.first < b.first; });
            return blocks;
        };

        vector<pair<uint32_t, uint32_t> > tx_exons;
        char ref_strand = '+';
        string rname;

        if (IX.is_jx_tid(tid))
        {
            const auto &J = IX.jx[IX.jx_index(tid)];
            rname = J.chrom;
            ref_strand = J.strand;
            uint32_t left_s = (uint32_t)(J.donor_pos) - IX.flank;
            uint32_t left_e = (uint32_t)(J.donor_pos);
            uint32_t right_s = (uint32_t)(J.acceptor_pos - 1);
            uint32_t right_e = right_s + IX.flank;
            if (ref_strand == '+')
            {
                tx_exons = {{left_s, left_e}, {right_s, right_e}};
            }
            else
            {
                tx_exons = {{right_s, right_e}, {left_s, left_e}};
            }
        }
        else
        {
            const auto &T = IX.tx[tid];
            rname = T.chrom;
            ref_strand = T.strand;
            if (T.genomic_exons.empty())
                return gp;
            if (ref_strand == '+')
            {
                tx_exons.assign(T.genomic_exons.begin(), T.genomic_exons.end());
            }
            else
            {
                tx_exons.reserve(T.genomic_exons.size());
                for (auto it = T.genomic_exons.rbegin(); it != T.genomic_exons.rend(); ++it)
                    tx_exons.push_back(*it);
            }
        }

        auto blocks = build_blocks(tx_exons, ref_strand, tx_pos, read_len);
        if (blocks.empty())
            return gp;

        gp.chrom = rname;
        gp.strand = ref_strand;
        gp.pos = blocks.front().first + 1;
        gp.cigar.clear();
        gp.ref_span = 0;

        for (size_t i = 0; i < blocks.size(); ++i)
        {
            if (i > 0)
            {
                uint32_t prev_end = blocks[i - 1].first + blocks[i - 1].second;
                uint32_t intron = blocks[i].first - prev_end;
                gp.cigar += to_string(intron) + "N";
                gp.ref_span += intron;
            }
            gp.cigar += to_string(blocks[i].second) + "M";
            gp.ref_span += blocks[i].second;
        }

        gp.valid = true;
        return gp;
    }

    static void write_sam_header(OutputBuffer &outbuf, const IndexVX &IX)
    {
        string header = "@HD\tVN:1.6\tSO:unsorted\n";
        if (!IX.contigs.empty())
        {
            for (const auto &[chr_name, chr_len] : IX.contigs)
            {
                header += "@SQ\tSN:" + chr_name + "\tLN:" + to_string(chr_len) + "\n";
            }
        }
        else
        {
            unordered_map<string, uint32_t> chr_max_pos;
            for (const auto &t : IX.tx)
            {
                if (!t.genomic_exons.empty())
                {
                    uint32_t max_end = 0;
                    for (const auto &[s, e] : t.genomic_exons)
                    {
                        if (e > max_end)
                            max_end = e;
                    }
                    uint32_t &current_max = chr_max_pos[t.chrom];
                    if (max_end > current_max)
                        current_max = max_end;
                }
            }
            for (const auto &j : IX.jx)
            {
                uint32_t jx_end = max(j.donor_pos, j.acceptor_pos) + IX.flank + 100;
                uint32_t &current_max = chr_max_pos[j.chrom];
                if (jx_end > current_max)
                    current_max = jx_end;
            }
            vector<pair<string, uint32_t> > sorted_chrs(chr_max_pos.begin(), chr_max_pos.end());
            sort(sorted_chrs.begin(), sorted_chrs.end());
            for (const auto &[chr_name, chr_len] : sorted_chrs)
            {
                header += "@SQ\tSN:" + chr_name + "\tLN:" + to_string(chr_len) + "\n";
            }
        }
        header += "@RG\tID:rg1\tSM:sample\tPL:ILLUMINA\n";
        header += "@PG\tID:mapper_v17_profile_detailed\tPN:mapper_v17_profile_detailed\tVN:17\n";
        outbuf.append(header);
    }

    static inline string normalize_qname(const string &s)
    {
        size_t len = s.size();
        if (len >= 2 && s[len - 2] == '/' && (s[len - 1] == '1' || s[len - 1] == '2'))
        {
            return s.substr(0, len - 2);
        }
        return s;
    }

    inline void append_nm_md_tags(string &line, const Alignment *aln, const ReadRecord &rec) const
    {
        if (!aln || !output_nm_md)
            return;
        const uint8_t *ref_seq = IX.is_jx_tid(aln->tid)
                                     ? IX.jx[IX.jx_index(aln->tid)].seq2b
                                     : IX.tx[aln->tid].seq2b;
        int pos0 = aln->pos;
        int L = (int)rec.seq.size();
        int nm = 0;
        int run = 0;
        string md;
        md.reserve(L * 2);
        const char *nt = "ACGT";
        for (int i = 0; i < L; ++i)
        {
            char read_base = rec.seq[i];
            if (read_base == 'a')
                read_base = 'A';
            else if (read_base == 'c')
                read_base = 'C';
            else if (read_base == 'g')
                read_base = 'G';
            else if (read_base == 't' || read_base == 'u' || read_base == 'U')
                read_base = 'T';
            uint8_t ref_2b = get2b(ref_seq, pos0 + i);
            char ref_base = nt[ref_2b];
            if (read_base == ref_base)
            {
                ++run;
            }
            else
            {
                md += to_string(run);
                md += ref_base;
                run = 0;
                ++nm;
            }
        }
        md += to_string(run);
        line += "\tNM:i:";
        line += to_string(nm);
        line += "\tMD:Z:";
        line += md;
    }

    // ===== ADD THIS NEW FUNCTION BEFORE write_sam_record() =====
    uint8_t calculate_mapq(
        const Alignment *best,
        const vector<Alignment> &all_hits,
        int read_len) const
    {
        if (!best)
            return 0;

        // Step 1: Base MAPQ from alignment quality
        int base_mapq = 0;
        if (best->mm == 0)
        {
            base_mapq = 60; // Perfect match
        }
        else
        {
            double err_rate = (double)best->mm / max(1, read_len);
            if (err_rate < 0.01)
                err_rate = 0.01; // Avoid log(0)
            base_mapq = (int)(-10.0 * log10(err_rate));
            base_mapq = min(60, max(0, base_mapq));
        }

        // Step 2: Check for alternative hits (ambiguity penalty)
        if (all_hits.size() >= 2)
        {
            int score_diff = best->score - all_hits[1].score;

            // Get chromosomes of best and 2nd-best
            GenomicPos gp1 = transcript_to_genomic(best->tid, best->pos, read_len);
            GenomicPos gp2 = transcript_to_genomic(all_hits[1].tid, all_hits[1].pos, read_len);

            bool diff_chr = (gp1.valid && gp2.valid && gp1.chrom != gp2.chrom);

            // Apply penalties based on ambiguity
            if (score_diff == 0 && diff_chr)
            {
                // Perfect tie on different chromosomes → AMBIGUOUS
                return 0;
            }
            else if (score_diff <= 2 && diff_chr)
            {
                // Near-tie on different chromosomes → LOW CONFIDENCE
                return min(base_mapq, 3);
            }
            else if (score_diff <= 5 && diff_chr)
            {
                // Some ambiguity across chromosomes
                return min(base_mapq, 15);
            }
            else if (score_diff <= 2)
            {
                // Near-tie on same chromosome (isoform ambiguity)
                return min(base_mapq, 20);
            }
            else if (score_diff <= 5)
            {
                // Minor ambiguity
                return min(base_mapq, 40);
            }
        }

        // Step 3: Single hit but marked ambiguous during discovery
        if (all_hits.size() == 1 && best->is_ambiguous)
        {
            return min(base_mapq, 10);
        }

        return base_mapq;
    }
    // ===== END NEW FUNCTION =====

    // ===== ADD PAIR SCORING LOGIC =====
    struct PairScore
    {
        int r1_score;
        int r2_score;
        int concordance_bonus;
        int insert_quality;
        int total;

        bool operator<(const PairScore &other) const
        {
            return total > other.total; // Higher is better
        }
    };

    PairScore score_pair(
        const Alignment &a1, const Alignment &a2,
        const ReadRecord &rec1, const ReadRecord &rec2) const
    {
        PairScore ps;
        ps.r1_score = a1.score;
        ps.r2_score = a2.score;
        ps.concordance_bonus = 0;
        ps.insert_quality = 0;

        GenomicPos gp1 = transcript_to_genomic(a1.tid, a1.pos, (int)rec1.seq.size());
        GenomicPos gp2 = transcript_to_genomic(a2.tid, a2.pos, (int)rec2.seq.size());

        if (gp1.valid && gp2.valid)
        {
            if (gp1.chrom == gp2.chrom)
            {
                // Same chromosome → good sign
                int ins = abs((int)gp1.pos - (int)gp2.pos);

                if (ins >= min_ins && ins <= max_ins)
                {
                    // Perfect concordant pair
                    ps.concordance_bonus = 15;

                    // Bonus for ideal insert size (300-600 bp for RNA-seq)
                    int ideal_ins = 400;
                    int ins_dev = abs(ins - ideal_ins);
                    if (ins_dev < 50)
                    {
                        ps.insert_quality = 20;
                    }
                    else if (ins_dev < 150)
                    {
                        ps.insert_quality = 10;
                    }
                    else
                    {
                        ps.insert_quality = 5;
                    }
                }
                else if (ins <= max_ins * 2)
                {
                    // Reasonable but not perfect
                    ps.concordance_bonus = 20;
                }
            }
            else
            {
                // Different chromosomes → penalty
                ps.concordance_bonus = -15;
            }
        }

        ps.total = ps.r1_score + ps.r2_score + ps.concordance_bonus + ps.insert_quality;
        return ps;
    }
    // ===== END PAIR SCORING =====

    // ===== ADD PAIR SELECTION FUNCTION =====
    struct BestPair
    {
        Alignment *best1;
        Alignment *best2;
        PairScore score;
        bool found;
    };

    BestPair select_best_pair(
        vector<Alignment> &H1, vector<Alignment> &H2,
        const ReadRecord &rec1, const ReadRecord &rec2) const
    {
        BestPair result;
        result.best1 = nullptr;
        result.best2 = nullptr;
        result.found = false;
        result.score.total = INT_MIN;

        if (H1.empty() || H2.empty())
        {
            return result;
        }

        // Score top N pairs (not all combinations - too slow)
        const int MAX_CHECK = 10;
        int n1 = min<int>(MAX_CHECK, (int)H1.size());
        int n2 = min<int>(MAX_CHECK, (int)H2.size());

        for (int i = 0; i < n1; i++)
        {
            for (int j = 0; j < n2; j++)
            {
                PairScore ps = score_pair(H1[i], H2[j], rec1, rec2);

                if (ps.total > result.score.total)
                {
                    result.score = ps;
                    result.best1 = &H1[i];
                    result.best2 = &H2[j];
                    result.found = true;
                }
            }
        }

        return result;
    }
    // ===== END PAIR SELECTION =====

    void write_sam_record(OutputBuffer &buf, const ReadRecord &rec, const Alignment *aln,
                          bool is_paired, bool is_r2, const Alignment *mate_aln,
                          const ReadRecord *mate_rec,
                          const vector<Alignment> *all_hits = nullptr) // ← ADD THIS)
    {
        string line;
        line.reserve(output_nm_md ? 384 : 256);
        line += normalize_qname(rec.id) + "\t";

        bool is_mapped = (aln != nullptr);
        bool mate_mapped = (mate_aln != nullptr);

        bool is_reverse = false, mate_reverse = false;
        if (aln)
        {
            is_reverse = aln->is_reverse;
        }
        if (mate_aln)
        {
            mate_reverse = mate_aln->is_reverse;
        }

        string rname = "*", mrname = "*";
        uint32_t pos = 0, mpos = 0;
        string cigar = "*";
        int read_len = (int)rec.seq.size();
        GenomicPos gpos, mate_gpos;
        if (is_mapped)
        {
            gpos = transcript_to_genomic(aln->tid, aln->pos, read_len);
            if (gpos.valid)
            {
                rname = gpos.chrom;
                pos = gpos.pos;
                cigar = gpos.cigar;
            }
            else
            {
                is_mapped = false;
                ST.f_genomic_conv++;
                ST.reads_with_tx_hit_but_bad_genomic++;
            }
        }
        if (mate_mapped)
        {
            int mate_len = mate_rec ? (int)mate_rec->seq.size() : read_len;
            mate_gpos = transcript_to_genomic(mate_aln->tid, mate_aln->pos, mate_len);
            if (mate_gpos.valid)
            {
                mrname = mate_gpos.chrom;
                mpos = mate_gpos.pos;
            }
            else
            {
                mate_mapped = false;
                ST.f_genomic_conv++;
                ST.reads_with_tx_hit_but_bad_genomic++;
            }
        }

        bool proper_pair = false;
        if (is_paired && is_mapped && mate_mapped && rname == mrname)
        {
            bool fr_orientation = (!is_reverse && mate_reverse && pos <= mpos) ||
                                  (is_reverse && !mate_reverse && mpos <= pos);
            if (fr_orientation)
            {
                int ins = (int)pos + read_len - (int)mpos;
                if (ins < 0)
                    ins = -ins;
                if (ins >= min_ins && ins <= max_ins)
                {
                    proper_pair = true;
                }
            }
        }

        uint16_t flag = 0;
        if (is_paired)
            flag |= 0x1;
        if (proper_pair)
            flag |= 0x2;
        if (!is_mapped)
            flag |= 0x4;
        if (is_paired && !mate_mapped)
            flag |= 0x8;
        if (is_reverse)
            flag |= 0x10;
        if (is_paired && mate_reverse)
            flag |= 0x20;
        if (is_r2)
            flag |= 0x80;
        else if (is_paired)
            flag |= 0x40;

        line += to_string(flag) + "\t";

        if (is_mapped)
        {
            line += rname + "\t";
            line += to_string(pos) + "\t";
        }
        else
        {
            line += "*\t0\t";
        }

        // ===== FINAL MAPQ CALCULATION WITH ALTERNATIVES =====
        uint8_t mapq = 0;
        if (is_mapped)
        {
            if (all_hits && !all_hits->empty())
            {
                // Use new MAPQ with alternative consideration
                mapq = calculate_mapq(aln, *all_hits, read_len);
            }
            else
            {
                // Fallback: basic MAPQ
                if (aln->mm == 0)
                    mapq = 60;
                else
                {
                    double err_rate = (double)aln->mm / max(1, read_len);
                    if (err_rate < 0.01)
                        err_rate = 0.01;
                    int q = (int)(-10.0 * log10(err_rate));
                    mapq = (uint8_t)min(60, max(0, q));
                }

                if (aln->is_ambiguous)
                {
                    mapq = min(mapq, (uint8_t)10);
                }
            }
        }
        // ===== END FINAL MAPQ =====

        line += to_string((int)mapq) + "\t";

        line += cigar + "\t";

        if (is_paired && mate_mapped)
        {
            if (rname == mrname)
            {
                line += "=";
            }
            else
            {
                line += mrname;
            }
            line.push_back('\t');
            line += to_string(mpos);
            line.push_back('\t');

            int tlen = 0;
            if (is_mapped && rname == mrname && mate_rec)
            {
                const int r1_len = read_len;
                const int r2_len = (int)mate_rec->seq.size();
                const int end1 = (int)pos + r1_len - 1;
                const int end2 = (int)mpos + r2_len - 1;
                const int left = min((int)pos, (int)mpos);
                const int right = max(end1, end2);
                tlen = right - left + 1;
                if ((int)pos != left)
                    tlen = -tlen;
            }
            line += to_string(tlen);
            line.push_back('\t');
        }
        else
        {
            line += "*\t0\t0\t";
        }

        line += rec.seq;
        line.push_back('\t');
        line += rec.qual;
        line.push_back('\t');
        line += "RG:Z:rg1";

        if (is_mapped)
        {
            line.push_back('\t');
            line += "XS:A:";
            line.push_back(
                IX.is_jx_tid(aln->tid)
                    ? IX.jx[IX.jx_index(aln->tid)].strand
                    : IX.tx[aln->tid].strand);
            line.push_back('\t');
            line += "NH:i:1";
        }

        append_nm_md_tags(line, aln, rec);

        line.push_back('\n');
        buf.append(line);
    }

    // ---------- Reprocess logic ----------
    struct ReprocessResult
    {
        bool recovered = false;
        Alignment aln;
        string tier_name;
        string notes;
    };

    // No-truth version: attempt reprocessing / rescue tiers. Accepts any valid genomic projection as recovery.
    ReprocessResult attempt_reprocess_single(const string &qname, const string &seq, const Alignment *orig_aln = nullptr)
    {
        ReprocessResult out;

        // ---- Snapshot mapper state to restore later ----
        int old_min_votes = min_votes;
        int old_lead_margin = lead_margin;
        bool old_skip_common = skip_common;
        int old_sieve_seeds = sieve_seeds;
        bool old_skip_rev = skip_rev_on_clear;
        int old_verify_topk = verify_topk;
        int old_verify_topk_lo = verify_topk_lo;
        uint32_t old_df1 = df1_bytes;
        uint32_t old_df2 = df2_bytes;
        uint32_t old_post_budget = post_budget;
        int old_stride = stride;
        int old_bin_shift = bin_shift;
        bool old_dedup_seeds = dedup_seeds;

        ST.reprocess_attempts++;

        // ---- Rescue tuning (non-CLI): dense seeding and relaxed controls for rescue only ----
        const uint32_t DF_DISABLE = 1000000u;
        min_votes = 1;
        lead_margin = 0;
        skip_common = false; // pass all seeds
        sieve_seeds = max(80, old_sieve_seeds * 4);
        df1_bytes = DF_DISABLE;
        df2_bytes = DF_DISABLE;
        post_budget = 0;
        skip_rev_on_clear = false;
        stride = 1;                // dense seeding for rescue
        bin_shift = old_bin_shift; // keep original
        verify_topk = max(old_verify_topk, 200);
        verify_topk_lo = min(verify_topk, old_verify_topk_lo + 40);
        dedup_seeds = false;

        auto restore_state = [&]()
        {
            min_votes = old_min_votes;
            lead_margin = old_lead_margin;
            skip_common = old_skip_common;
            sieve_seeds = old_sieve_seeds;
            skip_rev_on_clear = old_skip_rev;
            verify_topk = old_verify_topk;
            verify_topk_lo = old_verify_topk_lo;
            df1_bytes = old_df1;
            df2_bytes = old_df2;
            post_budget = old_post_budget;
            stride = old_stride;
            bin_shift = old_bin_shift;
            dedup_seeds = old_dedup_seeds;
        };

        // Prepare fwd/rev packed bytes for verify attempts later
        string fwd = upper(seq);
        string rev = revcomp(fwd);
        auto f2b = pack2(fwd);
        auto r2b = pack2(rev);

        // ---- Variables to collect summary columns ----
        bool final_recovered = false;
        size_t kmers_total = 0, bloom_pass = 0, key_found = 0, df_pass = 0, common_count = 0, after_skip_common = 0;
        size_t s_sieve = 0, postings_hits = 0;
        int shortlist_size = 0, budgetVerify = 0;
        int best_tid = -1, best_pos = -1, best_mm = -1, best_score = INT_MIN;

        // Discover forward/rev captured values
        int discover_fwd_had_any_seed = 0, discover_fwd_best_votes = 0, discover_fwd_best_runner = 0, discover_fwd_shortlist_sz = 0;
        int discover_rev_had_any_seed = 0, discover_rev_best_votes = 0, discover_rev_best_runner = 0, discover_rev_shortlist_sz = 0;

        // ---- Local helper: try discover+verify and capture discover stats ----
        auto try_discover_verify = [&](const string &s, const vector<uint8_t> &s2b, bool is_rev) -> optional<Alignment>
        {
            auto R = discover(s, is_rev);
            if (R.all_seeds_filtered)
                ST.f_all_seeds_filtered++;

            // capture discover stats
            if (!is_rev)
            {
                discover_fwd_had_any_seed = R.had_any_seed ? 1 : 0;
                discover_fwd_best_votes = R.best_votes;
                discover_fwd_best_runner = R.best_runner;
                discover_fwd_shortlist_sz = (int)R.shortlist.size();
            }
            else
            {
                discover_rev_had_any_seed = R.had_any_seed ? 1 : 0;
                discover_rev_best_votes = R.best_votes;
                discover_rev_best_runner = R.best_runner;
                discover_rev_shortlist_sz = (int)R.shortlist.size();
            }

            if (R.shortlist.empty())
                return nullopt;

            int budget = choose_verify_budget(R.shortlist, R.best_votes, R.best_runner);
            bool has_clear = (R.best_votes >= (uint16_t)min_votes) && ((int)R.best_votes - (int)R.best_runner >= lead_margin);
            if (!has_clear)
                budget = min<int>((int)R.shortlist.size(), verify_topk);
            budget = max(1, min((int)R.shortlist.size(), budget));
            budgetVerify = budget;

            // relaxed mm fractions for rescue
            const double mm_frac_tx = 0.10;
            const double mm_frac_jx = 0.06;

            optional<Alignment> bestA = nullopt;
            for (int i = 0; i < budget; i++)
            {
                auto a = verify_one(s, s2b, R.shortlist[i], is_rev, mm_frac_tx, mm_frac_jx);
                if (!a)
                {
                    std::ostringstream ss;
                    continue;
                }
                Alignment A = *a;
                std::ostringstream ss_ok;

                if (!bestA.has_value() || (A.score > bestA->score) || (A.score == bestA->score && A.mm < bestA->mm))
                {
                    bestA = A;
                }
                if (A.mm == 0)
                    return bestA; // immediate perfect
            }
            return bestA;
        };

        // Try forward then reverse quickly with current rescue settings
        if (auto Aopt = try_discover_verify(fwd, f2b, false); Aopt.has_value())
        {
            Alignment A = *Aopt;
            GenomicPos gp = transcript_to_genomic(A.tid, A.pos, (int)seq.size());
            out.recovered = true;
            out.aln = A;
            out.tier_name = "rescue_discover_forward";
            out.notes = "discover_verify";
            final_recovered = true;

            ST.reprocess_recovered++;

            restore_state();
            return out;
        }

        if (auto Aopt = try_discover_verify(rev, r2b, true); Aopt.has_value())
        {
            Alignment A = *Aopt;
            GenomicPos gp = transcript_to_genomic(A.tid, A.pos, (int)seq.size());
            out.recovered = true;
            out.aln = A;
            out.tier_name = "rescue_discover_reverse";
            out.notes = "discover_verify_rev";
            final_recovered = true;

            best_tid = (int)A.tid;
            best_pos = A.pos;
            best_mm = A.mm;
            best_score = A.score;

            ST.reprocess_recovered++;

            restore_state();
            return out;
        }

        // 1) Generate dense kmers and instrument filtering:
        struct SeedDiag
        {
            uint64_t kmer;
            int pos;
            bool bloom_ok;
            bool found_key;
            bool is_common;
            uint32_t key_idx;
            uint32_t byte_span;
            bool df_pass;
        };
        vector<SeedDiag> seeds;
        seeds.reserve(max(1, (int)fwd.size() - (int)IX.k + 1));
        const int K = (int)IX.k;
        uint64_t KMASK = (K == 32) ? 0xFFFFFFFFFFFFFFFFULL : ((1ULL << (2 * K)) - 1ULL);
        uint64_t roll = 0;
        int good = 0;
        unordered_set<uint64_t> seen_kmers;
        for (size_t i = 0; i < fwd.size(); ++i)
        {
            int b = nt2b(fwd[i]);
            if (b < 0)
            {
                roll = 0;
                good = 0;
                continue;
            }
            roll = ((roll << 2) | (uint64_t)b) & KMASK;
            if (++good >= K)
            {
                int start = (int)i - K + 1;
                if ((start % max(1, stride)) != 0)
                    continue;
                uint64_t kmer = roll;
                if (IX.uses_canonical)
                {
                    uint64_t rc = 0, kx = kmer;
                    for (int j = 0; j < K; j++)
                    {
                        uint64_t base = kx & 3ULL;
                        rc = (rc << 2) | (3ULL - base);
                        kx >>= 2;
                    }
                    if (rc < kmer)
                        kmer = rc;
                }
                if (dedup_seeds)
                {
                    if (!seen_kmers.insert(kmer).second)
                        continue;
                }
                SeedDiag sd;
                sd.kmer = kmer;
                sd.pos = start;
                sd.bloom_ok = IX.bloom_may_contain(kmer);
                auto [ok, idx] = IX.find_key(kmer);
                sd.found_key = ok;
                sd.key_idx = ok ? idx : UINT32_MAX;
                sd.byte_span = ok ? IX.postings_span_bytes(idx) : 0;
                sd.df_pass = ok ? (sd.byte_span <= df2_bytes) : false;
                sd.is_common = ok ? IX.is_common(sd.key_idx) : false;
                seeds.push_back(sd);
            }
        }

        // ---- Sieve stage: sample first n seeds or up to sieve_seeds ----
        s_sieve = (size_t)min<int>((int)seeds.size(), max(1, sieve_seeds));
        postings_hits = 0;
        unordered_set<uint32_t> scanned_key_idxs;
        for (size_t i = 0; i < s_sieve; i++)
        {
            const auto &sd = seeds[i];
            if (!sd.found_key)
                continue;
            uint32_t kidx = sd.key_idx;
            if (scanned_key_idxs.insert(kidx).second)
            {
                uint32_t cnt = 0;
                IX.for_each_posting_breakable(kidx, [&](uint32_t tid, uint32_t pos) -> bool
                                              {
                cnt++; return true; }, false);
                postings_hits += cnt;
            }
        }

        // Now perform postings-driven candidate aggregation
        next_read_epoch();
        unordered_set<uint32_t> used_keyidxs;
        for (const auto &sd : seeds)
        {
            if (!sd.found_key)
                continue;
            uint32_t kid = sd.key_idx;
            if (!used_keyidxs.insert(kid).second)
                continue;
            IX.for_each_posting_breakable(kid, [&](uint32_t tid, uint32_t pos) -> bool
                                          {
            int32_t off_q = (int32_t)pos - (int32_t)sd.pos;
            off_q = off_q >> bin_shift;
            uint8_t flags = 0;
            if(IX.is_jx_tid(tid)){
                uint32_t jxid = IX.jx_index(tid);
                const auto &J = IX.jx[jxid];
                uint32_t mid = J.len_bp >> 1;
                flags |= (pos<mid)?1:2;
            }
            uint16_t c = bump_hh(tid, off_q, flags);
            on_vote(tid, off_q, c);
            return true; }, prefetch);
        }

        // collect shortlist
        struct CandLocal
        {
            uint32_t tid;
            int32_t off_q;
            uint16_t votes;
            uint8_t flags;
        };
        vector<CandLocal> tmp;
        tmp.reserve(hh_touched.size());
        for (uint32_t tid : hh_touched)
        {
            auto &tb = per_tid[tid];
            if (tb.gen != per_tid_gen)
                continue;
            auto &H = hh[tid];
            int win = -1;
            for (int i = 0; i < 4; i++)
            {
                if (H.cnt[i] == 0)
                    continue;
                if (H.off_q[i] == tb.best_off_q && H.cnt[i] == tb.best)
                {
                    win = i;
                    break;
                }
            }
            if (win == -1)
            {
                for (int i = 0; i < 4; i++)
                    if (H.cnt[i] == tb.best)
                    {
                        win = i;
                        break;
                    }
            }
            if (win != -1)
                tmp.push_back({tid, H.off_q[win], H.cnt[win], H.flg[win]});
        }

        tmp.erase(remove_if(tmp.begin(), tmp.end(), [&](const CandLocal &c)
                            {
        if(!IX.is_jx_tid(c.tid)) return false;
        return ((c.flags & 3) != 3); }),
                  tmp.end());

        sort(tmp.begin(), tmp.end(), [](const CandLocal &a, const CandLocal &b)
             {
        if(a.votes != b.votes) return a.votes > b.votes;
        if(a.tid != b.tid) return a.tid < b.tid;
        return a.off_q < b.off_q; });

        shortlist_size = (int)tmp.size();

        // Convert to Cand and verify top candidates
        vector<Cand> shortlist;
        shortlist.reserve(tmp.size());
        for (const auto &c : tmp)
            shortlist.push_back({c.tid, c.off_q, c.votes, c.flags});

        budgetVerify = min<int>((int)shortlist.size(), verify_topk);
        const double mm_frac_tx = 0.10;
        const double mm_frac_jx = 0.06;
        optional<Alignment> bestA = nullopt;
        for (int i = 0; i < budgetVerify; i++)
        {
            auto a = verify_one(fwd, f2b, shortlist[i], false, mm_frac_tx, mm_frac_jx);
            if (!a)
            {
                std::ostringstream ss;
                continue;
            }
            Alignment A = *a;
            {
                std::ostringstream ss;
            }
            if (!bestA.has_value() || (A.score > bestA->score) || (A.score == bestA->score && A.mm < bestA->mm))
            {
                bestA = A;
            }
            if (A.mm == 0)
                break;
        }

        if (bestA.has_value())
        {
            Alignment A = *bestA;
            GenomicPos gp = transcript_to_genomic(A.tid, A.pos, (int)seq.size());
            out.recovered = true;
            out.aln = A;
            out.tier_name = "rescue_postings_verify";
            out.notes = "fallback_verify";
            final_recovered = true;

            best_tid = (int)A.tid;
            best_pos = A.pos;
            best_mm = A.mm;
            best_score = A.score;

            ST.reprocess_recovered++;

            restore_state();
            return out;
        }

        ST.reprocess_failed++;

        restore_state();
        return out;
    }

    // Add this inside the Mapper struct, before map_paired()

    bool attempt_mate_rescue(
        vector<Alignment> &H1, Alignment *&best1,
        vector<Alignment> &H2, Alignment *&best2,
        const ReadRecord &rec1, const ReadRecord &rec2)
    {
        if (!best1 || !best2)
            return false;

        GenomicPos gp1 = transcript_to_genomic(best1->tid, best1->pos, (int)rec1.seq.size());
        GenomicPos gp2 = transcript_to_genomic(best2->tid, best2->pos, (int)rec2.seq.size());

        if (!gp1.valid || !gp2.valid)
            return false;

        // If already concordant, nothing to rescue
        if (gp1.chrom == gp2.chrom)
        {
            int ins = abs((int)gp1.pos - (int)gp2.pos);
            if (ins >= min_ins && ins <= max_ins)
            {
                return false;
            }
        }

        // Discordant pair - find best rescue option
        struct RescueOption
        {
            Alignment *new_best1;
            Alignment *new_best2;
            int insert_size;
            int score_penalty; // How much score we lose by switching
            double quality;    // Overall quality score
        };

        vector<RescueOption> options;

        // Option 1: Keep R1, try alternatives for R2
        for (size_t i = 0; i < H2.size() && i < 10; i++)
        {
            GenomicPos alt_gp2 = transcript_to_genomic(H2[i].tid, H2[i].pos, (int)rec2.seq.size());
            if (alt_gp2.valid && alt_gp2.chrom == gp1.chrom)
            {
                int ins = abs((int)gp1.pos - (int)alt_gp2.pos);
                if (ins <= max_ins * 2)
                {
                    int penalty = best2->score - H2[i].score;
                    // Quality = high score - low penalty - bad insert size deviation
                    int ideal_ins = 400; // Expected insert size for RNA-seq
                    int ins_dev = abs(ins - ideal_ins);
                    double quality = 100.0 - penalty - (ins_dev / 10.0);

                    options.push_back({best1, &H2[i], ins, penalty, quality});
                }
            }
        }

        // Option 2: Keep R2, try alternatives for R1
        for (size_t i = 0; i < H1.size() && i < 10; i++)
        {
            GenomicPos alt_gp1 = transcript_to_genomic(H1[i].tid, H1[i].pos, (int)rec1.seq.size());
            if (alt_gp1.valid && alt_gp1.chrom == gp2.chrom)
            {
                int ins = abs((int)alt_gp1.pos - (int)gp2.pos);
                if (ins <= max_ins * 2)
                {
                    int penalty = best1->score - H1[i].score;
                    int ideal_ins = 400;
                    int ins_dev = abs(ins - ideal_ins);
                    double quality = 100.0 - penalty - (ins_dev / 10.0);

                    options.push_back({&H1[i], best2, ins, penalty, quality});
                }
            }
        }

        if (options.empty())
            return false;

        // Pick the best option based on quality score
        sort(options.begin(), options.end(),
             [](const RescueOption &a, const RescueOption &b)
             {
                 return a.quality > b.quality;
             });

        const auto &best_option = options[0];

        // Only accept rescue if:
        // 1. Insert size is reasonable (< 1500)
        // 2. Score penalty is acceptable (< 10)
        // 3. Quality is positive
        if (best_option.insert_size <= max_ins &&
            best_option.score_penalty <= 10 &&
            best_option.quality > 0)
        {
            best1 = best_option.new_best1;
            best2 = best_option.new_best2;
            return true;
        }

        return false;
    }

    // ---------- mapping loops (paired and single) ----------
    void map_paired(const string &r1, const string &r2, const string &out_path, uint64_t max_pairs = 0, bool do_rescue = false)
    {
        FastqReader fq1(r1), fq2(r2);
        OutputBuffer outbuf(out_path);
        write_sam_header(outbuf, IX);

        const size_t BATCH_SIZE = 50000;
        vector<ReadRecord> batch1, batch2;
        uint64_t total = 0;

        // temporary containers for recovered alignments
        vector<Alignment> temp_alns1, temp_alns2;
        temp_alns1.reserve(8);
        temp_alns2.reserve(8);

        while (true)
        {
            size_t n1 = fq1.next_batch(batch1, BATCH_SIZE);
            size_t n2 = fq2.next_batch(batch2, BATCH_SIZE);
            if (n1 == 0 || n2 == 0)
                break;
            if (n1 != n2)
            {
                fprintf(stderr, "ERROR: FASTQ mismatch\n");
                exit(1);
            }

            for (size_t idx = 0; idx < n1; ++idx)
            {
                total++;
                ST.total_pairs++;

                if (max_pairs && total > max_pairs)
                    goto finish;

                const auto &rec1 = batch1[idx];
                const auto &rec2 = batch2[idx];

                auto H1 = map_read_single(rec1.id, rec1.seq);
                auto H2 = map_read_single(rec2.id, rec2.seq);

                // ===== USE PAIR SCORING TO SELECT BEST PAIR =====
                Alignment *best1 = nullptr;
                Alignment *best2 = nullptr;

                if (!H1.empty() && !H2.empty())
                {
                    // Default to independent best
                    best1 = &H1[0];
                    best2 = &H2[0];

                    // Only use pair scoring if:
                    // 1. There are alternatives, AND
                    // 2. Top hits are close in score, OR
                    // 3. Mates are discordant

                    bool has_alternatives = (H1.size() > 1 || H2.size() > 1);

                    if (has_alternatives)
                    {
                        GenomicPos gp1 = transcript_to_genomic(best1->tid, best1->pos, (int)rec1.seq.size());
                        GenomicPos gp2 = transcript_to_genomic(best2->tid, best2->pos, (int)rec2.seq.size());

                        bool discordant = (gp1.valid && gp2.valid && gp1.chrom != gp2.chrom);

                        // Check if scores are close
                        bool r1_close = (H1.size() > 1 && H1[0].score - H1[1].score <= 10);
                        bool r2_close = (H2.size() > 1 && H2[0].score - H2[1].score <= 10);

                        // Only apply pair scoring if needed
                        if (discordant || r1_close || r2_close)
                        {
                            auto best_pair = select_best_pair(H1, H2, rec1, rec2);

                            if (best_pair.found)
                            {
                                // But only accept if pair score is CLEARLY better
                                int independent_score = H1[0].score + H2[0].score;
                                int pair_score = best_pair.score.total - best_pair.score.concordance_bonus; // Exclude bonus

                                // Accept pair if it's not much worse in pure alignment quality
                                if (pair_score >= independent_score - 10)
                                {
                                    best1 = best_pair.best1;
                                    best2 = best_pair.best2;
                                }
                            }
                        }
                    }
                }
                else
                {
                    // One or both unmapped
                    best1 = H1.empty() ? nullptr : &H1[0];
                    best2 = H2.empty() ? nullptr : &H2[0];
                }
                // ===== END PAIR SCORING =====
                // ===== TRACK PRE-RESCUE STATE =====
                bool rescue_attempted = false;
                bool rescue_succeeded = false;
                string r1_chr_pre = "*", r2_chr_pre = "*";
                int r1_pos_pre = -1, r2_pos_pre = -1;

                if (best1)
                {
                    GenomicPos gp1_pre = transcript_to_genomic(best1->tid, best1->pos, (int)rec1.seq.size());
                    if (gp1_pre.valid)
                    {
                        r1_chr_pre = gp1_pre.chrom;
                        r1_pos_pre = gp1_pre.pos;
                    }
                }
                if (best2)
                {
                    GenomicPos gp2_pre = transcript_to_genomic(best2->tid, best2->pos, (int)rec2.seq.size());
                    if (gp2_pre.valid)
                    {
                        r2_chr_pre = gp2_pre.chrom;
                        r2_pos_pre = gp2_pre.pos;
                    }
                }

                // ===== ATTEMPT MATE RESCUE =====
                if (best1 && best2 && r1_chr_pre != "*" && r2_chr_pre != "*")
                {
                    if (r1_chr_pre != r2_chr_pre)
                    {
                        rescue_attempted = true;
                        rescue_succeeded = attempt_mate_rescue(H1, best1, H2, best2, rec1, rec2);
                    }
                }
                // ===== END MATE RESCUE =====

                // Production (no-truth) behaviour for read1: do not consult any truth DB.
                // If user enabled rescue (do_rescue) and the read is unmapped, attempt rescue.
                // if (do_rescue && !best1)
                // Production (no-truth) behaviour for read1: do not consult any truth DB.
                // If user enabled rescue (do_rescue) and the read is unmapped, attempt rescue.
                if (do_rescue && !best1)
                {
                    ST.reprocess_attempts++;
                    auto rres = attempt_reprocess_single(rec1.id, rec1.seq);
                    if (rres.recovered)
                    {
                        ST.reprocess_recovered++;
                        temp_alns1.clear();
                        temp_alns1.push_back(rres.aln);
                        best1 = &temp_alns1[0];
                        GenomicPos gp = transcript_to_genomic(best1->tid, best1->pos, (int)rec1.seq.size());
                    }
                    else
                    {
                        ST.reprocess_failed++;
                    }
                }

                // Production (no-truth) behaviour for read2: same logic
                if (do_rescue && !best2)
                {
                    ST.reprocess_attempts++;
                    auto rres = attempt_reprocess_single(rec2.id, rec2.seq);
                    if (rres.recovered)
                    {
                        ST.reprocess_recovered++;
                        temp_alns2.clear();
                        temp_alns2.push_back(rres.aln);
                        best2 = &temp_alns2[0];
                        GenomicPos gp = transcript_to_genomic(best2->tid, best2->pos, (int)rec2.seq.size());
                    }
                    else
                    {
                        ST.reprocess_failed++;
                    }
                }

                bool conc = false;
                if (best1 && best2)
                {
                    if (tid_compatible_same_gene(best1->tid, best2->tid))
                    {
                        if (best1->tid == best2->tid)
                        {
                            int ins = (best1->pos + (int)rec1.seq.size()) - best2->pos;
                            if (ins < 0)
                                ins = -ins;
                            if (ins >= min_ins && ins <= max_ins)
                            {
                                conc = true;
                            }
                        }
                        else
                        {
                            conc = true;
                        }
                    }
                }

                if (!best1 && !best2)
                {
                    ST.f_pair_both_unmapped++;
                }
                else if (!best1 || !best2)
                {
                    ST.f_pair_one_unmapped++;
                }
                else
                {
                    if (conc)
                    {
                        ST.ok_pairs++;
                        ST.ok_mapped_reads += 2;
                    }
                    else
                    {
                        ST.f_pair_incompatible++;
                        ST.ok_mapped_reads += 2;
                    }
                }

                write_sam_record(outbuf, rec1, best1, true, false, best2, &rec2, &H1);
                write_sam_record(outbuf, rec2, best2, true, true, best1, &rec1, &H2);
            }
        }

    finish:
        ;
    }

    void map_single_end(const string &reads, const string &out_path, uint64_t max_reads = 0, bool do_rescue = false)
    {

        FastqReader fq(reads);
        OutputBuffer outbuf(out_path);
        write_sam_header(outbuf, IX);

        uint64_t total = 0;
        const size_t BATCH_SIZE = 100000;
        vector<ReadRecord> batch;
        vector<Alignment> temp_alns;
        temp_alns.reserve(8);

        while (true)
        {
            size_t n = fq.next_batch(batch, BATCH_SIZE);
            if (n == 0)
                break;
            for (size_t idx = 0; idx < n; ++idx)
            {
                total++;
                if (max_reads && total > max_reads)
                    break;

                const auto &rec = batch[idx];
                auto H = map_read_single(rec.id, rec.seq);
                Alignment *best = H.empty() ? nullptr : &H[0];

                // Optional rescue only.
                if (do_rescue && !best)
                {
                    ST.reprocess_attempts++;
                    auto rres = attempt_reprocess_single(rec.id, rec.seq);
                    if (rres.recovered)
                    {
                        ST.reprocess_recovered++;
                        temp_alns.clear();
                        temp_alns.push_back(rres.aln);
                        best = &temp_alns[0];
                        GenomicPos gp = transcript_to_genomic(best->tid, best->pos, (int)rec.seq.size());
                    }
                    else
                    {
                        ST.reprocess_failed++;
                    }
                }

                if (best)
                    ST.ok_mapped_reads++;
                write_sam_record(outbuf, rec, best, false, false, nullptr, nullptr, &H);
            }
            if (max_reads && total >= max_reads)
                break;
        }
    }
};

// ======================= MAIN ======================================
int main(int argc, char **argv)
{
    if (argc < 4)
    {
        fprintf(stderr, "Usage:\n");
        fprintf(stderr, "  Paired: %s <index.bin> <R1.fq[.gz]> <R2.fq[.gz]> <out.sam> [options]\n", argv[0]);
        fprintf(stderr, "  Single: %s <index.bin> <reads.fq[.gz]> <out.sam> [options]\n", argv[0]);
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  --max-pairs N        limit paired mapping to N pairs\n");
        fprintf(stderr, "  --max N              alias for --max-pairs\n");
        fprintf(stderr, "  --stride N           seed stride\n");
        fprintf(stderr, "  --df1-bytes N        early DF cutoff\n");
        fprintf(stderr, "  --df2-bytes N        late DF cutoff\n");
        fprintf(stderr, "  --min-votes N        seed votes to pass\n");
        fprintf(stderr, "  --lead-margin N      lead over runner\n");
        fprintf(stderr, "  --early-seeds N      min seeds before early exit\n");
        fprintf(stderr, "  --verify-topk N      max verifies\n");
        fprintf(stderr, "  --verify-topk-lo N   verify when clear\n");
        fprintf(stderr, "  --output-nm-md 1     emit NM/MD\n");
        fprintf(stderr, "  --rescue             enable rescue/reprocess tiers for unmapped reads\n");
        return 1;
    }

    auto looks_like_fastq = [](const string &s) -> bool
    {
        return ends_with(s, ".fq") || ends_with(s, ".fastq") ||
               ends_with(s, ".fq.gz") || ends_with(s, ".fastq.gz");
    };

    string idx_path = argv[1];
    string arg2 = argv[2];
    string arg3 = argv[3];

    bool paired = false;
    string reads1 = arg2, reads2, out_path;
    int opt_start = 0;

    if (argc >= 5 && looks_like_fastq(arg3))
    {
        paired = true;
        reads2 = arg3;
        out_path = argv[4];
        opt_start = 5;
    }
    else
    {
        paired = false;
        out_path = arg3;
        opt_start = 4;
    }

    uint64_t max_pairs = 0;
    uint64_t max_reads = 0;
    int min_ins = 100, max_ins = 800, stride = 1, min_votes = 3, lead_margin = 2, early_seeds = 6;
    int verify_topk = 12, verify_topk_lo = 5, sieve_seeds = 10, bin_shift = 3;
    uint32_t post_budget = 0, df1_bytes = 128, df2_bytes = 512;
    bool prefetch = true, dedup_seeds = false, skip_rev_on_clear = true;
    bool output_nm_md = false;
    bool rescue_unmapped = false;

    for (int i = opt_start; i < argc; i++)
    {
        if (strcmp(argv[i], "--max-pairs") == 0 && i + 1 < argc)
        {
            max_pairs = parse_u64(argv[++i], 0);
        }
        else if (strcmp(argv[i], "--max") == 0 && i + 1 < argc)
        {
            max_pairs = parse_u64(argv[++i], 0);
        }
        else if (strcmp(argv[i], "--max-reads") == 0 && i + 1 < argc)
        {
            max_reads = parse_u64(argv[++i], 0);
        }
        else if (strcmp(argv[i], "--min-ins") == 0 && i + 1 < argc)
        {
            min_ins = parse_i32(argv[++i], 100);
        }
        else if (strcmp(argv[i], "--max-ins") == 0 && i + 1 < argc)
        {
            max_ins = parse_i32(argv[++i], 800);
        }
        else if (strcmp(argv[i], "--stride") == 0 && i + 1 < argc)
        {
            stride = parse_i32(argv[++i], 1);
            if (stride < 1)
                stride = 1;
        }
        else if (strcmp(argv[i], "--df1-bytes") == 0 && i + 1 < argc)
        {
            df1_bytes = parse_u32(argv[++i], 128);
        }
        else if (strcmp(argv[i], "--df2-bytes") == 0 && i + 1 < argc)
        {
            df2_bytes = parse_u32(argv[++i], 512);
        }
        else if (strcmp(argv[i], "--post-budget") == 0 && i + 1 < argc)
        {
            post_budget = parse_u32(argv[++i], 0);
        }
        else if (strcmp(argv[i], "--min-votes") == 0 && i + 1 < argc)
        {
            min_votes = parse_i32(argv[++i], 3);
        }
        else if (strcmp(argv[i], "--lead-margin") == 0 && i + 1 < argc)
        {
            lead_margin = parse_i32(argv[++i], 2);
        }
        else if (strcmp(argv[i], "--early-seeds") == 0 && i + 1 < argc)
        {
            early_seeds = parse_i32(argv[++i], 6);
        }
        else if (strcmp(argv[i], "--verify-topk") == 0 && i + 1 < argc)
        {
            verify_topk = parse_i32(argv[++i], 12);
        }
        else if (strcmp(argv[i], "--verify-topk-lo") == 0 && i + 1 < argc)
        {
            verify_topk_lo = parse_i32(argv[++i], 5);
        }
        else if (strcmp(argv[i], "--prefetch") == 0 && i + 1 < argc)
        {
            prefetch = parse_i32(argv[++i], 1) != 0;
        }
        else if (strcmp(argv[i], "--dedup-seeds") == 0 && i + 1 < argc)
        {
            dedup_seeds = parse_i32(argv[++i], 0) != 0;
        }
        else if (strcmp(argv[i], "--sieve-seeds") == 0 && i + 1 < argc)
        {
            sieve_seeds = parse_i32(argv[++i], 10);
        }
        else if (strcmp(argv[i], "--skip-rev-on-clear") == 0 && i + 1 < argc)
        {
            skip_rev_on_clear = parse_i32(argv[++i], 1) != 0;
        }
        else if (strcmp(argv[i], "--bin-shift") == 0 && i + 1 < argc)
        {
            bin_shift = parse_i32(argv[++i], 3);
        }
        else if (strcmp(argv[i], "--output-nm-md") == 0 && i + 1 < argc)
        {
            output_nm_md = parse_i32(argv[++i], 0) != 0;
        }
        else if (strcmp(argv[i], "--rescue") == 0)
        {
            rescue_unmapped = true;
        }
    }

    IndexVX IX;
    IX.load(idx_path);

    FailureStats ST;
    Mapper mp(IX, ST);
    mp.min_ins = min_ins;
    mp.max_ins = max_ins;
    mp.stride = stride;
    mp.min_votes = min_votes;
    mp.lead_margin = lead_margin;
    mp.min_seeds_before_early = early_seeds;
    mp.verify_topk = verify_topk;
    mp.verify_topk_lo = verify_topk_lo;
    mp.df1_bytes = df1_bytes;
    mp.df2_bytes = df2_bytes;
    mp.post_budget = post_budget;
    mp.prefetch = prefetch;
    mp.dedup_seeds = dedup_seeds;
    mp.sieve_seeds = sieve_seeds;
    mp.skip_rev_on_clear = skip_rev_on_clear;
    mp.bin_shift = bin_shift;
    mp.output_nm_md = output_nm_md;

    if (!paired)
    {
        mp.map_single_end(reads1, out_path, max_reads, rescue_unmapped);
    }
    else
    {
        mp.map_paired(reads1, reads2, out_path, max_pairs, rescue_unmapped);
    }

    return 0;
}