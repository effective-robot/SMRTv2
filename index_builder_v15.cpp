// index_builder_v15.cpp
// v15 = v14 + CANONICAL K-MERS + SKIP EXONIC JX K-MERS + COMPACT DF=1 ENCODING
// Major improvements:
// 1. Canonical k-mer indexing (50% fewer k-mers)
// 2. Skip redundant exonic k-mers in junctions (cleaner junction evidence)
// 3. Compact encoding for singleton k-mers (better compression)

#include <bits/stdc++.h>
#include <zlib.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

// ---------- utility ----------
static inline bool ends_with(const string& s, const string& suf){
    return s.size()>=suf.size() && equal(suf.rbegin(),suf.rend(),s.rbegin());
}
static inline string upper(const string& s){
    string o(s);
    for(char& c: o){
        unsigned char u = toupper((unsigned char)c);
        c = (u=='U')?'T':u;
    }
    return o;
}
static inline string revcomp(const string& s){
    auto rc=[](char c)->char{
        switch(c){ case 'A':return 'T'; case 'C':return 'G'; case 'G':return 'C'; case 'T':return 'A'; default:return 'N'; }
    };
    string out(s.size(),'N');
    for(size_t i=0;i<s.size();++i) out[i]=rc(s[s.size()-1-i]);
    return out;
}
static inline int nt2b(char c){
    switch(c){
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}
static inline bool encode_kmer64(const char* s, int k, uint64_t& out){
    uint64_t v=0;
    for(int i=0;i<k;i++){
        int b=nt2b(s[i]);
        if(b<0) return false;
        v=(v<<2)|(uint64_t)b;
    }
    out=v; return true;
}

// NEW: Reverse complement k-mer
static inline uint64_t revcomp_kmer(uint64_t kmer, int k){
    uint64_t rc = 0;
    for(int i=0; i<k; ++i){
        uint64_t base = kmer & 3ULL;
        rc = (rc << 2) | (3ULL - base);  // A<->T, C<->G
        kmer >>= 2;
    }
    return rc;
}

// NEW: Canonical k-mer (lexicographically smaller)
static inline uint64_t canonical_kmer(uint64_t kmer, int k){
    uint64_t rc = revcomp_kmer(kmer, k);
    return min(kmer, rc);
}

static inline vector<uint8_t> pack2bit(const string& seq){
    vector<uint8_t> out((seq.size()+3)/4,0);
    size_t o=0; int sh=6;
    for(char c: seq){
        int b=nt2b(c); if(b<0) b=0;
        out[o] |= (uint8_t)(b<<sh);
        sh-=2;
        if(sh<0){ sh=6; ++o; }
    }
    return out;
}
template<typename T> static inline void write_pod(ofstream& out, const T& v){
    out.write((const char*)&v, sizeof(T));
}
static inline void write_str(ofstream& out, const string& s){
    uint32_t n=(uint32_t)s.size();
    write_pod(out, n);
    if(n) out.write(s.data(), n);
}

// ---------- file helpers ----------
struct LineReader{
    bool gz; FILE* fp; gzFile gzfp; vector<char> buf;
    explicit LineReader(const string& path, size_t bufsz=1<<18)
        :gz(ends_with(path,".gz")),fp(nullptr),gzfp(nullptr),buf(bufsz)
    {
        if(gz){
            gzfp=gzopen(path.c_str(),"rb");
            if(!gzfp){ fprintf(stderr,"ERROR: open %s\n", path.c_str()); exit(1); }
            gzbuffer(gzfp, (unsigned int)bufsz);
        }else{
            fp=fopen(path.c_str(),"rb");
            if(!fp){ fprintf(stderr,"ERROR: open %s\n", path.c_str()); exit(1); }
            setvbuf(fp, nullptr, _IOFBF, bufsz);
        }
    }
    bool getline(string& out){
        out.clear();
        if(gz){
            for(;;){
                char* r=gzgets(gzfp, buf.data(), (int)buf.size());
                if(!r) return !out.empty();
                size_t len=strlen(r);
                if(len && r[len-1]=='\n'){ out.append(r, len-1); return true; }
                out.append(r, len);
                if(len < buf.size()-1) return true;
            }
        }else{
            int ch;
            while((ch=fgetc(fp))!=EOF){
                if(ch=='\n') return true;
                out.push_back((char)ch);
            }
            return !out.empty();
        }
    }
    ~LineReader(){ if(gz&&gzfp)gzclose(gzfp); if(!gz&&fp)fclose(fp); }
};

static void ensure_dir(const string& d){
    if(d.empty()) return;
    struct stat st{};
    if(stat(d.c_str(), &st)==0){
        if(S_ISDIR(st.st_mode)) return;
    }
#if defined(_WIN32)
    _mkdir(d.c_str());
#else
    mkdir(d.c_str(), 0775);
#endif
}

// ---------- annotations ----------
struct Exon { uint32_t start, end; };
struct TranscriptInfo{
    string id, chrom;
    char strand='+';
    vector<Exon> exons;
    string seq;
};
struct JunctionInfo{
    string id, gene_id, chrom;
    char strand='+';
    uint32_t donor_pos=0, acceptor_pos=0;
    string seq;
    vector<string> source_transcripts;
};

// ---------- config ----------
struct Config{
    int k=31;
    int flank=25;
    int prefix_bits=20;
    size_t chunk_mb=256;
    string tmpdir="./tmp";
    uint32_t df_cap=0;
    uint32_t min_junction_size=20;
    uint32_t junction_flank_overlap=15;  // NEW: bp overlap to check for exonic k-mers
    string fasta, gtf, outprefix;
};

// ---------- metadata ----------
struct TargetMetaTX{
    string id, chrom;
    char strand;
    uint32_t len_bp;
    uint64_t seq_off;
    vector<pair<uint32_t, uint32_t>> genomic_exons;
};

struct TargetMetaJX{
    string id, gene_id, chrom;
    char strand;
    uint32_t donor_pos, acceptor_pos, len_bp;
    uint64_t seq_off;
};

// ---------- postings ----------
#pragma pack(push,1)
struct PostRec{ 
    uint64_t kmer;  // Already canonical
    uint32_t tid; 
    uint32_t pos; 
};
#pragma pack(pop)

struct RunWriter{
    FILE* f=nullptr;
    string path;
    vector<PostRec> buf;
    size_t bytes_limit;
    explicit RunWriter(const string& p, size_t lim):path(p),bytes_limit(lim){}
    void open(){
        ensure_dir(std::filesystem::path(path).parent_path().string());
        f=fopen(path.c_str(),"wb");
        if(!f){ fprintf(stderr,"ERROR: cannot create %s\n", path.c_str()); exit(1); }
        setvbuf(f, nullptr, _IOFBF, 1<<20);
    }
    void add(const PostRec& r){ buf.push_back(r); }
    size_t buffered_bytes() const { return buf.size()*sizeof(PostRec); }
    bool has_data() const { return !buf.empty(); }
    void flush(){
        if(buf.empty()) return;
        sort(buf.begin(), buf.end(), [](const PostRec& a, const PostRec& b){
            if(a.kmer!=b.kmer) return a.kmer < b.kmer;
            if(a.tid != b.tid) return a.tid < b.tid;
            return a.pos < b.pos;
        });
        if(!f) open();
        fwrite(buf.data(), sizeof(PostRec), buf.size(), f);
        buf.clear();
    }
    void close(){
        flush();
        if(f){ fclose(f); f=nullptr; }
    }
};

struct RunReader{
    FILE* f=nullptr;
    string path;
    vector<PostRec> buf;
    size_t i=0, n=0;
    explicit RunReader(const string& p):path(p){}
    void open(){
        f=fopen(path.c_str(),"rb");
        if(!f){ fprintf(stderr,"ERROR: cannot open %s\n", path.c_str()); exit(1); }
        setvbuf(f, nullptr, _IOFBF, 1<<20);
        refill();
    }
    void refill(){
        buf.resize(1<<16);
        n=fread(buf.data(), sizeof(PostRec), buf.size(), f);
        i=0;
    }
    bool next(PostRec& out){
        if(i>=n){
            refill();
            if(n==0) return false;
        }
        out=buf[i++];
        return true;
    }
    ~RunReader(){ if(f) fclose(f); }
};

// ---------- builder ----------
class IndexBuilderV15{
public:
    Config cfg;

    unordered_map<string,string> genome;
    vector<TranscriptInfo> transcripts;
    vector<JunctionInfo> junctions;

    vector<TargetMetaTX> tx_meta;
    vector<TargetMetaJX> jx_meta;
    vector<uint8_t> seq_blob;

    vector<string> run_files;
    
    // NEW: Track transcript k-mers for junction filtering
    unordered_map<uint64_t, vector<uint32_t>> transcript_kmer_map;  // canonical_kmer -> [tid1, tid2, ...]

    void parse_fasta(const string& path){
        fprintf(stderr,"Reading genome from %s...\n", path.c_str());
        LineReader lr(path);
        string line, cur_id, seqbuf;
        size_t n=0;
        while(lr.getline(line)){
            if(!line.empty() && line[0]=='>'){
                if(!cur_id.empty()){
                    genome[cur_id]=upper(seqbuf);
                    seqbuf.clear(); n++;
                }
                size_t sp=line.find_first_of(" \t");
                cur_id = (sp==string::npos)? line.substr(1) : line.substr(1, sp-1);
            }else if(!line.empty()){
                seqbuf.append(upper(line));
            }
        }
        if(!cur_id.empty()){
            genome[cur_id]=upper(seqbuf);
            n++;
        }
        fprintf(stderr,"  Loaded %zu chromosomes/contigs\n", n);
    }

    void parse_gtf(const string& path){
        fprintf(stderr,"Parsing annotations from %s...\n", path.c_str());
        LineReader lr(path);
        string line;

        struct GTmp{
            vector<Exon> exons;
            vector<Exon> cds;
            char strand='+';
            string chrom;
        };

        unordered_map<string,GTmp> tmp;

        while(lr.getline(line)){
            if(line.empty() || line[0]=='#') continue;
            istringstream iss(line);
            string chrom, src, feat, sstr, send, sc, strand, frame, attrs;
            if(!(iss>>chrom>>src>>feat>>sstr>>send>>sc>>strand>>frame)) continue;
            getline(iss, attrs);
            if(feat!="CDS" && feat!="exon") continue;

            string tid;
            size_t p=attrs.find("transcript_id");
            if(p!=string::npos){
                size_t q1=attrs.find('"',p);
                size_t q2=(q1==string::npos)?string::npos:attrs.find('"',q1+1);
                if(q1!=string::npos && q2!=string::npos) tid=attrs.substr(q1+1, q2-q1-1);
            }
            if(tid.empty()) continue;

            uint32_t st = (uint32_t)max(0, stoi(sstr)-1);
            uint32_t en = (uint32_t)stoi(send);
            char stn = strand.empty()?'+':strand[0];

            auto& g=tmp[tid];
            if(feat=="exon") g.exons.push_back({st,en});
            else             g.cds.push_back({st,en});
            g.strand=stn; g.chrom=chrom;
        }

        transcripts.clear();
        transcripts.reserve(tmp.size());
        for(auto& kv: tmp){
            auto& g=kv.second;
            const vector<Exon>* spans=nullptr;
            if(!g.exons.empty()) spans=&g.exons;
            else if(!g.cds.empty()) spans=&g.cds;
            if(!spans) continue;
            if(!genome.count(g.chrom)) continue;

            vector<Exon> segs=*spans;
            sort(segs.begin(), segs.end(), [](const Exon& a, const Exon& b){ return a.start<b.start; });

            const string& chr=genome[g.chrom];
            string cat; cat.reserve(1024);
            bool ok=true;
            for(const auto& e: segs){
                if(e.end>chr.size()){ ok=false; break; }
                cat.append(chr.data()+e.start, e.end-e.start);
            }
            if(!ok || (int)cat.size()<cfg.k) continue;

            string final_seq = (g.strand=='-') ? revcomp(cat) : cat;
            if((int)final_seq.size()<cfg.k) continue;

            TranscriptInfo t;
            t.id=kv.first;
            t.chrom=g.chrom;
            t.strand=g.strand;
            t.exons=move(segs);
            t.seq=move(final_seq);
            transcripts.push_back(move(t));
        }
        fprintf(stderr,"  Extracted %zu transcripts\n", transcripts.size());
    }

    void build_junctions(){
        fprintf(stderr,"Building splice junctions (flank=%d, deduplication ON)...\n", cfg.flank);
        
        struct JunctionKey {
            string chrom;
            uint32_t donor, acceptor;
            char strand;
            
            bool operator==(const JunctionKey& o) const {
                return chrom==o.chrom && donor==o.donor && 
                       acceptor==o.acceptor && strand==o.strand;
            }
        };
        
        struct JunctionKeyHash {
            size_t operator()(const JunctionKey& k) const {
                return hash<string>{}(k.chrom) ^ 
                       (hash<uint32_t>{}(k.donor) << 1) ^
                       (hash<uint32_t>{}(k.acceptor) << 2) ^
                       (hash<char>{}(k.strand) << 3);
            }
        };
        
        unordered_map<JunctionKey, vector<string>, JunctionKeyHash> junction_map;
        
        size_t total_edges = 0;
        size_t skipped_tiny = 0;
        size_t skipped_invalid = 0;
        
        for(const auto& t: transcripts){
            if(t.exons.size() < 2) continue;
            
            for(size_t i=0; i+1 < t.exons.size(); ++i){
                auto e1 = t.exons[i];
                auto e2 = t.exons[i+1];
                
                total_edges++;
                
                uint32_t intron_size = (e2.start > e1.end) ? (e2.start - e1.end) : 0;
                if(intron_size < cfg.min_junction_size){
                    skipped_tiny++;
                    continue;
                }
                
                if(e1.end >= e2.start){
                    skipped_invalid++;
                    continue;
                }
                
                JunctionKey key{t.chrom, e1.end, e2.start, t.strand};
                junction_map[key].push_back(t.id);
            }
        }
        
        fprintf(stderr,"  Total transcript edges: %zu\n", total_edges);
        fprintf(stderr,"  Skipped (tiny intron < %u bp): %zu\n", cfg.min_junction_size, skipped_tiny);
        fprintf(stderr,"  Skipped (invalid geometry): %zu\n", skipped_invalid);
        fprintf(stderr,"  Unique junction sites: %zu\n", junction_map.size());
        
        junctions.clear();
        junctions.reserve(junction_map.size());
        
        size_t skipped_seq = 0;
        
        for(const auto& [key, tx_list]: junction_map){
            if(!genome.count(key.chrom)) continue;
            
            const string& chr = genome[key.chrom];
            uint32_t donor = key.donor;
            uint32_t acceptor = key.acceptor;
            
            uint32_t left_start = (donor > (uint32_t)cfg.flank) ? (donor - cfg.flank) : 0;
            
            if(left_start >= donor || left_start >= chr.size()){
                skipped_seq++;
                continue;
            }
            
            string left(chr.data() + left_start, donor - left_start);
            
            uint32_t right_end = min<uint32_t>((uint32_t)chr.size(), acceptor + cfg.flank);
            if(acceptor >= chr.size() || acceptor >= right_end){
                skipped_seq++;
                continue;
            }
            
            string right(chr.data() + acceptor, right_end - acceptor);
            string junction_seq = left + right;
            
            if(key.strand == '-'){
                junction_seq = revcomp(junction_seq);
            }
            
            if((int)junction_seq.size() < cfg.k){
                skipped_seq++;
                continue;
            }
            
            JunctionInfo j;
            j.id = key.chrom + ":" + to_string(donor) + "-" + to_string(acceptor) + ":" + key.strand;
            j.gene_id = tx_list[0];
            j.chrom = key.chrom;
            j.strand = key.strand;
            j.donor_pos = donor;
            j.acceptor_pos = acceptor;
            j.seq = move(junction_seq);
            j.source_transcripts = tx_list;
            
            junctions.push_back(move(j));
        }
        
        if(skipped_seq > 0){
            fprintf(stderr,"  Skipped (sequence issues): %zu\n", skipped_seq);
        }
        
        fprintf(stderr,"  Created %zu unique junctions (%.1fx reduction)\n", 
                junctions.size(), 
                total_edges > 0 ? (double)total_edges / junctions.size() : 1.0);
    }

    uint64_t append_seq(const string& seq){
        auto packed = pack2bit(seq);
        uint64_t off=seq_blob.size();
        seq_blob.insert(seq_blob.end(), packed.begin(), packed.end());
        return off;
    }

    void build_targets(){
        fprintf(stderr,"Packing sequences (2-bit) and assigning target IDs...\n");
        tx_meta.clear(); jx_meta.clear(); seq_blob.clear();
        tx_meta.reserve(transcripts.size());
        jx_meta.reserve(junctions.size());

        for(const auto& t: transcripts){
            TargetMetaTX m;
            m.id=t.id;
            m.chrom=t.chrom;
            m.strand=t.strand;
            m.len_bp=(uint32_t)t.seq.size();
            m.seq_off=append_seq(t.seq);

            m.genomic_exons.reserve(t.exons.size());
            for(const auto& e: t.exons){
                m.genomic_exons.emplace_back(e.start, e.end);
            }
            tx_meta.push_back(move(m));
        }

        for(const auto& j: junctions){
            TargetMetaJX m;
            m.id=j.id;
            m.gene_id=j.gene_id;
            m.chrom=j.chrom;
            m.strand=j.strand;
            m.donor_pos=j.donor_pos;
            m.acceptor_pos=j.acceptor_pos;
            m.len_bp=(uint32_t)j.seq.size();
            m.seq_off=append_seq(j.seq);
            jx_meta.push_back(move(m));
        }

        fprintf(stderr,"  Packed sequences: %.1f MB\n", seq_blob.size()/1024.0/1024.0);

        size_t exon_bytes = 0;
        for(const auto& m: tx_meta){
            exon_bytes += m.genomic_exons.size() * sizeof(pair<uint32_t,uint32_t>);
        }
        fprintf(stderr,"  Exon metadata: %.1f MB (for SAM generation)\n", exon_bytes/1024.0/1024.0);
    }

    // NEW: Build transcript k-mer map for junction filtering
    void build_transcript_kmer_map(){
    fprintf(stderr,"Building transcript k-mer map (for junction filtering)...\n");
    transcript_kmer_map.clear();
    transcript_kmer_map.reserve(10000000);
    
    size_t total_kmers = 0;
    
    for(uint32_t tid=0; tid<transcripts.size(); ++tid){
        const auto& t = transcripts[tid];
        int L = (int)t.seq.size();
        
        // FIX: Index k-mers at BOTH ends of the transcript sequence
        // These correspond to exon boundaries in the spliced transcript
        
        uint32_t overlap = min((uint32_t)cfg.junction_flank_overlap, (uint32_t)L);
        
        // Start of transcript (first exon start)
        for(int pos = 0; pos + cfg.k <= (int)overlap && pos + cfg.k <= L; ++pos){
            uint64_t kmer;
            if(encode_kmer64(t.seq.data() + pos, cfg.k, kmer)){
                uint64_t canonical = canonical_kmer(kmer, cfg.k);
                transcript_kmer_map[canonical].push_back(tid);
                total_kmers++;
            }
        }
        
        // End of transcript (last exon end)
        int end_start = max(0, L - (int)overlap);
        for(int pos = end_start; pos + cfg.k <= L; ++pos){
            uint64_t kmer;
            if(encode_kmer64(t.seq.data() + pos, cfg.k, kmer)){
                uint64_t canonical = canonical_kmer(kmer, cfg.k);
                transcript_kmer_map[canonical].push_back(tid);
                total_kmers++;
            }
        }
        
        // ALSO: For multi-exon transcripts, index k-mers at INTERNAL exon boundaries
        // This requires computing cumulative exon positions in the spliced transcript
        if(t.exons.size() > 1){
            int tx_pos = 0;
            for(size_t i=0; i<t.exons.size(); ++i){
                int exon_len = t.exons[i].end - t.exons[i].start;
                
                // End of this exon (before junction)
                int boundary_start = max(tx_pos + exon_len - (int)overlap, tx_pos);
                int boundary_end = tx_pos + exon_len;
                
                for(int pos = boundary_start; pos + cfg.k <= boundary_end && pos + cfg.k <= L; ++pos){
                    uint64_t kmer;
                    if(encode_kmer64(t.seq.data() + pos, cfg.k, kmer)){
                        uint64_t canonical = canonical_kmer(kmer, cfg.k);
                        transcript_kmer_map[canonical].push_back(tid);
                        total_kmers++;
                    }
                }
                
                tx_pos += exon_len;
                
                // Start of next exon (after junction)
                if(i + 1 < t.exons.size()){
                    int next_boundary_end = min(tx_pos + (int)overlap, L);
                    for(int pos = tx_pos; pos + cfg.k <= next_boundary_end && pos + cfg.k <= L; ++pos){
                        uint64_t kmer;
                        if(encode_kmer64(t.seq.data() + pos, cfg.k, kmer)){
                            uint64_t canonical = canonical_kmer(kmer, cfg.k);
                            transcript_kmer_map[canonical].push_back(tid);
                            total_kmers++;
                        }
                    }
                }
            }
        }
    }
    
    fprintf(stderr,"  Indexed %zu boundary k-mers from %zu transcripts\n", 
            total_kmers, transcripts.size());
    fprintf(stderr,"  Unique boundary k-mers: %zu\n", transcript_kmer_map.size());
}
    void emit_postings(){
        ensure_dir(cfg.tmpdir);
        const size_t run_bytes = max<size_t>(1, cfg.chunk_mb) * (1<<20);
        size_t run_id=0;
        auto make_path=[&](size_t rid){ return cfg.tmpdir + "/run_" + to_string(rid) + ".bin"; };
        RunWriter writer(make_path(run_id), run_bytes);

        auto rotate_run = [&](){
            if(!writer.has_data()) return;
            writer.close();
            run_files.push_back(writer.path);
            ++run_id;
            writer = RunWriter(make_path(run_id), run_bytes);
        };

        size_t emitted=0;
        size_t tx_kmers_emitted = 0;
        size_t jx_kmers_emitted = 0;
        size_t jx_kmers_skipped = 0;

        fprintf(stderr,"Streaming postings into runs (chunk=%zu MB, canonical k-mers)...\n", cfg.chunk_mb);

        // Emit transcript k-mers (canonical)
        for(uint32_t tid=0; tid<tx_meta.size(); ++tid){
            const auto& t=transcripts[tid];
            int L=(int)t.seq.size();
            for(int i=0;i+cfg.k<=L;++i){
                uint64_t kmer;
                if(!encode_kmer64(t.seq.data()+i, cfg.k, kmer)) continue;
                
                uint64_t canonical = canonical_kmer(kmer, cfg.k);
                writer.add(PostRec{canonical, tid, (uint32_t)i});
                ++emitted;
                ++tx_kmers_emitted;
                if(writer.buffered_bytes() >= run_bytes) rotate_run();
            }
        }
        
        fprintf(stderr,"  Emitted %zu transcript k-mers\n", tx_kmers_emitted);

        // Emit junction k-mers (canonical, skip exonic k-mers)
        uint32_t jx_base = (uint32_t)tx_meta.size();
        for(uint32_t j=0;j<jx_meta.size();++j){
            const auto& x=junctions[j];
            int L=(int)x.seq.size();
            for(int i=0;i+cfg.k<=L;++i){
                uint64_t kmer;
                if(!encode_kmer64(x.seq.data()+i, cfg.k, kmer)) continue;
                
                uint64_t canonical = canonical_kmer(kmer, cfg.k);
                
                // NEW: Skip if this k-mer is in transcript boundaries
                if(transcript_kmer_map.count(canonical)){
                    jx_kmers_skipped++;
                    continue;
                }
                
                writer.add(PostRec{canonical, jx_base + j, (uint32_t)i});
                ++emitted;
                ++jx_kmers_emitted;
                if(writer.buffered_bytes() >= run_bytes) rotate_run();
            }
        }
        
        if(writer.has_data()){
            writer.close();
            run_files.push_back(writer.path);
        }
        
        fprintf(stderr,"  Emitted %zu junction k-mers (skipped %zu exonic)\n", 
                jx_kmers_emitted, jx_kmers_skipped);
        fprintf(stderr,"  Total: ~%zu postings into %zu run(s)\n", emitted, run_files.size());
        fprintf(stderr,"  Junction filtering: %.1f%% of junction k-mers were exonic (skipped)\n",
                (jx_kmers_skipped + jx_kmers_emitted > 0) ? 
                100.0 * jx_kmers_skipped / (jx_kmers_skipped + jx_kmers_emitted) : 0.0);
    }

    void merge_and_write(){
        string tmp_post_path = cfg.tmpdir + "/postings.tmp";

        vector<unique_ptr<RunReader>> readers;
        readers.reserve(run_files.size());
        for(const auto& p: run_files){
            auto rr = make_unique<RunReader>(p);
            rr->open();
            readers.push_back(move(rr));
        }

        struct Item{
            PostRec rec;
            size_t rid;
        };
        struct Cmp{
            bool operator()(const Item& a, const Item& b) const{
                if(a.rec.kmer!=b.rec.kmer) return a.rec.kmer > b.rec.kmer;
                if(a.rec.tid != b.rec.tid) return a.rec.tid > b.rec.tid;
                return a.rec.pos > b.rec.pos;
            }
        };
        priority_queue<Item, vector<Item>, Cmp> pq;
        for(size_t i=0;i<readers.size();++i){
            PostRec r;
            if(readers[i]->next(r)) pq.push(Item{r,i});
        }

        FILE* postf=fopen(tmp_post_path.c_str(),"wb");
        if(!postf){ fprintf(stderr,"ERROR: cannot create %s\n", tmp_post_path.c_str()); exit(1); }
        setvbuf(postf, nullptr, _IOFBF, 1<<20);
        uint32_t cur_post_off=0;

        const uint32_t num_buckets = 1u << cfg.prefix_bits;
        vector<uint32_t> bucket_counts(num_buckets, 0);
        vector<uint8_t> suffix6;
        vector<uint32_t> post_offs;
        vector<uint8_t> common_bits;

        suffix6.reserve(6u * 1000000);
        post_offs.reserve(1000000);

        auto set_common=[&](uint32_t idx){
            size_t byte=idx>>3, bit=idx&7;
            if(byte >= common_bits.size()) common_bits.resize(byte+1,0);
            common_bits[byte] |= (1u<<bit);
        };

        auto zigzag_encode=[](int64_t n)->uint64_t{
            return (uint64_t)((n << 1) ^ (n >> 63));
        };

        vector<pair<uint32_t, vector<uint32_t>>> per_tid;
        per_tid.reserve(64);

        const int suffix_bits = 62 - cfg.prefix_bits;
        const uint64_t suffix_mask = (suffix_bits==64)?~0ULL : ((suffix_bits==0)?0ULL : ((1ULL<<suffix_bits)-1ULL));

        uint64_t prev_kmer=(uint64_t)(-1);
        uint32_t num_keys=0;
        uint64_t skipped_high_df = 0;
        uint64_t singleton_count = 0;  // NEW: Count DF=1 k-mers

        // NEW: Compact encoding for DF=1 k-mers
        auto flush_group = [&](uint64_t kmer){
            if(per_tid.empty()) return;
            
            uint32_t total_count = 0;
            for(const auto& g: per_tid) total_count += (uint32_t)g.second.size();
            
            if(cfg.df_cap > 0 && total_count > cfg.df_cap){
                skipped_high_df++;
                per_tid.clear();
                return;
            }
            
            uint32_t bucket = (uint32_t)(kmer >> (62 - cfg.prefix_bits));
            bucket_counts[bucket]++;

            uint64_t suf = kmer & suffix_mask;
            suffix6.push_back((uint8_t)((suf>>0)&0xFF));
            suffix6.push_back((uint8_t)((suf>>8)&0xFF));
            suffix6.push_back((uint8_t)((suf>>16)&0xFF));
            suffix6.push_back((uint8_t)((suf>>24)&0xFF));
            suffix6.push_back((uint8_t)((suf>>32)&0xFF));
            suffix6.push_back((uint8_t)((suf>>40)&0xFF));

            post_offs.push_back(cur_post_off);

            sort(per_tid.begin(), per_tid.end(),
                 [](const auto& a, const auto& b){ return a.first<b.first; });

            if(total_count > 100){
                set_common(num_keys);
            }

            vector<uint8_t> enc;
            
            // NEW: COMPACT ENCODING FOR DF=1 (single target, single position)
            if(per_tid.size() == 1 && per_tid[0].second.size() == 1){
                singleton_count++;
                // Format: 0xFF (marker) + tid(24 bits) + pos(24 bits) = 7 bytes
                enc.push_back(0xFF);  // Singleton marker
                uint32_t tid = per_tid[0].first;
                uint32_t pos = per_tid[0].second[0];
                enc.push_back((uint8_t)(tid & 0xFF));
                enc.push_back((uint8_t)((tid >> 8) & 0xFF));
                enc.push_back((uint8_t)((tid >> 16) & 0xFF));
                enc.push_back((uint8_t)(pos & 0xFF));
                enc.push_back((uint8_t)((pos >> 8) & 0xFF));
                enc.push_back((uint8_t)((pos >> 16) & 0xFF));
            }
            else {
                // STANDARD ENCODING (varint)
                enc.reserve(per_tid.size()*4);
                auto put_var=[&](uint64_t v){
                    while(v>=0x80){ enc.push_back((uint8_t)((v&0x7F)|0x80)); v>>=7; }
                    enc.push_back((uint8_t)v);
                };

                put_var(per_tid.size());
                uint32_t prev_tid=0;
                for(const auto& g: per_tid){
                    put_var((uint64_t)(g.first - prev_tid));
                    prev_tid = g.first;
                    put_var((uint64_t)g.second.size());

                    if(g.second.size() == 1){
                        put_var((uint64_t)g.second[0]);
                    }else{
                        put_var((uint64_t)g.second[0]);
                        int64_t prev_delta = 0;
                        for(size_t i=1; i<g.second.size(); ++i){
                            int64_t delta = (int64_t)g.second[i] - (int64_t)g.second[i-1];
                            int64_t delta_of_delta = delta - prev_delta;
                            put_var(zigzag_encode(delta_of_delta));
                            prev_delta = delta;
                        }
                    }
                }
            }
            
            if(!enc.empty()){
                fwrite(enc.data(), 1, enc.size(), postf);
                cur_post_off += (uint32_t)enc.size();
            }

            per_tid.clear();
            ++num_keys;
        };

        size_t merged=0;
        while(!pq.empty()){
            Item it=pq.top(); pq.pop();
            PostRec r=it.rec;
            PostRec nxt;
            if(readers[it.rid]->next(nxt)) pq.push(Item{nxt, it.rid});

            if(prev_kmer==(uint64_t)(-1)) prev_kmer=r.kmer;
            if(r.kmer != prev_kmer){
                flush_group(prev_kmer);
                prev_kmer=r.kmer;
            }
            if(per_tid.empty() || per_tid.back().first != r.tid){
                per_tid.emplace_back(r.tid, vector<uint32_t>{r.pos});
            }else{
                auto& v = per_tid.back().second;
                if(v.empty() || v.back()!=r.pos) v.push_back(r.pos);
            }

            if(++merged % 5000000 == 0){
                fprintf(stderr,"  merged %zu postings...\r", merged);
            }
        }
        if(prev_kmer!=(uint64_t)(-1)) flush_group(prev_kmer);

        if(skipped_high_df > 0){
            fprintf(stderr,"\n  Skipped %llu high-DF k-mers (DF > %u)\n", 
                    (unsigned long long)skipped_high_df, cfg.df_cap);
        }
        fprintf(stderr,"  Compact encoding: %llu singletons (%.1f%% of k-mers)\n",
                (unsigned long long)singleton_count,
                num_keys > 0 ? 100.0 * singleton_count / num_keys : 0.0);

        vector<uint32_t> bucket_cum(num_buckets+1, 0);
        uint64_t acc=0;
        for(uint32_t b=0;b<num_buckets;b++){
            bucket_cum[b]=(uint32_t)acc;
            acc += bucket_counts[b];
        }
        bucket_cum[num_buckets]=(uint32_t)acc;

        post_offs.push_back(cur_post_off);
        fclose(postf);

        size_t need_bytes = (num_keys+7)/8;
        if(common_bits.size()<need_bytes) common_bits.resize(need_bytes,0);

        // ---------- WRITE INDEX (v15) ----------
        const string outpath = cfg.outprefix + ".index.v15.bin";
        ofstream out(outpath, ios::binary);
        if(!out){ fprintf(stderr,"ERROR: cannot open %s for write\n", outpath.c_str()); exit(1); }

        const string magic = "RMIDX15\n";
        out.write(magic.data(), magic.size());

        write_pod(out, (uint32_t)cfg.k);
        write_pod(out, (uint32_t)cfg.flank);
        write_pod(out, (uint32_t)cfg.prefix_bits);

        // Contig table
        {
            const char sqtb[4] = {'S','Q','T','B'};
            out.write(sqtb, 4);

            vector<pair<string,uint32_t>> contigs;
            contigs.reserve(genome.size());
            for (const auto& kv : genome) {
                contigs.emplace_back(kv.first, (uint32_t)kv.second.size());
            }
            sort(contigs.begin(), contigs.end(),
                 [](const auto& a, const auto& b){ return a.first < b.first; });

            uint32_t n_contigs = (uint32_t)contigs.size();
            write_pod(out, n_contigs);
            for (const auto& pr : contigs) {
                write_str(out, pr.first);
                write_pod(out, pr.second);
            }
        }

        uint32_t n_tx = (uint32_t)tx_meta.size();
        uint32_t n_jx = (uint32_t)jx_meta.size();
        uint32_t n_targets = n_tx + n_jx;
        write_pod(out, n_tx);
        write_pod(out, n_jx);
        write_pod(out, n_targets);

        for(const auto& m: tx_meta){
            write_str(out, m.id);
            write_str(out, m.chrom);
            out.write((const char*)&m.strand, 1);
            write_pod(out, m.len_bp);
            write_pod(out, m.seq_off);

            uint32_t n_exons = (uint32_t)m.genomic_exons.size();
            write_pod(out, n_exons);
            for(const auto& [start, end]: m.genomic_exons){
                write_pod(out, start);
                write_pod(out, end);
            }
        }

        for(const auto& m: jx_meta){
            write_str(out, m.id);
            write_str(out, m.gene_id);
            write_str(out, m.chrom);
            out.write((const char*)&m.strand, 1);
            write_pod(out, m.donor_pos);
            write_pod(out, m.acceptor_pos);
            write_pod(out, m.len_bp);
            write_pod(out, m.seq_off);
        }

        uint64_t seq_size = (uint64_t)seq_blob.size();
        write_pod(out, seq_size);
        if(seq_size) out.write((const char*)seq_blob.data(), (streamsize)seq_blob.size());

        write_pod(out, (uint32_t)num_buckets);
        write_pod(out, (uint32_t)num_keys);

        for(uint32_t v: bucket_cum) write_pod(out, v);
        if(!suffix6.empty()) out.write((const char*)suffix6.data(), (streamsize)suffix6.size());
        for(uint32_t v: post_offs) write_pod(out, v);

        uint32_t cb_size=(uint32_t)common_bits.size();
        write_pod(out, cb_size);
        if(cb_size) out.write((const char*)common_bits.data(), cb_size);

        {
            ifstream in(tmp_post_path, ios::binary);
            in.seekg(0, ios::end);
            uint64_t ps=(uint64_t)in.tellg();
            in.seekg(0, ios::beg);
            write_pod(out, ps);
            const size_t BUFSZ=1<<20;
            vector<char> buf(BUFSZ);
            while(in){
                in.read(buf.data(), BUFSZ);
                streamsize got=in.gcount();
                if(got>0) out.write(buf.data(), got);
            }
        }
        out.close();

        remove(tmp_post_path.c_str());
        for(const auto& p: run_files) remove(p.c_str());

        struct stat st{};
        stat(outpath.c_str(), &st);
        fprintf(stderr,"\n=== Index Statistics (v15 - Canonical + Filtered + Compact) ===\n");
        fprintf(stderr,"K = %d, flank = %d, prefix bits = %d (buckets=%u)\n",
                cfg.k, cfg.flank, cfg.prefix_bits, (1u<<cfg.prefix_bits));
        fprintf(stderr,"Transcripts: %u, Junctions: %u (targets=%u)\n", n_tx, n_jx, n_targets);
        fprintf(stderr,"Unique k-mers: %u (canonical)\n", num_keys);
        fprintf(stderr,"Sequences blob: %.1f MB\n", seq_blob.size()/1024.0/1024.0);
        fprintf(stderr,"Index file size: %.1f MB\n", st.st_size/1024.0/1024.0);
        fprintf(stderr,"Features: Canonical k-mers, junction filtering, compact DF=1 encoding\n");
        if(cfg.df_cap > 0){
            fprintf(stderr,"Physical DF cap: %u (skipped %llu k-mers)\n", 
                    cfg.df_cap, (unsigned long long)skipped_high_df);
        }
        fprintf(stderr,"Compression: %llu singletons use compact encoding (7 bytes)\n",
                (unsigned long long)singleton_count);
        fprintf(stderr,"=====================================================\n");
    }

    void run(){
        parse_fasta(cfg.fasta);
        parse_gtf(cfg.gtf);
        build_junctions();
        build_targets();
        build_transcript_kmer_map();  // NEW
        emit_postings();
        merge_and_write();
        fprintf(stderr,"Index build complete!\n");
    }
};

// ---------- CLI ----------
static Config parse_cli(int argc, char** argv){
    if(argc<4){
        fprintf(stderr,
            "Usage: %s <genome.fa[.gz]> <annotations.gtf[.gz]> <output_prefix>\n"
            "       [--k 31] [--flank 25] [--prefix-bits 20]\n"
            "       [--chunk-mb 256] [--tmpdir ./tmp]\n"
            "       [--df-cap 0] [--min-junction-size 20]\n"
            "       [--junction-flank-overlap 15]\n"
            "\n"
            "v15 features:\n"
            "  - Canonical k-mer indexing (50%% fewer k-mers)\n"
            "  - Skip exonic k-mers in junctions (cleaner evidence)\n"
            "  - Compact encoding for DF=1 k-mers (better compression)\n"
            "  - Physical DF capping\n"
            "\n"
            "Options:\n"
            "  --df-cap N               Skip k-mers with DF > N (0 = no cap, 100-200 recommended)\n"
            "  --min-junction-size      Minimum intron size (default: 20bp)\n"
            "  --junction-flank-overlap Check N bp at exon boundaries for filtering (default: 15bp)\n"
            , argv[0]);
        exit(1);
    }
    Config c;
    c.fasta=argv[1];
    c.gtf=argv[2];
    c.outprefix=argv[3];
    for(int i=4;i<argc;i++){
        string a=argv[i];
        auto need=[&](int more){
            if(i+more>=argc){
                fprintf(stderr,"Missing argument after %s\n", a.c_str());
                exit(1);
            }
        };
        if(a=="--k"){ need(1); c.k=stoi(argv[++i]); }
        else if(a=="--flank"){ need(1); c.flank=stoi(argv[++i]); }
        else if(a=="--prefix-bits"){ need(1); c.prefix_bits=stoi(argv[++i]);
            if(c.prefix_bits<10 || c.prefix_bits>26){
                fprintf(stderr,"prefix-bits out of range [10..26]\n"); exit(1);
            }
        }
        else if(a=="--chunk-mb"){ need(1); c.chunk_mb=(size_t)stoull(argv[++i]); }
        else if(a=="--tmpdir"){ need(1); c.tmpdir=argv[++i]; }
        else if(a=="--df-cap"){ need(1); c.df_cap=(uint32_t)stoul(argv[++i]); }
        else if(a=="--min-junction-size"){ need(1); c.min_junction_size=(uint32_t)stoul(argv[++i]); }
        else if(a=="--junction-flank-overlap"){ need(1); c.junction_flank_overlap=(uint32_t)stoul(argv[++i]); }
        else{
            fprintf(stderr,"Unknown argument: %s\n", a.c_str());
            exit(1);
        }
    }
    return c;
}

int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    IndexBuilderV15 builder;
    builder.cfg = parse_cli(argc, argv);

    auto t0 = chrono::steady_clock::now();
    ensure_dir(builder.cfg.tmpdir);
    builder.run();
    auto t1 = chrono::steady_clock::now();
    double sec = chrono::duration<double>(t1-t0).count();
    fprintf(stderr,"Total build time: %.2f s\n", sec);
    return 0;
}
