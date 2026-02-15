#include <bits/stdc++.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#endif

#ifdef __SSE4_2__
#include <nmmintrin.h>
#endif

using namespace std;

// --------------------- Configuration ---------------------
static constexpr uint64_t HASH_EMPTY = 0;

// --------------------- Common structs ---------------------
struct Record { 
    string header; 
    string seq; 
    
    Record() = default;
    Record(string h, string s) : header(std::move(h)), seq(std::move(s)) {}
};

// --------------------- Fast string utilities---------------------
static inline string trim(const string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n");
    if (a == string::npos) return "";
    const size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// Fast character conversion using lookup table
static const array<char, 256> TOUPPER_TABLE = []() {
    array<char, 256> table{};
    for (int i = 0; i < 256; ++i) {
        table[i] = static_cast<char>(toupper(i));
    }
    return table;
}();

static inline char fast_toupper(char c) {
    return TOUPPER_TABLE[static_cast<unsigned char>(c)];
}

// --------------------- Fast hash functions ---------------------
// xxHash-inspired fast hash
static inline uint64_t hash_bytes_fast(const char* data, size_t len) {
    constexpr uint64_t PRIME1 = 11400714785074694791ULL;
    constexpr uint64_t PRIME2 = 14029467366897019727ULL;
    constexpr uint64_t PRIME3 = 1609587929392839161ULL;
    constexpr uint64_t PRIME5 = 2870177450012600261ULL;
    
    uint64_t hash = PRIME5;
    // Process 8 bytes at a time
    while (len >= 8) {
        constexpr uint64_t PRIME4 = 9650029242287828579ULL;
        uint64_t val;
        memcpy(&val, data, 8);
        val *= PRIME2;
        val = ((val << 31) | (val >> 33)) * PRIME1;
        hash ^= val;
        hash = ((hash << 27) | (hash >> 37)) * PRIME1 + PRIME4;
        data += 8;
        len -= 8;
    }
    
    // Process remaining bytes
    while (len > 0) {
        hash ^= (*data++) * PRIME5;
        hash = ((hash << 11) | (hash >> 53)) * PRIME1;
        len--;
    }
    
    // Final avalanche
    hash ^= hash >> 33;
    hash *= PRIME2;
    hash ^= hash >> 29;
    hash *= PRIME3;
    hash ^= hash >> 32;
    
    return hash == HASH_EMPTY ? 1 : hash;
}

static inline uint64_t hash_sequence_fast(const string& s) {
    return hash_bytes_fast(s.data(), s.size());
}

// --------------------- Memory-mapped FASTA reader ---------------------
class MemoryMappedFASTA {
private:
#ifdef _WIN32
    HANDLE hFile;
    HANDLE hMapping;
#else
    int fd;
#endif
    char* data;
    size_t file_size;

public:
    explicit MemoryMappedFASTA(const string& filename) : data(nullptr), file_size(0) {
#ifdef _WIN32
        hFile = CreateFileA(filename.c_str(), GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
        if (hFile == INVALID_HANDLE_VALUE) {
            throw runtime_error("Cannot open file: " + filename);
        }

        LARGE_INTEGER size_li;
        if (!GetFileSizeEx(hFile, &size_li)) {
            CloseHandle(hFile);
            throw runtime_error("Cannot stat file: " + filename);
        }
        file_size = size_li.QuadPart;

        if (file_size == 0) {
            CloseHandle(hFile);
            throw runtime_error("Empty file: " + filename);
        }

        hMapping = CreateFileMappingA(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
        if (hMapping == NULL) {
            CloseHandle(hFile);
            throw runtime_error("Cannot create file mapping: " + filename);
        }

        data = static_cast<char*>(MapViewOfFile(hMapping, FILE_MAP_READ, 0, 0, 0));
        if (data == NULL) {
            CloseHandle(hMapping);
            CloseHandle(hFile);
            throw runtime_error("Cannot map file view: " + filename);
        }
#else
        fd = open(filename.c_str(), O_RDONLY);
        if (fd == -1) throw runtime_error("Cannot open file: " + filename);
        
        struct stat sb{};
        if (fstat(fd, &sb) == -1) {
            close(fd);
            throw runtime_error("Cannot stat file: " + filename);
        }
        
        file_size = sb.st_size;
        if (file_size == 0) {
            close(fd);
            throw runtime_error("Empty file: " + filename);
        }
        
        data = static_cast<char*>(mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0));
        if (data == MAP_FAILED) {
            close(fd);
            throw runtime_error("Cannot mmap file: " + filename);
        }
        
        // Advise OS about access pattern
        madvise(data, file_size, MADV_SEQUENTIAL);
#endif
    }

    ~MemoryMappedFASTA() {
#ifdef _WIN32
        if (data) UnmapViewOfFile(data);
        if (hMapping) CloseHandle(hMapping);
        if (hFile != INVALID_HANDLE_VALUE) CloseHandle(hFile);
#else
        if (data && data != MAP_FAILED) {
            munmap(data, file_size);
        }
        if (fd != -1) close(fd);
#endif
    }

    [[nodiscard]] vector<Record> parse() const {
        vector<Record> records;
        
        // Quick scan to estimate record count
        size_t header_count = 0;
        for (size_t i = 0; i < file_size; ++i) {
            if (data[i] == '>') header_count++;
        }
        records.reserve(header_count);
        
        const char* ptr = data;
        const char* end = data + file_size;
        
        string header, sequence;
        header.reserve(256);
        sequence.reserve(10000);
        
        while (ptr < end) {
            if (*ptr == '>') {
                // Save previous record
                if (!header.empty()) {
                    records.emplace_back(std::move(header), std::move(sequence));
                    header.clear();
                    sequence.clear();
                }
                
                // Read header line
                const char* line_start = ptr;
                while (ptr < end && *ptr != '\n' && *ptr != '\r') ptr++;
                header.assign(line_start, ptr);
                
                // Skip newline characters
                while (ptr < end && (*ptr == '\n' || *ptr == '\r')) ptr++;
            } else {
                // Read sequence data, converting to uppercase
                while (ptr < end && *ptr != '>') {
                    char c = *ptr;
                    if (c != '\n' && c != '\r' && c != ' ' && c != '\t') {
                        sequence.push_back(fast_toupper(c));
                    }
                    ptr++;
                }
            }
        }
        
        // Don't forget last record
        if (!header.empty()) {
            records.emplace_back(std::move(header), std::move(sequence));
        }
        
        return records;
    }
};

// Fallback stream-based reader for when mmap fails
vector<Record> read_fasta_stream(const string& filename) {
    ifstream fin(filename);
    if (!fin) throw runtime_error("Cannot open file: " + filename);
    
    vector<Record> recs;
    string line, header, seq;
    
    while (getline(fin, line)) {
        if (!line.empty() && line[0] == '>') {
            if (!header.empty()) { 
                recs.emplace_back(std::move(header), std::move(seq)); 
                seq.clear(); 
            }
            header = trim(line);
        } else {
            for (char c : line) {
                if (!isspace(static_cast<unsigned char>(c))) {
                    seq.push_back(fast_toupper(c));
                }
            }
        }
    }
    if (!header.empty()) recs.emplace_back(std::move(header), std::move(seq));
    return recs;
}

// --------------------- Optimized edit distance ---------------------
class OptimizedEditDistance {
public:
    // Quick rejection based on length difference and character frequency.
    static bool quick_reject_by_edits(const string& a, const string& b, const int max_edits) {
        if (abs((int)a.size() - (int)b.size()) > max_edits) {
            return true;
        }

        const int max_len = max(a.size(), b.size());
        if (max_len > 50) {
            array<int, 256> freq_diff = {0};

            for (char c : a) freq_diff[static_cast<unsigned char>(c)]++;
            for (char c : b) freq_diff[static_cast<unsigned char>(c)]--;

            int total_diff = 0;
            for (int f : freq_diff) {
                total_diff += abs(f);
            }

            if (total_diff / 2 > max_edits) {
                return true;
            }
        }

        return false;
    }
    
#ifdef __SSE4_2__
    // SIMD-optimized exact match check
    static bool sequences_equal_simd(const string& a, const string& b) {
        if (a.size() != b.size()) return false;
        
        const char* p1 = a.data();
        const char* p2 = b.data();
        size_t len = a.size();
        
        // Process 16 bytes at a time
        while (len >= 16) {
            const __m128i v1 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(p1));
            const __m128i v2 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(p2));

            if (__m128i cmp = _mm_cmpeq_epi8(v1, v2); _mm_movemask_epi8(cmp) != 0xFFFF) {
                return false;
            }
            
            p1 += 16;
            p2 += 16;
            len -= 16;
        }
        
        // Handle remainder
        return memcmp(p1, p2, len) == 0;
    }
#endif
    
    // Optimized banded edit distance with early termination
    static int compute_banded(const string& a, const string& b, const int max_edits) {
        if (quick_reject_by_edits(a, b, max_edits)) {
            return max_edits + 1;
        }

        const int n = (int)a.size();
        const int m = (int)b.size();
        if (abs(n - m) > max_edits) return max_edits + 1;
        
        const int INF = max_edits + 1;
        const int band = max_edits;
        
        vector<int> prev, cur;
        int jmin_prev = max(0, 0 - band);
        int jmax_prev = min(m, 0 + band);
        prev.assign(jmax_prev - jmin_prev + 1, INF);
        
        for (int j = jmin_prev; j <= jmax_prev; ++j) {
            prev[j - jmin_prev] = j;
            if (prev[j - jmin_prev] > max_edits) prev[j - jmin_prev] = INF;
        }
        
        for (int i = 1; i <= n; ++i) {
           const int jmin_cur = max(0, i - band);
            const int jmax_cur = min(m, i + band);
            cur.assign(jmax_cur - jmin_cur + 1, INF);
            int best = INF;
            
            for (int j = jmin_cur; j <= jmax_cur; ++j) {
                int cost_sub = INF, cost_del = INF, cost_ins = INF;
                
                // Substitution
                if (j - 1 >= jmin_prev && j - 1 <= jmax_prev) {
                    if (int pv = prev[(j - 1) - jmin_prev]; pv != INF) cost_sub = pv + (a[i-1] == b[j-1] ? 0 : 1);
                }
                
                // Deletion
                if (j >= jmin_prev && j <= jmax_prev) {
                    if (const int pv = prev[j - jmin_prev]; pv != INF) cost_del = pv + 1;
                }
                
                // Insertion
                if (j - 1 >= jmin_cur && j - 1 <= jmax_cur) {
                    if (const int cv = cur[(j - 1) - jmin_cur]; cv != INF) cost_ins = cv + 1;
                }
                
                int v = min({cost_sub, cost_del, cost_ins});
                if (v > max_edits) v = INF;
                cur[j - jmin_cur] = v;
                if (v < best) best = v;
            }
            
            if (best == INF) return max_edits + 1;
            prev.swap(cur);
            jmin_prev = jmin_cur; 
            jmax_prev = jmax_cur;
        }
        
        int res = prev[m - jmin_prev];
        return (res == INF) ? (max_edits + 1) : res;
    }
};

static inline bool identity_at_least(const string& a, const string& b, const int pct_threshold) {
    const int max_len = max(a.size(), b.size());
    if (max_len == 0) return true;
    const int allowed = (int)ceil((100.0 - pct_threshold) * max_len / 100.0);
    
    if (abs((int)a.size() - (int)b.size()) > allowed) return false;
    
#ifdef __SSE4_2__
    if (pct_threshold == 100) {
        return OptimizedEditDistance::sequences_equal_simd(a, b);
    }
#endif

    const int ed = OptimizedEditDistance::compute_banded(a, b, allowed);
    if (ed > allowed) return false;

    const double pid = 100.0 * (1.0 - (double)ed / (double)max_len);
    return pid + 1e-9 >= pct_threshold;
}

static inline bool sequences_equal_fast(const string& a, const string& b) {
#ifdef __SSE4_2__
    return OptimizedEditDistance::sequences_equal_simd(a, b);
#else
    return a == b;
#endif
}

static vector<uint64_t> compute_hashes_parallel(const vector<Record>& recs, int threads) {
    vector<uint64_t> hashes(recs.size(), 1);
    if (recs.empty()) return hashes;

    const int worker_count = min<int>(max(1, threads), recs.size());
    if (worker_count == 1) {
        for (size_t i = 0; i < recs.size(); ++i) {
            hashes[i] = hash_sequence_fast(recs[i].seq);
        }
        return hashes;
    }

    vector<thread> workers;
    workers.reserve(worker_count);
    const size_t block = (recs.size() + worker_count - 1) / worker_count;
    for (int w = 0; w < worker_count; ++w) {
        const size_t start = w * block;
        const size_t stop = min(recs.size(), start + block);
        if (start >= stop) break;
        workers.emplace_back([&, start, stop] {
            for (size_t i = start; i < stop; ++i) {
                hashes[i] = hash_sequence_fast(recs[i].seq);
            }
        });
    }
    for (auto& worker : workers) worker.join();
    return hashes;
}

static pair<int, int> feasible_length_bounds(const int len, const int pct_threshold) {
    if (len == 0) return {0, 0};
    const double ratio = pct_threshold / 100.0;
    const int min_len = max(0, (int)ceil(len * ratio - 1e-12));
    const int max_len = (int)floor(len / ratio + 1e-12);
    return {min_len, max_len};
}

static vector<uint64_t> sequence_seed_hashes(const string& seq, const int k) {
    vector<uint64_t> seeds;
    if (seq.empty()) {
        seeds.push_back(1);
        return seeds;
    }

    const int n = (int)seq.size();
    if (n <= k) {
        seeds.push_back(hash_sequence_fast(seq));
        return seeds;
    }

    const int stride = max(1, k / 2);
    for (int pos = 0; pos + k <= n; pos += stride) {
        seeds.push_back(hash_bytes_fast(seq.data() + pos, k));
    }
    const int tail = n - k;
    if (tail % stride != 0) {
        seeds.push_back(hash_bytes_fast(seq.data() + tail, k));
    }

    sort(seeds.begin(), seeds.end());
    seeds.erase(unique(seeds.begin(), seeds.end()), seeds.end());
    return seeds;
}

// --------------------- Arguments ---------------------
struct Args {
    string in_path, out_path;
    int pct = 100;
    int kmer = 12;
    int threads = thread::hardware_concurrency();
    bool use_mmap = true;
};

void usage(const char* prog) {
    cerr << "Usage: " << prog << " -i IN.fasta -o OUT.fasta [-p 80|85|90|95|100] [-k KMER] [-t THREADS]\n";
    cerr << "Options:\n";
    cerr << "  -i FILE     Input FASTA file\n";
    cerr << "  -o FILE     Output FASTA file\n";
    cerr << "  -p INT      Identity percentage threshold (80,85,90,95,100) [100]\n";
    cerr << "  -k INT      K-mer size for indexing [12]\n";
    cerr << "  -t INT      Number of threads [" << thread::hardware_concurrency() << "]\n";
    cerr << "  -h          Show this help\n";
}

Args parse_args(const int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; i++) {
        string s = argv[i];
        if (s == "-i" && i + 1 < argc) {
            a.in_path = argv[++i];
        } else if (s == "-o" && i + 1 < argc) {
            a.out_path = argv[++i];
        } else if ((s == "-p" || s == "-k" || s == "-t") && i + 1 < argc) {
            const char* value = argv[++i];
            try {
                const int parsed = stoi(value);
                if (s == "-p") a.pct = parsed;
                else if (s == "-k") a.kmer = parsed;
                else a.threads = parsed;
            } catch (const exception&) {
                cerr << "Invalid numeric value for " << s << ": " << value << "\n";
                exit(1);
            }
        }
        else if (s == "-h" || s == "--help") { usage(argv[0]); exit(0); }
        else { 
            cerr << "Unknown or incomplete option: " << s << "\n"; 
            usage(argv[0]); 
            exit(1); 
        }
    }
    
    if (a.in_path.empty() || a.out_path.empty()) { 
        usage(argv[0]); 
        exit(1); 
    }
    
    if (!(a.pct == 80 || a.pct == 85 || a.pct == 90 || a.pct == 95 || a.pct == 100)) {
        cerr << "Error: -p must be 80, 85, 90, 95, or 100\n"; 
        exit(1);
    }
    
    if (a.kmer < 4) a.kmer = 4;
    if (a.threads < 1) a.threads = 1;
    return a;
}

// --------------------- Output writing ---------------------
void write_fasta(const string& filename, const vector<Record>& recs, size_t width = 80) {
    ofstream out(filename);
    if (!out) throw runtime_error("Cannot open output file: " + filename);
    
    for (const auto& r : recs) {
        out << r.header << '\n';
        const string& s = r.seq;
        for (size_t i = 0; i < s.size(); i += width) {
            out << s.substr(i, min(width, s.size() - i)) << '\n';
        }
    }
}

// --------------------- Main logic ---------------------
int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    Args args = parse_args(argc, argv);
    
    // Read input file
    vector<Record> recs;
    try {
        if (args.use_mmap) {
            MemoryMappedFASTA reader(args.in_path);
            recs = reader.parse();
        } else {
            recs = read_fasta_stream(args.in_path);
        }
    } catch (const exception& e) {
        cerr << "Failed to read input file: " << e.what() << "\n";
        // Fallback to stream reading
        try {
            recs = read_fasta_stream(args.in_path);
        } catch (const exception& e2) {
            cerr << "Failed to read input file with fallback: " << e2.what() << "\n";
            return 1;
        }
    }
    
    if (recs.empty()) {
        cerr << "No sequences found in input file\n";
        return 1;
    }
    
    cerr << "Loaded " << recs.size() << " sequences\n";
    
    vector<int> kept_indices;
    kept_indices.reserve(recs.size());
    size_t removed_count = 0;

    if (args.pct == 100) {
        const vector<uint64_t> seq_hashes = compute_hashes_parallel(recs, args.threads);
        unordered_map<uint64_t, vector<int>> hash_to_kept;
        hash_to_kept.reserve(recs.size() * 2);

        for (int i = 0; i < (int)recs.size(); ++i) {
            bool is_dup = false;
            auto& bucket = hash_to_kept[seq_hashes[i]];
            for (const int candidate_idx : bucket) {
                if (sequences_equal_fast(recs[i].seq, recs[candidate_idx].seq)) {
                    is_dup = true;
                    break;
                }
            }

            if (is_dup) {
                ++removed_count;
            } else {
                bucket.push_back(i);
                kept_indices.push_back(i);
            }
        }
    } else {
        unordered_map<uint64_t, vector<int>> seed_index;
        seed_index.reserve(recs.size() * 2);
        map<int, vector<int>> kept_by_length;
        vector<int> seen_stamp(recs.size(), 0);
        int stamp = 1;
        vector<int> candidates;
        candidates.reserve(256);

        for (int i = 0; i < (int)recs.size(); ++i) {
            const string& seq = recs[i].seq;
            const int seq_len = (int)seq.size();
            bool is_dup = false;
            ++stamp;
            if (stamp == INT_MAX) {
                fill(seen_stamp.begin(), seen_stamp.end(), 0);
                stamp = 1;
            }
            candidates.clear();

            const vector<uint64_t> seeds = sequence_seed_hashes(seq, args.kmer);
            for (const uint64_t seed : seeds) {
                auto it = seed_index.find(seed);
                if (it == seed_index.end()) continue;
                for (const int candidate_idx : it->second) {
                    if (seen_stamp[candidate_idx] == stamp) continue;
                    seen_stamp[candidate_idx] = stamp;
                    candidates.push_back(candidate_idx);
                }
            }

            for (const int candidate_idx : candidates) {
                if (identity_at_least(seq, recs[candidate_idx].seq, args.pct)) {
                    is_dup = true;
                    break;
                }
            }

            if (!is_dup) {
                const auto [min_len, max_len] = feasible_length_bounds(seq_len, args.pct);
                for (auto it = kept_by_length.lower_bound(min_len); it != kept_by_length.end() && it->first <= max_len; ++it) {
                    for (const int candidate_idx : it->second) {
                        if (seen_stamp[candidate_idx] == stamp) continue;
                        if (identity_at_least(seq, recs[candidate_idx].seq, args.pct)) {
                            is_dup = true;
                            break;
                        }
                    }
                    if (is_dup) break;
                }
            }

            if (is_dup) {
                ++removed_count;
            } else {
                kept_indices.push_back(i);
                kept_by_length[seq_len].push_back(i);
                for (const uint64_t seed : seeds) {
                    seed_index[seed].push_back(i);
                }
            }
        }
    }

    // Build final result in original input order.
    vector<Record> kept_records;
    kept_records.reserve(kept_indices.size());
    for (const int idx : kept_indices) {
        kept_records.push_back(std::move(recs[idx]));
    }
    
    // Write output
    try {
        write_fasta(args.out_path, kept_records);
    } catch (const exception& e) {
        cerr << "Failed to write output: " << e.what() << "\n";
        return 1;
    }
    
    // Report statistics
    size_t kept_count = kept_records.size();
    
    cerr << "Input sequences:  " << recs.size() << "\n";
    cerr << "Kept sequences:   " << kept_count << "\n";
    cerr << "Removed sequences:" << removed_count << "\n";
    cerr << "Threshold:        " << args.pct << "% identity\n";
    if (args.pct < 100) cerr << "Seed k-mer:       " << args.kmer << "\n";
    cerr << "Threads:          " << args.threads << "\n";
    
    return 0;
}
