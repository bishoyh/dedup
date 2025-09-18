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
static constexpr size_t CACHE_LINE_SIZE = 64;
static constexpr size_t HASH_TABLE_SIZE = 2097152; // 2M entries, power of 2
static constexpr uint64_t HASH_EMPTY = 0;
static constexpr int MAX_BATCH_SIZE = 1000;

// --------------------- Common structs ---------------------
struct Record { 
    string header; 
    string seq; 
    
    Record() = default;
    Record(string h, string s) : header(std::move(h)), seq(std::move(s)) {}
};

// --------------------- Fast string utilities ---------------------
static inline string trim(const string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

static inline string kslice(const string& s, int k) {
    if ((int)s.size() <= k) return s;
    return s.substr(0, k);
}

// Fast character conversion using lookup table
static const array<char, 256> TOUPPER_TABLE = []() {
    array<char, 256> table;
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
static inline uint64_t hash_sequence_fast(const string& s) {
    const uint64_t PRIME1 = 11400714785074694791ULL;
    const uint64_t PRIME2 = 14029467366897019727ULL;
    const uint64_t PRIME3 = 1609587929392839161ULL;
    const uint64_t PRIME4 = 9650029242287828579ULL;
    const uint64_t PRIME5 = 2870177450012600261ULL;
    
    uint64_t hash = PRIME5;
    const char* data = s.data();
    size_t len = s.size();
    
    // Process 8 bytes at a time
    while (len >= 8) {
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

// Rolling hash for k-mers
class RollingHasher {
private:
    static constexpr uint64_t BASE = 257;
    static constexpr uint64_t MOD = (1ULL << 61) - 1;
    uint64_t base_pow;
    int k;
    
public:
    RollingHasher(int kmer_size) : k(kmer_size), base_pow(1) {
        for (int i = 0; i < k - 1; ++i) {
            base_pow = (base_pow * BASE) % MOD;
        }
    }
    
    uint64_t hash_kmer(const string& kmer) {
        uint64_t hash = 0;
        for (char c : kmer) {
            hash = (hash * BASE + static_cast<unsigned char>(c)) % MOD;
        }
        return hash == HASH_EMPTY ? 1 : hash;
    }
};

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
    MemoryMappedFASTA(const string& filename) : data(nullptr), file_size(0) {
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
        
        struct stat sb;
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

    vector<Record> parse() {
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

// --------------------- Lock-free hash set ---------------------
class alignas(CACHE_LINE_SIZE) LockFreeHashSet {
private:
    vector<atomic<uint64_t>> table;
    size_t mask;
    
public:
    LockFreeHashSet() : table(HASH_TABLE_SIZE), mask(HASH_TABLE_SIZE - 1) {
        for (auto& slot : table) {
            slot.store(HASH_EMPTY, memory_order_relaxed);
        }
    }
    
    // Returns true if hash was newly inserted, false if it was already present
    bool insert_if_new(uint64_t hash) {
        if (hash == HASH_EMPTY) hash = 1;
        
        size_t pos = hash & mask;
        
        for (size_t i = 0; i < HASH_TABLE_SIZE; ++i) {
            uint64_t expected = HASH_EMPTY;
            
            // Try to insert into empty slot
            if (table[pos].compare_exchange_weak(expected, hash, memory_order_acq_rel)) {
                return true; // Successfully inserted
            }
            
            // Check if already exists
            if (table[pos].load(memory_order_acquire) == hash) {
                return false; // Already exists
            }
            
            // Linear probing
            pos = (pos + 1) & mask;
        }
        
        // Hash table full (very unlikely) - assume it's new to avoid data loss
        return true;
    }
};

// --------------------- Optimized edit distance ---------------------
class OptimizedEditDistance {
public:
    // Quick rejection based on length difference and character frequency
    static bool quick_reject(const string& a, const string& b, int threshold_pct) {
        int max_len = max(a.size(), b.size());
        int allowed_edits = (int)ceil((100.0 - threshold_pct) * max_len / 100.0);
        
        // Length difference check
        if (abs((int)a.size() - (int)b.size()) > allowed_edits) {
            return true;
        }
        
        // For high similarity thresholds, do character frequency check
        if (threshold_pct >= 90 && max_len > 50) {
            array<int, 256> freq_diff = {0};
            
            for (char c : a) freq_diff[static_cast<unsigned char>(c)]++;
            for (char c : b) freq_diff[static_cast<unsigned char>(c)]--;
            
            int total_diff = 0;
            for (int f : freq_diff) {
                total_diff += abs(f);
            }
            
            if (total_diff / 2 > allowed_edits) {
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
            __m128i v1 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(p1));
            __m128i v2 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(p2));
            __m128i cmp = _mm_cmpeq_epi8(v1, v2);
            
            if (_mm_movemask_epi8(cmp) != 0xFFFF) {
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
    static int compute_banded(const string& a, const string& b, int max_edits) {
        if (quick_reject(a, b, 100 * (1.0 - (double)max_edits / max(a.size(), b.size())))) {
            return max_edits + 1;
        }
        
        int n = a.size(), m = b.size();
        if (abs(n - m) > max_edits) return max_edits + 1;
        
        const int INF = max_edits + 1;
        int band = max_edits;
        
        vector<int> prev, cur;
        int jmin_prev = max(0, 0 - band);
        int jmax_prev = min(m, 0 + band);
        prev.assign(jmax_prev - jmin_prev + 1, INF);
        
        for (int j = jmin_prev; j <= jmax_prev; ++j) {
            prev[j - jmin_prev] = j;
            if (prev[j - jmin_prev] > max_edits) prev[j - jmin_prev] = INF;
        }
        
        for (int i = 1; i <= n; ++i) {
            int jmin_cur = max(0, i - band);
            int jmax_cur = min(m, i + band);
            cur.assign(jmax_cur - jmin_cur + 1, INF);
            int best = INF;
            
            for (int j = jmin_cur; j <= jmax_cur; ++j) {
                int cost_sub = INF, cost_del = INF, cost_ins = INF;
                
                // Substitution
                if (j - 1 >= jmin_prev && j - 1 <= jmax_prev) {
                    int pv = prev[(j - 1) - jmin_prev];
                    if (pv != INF) cost_sub = pv + (a[i-1] == b[j-1] ? 0 : 1);
                }
                
                // Deletion
                if (j >= jmin_prev && j <= jmax_prev) {
                    int pv = prev[j - jmin_prev];
                    if (pv != INF) cost_del = pv + 1;
                }
                
                // Insertion
                if (j - 1 >= jmin_cur && j - 1 <= jmax_cur) {
                    int cv = cur[(j - 1) - jmin_cur];
                    if (cv != INF) cost_ins = cv + 1;
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

static inline bool identity_at_least(const string& a, const string& b, int pct_threshold) {
    int max_len = max(a.size(), b.size());
    int allowed = (int)ceil((100.0 - pct_threshold) * max_len / 100.0);
    
    if (abs((int)a.size() - (int)b.size()) > allowed) return false;
    
#ifdef __SSE4_2__
    if (pct_threshold == 100) {
        return OptimizedEditDistance::sequences_equal_simd(a, b);
    }
#endif
    
    int ed = OptimizedEditDistance::compute_banded(a, b, allowed);
    if (ed > allowed) return false;
    
    double pid = 100.0 * (1.0 - (double)ed / (double)max_len);
    return pid + 1e-9 >= pct_threshold;
}

// --------------------- Thread-safe queue with batching ---------------------
template<typename T>
class BatchingQueue {
private:
    deque<T> dq;
    mutex mx;
    condition_variable cv;
    bool done = false;
    size_t batch_size;
    
public:
    BatchingQueue(size_t batch_sz = MAX_BATCH_SIZE) : batch_size(batch_sz) {}
    
    void push(T v) {
        {
            lock_guard<mutex> lk(mx);
            dq.push_back(std::move(v));
        }
        cv.notify_one();
    }
    
    bool pop_batch(vector<T>& batch) {
        unique_lock<mutex> lk(mx);
        cv.wait(lk, [&]{ return done || !dq.empty(); });
        
        if (dq.empty()) return false;
        
        batch.clear();
        batch.reserve(batch_size);
        
        size_t count = min(batch_size, dq.size());
        for (size_t i = 0; i < count; ++i) {
            batch.push_back(std::move(dq.front()));
            dq.pop_front();
        }
        
        return true;
    }
    
    void close() {
        {
            lock_guard<mutex> lk(mx);
            done = true;
        }
        cv.notify_all();
    }
};

// --------------------- Processing state ---------------------
struct ProcessingState {
    vector<atomic<bool>> is_kept;
    vector<atomic<int>> kept_order;
    atomic<int> next_order{0};
    atomic<size_t> removed{0};
    
    ProcessingState(size_t n) : is_kept(n), kept_order(n) {
        for (size_t i = 0; i < n; ++i) {
            is_kept[i].store(false, memory_order_relaxed);
            kept_order[i].store(-1, memory_order_relaxed);
        }
    }
};

// --------------------- Sharded index for <100% identity ---------------------
struct alignas(CACHE_LINE_SIZE) Shard {
    unordered_map<string, vector<int>> prefix_index;
    shared_mutex mtx;
    
    Shard() {
        prefix_index.reserve(1024);
    }
};

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

Args parse_args(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; i++) {
        string s = argv[i];
        if (s == "-i" && i + 1 < argc) a.in_path = argv[++i];
        else if (s == "-o" && i + 1 < argc) a.out_path = argv[++i];
        else if (s == "-p" && i + 1 < argc) a.pct = stoi(argv[++i]);
        else if (s == "-k" && i + 1 < argc) a.kmer = stoi(argv[++i]);
        else if (s == "-t" && i + 1 < argc) a.threads = stoi(argv[++i]);
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
    
    // Initialize processing state
    ProcessingState state(recs.size());
    
    // Work queue for batched processing
    BatchingQueue<int> work_queue;
    
    if (args.pct == 100) {
        // 100% identity case - use lock-free hash set
        LockFreeHashSet global_hash_set;
        
        auto worker = [&](int tid) {
            vector<int> batch;
            while (work_queue.pop_batch(batch)) {
                for (int idx : batch) {
                    uint64_t hash = hash_sequence_fast(recs[idx].seq);
                    
                    if (global_hash_set.insert_if_new(hash)) {
                        // New sequence
                        state.is_kept[idx].store(true, memory_order_relaxed);
                        state.kept_order[idx].store(
                            state.next_order.fetch_add(1, memory_order_relaxed), 
                            memory_order_relaxed
                        );
                    } else {
                        // Duplicate
                        state.removed.fetch_add(1, memory_order_relaxed);
                    }
                }
            }
        };
        
        // Launch workers
        vector<thread> workers;
        workers.reserve(args.threads);
        for (int t = 0; t < args.threads; ++t) {
            workers.emplace_back(worker, t);
        }
        
        // Enqueue work
        for (int i = 0; i < (int)recs.size(); ++i) {
            work_queue.push(i);
        }
        work_queue.close();
        
        // Wait for completion
        for (auto& w : workers) {
            w.join();
        }
        
    } else {
        // <100% identity case - use sharded indexing
        int shard_count = max(args.threads * 4, 16);
        vector<Shard> shards(shard_count);
        
        auto shard_for_key = [&](const string& key) -> size_t {
            return hash<string>{}(key) % shard_count;
        };
        
        auto worker = [&](int tid) {
            RollingHasher hasher(args.kmer);
            vector<int> batch;
            
            while (work_queue.pop_batch(batch)) {
                for (int idx : batch) {
                    const string& s = recs[idx].seq;
                    string key = kslice(s, args.kmer);
                    size_t shard_id = shard_for_key(key);
                    Shard& shard = shards[shard_id];
                    
                    bool is_dup = false;
                    
                    // Check against candidates in this shard
                    {
                        shared_lock<shared_mutex> rl(shard.mtx);
                        auto it = shard.prefix_index.find(key);
                        if (it != shard.prefix_index.end()) {
                            for (int candidate_idx : it->second) {
                                if (state.is_kept[candidate_idx].load(memory_order_acquire)) {
                                    if (identity_at_least(s, recs[candidate_idx].seq, args.pct)) {
                                        is_dup = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    
                    if (is_dup) {
                        state.removed.fetch_add(1, memory_order_relaxed);
                    } else {
                        // Add to index and mark as kept
                        {
                            unique_lock<shared_mutex> wl(shard.mtx);
                            shard.prefix_index[key].push_back(idx);
                        }
                        
                        state.is_kept[idx].store(true, memory_order_relaxed);
                        state.kept_order[idx].store(
                            state.next_order.fetch_add(1, memory_order_relaxed),
                            memory_order_relaxed
                        );
                    }
                }
            }
        };
        
        // Launch workers
        vector<thread> workers;
        workers.reserve(args.threads);
        for (int t = 0; t < args.threads; ++t) {
            workers.emplace_back(worker, t);
        }
        
        // Enqueue work
        for (int i = 0; i < (int)recs.size(); ++i) {
            work_queue.push(i);
        }
        work_queue.close();
        
        // Wait for completion
        for (auto& w : workers) {
            w.join();
        }
    }
    
    // Collect results in order
    vector<pair<int, int>> kept_with_order; // (order, index)
    kept_with_order.reserve(recs.size());
    
    for (size_t i = 0; i < recs.size(); ++i) {
        if (state.is_kept[i].load(memory_order_relaxed)) {
            int order = state.kept_order[i].load(memory_order_relaxed);
            kept_with_order.emplace_back(order, i);
        }
    }
    
    // Sort by order to maintain original sequence
    sort(kept_with_order.begin(), kept_with_order.end());
    
    // Build final result
    vector<Record> kept_records;
    kept_records.reserve(kept_with_order.size());
    for (const auto& [order, idx] : kept_with_order) {
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
    size_t removed_count = state.removed.load(memory_order_relaxed);
    
    cerr << "Input sequences:  " << recs.size() << "\n";
    cerr << "Kept sequences:   " << kept_count << "\n";
    cerr << "Removed sequences:" << removed_count << "\n";
    cerr << "Threshold:        " << args.pct << "% identity\n";
    if (args.pct < 100) cerr << "Prefix k-mer:     " << args.kmer << "\n";
    cerr << "Threads:          " << args.threads << "\n";
    
    return 0;
}
