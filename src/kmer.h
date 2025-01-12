#ifndef TROW_KMER_H
#define TROW_KMER_H

#ifdef _MSC_VER
#define NOMINMAX
#endif

#define LENGTH (1 << 17)
#define WINSIZE (1U << 15)  // sliding window size
#define CHUNK (1 << 15)

#define SPAN (1L << 20)       // desired distance between access points

#define MAX_SEQ 1000

#define ABS_MAX_ANS_NUM 10

#define ABS_MIN_DNA_COUNT 1
#define ABS_MIN_PRINT_COUNT 10
#define ABS_MIN_ANS_COUNT 20

#define NUM_FOR_MAX_COUNT 4
#define NUM_TOT_MAX_COUNT 4
#define NUM_RAT_MAX_COUNT 4
#define NUM_RAT_CAND 20

#define ABS_MIN_MER 3
#define ABS_TABLE_MAX_MER 15
#define ABS_UINT64_MAX_MER 32
#define ABS_MAX_MER 64

// k_mer_data macro
#define K_MER_DATA_COUNT 0
#define K_MER_DATA_MAX 1
#define K_MER_DATA_MAX_SEQ 2

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#include <absl/numeric/int128.h>
#include <absl/container/flat_hash_map.h>

#include <tuple>
#include <cstdint>
#include <cinttypes>
#include <utility>

#include <tbb/concurrent_queue.h>
#include <tbb/tbb.h>

#include "zran.h"
#include <zlib.h>

extern int MIN_MER;
extern int MAX_MER;
extern int TABLE_MAX_MER;
extern int NUM_THREAD;
extern int SLICE_LENGTH;
extern int QUEUE_SIZE;
extern double LOW_BASELINE;
extern double HIGH_BASELINE;
extern bool INDEX;

template <typename T>
struct FinalData {
    T forward;
    T backward;
    T both;
};

typedef std::pair<int, int> KmerData;
typedef std::vector<std::pair<int, int>> LocationVector;
typedef absl::uint128 uint128_t;

struct FastqLocData {
    std::pair<int64_t, int64_t> pos;
    KmerData rht_kmer;
    KmerData lef_kmer;
    std::pair<uint128_t, uint128_t> rht_seq;
    std::pair<uint128_t, uint128_t> lef_seq;

    FastqLocData(const std::pair<int64_t, int64_t>& pos_,
                 KmerData  rht_kmer_,
                 KmerData  lef_kmer_,
                 const std::pair<uint128_t, uint128_t>& rht_seq_,
                 const std::pair<uint128_t, uint128_t>& lef_seq_)
            : pos(pos_),
              rht_kmer(std::move(rht_kmer_)),
              lef_kmer(std::move(lef_kmer_)),
              rht_seq(rht_seq_),
              lef_seq(lef_seq_) {}
};

typedef std::tuple<uint16_t**, uint64_t**, uint128_t**, uint32_t**> ThreadDataTuple;
typedef std::pair<int, uint128_t> KmerSeq;

struct KmerAnsData {
    KmerSeq first;
    FinalData<int64_t> second;
    bool dir;

    KmerAnsData(KmerSeq  first_,
                const FinalData<int64_t>& second_,
                bool dir_)
            : first(std::move(first_)),
              second(second_),
              dir(dir_) {}
};

typedef absl::flat_hash_map<KmerSeq, uint32_t> ResultMap;
typedef std::pair<ResultMap*, ResultMap*> ResultMapPair;
typedef FinalData<ResultMapPair> ResultMapData;
typedef std::pair<ResultMapData, std::vector<FastqLocData>*> ResultMapPairData;

typedef absl::flat_hash_map<uint64_t, uint16_t> CounterMap;
typedef absl::flat_hash_map<uint128_t, uint16_t> CounterMap_128;

typedef absl::flat_hash_map<KmerSeq, FinalData<int64_t>> FinalFastqData;
typedef std::vector<KmerAnsData> FinalFastqVector;
typedef std::vector<std::pair<KmerSeq, FinalData<int64_t>>> LeagcyFinalFastqVector;
typedef std::pair<FinalFastqVector*, FinalFastqVector*> FinalFastqVectorPair;

struct QueueData {
    char* buffer;
    LocationVector* loc_vector;
    int64_t buf_off;
};

struct FinalFastqOutput {
    FinalFastqVector* high;
    FinalFastqVector* low;
    std::vector<FastqLocData>* k_mer_loc_vector;
};

typedef tbb::concurrent_bounded_queue<QueueData> TBBQueue;
typedef tbb::task_group TaskGroup;

uint8_t** set_repeat_check_table();
uint32_t** set_rotation_table(uint8_t** repeat_check_table);

uint64_t* set_extract_k_mer();
uint128_t* set_extract_k_mer_128();
uint128_t* set_extract_k_mer_ans();

void fill_rotation_table(uint32_t* rot, int k, uint32_t seq, uint8_t** repeat_check_table);

uint16_t** set_k_mer_counter();

uint64_t** set_k_mer_data();
uint128_t** set_k_mer_data_128();

uint32_t** set_k_mer_counter_list();

class ThreadData {
private:
    uint16_t **k_mer_counter = nullptr;
    uint64_t **k_mer_data = nullptr;
    uint128_t **k_mer_data_128 = nullptr;
    uint32_t **k_mer_counter_list = nullptr;
public:
    ThreadDataTuple init_check() {
        if (k_mer_counter == nullptr) {
            if (MAX_MER <= ABS_UINT64_MAX_MER) {
                k_mer_data = set_k_mer_data();
            } else {
                k_mer_data_128 = set_k_mer_data_128();
            }

            if (MIN_MER <= TABLE_MAX_MER) {
                k_mer_counter = set_k_mer_counter();
                k_mer_counter_list = set_k_mer_counter_list();
            }
        }
        return ThreadDataTuple {k_mer_counter, k_mer_data, k_mer_data_128, k_mer_counter_list};
    }
};

ResultMapPairData buffer_task(TBBQueue* task_queue, ThreadData* thread_data, uint32_t** rot_table, const uint64_t *extract_k_mer, const uint128_t *extract_k_mer_128, uint8_t** repeat_check_table);

void read_fastq_thread(FILE* fp, TBBQueue* buffer_task_queue);
void read_fastq_gz_thread(FILE* fp, gz_index **built, TBBQueue* buffer_task_queue);

ResultMapPairData buffer_task_long(TBBQueue* task_queue, ThreadData* thread_data, uint32_t** rot_table, const uint64_t *extract_k_mer, const uint128_t *extract_k_mer_128, uint8_t** repeat_check_table);

void read_fastq_thread_long(FILE* fp, TBBQueue* buffer_task_queue);
void read_fastq_gz_thread_long(FILE* fp, gz_index **built, TBBQueue* buffer_task_queue);

void int_to_four(char* buffer, uint128_t seq, int n);

FinalFastqOutput process_kmer(const char* file_name, uint8_t **repeat_check_table, uint32_t **rot_table,
                                  uint64_t *extract_k_mer, uint128_t *extract_k_mer_128, uint128_t *extract_k_mer_ans,
                                  ThreadData* thread_data_list, bool is_gz, gz_index **built);

FinalFastqOutput process_kmer_long(const char* file_name, uint8_t **repeat_check_table, uint32_t **rot_table,
                                    uint64_t *extract_k_mer, uint128_t *extract_k_mer_128, uint128_t *extract_k_mer_ans,
                                    ThreadData* thread_data_list, bool is_gz, gz_index **built);

FinalFastqOutput process_output(const char* file_name, ResultMapPairData* result_list, uint32_t **rot_table, uint128_t *extract_k_mer_ans);

KmerData k_mer_check(const char* seq, int st, int nd, uint32_t** rot_table, const uint64_t *extract_k_mer,
                     uint16_t** k_mer_counter, CounterMap* k_mer_counter_map,
                     uint64_t** k_mer_data, uint32_t** k_mer_counter_list, uint8_t** repeat_check_table,
                     ResultMapPair result_pair,
                     int16_t* k_mer_total_cnt, int min_mer, int max_mer, std::pair<uint64_t, uint64_t>* repeat_seq = nullptr);

KmerData k_mer_check_128(const char* seq, int st, int nd, uint32_t** rot_table, const uint128_t *extract_k_mer,
                         uint16_t** k_mer_counter, CounterMap_128* k_mer_counter_map,
                         uint128_t** k_mer_data, uint32_t** k_mer_counter_list, uint8_t** repeat_check_table,
                         ResultMapPair result_pair,
                         int16_t* k_mer_total_cnt, int min_mer, int max_mer,  std::pair<uint128_t, uint128_t>* repeat_seq = nullptr);

void k_mer_target(const char* seq, int st, int nd, uint32_t** rot_table, const uint64_t *extract_k_mer,
                  uint16_t** k_mer_counter, CounterMap* k_mer_counter_map,
                  uint64_t** k_mer_data, uint32_t** k_mer_counter_list, uint8_t** repeat_check_table,
                  ResultMapPair result_pair,
                  int16_t* k_mer_total_cnt,
                  int k);

void k_mer_target_128(const char* seq, int st, int nd, uint32_t** rot_table, const uint128_t *extract_k_mer,
                      uint16_t** k_mer_counter, CounterMap_128* k_mer_counter_map,
                      uint128_t** k_mer_data, uint32_t** k_mer_counter_list, uint8_t** repeat_check_table,
                      ResultMapPair result_pair,
                      int16_t* k_mer_total_cnt,
                      int k);

bool check_ans_seq(const KmerSeq& seq, uint128_t* extract_k_mer_ans, uint32_t** rot_table);

uint64_t get_rot_seq(uint64_t seq, int k);
uint128_t get_rot_seq_128(const uint128_t& seq, int k);

uint8_t get_repeat_check(uint64_t seq, int k);
uint8_t get_repeat_check(uint128_t seq, int k);
int get_dna_count(uint128_t seq, int k);

uint32_t reverse_complement_32(uint32_t x);
uint64_t reverse_complement_64(uint64_t x);
uint128_t reverse_complement_128(uint128_t x);

FinalData<int64_t> add_data(FinalData<int64_t> a, FinalData<int64_t> b);

std::vector<std::pair<KmerSeq, int>>* final_process_output(FinalFastqData* total_result_high, FinalFastqData* total_result_low);
ResultMap* get_score_map(FinalFastqData* total_result);

ptrdiff_t deflate_index_extract(FILE *in, gz_index *index,
                                off_t offset, unsigned char *buf, size_t len);

ptrdiff_t file_extract(FILE *in, off_t offset, unsigned char *buf, size_t len);

void deflate_index_free(gz_index *index);

void get_trm_read(const std::filesystem::path &fastq_path, std::vector<std::pair<KmerSeq, int>>* put_trm,
                   FinalFastqOutput fastq_file_data, gz_index* index, size_t st, size_t nd, char* temp_path);

gz_index *get_thread_safe_index(gz_index* index);
#endif //TROW_KMER_H