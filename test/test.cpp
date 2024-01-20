#define BOOST_TEST_MODULE TREW_TEST
#include "../src/kmer.h"
#include <string>
#include <filesystem>
#include <boost/test/included/unit_test.hpp>
namespace fs = std::filesystem;


extern char trans_arr[4];
#define O (-1)
#define A 3
#define C 2
#define G 1
#define T 0

static const int codes[256] = {
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, A, O, C, O, O, O, G, O, O, O, O, O, O, O, O,
        O, O, O, O, T, O, O, O, O, O, O, O, O, O, O, O,
        O, A, O, C, O, O, O, G, O, O, O, O, O, O, O, O,
        O, O, O, O, T, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O,
        O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, O
};

#undef O
#undef A
#undef C
#undef G
#undef T

int MAX_MER;
int MIN_MER;
int TABLE_MAX_MER;
int NUM_THREAD;
int SLICE_LENGTH;
int QUEUE_SIZE = -1;
double BASELINE;

// test function
std::string operator*(std::string& s, size_t n) {
    std::string result;
    result.reserve(s.size() * n);
    for(size_t i = 0; i < n; ++i) {
        result += s;
    }
    return result;
}

uint64_t four_to_int(std::string seq) {
    uint64_t tmp = 0;
    for (int i = 0; i <= seq.length() - 1; i++) {
        tmp *= 4;
        tmp += codes[seq[i]];
    }

    return tmp;
}

uint128_t four_to_int_128(std::string seq) {
    uint128_t tmp = 0;
    for (int i = 0; i <= seq.length() - 1; i++) {
        tmp *= 4;
        tmp += codes[seq[i]];
    }

    return tmp;
}

BOOST_AUTO_TEST_CASE(get_rot_seq_test) {
    std::string bef, aft;

    bef = "ATATATTTT";
    aft = "TTTTATATA";
    BOOST_CHECK_EQUAL(get_rot_seq(four_to_int(bef), bef.length()), four_to_int(aft));

    bef = "GCGACTTGACGC";
    aft = "TTGACGCGCGAC";
    BOOST_CHECK_EQUAL(get_rot_seq(four_to_int(bef), bef.length()), four_to_int(aft));

    bef = "GGGGGGGTGGG";
    aft = "TGGGGGGGGGG";
    BOOST_CHECK_EQUAL(get_rot_seq(four_to_int(bef), bef.length()), four_to_int(aft));
}

BOOST_AUTO_TEST_CASE(get_repeat_check_test) {
    std::string buffer_s;
    buffer_s = "ATTTTTTT";
    BOOST_CHECK_EQUAL(get_repeat_check(four_to_int(buffer_s), buffer_s.length()), 1);

    buffer_s = "ATTTTTTTGC";
    BOOST_CHECK_EQUAL(get_repeat_check(four_to_int(buffer_s), buffer_s.length()), 0);

    buffer_s = "ATTATAGCGATCGTCACCATTGC";
    BOOST_CHECK_EQUAL(get_repeat_check(four_to_int(buffer_s), buffer_s.length()), 0);
}

BOOST_AUTO_TEST_CASE(k_mer_total_cnt_test) {
    MAX_MER = 21;
    MIN_MER = 5;
    TABLE_MAX_MER = 13;
    BASELINE = 0.8;
    NUM_THREAD = 1;

    std::string buffer_s = "ATGCATCACACTCGCCGATGCATCACNNNNNNNNNGCCGATGCATCACACTCGCCGNTGCATCACACTCGCCGATGCATC"
                           "ACACTCGCCGATGCATCACANNNGCCGATGCATCACACNNGCCGATGCATCACACTCNNCCGATGCATCACACTCGCCGA";
    const char* buffer = buffer_s.c_str();

    int16_t* k_mer_total_cnt = (int16_t*) calloc(MAX_MER - MIN_MER + 2, sizeof(int16_t));
    if (k_mer_total_cnt == nullptr) {
        fprintf(stderr, "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    int st = 0;
    int nd = (int) buffer_s.length() - 1;
    int len;
    int chr_code;
    int err_loc;

    err_loc = st - 1;
    for (int i = st ; i <= nd; i++) {
        chr_code = codes[buffer[i]];
        if (chr_code >= 0) {
            // i - err_loc가 i로 끝나는 길이
            len = MIN(i - err_loc, MAX_MER);
            // MIN_MER ~ len에 k-total count가 증가한다
            if(len >= MIN_MER){
                k_mer_total_cnt[0]++;
                k_mer_total_cnt[len - MIN_MER + 1]--;
            }
        } else {
            err_loc = i;
        }
    }

    for(int i = 1; i <= MAX_MER - MIN_MER; i++) {
        k_mer_total_cnt[i] = (int16_t)(k_mer_total_cnt[i] + k_mer_total_cnt[i - 1]);
    }

    for (int k = MIN_MER; k <= MAX_MER; k++) {
        err_loc = st - 1;
        for (int i = st ; i <= nd; i++) {
            chr_code = codes[buffer[i]];
            if (chr_code >= 0) {
                if(i - err_loc >= k) {
                    k_mer_total_cnt[k - MIN_MER]--;
                }
            } else {
                err_loc = i;
            }
        }
    }

    for (int i = 0; i <= MAX_MER - MIN_MER; i ++) {
        BOOST_CHECK_EQUAL(k_mer_total_cnt[i], 0);
    }
}

BOOST_AUTO_TEST_CASE(k_mer_test) {
    MAX_MER = 32;
    MIN_MER = 5;
    TABLE_MAX_MER = 13;
    BASELINE = 0.8;
    NUM_THREAD = 1;

    uint8_t **repeat_check_table = set_repeat_check_table();
    uint32_t **rot_table = set_rotation_table(repeat_check_table);
    uint64_t *extract_k_mer = set_extract_k_mer();
    uint16_t **k_mer_counter = set_k_mer_counter();
    uint64_t **k_mer_data = set_k_mer_data();
    uint32_t **k_mer_counter_list = set_k_mer_counter_list();

    CounterMap* k_mer_counter_map = new CounterMap[MAX_MER - TABLE_MAX_MER];
    ResultMap* result = new ResultMap {};
    int16_t* k_mer_total_cnt = (int16_t*) malloc(sizeof(int16_t) * (MAX_MER - MIN_MER + 2));

    std::vector<std::string> repeat_list = {"TTGCATCACACCCTCGCCG", "TTAGGG", "TTAGAGCCCACA", "TTTTGCCCTCATCACACCCTCGCCTCCTTCGC"};
    for(auto& repeat : repeat_list) {
        int repeat_len = (int) repeat.length();
        std::string buffer = repeat * 20;
        uint64_t min_repeat_seq = four_to_int(repeat);
        min_repeat_seq = MIN(min_repeat_seq, get_rot_seq(reverse_complement_64(min_repeat_seq >> (2 * (32 - repeat_len))), repeat_len));

        k_mer(buffer.c_str(), 0, (int) buffer.length() - 1, rot_table, extract_k_mer, k_mer_counter, k_mer_counter_map,
              k_mer_data, k_mer_counter_list, repeat_check_table, result, k_mer_total_cnt);

        BOOST_CHECK_EQUAL(result -> size(), 1);
        for(auto& [k, v] : (*result)) {
            BOOST_CHECK_EQUAL(repeat_len, k.first);
            BOOST_CHECK_EQUAL(min_repeat_seq, k.second);
            BOOST_CHECK_EQUAL(repeat_len * (20 - 1) + 1, v);
        }
        result -> clear();
    }
}

BOOST_AUTO_TEST_CASE(k_mer_128_test) {
    MAX_MER = 64;
    MIN_MER = 5;
    TABLE_MAX_MER = 13;
    BASELINE = 0.8;
    NUM_THREAD = 1;

    uint8_t **repeat_check_table = set_repeat_check_table();
    uint32_t **rot_table = set_rotation_table(repeat_check_table);
    uint128_t *extract_k_mer = set_extract_k_mer_128();
    uint16_t **k_mer_counter = set_k_mer_counter();
    uint128_t **k_mer_data = set_k_mer_data_128();
    uint32_t **k_mer_counter_list = set_k_mer_counter_list();

    CounterMap_128* k_mer_counter_map = new CounterMap_128[MAX_MER - TABLE_MAX_MER];
    ResultMap* result = new ResultMap {};
    int16_t* k_mer_total_cnt = (int16_t*) malloc(sizeof(int16_t) * (MAX_MER - MIN_MER + 2));

    std::vector<std::string> repeat_list = {"TGCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "TTAGGG", "TTAGAGCCCACA", "TTTTGCCCTCATCACACCCTCGCCTCCTTCGC", "TTTTGCCCTCATCACACCCTCGCCTCCTTCGTGCTTGCCCCCACACTGACTGACGTGCAGTCTG"};
    for(auto& repeat : repeat_list) {
        int repeat_len = (int) repeat.length();
        std::string buffer = repeat * 10;
        uint128_t min_repeat_seq = four_to_int_128(repeat);
        min_repeat_seq = MIN(min_repeat_seq, get_rot_seq_128(reverse_complement_128(min_repeat_seq) >> (2 * (64 - repeat_len)), repeat_len));

        k_mer_128(buffer.c_str(), 0, (int) buffer.length() - 1, rot_table, extract_k_mer, k_mer_counter, k_mer_counter_map,
                  k_mer_data, k_mer_counter_list, repeat_check_table, result, k_mer_total_cnt);

        BOOST_CHECK_EQUAL(result -> size(), 1);
        for(auto& [k, v] : (*result)) {
            BOOST_CHECK_EQUAL(repeat_len, k.first);
            BOOST_CHECK_EQUAL(min_repeat_seq, k.second);
            BOOST_CHECK_EQUAL(repeat_len * (10 - 1) + 1, v);
        }
        result -> clear();
    }
}

BOOST_AUTO_TEST_CASE(main_test_32) {
    MAX_MER = 32;
    MIN_MER = 5;
    TABLE_MAX_MER = 1;
    BASELINE = 0.8;

    uint8_t **repeat_check_table = nullptr;
    uint32_t **rot_table = nullptr;

    uint64_t *extract_k_mer = nullptr;
    uint128_t *extract_k_mer_128 = nullptr;
    if (MAX_MER <= ABS_UINT64_MAX_MER) {
        extract_k_mer = set_extract_k_mer();
    }
    else {
        extract_k_mer_128 = set_extract_k_mer_128();
    }

    std::vector<std::pair<fs::path, bool>> FILE_LOC_VECTOR;
    if (boost::unit_test::framework::master_test_suite().argc == 3) {
        FILE_LOC_VECTOR = {{fs::path(boost::unit_test::framework::master_test_suite().argv[1]), true},
                           {fs::path(boost::unit_test::framework::master_test_suite().argv[2]), false}};
    }
    else {
        FILE_LOC_VECTOR = {{fs::current_path().parent_path() / fs::path("test") / fs::path("test.fastq.gz"), true},
                           {fs::current_path().parent_path() / fs::path("test") / fs::path("test.fastq"), false}};
    }

    std::vector<int> NUM_THREAD_VECTOR = {1, 2};

    ThreadData* thread_data_list = new ThreadData[*std::max_element(NUM_THREAD_VECTOR.begin(), NUM_THREAD_VECTOR.end())];
    for (auto& [fastq_path, is_gz] : FILE_LOC_VECTOR) {
        for(auto& num_thread : NUM_THREAD_VECTOR) {
            NUM_THREAD = num_thread;
            BOOST_CHECK_NO_THROW(process_kmer(fastq_path.string().c_str(), repeat_check_table, rot_table,
                                 extract_k_mer, extract_k_mer_128,
                                 thread_data_list, is_gz));
        }
    }
    delete[] thread_data_list;
}

BOOST_AUTO_TEST_CASE(main_test_64) {
    MAX_MER = 64;
    MIN_MER = 5;
    TABLE_MAX_MER = 1;
    BASELINE = 0.8;

    uint8_t **repeat_check_table = nullptr;
    uint32_t **rot_table = nullptr;

    uint64_t *extract_k_mer = nullptr;
    uint128_t *extract_k_mer_128 = nullptr;
    if (MAX_MER <= ABS_UINT64_MAX_MER) {
        extract_k_mer = set_extract_k_mer();
    }
    else {
        extract_k_mer_128 = set_extract_k_mer_128();
    }

    std::vector<std::pair<fs::path, bool>> FILE_LOC_VECTOR;
    if (boost::unit_test::framework::master_test_suite().argc == 3) {
        FILE_LOC_VECTOR = {{fs::path(boost::unit_test::framework::master_test_suite().argv[1]), true},
                           {fs::path(boost::unit_test::framework::master_test_suite().argv[2]), false}};
    }
    else {
        FILE_LOC_VECTOR = {{fs::current_path().parent_path() / fs::path("test") / fs::path("test.fastq.gz"), true},
                           {fs::current_path().parent_path() / fs::path("test") / fs::path("test.fastq"), false}};
    }

    std::vector<int> NUM_THREAD_VECTOR = {1, 2};

    ThreadData* thread_data_list = new ThreadData[*std::max_element(NUM_THREAD_VECTOR.begin(), NUM_THREAD_VECTOR.end())];
    for (auto& [fastq_path, is_gz] : FILE_LOC_VECTOR) {
        for(auto& num_thread : NUM_THREAD_VECTOR) {
            NUM_THREAD = num_thread;
            BOOST_CHECK_NO_THROW(process_kmer(fastq_path.string().c_str(), repeat_check_table, rot_table,
                                              extract_k_mer, extract_k_mer_128,
                                              thread_data_list, is_gz));
        }
    }
    delete[] thread_data_list;
}

BOOST_AUTO_TEST_CASE(main_test_long_32) {
    MAX_MER = 32;
    MIN_MER = 5;
    TABLE_MAX_MER = 1;
    SLICE_LENGTH = 150;
    BASELINE = 0.8;

    uint8_t **repeat_check_table = nullptr;
    uint32_t **rot_table = nullptr;

    uint64_t *extract_k_mer = nullptr;
    uint128_t *extract_k_mer_128 = nullptr;
    if (MAX_MER <= ABS_UINT64_MAX_MER) {
        extract_k_mer = set_extract_k_mer();
    }
    else {
        extract_k_mer_128 = set_extract_k_mer_128();
    }

    std::vector<std::pair<fs::path, bool>> FILE_LOC_VECTOR;
    if (boost::unit_test::framework::master_test_suite().argc == 3) {
        FILE_LOC_VECTOR = {{fs::path(boost::unit_test::framework::master_test_suite().argv[3]), true},
                           {fs::path(boost::unit_test::framework::master_test_suite().argv[4]), false}};
    }
    else {
        FILE_LOC_VECTOR = {{fs::current_path().parent_path() / fs::path("test") / fs::path("test_long.fastq.gz"), true},
                           {fs::current_path().parent_path() / fs::path("test") / fs::path("test_long.fastq"), false}};
    }

    std::vector<int> NUM_THREAD_VECTOR = {1, 2};

    ThreadData* thread_data_list = new ThreadData[*std::max_element(NUM_THREAD_VECTOR.begin(), NUM_THREAD_VECTOR.end())];
    for (auto& [fastq_path, is_gz] : FILE_LOC_VECTOR) {
        for(auto& num_thread : NUM_THREAD_VECTOR) {
            NUM_THREAD = num_thread;
            BOOST_CHECK_NO_THROW(process_kmer_long(fastq_path.string().c_str(), repeat_check_table, rot_table,
                                                   extract_k_mer, extract_k_mer_128,
                                                   thread_data_list, is_gz));
        }
    }
    delete[] thread_data_list;
}

BOOST_AUTO_TEST_CASE(main_test_long_64) {
    MAX_MER = 64;
    MIN_MER = 5;
    TABLE_MAX_MER = 1;
    SLICE_LENGTH = 150;
    BASELINE = 0.8;

    uint8_t **repeat_check_table = nullptr;
    uint32_t **rot_table = nullptr;

    uint64_t *extract_k_mer = nullptr;
    uint128_t *extract_k_mer_128 = nullptr;
    if (MAX_MER <= ABS_UINT64_MAX_MER) {
        extract_k_mer = set_extract_k_mer();
    }
    else {
        extract_k_mer_128 = set_extract_k_mer_128();
    }

    std::vector<std::pair<fs::path, bool>> FILE_LOC_VECTOR;
    if (boost::unit_test::framework::master_test_suite().argc == 5) {
        FILE_LOC_VECTOR = {{fs::path(boost::unit_test::framework::master_test_suite().argv[3]), true},
                           {fs::path(boost::unit_test::framework::master_test_suite().argv[4]), false}};
    }
    else {
        FILE_LOC_VECTOR = {{fs::current_path().parent_path() / fs::path("test") / fs::path("test_long.fastq.gz"), true},
                           {fs::current_path().parent_path() / fs::path("test") / fs::path("test_long.fastq"), false}};
    }

    std::vector<int> NUM_THREAD_VECTOR = {1, 2};

    ThreadData* thread_data_list = new ThreadData[*std::max_element(NUM_THREAD_VECTOR.begin(), NUM_THREAD_VECTOR.end())];
    for (auto& [fastq_path, is_gz] : FILE_LOC_VECTOR) {
        for(auto& num_thread : NUM_THREAD_VECTOR) {
            NUM_THREAD = num_thread;
            BOOST_CHECK_NO_THROW(process_kmer_long(fastq_path.string().c_str(), repeat_check_table, rot_table,
                                                   extract_k_mer, extract_k_mer_128,
                                                   thread_data_list, is_gz));
        }
    }
    delete[] thread_data_list;
}