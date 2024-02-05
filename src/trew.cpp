#include "kmer.h"
#include "argparse/argparse.hpp"

#include <algorithm>
#include <filesystem>

int MAX_MER;
int MIN_MER;
int TABLE_MAX_MER;
int NUM_THREAD;
int SLICE_LENGTH;
int QUEUE_SIZE;
double BASELINE;

int main(int argc, char** argv) {
    argparse::ArgumentParser program("trew", "0.2.0");

    argparse::ArgumentParser long_command("long");
    argparse::ArgumentParser short_command("short");

    long_command.add_description("Estimate TRM from long-read sequencing data.");
    short_command.add_description("Estimate TRM from short-read sequencing data.");

    long_command.add_argument("MIN_MER")
            .help("minimum length of sequence to find telomere [MIN_MER >= " + std::to_string(ABS_MIN_MER) + "]")
            .nargs(1)
            .scan<'d', int>();

    long_command.add_argument("MAX_MER")
            .help("maximum length of sequence to find telomere [MAX_MER <= " + std::to_string(ABS_MAX_MER) + "]")
            .nargs(1)
            .scan<'d', int>();

    long_command.add_argument("LONG_FASTQ_LOC")
            .help("locations of FASTQ file")
            .nargs(argparse::nargs_pattern::at_least_one);

    long_command.add_argument("-t", "--thread")
            .help("number of threads")
            .default_value(1)
            .scan<'d', int>()
            .metavar("THREAD");

    long_command.add_argument("-m", "--table_max_mer")
            .help("maximum length of sequence to use table [TABLE_MAX_MER <= " + std::to_string(ABS_TABLE_MAX_MER) + "]")
            .default_value(12)
            .scan<'d', int>()
            .metavar("TABLE_MAX_MER");

    long_command.add_argument("-s", "--slice_length")
            .help("length of sequence to slice each side of read [SLICE_LENGTH >= 2 * MAX_MER]")
            .default_value(150)
            .scan<'d', int>()
            .metavar("SLICE_LENGTH");

    long_command.add_argument("-q", "--queue_size")
            .help("size of buffer queue in MiB [QUEUE_SIZE >= 4, unlimited : -1]")
            .default_value(-1)
            .scan<'d', int>()
            .metavar("QUEUE_SIZE");

    long_command.add_argument("-b", "--baseline")
            .help("BASELINE for repeat read [" + std::to_string(LOW_BASELINE) + " <= BASELINE <= 1]")
            .default_value(0.8)
            .scan<'g', double>()
            .metavar("BASELINE");

    short_command.add_argument("MIN_MER")
            .help("minimum length of sequence to find telomere [MIN_MER >= " + std::to_string(ABS_MIN_MER) + "]")
            .nargs(1)
            .scan<'d', int>();

    short_command.add_argument("MAX_MER")
            .help("maximum length of sequence to find telomere [MAX_MER <= " + std::to_string(ABS_MAX_MER) + "]")
            .nargs(1)
            .scan<'d', int>();

    short_command.add_argument("FASTQ_LOC")
            .help("locations of FASTQ file")
            .nargs(argparse::nargs_pattern::at_least_one);

    short_command.add_argument("-t", "--thread")
            .help("number of threads")
            .default_value(1)
            .scan<'d', int>()
            .metavar("THREAD");

    short_command.add_argument("-m", "--table_max_mer")
            .help("maximum length of sequence to use table (reduce this option if memory usage is high) [TABLE_MAX_MER <= " + std::to_string(ABS_TABLE_MAX_MER) + "]")
            .default_value(12)
            .scan<'d', int>()
            .metavar("TABLE_MAX_MER");

    short_command.add_argument("-q", "--queue_size")
            .help("size of buffer queue in MiB [QUEUE_SIZE >= 4, unlimited : -1]")
            .default_value(-1)
            .scan<'d', int>()
            .metavar("QUEUE_SIZE");

    short_command.add_argument("-b", "--baseline")
            .help("BASELINE for repeat read [" + std::to_string(LOW_BASELINE) + " <= BASELINE <= 1]")
            .default_value(0.8)
            .scan<'g', double>()
            .metavar("BASELINE");

    program.add_subparser(long_command);
    program.add_subparser(short_command);

    try {
        program.parse_args(argc, argv);
    }
    catch (...) {
        if (program.is_subcommand_used("long")) {
            std::cerr << long_command;
        }
        else if (program.is_subcommand_used("short")) {
            std::cerr << short_command;
        }
        else {
            std::cerr << program;
        }

        return 1;
    }

    if (program.is_subcommand_used("long")) {
        try {
            MIN_MER = long_command.get<int>("MIN_MER");
            MAX_MER = long_command.get<int>("MAX_MER");
            NUM_THREAD = long_command.get<int>("--thread");
            TABLE_MAX_MER = long_command.get<int>("--table_max_mer");
            SLICE_LENGTH = long_command.get<int>("--slice_length");
            QUEUE_SIZE = long_command.get<int>("--queue_size");
            BASELINE = long_command.get<double>("--baseline");

            // argument check
            if (MIN_MER > MAX_MER) {
                fprintf(stderr, "MIN_MER must not be greater than MAX_MER.\n");
                throw std::exception();
            }

            if (MIN_MER < ABS_MIN_MER) {
                fprintf(stderr, "MIN_MER must be greater than or equal to %d.\n", ABS_MIN_MER);
                throw std::exception();
            }

            if (MAX_MER > ABS_MAX_MER) {
                fprintf(stderr, "MAX_MER must be less than or equal to %d.\n", ABS_MAX_MER);
                throw std::exception();
            }

            if (TABLE_MAX_MER > ABS_TABLE_MAX_MER) {
                fprintf(stderr, "TABLE_MAX_MER must be less than or equal to %d.\n", ABS_TABLE_MAX_MER);
                throw std::exception();
            }

            if (SLICE_LENGTH < 2 * MAX_MER) {
                fprintf(stderr, "SLICE_LENGTH must be greater than or equal to twice of MAX_MER.\n");
                throw std::exception();
            }

            if (QUEUE_SIZE != -1 && QUEUE_SIZE < 4) {
                fprintf(stderr, "QUEUE_SIZE must be -1 (unlimited) or greater than or equal to 4.\n");
                throw std::exception();
            }

            if (LOW_BASELINE > BASELINE || BASELINE > 1) {
                fprintf(stderr, "BASELINE must be at least %.1f and no more than 1.\n", LOW_BASELINE);
                throw std::exception();
            }

            if (TABLE_MAX_MER <= 0) {
                fprintf(stderr, "TABLE_MAX_MER must be positive.\n");
                throw std::exception();
            }

            if (NUM_THREAD <= 0) {
                fprintf(stderr, "number of threads must be positive.\n");
                throw std::exception();
            }
        }
        catch (...) {
            std::cerr << long_command;
            return 1;
        }

        std::vector<std::string> fastq_loc_list = long_command.get<std::vector<std::string>>("LONG_FASTQ_LOC");
        std::vector<std::filesystem::path> fastq_path_list {};

        for (const auto& fastq_loc : fastq_loc_list) {
            std::filesystem::path fastq_path {fastq_loc};
            if (! std::filesystem::is_regular_file(fastq_loc)) {
                fprintf(stderr, "%s : file not found\n", fastq_loc.c_str());
                return 1;
            }
            fastq_path_list.push_back(fastq_path);
        }

        uint8_t **repeat_check_table = nullptr;
        uint32_t **rot_table = nullptr;
        if (MIN_MER <= TABLE_MAX_MER) {
            repeat_check_table = set_repeat_check_table();
            rot_table = set_rotation_table(repeat_check_table);
        }

        uint64_t *extract_k_mer = nullptr;
        uint128_t *extract_k_mer_128 = nullptr;
        if (MAX_MER <= ABS_UINT64_MAX_MER) {
            extract_k_mer = set_extract_k_mer();
        }
        else {
            extract_k_mer_128 = set_extract_k_mer_128();
        }

        uint128_t *extract_k_mer_ans = nullptr;
        if (MIN_MER > ABS_MIN_MER) {
            extract_k_mer_ans = set_extract_k_mer_ans();
        }

        ThreadData* thread_data_list = new ThreadData[NUM_THREAD];
        FinalFastqData* total_result = new FinalFastqData {};
        FinalFastqVector* temp_fastq_vector;

        std::vector<std::string> gz_extension_list = std::vector<std::string> {".gz", ".bgz"};
        for (const auto& fastq_path : fastq_path_list) {
            fprintf(stdout, ">%s\n", std::filesystem::canonical(fastq_path).string().c_str());

            bool is_gz = false;
            std::string fastq_ext = fastq_path.extension().string();
            for (const auto& ext : gz_extension_list) {
                if (ext == fastq_ext) {
                    is_gz = true;
                    break;
                }
            }

            temp_fastq_vector = process_kmer_long(fastq_path.string().c_str(), repeat_check_table, rot_table,
                                                  extract_k_mer, extract_k_mer_128, extract_k_mer_ans,
                                                  thread_data_list, is_gz);
            
            for (auto& [k, v] : *temp_fastq_vector) {
                if (total_result -> contains(k)) {
                    (*total_result)[k] = add_data((*total_result)[k], v);
                } else {
                    (*total_result)[k] = v;
                }
            }
            delete temp_fastq_vector;
            
        }

        int64_t tmp_tot_cnt;
        int64_t max_tot_cnt = -1;
        for (auto& [k, v] : *total_result) {
            tmp_tot_cnt = v.forward + v.backward + v.both;
            if (tmp_tot_cnt > max_tot_cnt) {
                max_tot_cnt = tmp_tot_cnt;
            }
        }

        fprintf(stdout, ">Putative_TRM\n");
        if (max_tot_cnt >= ANS_COUNT) {
            int64_t _tcnt;
            FinalFastqVector total_result_vector {};
            for (auto& [k, v] : *total_result) {
                if (v.forward + v.backward + v.both > PRINT_COUNT) {
                    if (v.backward > v.forward) {
                        _tcnt = v.forward;
                        v.forward = v.backward;
                        v.backward = _tcnt;
                    }

                    total_result_vector.emplace_back(k ,v);
                }
            }
            FinalFastqData ratio_result {};
            ResultMap score_result_map {};

            std::sort(total_result_vector.begin(), total_result_vector.end(), [](auto &a, auto &b) {
                return a.second.forward > b.second.forward;
            });

            for (int i = 0; i < MIN(10, total_result_vector.size()); i++) {
                if (total_result_vector[i].second.forward == 0) {
                    break;
                }
                if (i < 3) {
                    score_result_map[total_result_vector[i].first] += 1;
                }
                if (total_result_vector[i].second.backward >= 0) {
                    ratio_result[total_result_vector[i].first] = total_result_vector[i].second;
                }
            }

            std::sort(total_result_vector.begin(), total_result_vector.end(), [](auto &a, auto &b) {
                return a.second.forward + a.second.backward + a.second.both > b.second.forward + b.second.backward + b.second.both;
            });

            for (int i = 0; i < MIN(10, total_result_vector.size()); i++) {
                if (i < 3) {
                    score_result_map[total_result_vector[i].first] += 1;
                }
                if (total_result_vector[i].second.forward > 0 && total_result_vector[i].second.backward >= 0) {
                    ratio_result[total_result_vector[i].first] = total_result_vector[i].second;
                }
            }

            FinalFastqVector ratio_result_vector(ratio_result.begin(), ratio_result.end());
            std::sort(ratio_result_vector.begin(), ratio_result_vector.end(), [](auto &a, auto &b) {
                return (double) a.second.backward / a.second.forward < (double) b.second.backward / b.second.forward;
            });

            for (int i = 0; i < MIN(3, ratio_result_vector.size()); i++) {
                score_result_map[ratio_result_vector[i].first] += 1;
            }

            std::vector<std::pair<KmerSeq, uint32_t>> score_result_vector(score_result_map.begin(), score_result_map.end());
            std::sort(score_result_vector.begin(), score_result_vector.end(), [](auto &a, auto &b) {
                return a.second > b.second;
            });

            char buffer[ABS_MAX_MER + 1];
            for (auto& [k, v] : score_result_vector) {
                int_to_four(buffer, k.second, k.first);
                fprintf(stdout, "%s,%" PRIu32"\n", buffer, v);
            }
        } else {
            fprintf(stdout, "NO_PUTATIVE_TRM,-1\n");
        }

        delete total_result;
    }
    else if (program.is_subcommand_used("short")) {
        try {
            MIN_MER = short_command.get<int>("MIN_MER");
            MAX_MER = short_command.get<int>("MAX_MER");
            NUM_THREAD = short_command.get<int>("--thread");
            TABLE_MAX_MER = short_command.get<int>("--table_max_mer");
            QUEUE_SIZE = short_command.get<int>("--queue_size");
            BASELINE = short_command.get<double>("--baseline");

            // argument check
            if (MIN_MER > MAX_MER) {
                fprintf(stderr, "MIN_MER must not be greater than MAX_MER.\n");
                throw std::exception();
            }

            if (MIN_MER < ABS_MIN_MER) {
                fprintf(stderr, "MIN_MER must be greater than or equal to %d.\n", ABS_MIN_MER);
                throw std::exception();
            }

            if (MAX_MER > ABS_MAX_MER) {
                fprintf(stderr, "MAX_MER must be less than or equal to %d.\n", ABS_MAX_MER);
                throw std::exception();
            }

            if (TABLE_MAX_MER > ABS_TABLE_MAX_MER) {
                fprintf(stderr, "TABLE_MAX_MER must be less than or equal to %d.\n", ABS_TABLE_MAX_MER);
                throw std::exception();
            }

            if (QUEUE_SIZE != -1 && QUEUE_SIZE < 4) {
                fprintf(stderr, "QUEUE_SIZE must be -1 (unlimited) or greater than or equal to 4.\n");
                throw std::exception();
            }

            if (LOW_BASELINE > BASELINE || BASELINE > 1) {
                fprintf(stderr, "BASELINE must be at least %.1f and no more than 1.\n", LOW_BASELINE);
                throw std::exception();
            }

            if (TABLE_MAX_MER <= 0) {
                fprintf(stderr, "TABLE_MAX_MER must be positive.\n");
                throw std::exception();
            }

            if (NUM_THREAD <= 0) {
                fprintf(stderr, "number of threads must be positive.\n");
                throw std::exception();
            }
        }
        catch (...) {
            std::cerr << short_command;
            return 1;
        }

        std::vector<std::string> fastq_loc_list = short_command.get<std::vector<std::string>>("FASTQ_LOC");
        std::vector<std::filesystem::path> fastq_path_list {};

        for (const auto& fastq_loc : fastq_loc_list) {
            std::filesystem::path fastq_path {fastq_loc};
            if (! std::filesystem::is_regular_file(fastq_loc)) {
                fprintf(stderr, "%s : file not found\n", fastq_loc.c_str());
                return 1;
            }
            fastq_path_list.push_back(fastq_path);
        }

        uint8_t **repeat_check_table = nullptr;
        uint32_t **rot_table = nullptr;
        if (MIN_MER <= TABLE_MAX_MER) {
            repeat_check_table = set_repeat_check_table();
            rot_table = set_rotation_table(repeat_check_table);
        }

        uint64_t *extract_k_mer = nullptr;
        uint128_t *extract_k_mer_128 = nullptr;
        if (MAX_MER <= ABS_UINT64_MAX_MER) {
            extract_k_mer = set_extract_k_mer();
        }
        else {
            extract_k_mer_128 = set_extract_k_mer_128();
        }

        uint128_t *extract_k_mer_ans = nullptr;
        if (MIN_MER > ABS_MIN_MER) {
            extract_k_mer_ans = set_extract_k_mer_ans();
        }

        ThreadData* thread_data_list = new ThreadData[NUM_THREAD];
        FinalFastqData* total_result = new FinalFastqData {};
        FinalFastqVector* temp_fastq_vector;

        std::vector<std::string> gz_extension_list = std::vector<std::string> {".gz", ".bgz"};
        for (const auto& fastq_path : fastq_path_list) {
            fprintf(stdout, ">%s\n", std::filesystem::canonical(fastq_path).string().c_str());

            bool is_gz = false;
            std::string fastq_ext = fastq_path.extension().string();
            for (const auto& ext : gz_extension_list) {
                if (ext == fastq_ext) {
                    is_gz = true;
                    break;
                }
            }

            temp_fastq_vector = process_kmer(fastq_path.string().c_str(), repeat_check_table, rot_table,
                                             extract_k_mer, extract_k_mer_128, extract_k_mer_ans,
                                             thread_data_list, is_gz);

            for (auto& [k, v] : *temp_fastq_vector) {
                if (total_result -> contains(k)) {
                    (*total_result)[k] = add_data((*total_result)[k], v);
                } else {
                    (*total_result)[k] = v;
                }
            }
            delete temp_fastq_vector;
        }
        
        int64_t tmp_tot_cnt;
        int64_t max_tot_cnt = -1;
        for (auto& [k, v] : *total_result) {
            tmp_tot_cnt = v.forward + v.backward + v.both;
            if (tmp_tot_cnt > max_tot_cnt) {
                max_tot_cnt = tmp_tot_cnt;
            }
        }

        fprintf(stdout, ">Putative_TRM\n");
        if (max_tot_cnt >= ANS_COUNT) {
            int64_t _tcnt;
            for (auto& [_, v] : *total_result) {
                if (v.backward > v.forward) {
                    _tcnt = v.forward;
                    v.forward = v.backward;
                    v.backward = _tcnt;
                }
            }

            FinalFastqVector total_result_vector(total_result -> begin(), total_result -> end());
            FinalFastqData ratio_result {};
            ResultMap score_result_map {};

            std::sort(total_result_vector.begin(), total_result_vector.end(), [](auto &a, auto &b) {
                return a.second.forward > b.second.forward;
            });

            for (int i = 0; i < MIN(10, total_result_vector.size()); i++) {
                if (total_result_vector[i].second.forward == 0) {
                    break;
                }
                if (i < 3) {
                    score_result_map[total_result_vector[i].first] += 1;
                }
                if (total_result_vector[i].second.backward >= 0) {
                    ratio_result[total_result_vector[i].first] = total_result_vector[i].second;
                }
            }

            std::sort(total_result_vector.begin(), total_result_vector.end(), [](auto &a, auto &b) {
                return a.second.forward + a.second.backward + a.second.both > b.second.forward + b.second.backward + b.second.both;
            });

            for (int i = 0; i < MIN(10, total_result_vector.size()); i++) {
                if (i < 3) {
                    score_result_map[total_result_vector[i].first] += 1;
                }
                if (total_result_vector[i].second.forward > 0 && total_result_vector[i].second.backward >= 0) {
                    ratio_result[total_result_vector[i].first] = total_result_vector[i].second;
                }
            }

            FinalFastqVector ratio_result_vector(ratio_result.begin(), ratio_result.end());
            std::sort(ratio_result_vector.begin(), ratio_result_vector.end(), [](auto &a, auto &b) {
                return (double) a.second.backward / a.second.forward < (double) b.second.backward / b.second.forward;
            });

            for (int i = 0; i < MIN(3, ratio_result_vector.size()); i++) {
                score_result_map[ratio_result_vector[i].first] += 1;
            }

            std::vector<std::pair<KmerSeq, uint32_t>> score_result_vector(score_result_map.begin(), score_result_map.end());
            std::sort(score_result_vector.begin(), score_result_vector.end(), [](auto &a, auto &b) {
                return a.second > b.second;
            });

            char buffer[ABS_MAX_MER + 1];
            for (auto& [k, v] : score_result_vector) {
                int_to_four(buffer, k.second, k.first);
                fprintf(stdout, "%s,%" PRIu32"\n", buffer, v);
            }
        } else {
            fprintf(stdout, "NO_PUTATIVE_TRM,-1\n");
        }

        delete total_result;
        delete[] thread_data_list;
    }
    else {
        std::cerr << program;
        return 1;
    }

    return 0;
}