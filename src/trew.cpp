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

double LOW_BASELINE;
double HIGH_BASELINE;

int main(int argc, char** argv) {
    argparse::ArgumentParser program("trew", "0.5.0");

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

    long_command.add_argument("-l", "--low_baseline")
            .help("Low baseline for k-mer telomere counting")
            .default_value(0.5)
            .scan<'f', double>()
            .metavar("LOW_BASELINE");

    long_command.add_argument("-h", "--high_baseline")
            .help("High baseline for k-mer telomere counting")
            .default_value(0.8)
            .scan<'f', double>()
            .metavar("HIGH_BASELINE");

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

    short_command.add_argument("MIN_MER")
            .help("minimum length of sequence to find telomere [MIN_MER >= " + std::to_string(ABS_MIN_MER) + "]")
            .nargs(1)
            .scan<'d', int>();

    short_command.add_argument("MAX_MER")
            .help("maximum length of sequence to find telomere [MAX_MER <= " + std::to_string(ABS_MAX_MER) + "]")
            .nargs(1)
            .scan<'d', int>();

    short_command.add_argument("SHORT_FASTQ_LOC")
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

    short_command.add_argument("-l", "--low_baseline")
            .help("Low baseline for k-mer telomere counting")
            .default_value(0.5)
            .scan<'f', double>()
            .metavar("LOW_BASELINE");

    short_command.add_argument("-h", "--high_baseline")
            .help("High baseline for k-mer telomere counting")
            .default_value(0.8)
            .scan<'f', double>()
            .metavar("HIGH_BASELINE");

    short_command.add_argument("-q", "--queue_size")
            .help("size of buffer queue in MiB [QUEUE_SIZE >= 4, unlimited : -1]")
            .default_value(-1)
            .scan<'d', int>()
            .metavar("QUEUE_SIZE");

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
            LOW_BASELINE = long_command.get<double>("--low_baseline");
            HIGH_BASELINE = long_command.get<double>("--high_baseline");

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

            if (TABLE_MAX_MER <= 0) {
                fprintf(stderr, "TABLE_MAX_MER must be positive.\n");
                throw std::exception();
            }

            if (NUM_THREAD <= 0) {
                fprintf(stderr, "number of threads must be positive.\n");
                throw std::exception();
            }

            if (!(0 < LOW_BASELINE && LOW_BASELINE <= 1) || !(0 < HIGH_BASELINE && HIGH_BASELINE <= 1)) {
                fprintf(stderr, "Baseline must be in range 0 to 1.\n");
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

        FinalFastqData* total_result_low = new FinalFastqData {};
        FinalFastqData* total_result_high = new FinalFastqData {};

        ThreadData* thread_data_list = new ThreadData[NUM_THREAD];
        FinalFastqVectorPair temp_fastq_vector_pair;

        std::vector<std::string> gz_extension_list = std::vector<std::string> {".gz", ".bgz"};
        for (const auto& fastq_path : fastq_path_list) {
            bool is_gz = false;
            std::string fastq_ext = fastq_path.extension().string();
            for (const auto& ext : gz_extension_list) {
                if (ext == fastq_ext) {
                    is_gz = true;
                    break;
                }
            }

            temp_fastq_vector_pair = process_kmer_long(std::filesystem::canonical(fastq_path).string().c_str(), repeat_check_table, rot_table,
                                                       extract_k_mer, extract_k_mer_128, extract_k_mer_ans,
                                                       thread_data_list, is_gz);

            for (auto& [k, v] : *temp_fastq_vector_pair.first) {
                if (total_result_high -> contains(k)) {
                    (*total_result_high)[k] = add_data((*total_result_high)[k], v);
                } else {
                    (*total_result_high)[k] = v;
                }
            }
            for (auto& [k, v] : *temp_fastq_vector_pair.second) {
                if (total_result_low -> contains(k)) {
                    (*total_result_low)[k] = add_data((*total_result_low)[k], v);
                } else {
                    (*total_result_low)[k] = v;
                }
            }

            delete temp_fastq_vector_pair.first;
            delete temp_fastq_vector_pair.second;
        }

        delete[] thread_data_list;
        final_process_output(total_result_high, total_result_low);
    }
    else if (program.is_subcommand_used("short")) {
        try {
            MIN_MER = short_command.get<int>("MIN_MER");
            MAX_MER = short_command.get<int>("MAX_MER");
            NUM_THREAD = short_command.get<int>("--thread");
            TABLE_MAX_MER = short_command.get<int>("--table_max_mer");
            QUEUE_SIZE = short_command.get<int>("--queue_size");
            LOW_BASELINE = short_command.get<double>("--low_baseline");
            HIGH_BASELINE = short_command.get<double>("--high_baseline");

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

            if (TABLE_MAX_MER <= 0) {
                fprintf(stderr, "TABLE_MAX_MER must be positive.\n");
                throw std::exception();
            }

            if (NUM_THREAD <= 0) {
                fprintf(stderr, "number of threads must be positive.\n");
                throw std::exception();
            }

            if (!(0 < LOW_BASELINE && LOW_BASELINE <= 1) || !(0 < HIGH_BASELINE && HIGH_BASELINE <= 1)) {
                fprintf(stderr, "Baseline must be in range 0 to 1.\n");
                throw std::exception();
            }
        }
        catch (...) {
            std::cerr << short_command;
            return 1;
        }

        std::vector<std::string> fastq_loc_list = short_command.get<std::vector<std::string>>("SHORT_FASTQ_LOC");
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

        FinalFastqData* total_result_low = new FinalFastqData {};
        FinalFastqData* total_result_high = new FinalFastqData {};

        ThreadData* thread_data_list = new ThreadData[NUM_THREAD];
        FinalFastqVectorPair temp_fastq_vector_pair;

        std::vector<std::string> gz_extension_list = std::vector<std::string> {".gz", ".bgz"};
        for (const auto& fastq_path : fastq_path_list) {
            bool is_gz = false;
            std::string fastq_ext = fastq_path.extension().string();
            for (const auto& ext : gz_extension_list) {
                if (ext == fastq_ext) {
                    is_gz = true;
                    break;
                }
            }

            temp_fastq_vector_pair = process_kmer(std::filesystem::canonical(fastq_path).string().c_str(), repeat_check_table, rot_table,
                                                  extract_k_mer, extract_k_mer_128, extract_k_mer_ans,
                                                  thread_data_list, is_gz);

            for (auto& [k, v] : *temp_fastq_vector_pair.first) {
                if (total_result_high -> contains(k)) {
                    (*total_result_high)[k] = add_data((*total_result_high)[k], v);
                } else {
                    (*total_result_high)[k] = v;
                }
            }
            for (auto& [k, v] : *temp_fastq_vector_pair.second) {
                if (total_result_low -> contains(k)) {
                    (*total_result_low)[k] = add_data((*total_result_low)[k], v);
                } else {
                    (*total_result_low)[k] = v;
                }
            }

            delete temp_fastq_vector_pair.first;
            delete temp_fastq_vector_pair.second;
        }

        delete[] thread_data_list;
        final_process_output(total_result_high, total_result_low);
    }
    else {
        std::cerr << program;
        return 1;
    }

    return 0;
}