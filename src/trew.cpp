#include "kmer.h"
#include "argparse/argparse.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>


int MAX_MER;
int MIN_MER;
int TABLE_MAX_MER;
int NUM_THREAD;
int SLICE_LENGTH;
int QUEUE_SIZE;

double LOW_BASELINE;
double HIGH_BASELINE;

bool INDEX = true;

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

    long_command.add_argument("LONG_FASTQ")
            .help("locations of FASTQ file")
            .nargs(argparse::nargs_pattern::at_least_one);

    long_command.add_argument("-t", "--thread")
            .help("number of threads")
            .default_value(2)
            .scan<'d', int>()
            .metavar("THREAD");

    long_command.add_argument("-m", "--table_max_mer")
            .help("maximum length of sequence to use table [TABLE_MAX_MER <= " + std::to_string(ABS_TABLE_MAX_MER) + "]")
            .default_value(12)
            .scan<'d', int>()
            .metavar("TABLE_MAX_MER");

    long_command.add_argument("-L", "--low_baseline")
            .help("low baseline for k-mer telomere counting")
            .default_value(0.5)
            .scan<'f', double>()
            .metavar("LOW_BASELINE");

    long_command.add_argument("-H", "--high_baseline")
            .help("high baseline for k-mer telomere counting")
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

    short_command.add_argument("SHORT_FASTQ")
            .help("Paths to FASTQ file (required for single-end mode)")
            .nargs(argparse::nargs_pattern::any);

    short_command.add_argument("-t", "--thread")
            .help("number of threads")
            .default_value(2)
            .scan<'d', int>()
            .metavar("THREAD");

    short_command.add_argument("--paired_end")
            .help("use paired-end sequencing data")
            .default_value(false)
            .implicit_value(true);

    short_command.add_argument("--fq1")
            .help("path to front FASTQ file (required for paired-end mode)")
            .metavar("FASTQ_FRONT")
            .nargs(argparse::nargs_pattern::at_least_one);

    short_command.add_argument("--fq2")
            .help("Path to reverse FASTQ file (required for paired-end mode)")
            .metavar("FASTQ_REVERSE")
            .nargs(argparse::nargs_pattern::at_least_one);

    short_command.add_argument("-m", "--table_max_mer")
            .help("maximum length of sequence to use table (reduce this option if memory usage is high) [TABLE_MAX_MER <= " + std::to_string(ABS_TABLE_MAX_MER) + "]")
            .default_value(12)
            .scan<'d', int>()
            .metavar("TABLE_MAX_MER");

    short_command.add_argument("-L", "--low_baseline")
            .help("low baseline for k-mer telomere counting")
            .default_value(0.5)
            .scan<'f', double>()
            .metavar("LOW_BASELINE");

    short_command.add_argument("-H", "--high_baseline")
            .help("high baseline for k-mer telomere counting")
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

    bool IS_PAIRED_END = false;
    std::vector<int> paired_end_vector {};
    std::vector<std::filesystem::path> fastq_path_list {};
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

            if (LOW_BASELINE > HIGH_BASELINE) {
                fprintf(stderr, "Low baseline must be smaller than high baseline.\n");
                throw std::exception();
            }

            if (NUM_THREAD < 2) {
                fprintf(stderr, "You must use at least two threads.\n");
                throw std::exception();
            }

            std::vector<std::string> fastq_loc_list = long_command.get<std::vector<std::string>>("LONG_FASTQ");
            for (const auto& fastq_loc : fastq_loc_list) {
                std::filesystem::path fastq_path {fastq_loc};
                if (! std::filesystem::is_regular_file(fastq_loc)) {
                    fprintf(stderr, "%s : file not found\n", fastq_loc.c_str());
                    return 1;
                }
                fastq_path_list.push_back(fastq_path);
            }
        }
        catch (...) {
            std::cerr << long_command;
            return 1;
        }
    } else if (program.is_subcommand_used("short")) {
        try {
            MIN_MER = short_command.get<int>("MIN_MER");
            MAX_MER = short_command.get<int>("MAX_MER");
            NUM_THREAD = short_command.get<int>("--thread");
            TABLE_MAX_MER = short_command.get<int>("--table_max_mer");
            QUEUE_SIZE = short_command.get<int>("--queue_size");
            LOW_BASELINE = short_command.get<double>("--low_baseline");
            HIGH_BASELINE = short_command.get<double>("--high_baseline");
            IS_PAIRED_END = short_command.get<bool>("--paired_end");

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

            if (LOW_BASELINE > HIGH_BASELINE) {
                fprintf(stderr, "Low baseline must be smaller than high baseline.\n");
                throw std::exception();
            }

            if (NUM_THREAD < 2) {
                fprintf(stderr, "You must use at least two threads.\n");
                throw std::exception();
            }

            if (IS_PAIRED_END) {
                // Ensure SHORT_FASTQ_LOC has zero arguments
                std::vector<std::string> short_fastq_loc = short_command.get<std::vector<std::string>>("SHORT_FASTQ");
                if (!short_fastq_loc.empty()) {
                    std::cerr << "SHORT_FASTQ must not be provided when --IS_PAIRED_END is used.\n";
                    throw std::exception();
                }

                // Ensure --fq1 and --fq2 are used
                if (!short_command.is_used("--fq1") || !short_command.is_used("--fq2")) {
                    std::cerr << "--fq1 and --fq2 are required in paired-end mode.\n";
                    throw std::exception();
                }

                // Retrieve fq1 and fq2 lists
                std::vector<std::string> fq1_loc_list = short_command.get<std::vector<std::string>>("--fq1");
                std::vector<std::string> fq2_loc_list = short_command.get<std::vector<std::string>>("--fq2");

                // Check if the number of fq1 and fq2 files match
                if (fq1_loc_list.size() != fq2_loc_list.size()) {
                    std::cerr << "--fq1 and --fq2 must have the same number of files.\n";
                    throw std::exception();
                }

                for (size_t i = 0; i < fq1_loc_list.size(); i++) {
                    const auto& fq1_loc = fq1_loc_list[i];
                    if (!std::filesystem::is_regular_file(fq1_loc)) {
                        std::cerr << fq1_loc << " : file not found\n";
                        throw std::exception();
                    }

                    const auto& fq2_loc = fq2_loc_list[i];
                    if (!std::filesystem::is_regular_file(fq2_loc)) {
                        std::cerr << fq2_loc << " : file not found\n";
                        throw std::exception();
                    }

                    fastq_path_list.emplace_back(fq1_loc);
                    paired_end_vector.push_back(FOR_READ);

                    fastq_path_list.emplace_back(fq2_loc);
                    paired_end_vector.push_back(REV_READ);
                }
            } else {
                // Single-end mode, ensure SHORT_FASTQ_LOC has at least one argument
                std::vector<std::string> fastq_loc_list = short_command.get<std::vector<std::string>>("SHORT_FASTQ");
                if (fastq_loc_list.empty()) {
                    std::cerr << "SHORT_FASTQ_LOC is required in single-end mode.\n";
                    throw std::exception();
                }

                // Ensure --fq1 and --fq2 are not used
                if (short_command.is_used("--fq1") || short_command.is_used("--fq2")) {
                    std::cerr << "--fq1 and --fq2 should not be used in single-end mode.\n";
                    throw std::exception();
                }

                // Check if all SHORT_FASTQ_LOC files exist
                for (const auto& fastq_loc : fastq_loc_list) {
                    if (!std::filesystem::is_regular_file(fastq_loc)) {
                        std::cerr << fastq_loc << " : file not found\n";
                        throw std::exception();
                    }
                    fastq_path_list.emplace_back(fastq_loc);
                }
            }
        }
        catch (...) {
            std::cerr << short_command;
            return 1;
        }
    } else {
        std::cerr << program;
        return 1;
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
    std::vector<std::string> gz_extension_list = std::vector<std::string> {".gz", ".bgz"};

    bool IS_SHORT = program.is_subcommand_used("short");
    bool is_pair = IS_SHORT and IS_PAIRED_END;

    FinalFastqOutput fastq_output;
    for (size_t i = 0; i < fastq_path_list.size() / (is_pair ? 2 : 1); ++i) {
        std::vector<std::filesystem::path> fastq_tmp_path = {};
        if (is_pair) {
            fastq_tmp_path.emplace_back(fastq_path_list[2 * i]);
            fastq_tmp_path.emplace_back(fastq_path_list[2 * i + 1]);
        } else {
            fastq_tmp_path.emplace_back(fastq_path_list[i]);
        }

        std::vector<bool> is_gz_vec = {};
        for (auto& path : fastq_tmp_path) {
            std::string fastq_ext = path.extension().string();
            bool t = false;
            for (const auto& ext : gz_extension_list) {
                if (ext == fastq_ext) {
                    t = true;
                    break;
                }
            }
            is_gz_vec.emplace_back(t);
        }



        if (IS_SHORT) {
            if (IS_PAIRED_END) {
                fastq_output = process_kmer_pair(std::filesystem::canonical(fastq_tmp_path[0]).string().c_str(),
                                        std::filesystem::canonical(fastq_tmp_path[1]).string().c_str(), repeat_check_table, rot_table,
                                            extract_k_mer, extract_k_mer_128, extract_k_mer_ans,
                                            thread_data_list, is_gz_vec[0], is_gz_vec[1]);
            } else {
                fastq_output = process_kmer(std::filesystem::canonical(fastq_tmp_path[0]).string().c_str(), repeat_check_table, rot_table,
                                            extract_k_mer, extract_k_mer_128, extract_k_mer_ans,
                                            thread_data_list, is_gz_vec[0]);
            }
        } else {
            fastq_output = process_kmer_long(std::filesystem::canonical(fastq_tmp_path[0]).string().c_str(), repeat_check_table, rot_table,
                                              extract_k_mer, extract_k_mer_128, extract_k_mer_ans,
                                              thread_data_list, is_gz_vec[1]);
        }

        for (auto& [k, v] : *fastq_output.high) {
            if (total_result_high -> contains(k)) {
                (*total_result_high)[k] = add_data((*total_result_high)[k], v);
            } else {
                (*total_result_high)[k] = v;
            }
        }
        for (auto& [k, v] : *fastq_output.low) {
            if (total_result_low -> contains(k)) {
                (*total_result_low)[k] = add_data((*total_result_low)[k], v);
            } else {
                (*total_result_low)[k] = v;
            }
        }

        delete fastq_output.high;
        delete fastq_output.low;
    }

    delete[] thread_data_list;


    final_process_output(total_result_high, total_result_low);
    return 0;
}