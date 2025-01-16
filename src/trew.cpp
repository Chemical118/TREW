#include "kmer.h"
#include "argparse/argparse.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>

#ifdef _WIN32
#include <io.h>
#endif

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

    std::vector<gz_index *> gz_index_vector = {};
    std::vector<FinalFastqOutput> fastq_file_data_vector = {};
    TRMDirVector* put_trm;

    bool IS_SHORT = program.is_subcommand_used("short");

    FinalFastqOutput fastq_output;
    for (size_t i = 0; i < fastq_path_list.size(); ++i) {
        const auto& fastq_path = fastq_path_list[i];

        bool is_gz = false;
        std::string fastq_ext = fastq_path.extension().string();
        for (const auto& ext : gz_extension_list) {
            if (ext == fastq_ext) {
                is_gz = true;
                break;
            }
        }

        gz_index* index = NULL;
        if (IS_SHORT) {
            if (IS_PAIRED_END) {
                fastq_output = process_kmer(std::filesystem::canonical(fastq_path).string().c_str(), repeat_check_table, rot_table,
                                            extract_k_mer, extract_k_mer_128, extract_k_mer_ans,
                                            thread_data_list, is_gz, &index, paired_end_vector[i]);
            } else {
                fastq_output = process_kmer(std::filesystem::canonical(fastq_path).string().c_str(), repeat_check_table, rot_table,
                                            extract_k_mer, extract_k_mer_128, extract_k_mer_ans,
                                            thread_data_list, is_gz, &index, SIN_READ);
            }
        } else {
            fastq_output = process_kmer_long(std::filesystem::canonical(fastq_path).string().c_str(), repeat_check_table, rot_table,
                                              extract_k_mer, extract_k_mer_128, extract_k_mer_ans,
                                              thread_data_list, is_gz, &index);
        }
        gz_index_vector.push_back(index);
        fastq_file_data_vector.push_back(fastq_output);

        if (INDEX and ((is_gz and index == NULL) or (not is_gz and index != NULL))) {
            fprintf(stderr, "ZLIB index error\n");
            exit(-1);
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

    if (IS_PAIRED_END and INDEX) {
        ResultMapData result_data = ResultMapData {{new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}, {nullptr, nullptr}};

        auto extract_k_mer_128_ = set_extract_k_mer_128();
        auto k_mer_data_128 = set_k_mer_data_128();
        auto k_mer_counter = set_k_mer_counter();
        auto k_mer_counter_list = set_k_mer_counter_list();
        int16_t* k_mer_total_cnt = (int16_t*) malloc(sizeof(int16_t) * (MAX_MER - MIN_MER + 2));

        CounterMap_128* k_mer_counter_map = nullptr;
        if (TABLE_MAX_MER < MAX_MER) {
            k_mer_counter_map = new CounterMap_128[MAX_MER - TABLE_MAX_MER];
        }

        for (size_t i = 0; i < fastq_path_list.size() / 2; ++i) {
            FILE* fp1;
            if (gz_index_vector[2 * i] == NULL) {
                fp1 = fopen(fastq_path_list[2 * i].string().c_str(), "r");
            } else {
                fp1 = fopen(fastq_path_list[2 * i].string().c_str(), "rb");
            }

            FILE* fp2;
            if (gz_index_vector[2 * i + 1] == NULL) {
                fp2 = fopen(fastq_path_list[2 * i + 1].string().c_str(), "r");
            } else {
                fp2 = fopen(fastq_path_list[2 * i + 1].string().c_str(), "rb");
            }

            paired_end_bonus_result(result_data, rot_table, repeat_check_table, extract_k_mer_128_,
                                    k_mer_counter, k_mer_counter_map,
                                    k_mer_data_128, k_mer_counter_list, k_mer_total_cnt,
                                    fastq_file_data_vector[2 * i].k_mer_loc_vector, fastq_file_data_vector[2 * i + 1].k_mer_loc_vector,
                                    fp1, fp2,
                                    gz_index_vector[2 * i], gz_index_vector[2 * i + 1]);
        }

        for (auto& [seq, cnt] : *(result_data.backward.first)) {
            (*(result_data.forward.first))[KmerSeq {seq.first, get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first)}] += cnt;
        }
        for (auto& [seq, cnt] : *(result_data.backward.second)) {
            (*(result_data.forward.second))[KmerSeq {seq.first, get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first)}] += cnt;
        }

        uint128_t _t;
        uint128_t kseq;

        for (auto& [seq, cnt] : *(result_data.forward.second)) {
            _t = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
            kseq = MIN(_t, seq.second);

            if (kseq == seq.second) {
                (*total_result_low)[KmerSeq {seq.first, kseq}].forward += cnt;
            }
            else {
                (*total_result_low)[KmerSeq {seq.first, kseq}].backward += cnt;
            }
        }
        for (auto& [seq, cnt] : *(result_data.forward.first)) {
            _t = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
            kseq = MIN(_t, seq.second);

            if (kseq == seq.second) {
                (*total_result_high)[KmerSeq {seq.first, kseq}].forward += cnt;
            }
            else {
                (*total_result_high)[KmerSeq {seq.first, kseq}].backward += cnt;
            }
        }
    }

    put_trm = final_process_output(total_result_high, total_result_low);

    if (INDEX) {
        std::vector<std::vector<char*>> temp_file_loc_vector = {};

        for (int i = 0; i < fastq_file_data_vector.size(); i++) {
            size_t data_size = fastq_file_data_vector[i].k_mer_loc_vector -> size();
            size_t chunk_size = data_size / NUM_THREAD;

            std::vector<char*> temp_file_list;

            for (size_t j = 0; j < NUM_THREAD; ++j) {
				std::string template_loc = (std::filesystem::temp_directory_path() / "tmp_trew_XXXXXX").string();
				char* template_loc_buf = new char[template_loc.size() + 1];
				std::strcpy(template_loc_buf, template_loc.c_str());
#ifdef _WIN32
			  	template_loc_buf = _mktemp(template_loc_buf);
#else
				int fd = mkstemp(template_loc_buf);
				close(fd);
#endif
				temp_file_list.emplace_back(template_loc_buf);
            }

            tbb::task_arena arena(NUM_THREAD);
            arena.execute([&]{
                tbb::parallel_for(0, NUM_THREAD, 1, [&](size_t thread_idx) {
                    size_t st = thread_idx * chunk_size;
                    size_t nd = (thread_idx == NUM_THREAD -1) ? data_size : st + chunk_size;
                    get_trm_read(fastq_path_list[i], put_trm, fastq_file_data_vector[i],
								 get_thread_safe_index(gz_index_vector[i]), st, nd, temp_file_list[thread_idx]);
                });
            });

            temp_file_loc_vector.push_back(temp_file_list);
        }

        auto trm_out_loc = fastq_path_list[0].string() + ".trm_read.fa";
        std::remove(trm_out_loc.c_str());

        for (auto &tmp_list : temp_file_loc_vector) {
            for (auto &tmp_file : tmp_list) {
                std::ofstream final_fasta_out(trm_out_loc, std::ios_base::binary | std::ios_base::app);
                final_fasta_out.seekp(0, std::ios_base::end);

                std::ifstream tmp_fa(tmp_file, std::ios_base::binary);
                final_fasta_out << tmp_fa.rdbuf();

                tmp_fa.close();
                std::remove(tmp_file);
                delete[] tmp_file;
            }
        }
    }

    return 0;
}