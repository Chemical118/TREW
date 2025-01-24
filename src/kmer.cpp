#include "kmer.h"

#include <cstdio>
#include <cstdlib>
#include <vector>

char trans_arr[4] = {'T', 'G', 'C', 'A'};
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

uint32_t reverse_complement_32(uint32_t x) {
    x = (x >> 16) | (x << 16);
    x = ((x >> 8) & 0x00ff00ff) | ((x & 0x00ff00ff) << 8);
    x = ((x >> 4) & 0x0f0f0f0f) | ((x & 0x0f0f0f0f) << 4);
    x = ((x >> 2) & 0x33333333) | ((x & 0x33333333) << 2);
    return ~x;
}

uint64_t reverse_complement_64(uint64_t x) {
    x = (x >> 32) | (x << 32);
    x = ((x >> 16) & 0x0000ffff0000ffffULL) | ((x & 0x0000ffff0000ffffULL) << 16);
    x = ((x >> 8)  & 0x00ff00ff00ff00ffULL) | ((x & 0x00ff00ff00ff00ffULL) <<  8);
    x = ((x >> 4)  & 0x0f0f0f0f0f0f0f0fULL) | ((x & 0x0f0f0f0f0f0f0f0fULL) <<  4);
    x = ((x >> 2)  & 0x3333333333333333ULL) | ((x & 0x3333333333333333ULL) <<  2);
    return ~x;
}

static const uint128_t rev_cmp[5] = {absl::MakeUint128(0x00000000ffffffff, 0x00000000ffffffff),
                                     absl::MakeUint128(0x0000ffff0000ffff, 0x0000ffff0000ffff),
                                     absl::MakeUint128(0x00ff00ff00ff00ff, 0x00ff00ff00ff00ff),
                                     absl::MakeUint128(0x0f0f0f0f0f0f0f0f, 0x0f0f0f0f0f0f0f0f),
                                     absl::MakeUint128(0x3333333333333333, 0x3333333333333333)};

uint128_t reverse_complement_128(uint128_t x) {
    x = (x >> 64) | (x << 64);
    x = ((x >> 32) & rev_cmp[0]) | ((x & rev_cmp[0]) << 32);
    x = ((x >> 16) & rev_cmp[1]) | ((x & rev_cmp[1]) << 16);
    x = ((x >> 8)  & rev_cmp[2]) | ((x & rev_cmp[2]) <<  8);
    x = ((x >> 4)  & rev_cmp[3]) | ((x & rev_cmp[3]) <<  4);
    x = ((x >> 2)  & rev_cmp[4]) | ((x & rev_cmp[4]) <<  2);
    return ~x;
}

KmerSeq rot_reverse_complement(KmerSeq seq) {
    return KmerSeq{seq.first, get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first)};
}

FinalData<int64_t> add_data(FinalData<int64_t> a, FinalData<int64_t> b) {
    return FinalData<int64_t> {a.forward + b.forward, a.backward + b.backward, a.both + b.both};
}

ResultMapData buffer_task(TBBQueue* task_queue, ThreadData* thread_data, uint32_t** rot_table, const uint64_t *extract_k_mer, const uint128_t *extract_k_mer_128, uint8_t** repeat_check_table) {
    auto [k_mer_counter, k_mer_data, k_mer_data_128, k_mer_counter_list] = thread_data -> init_check();

    int16_t* k_mer_total_cnt = (int16_t*) malloc(sizeof(int16_t) * (MAX_MER - MIN_MER + 2));
    if (k_mer_total_cnt == nullptr) {
        fprintf(stderr, "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    ResultMapData result = ResultMapData {{new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}};

    ResultMapPair temp_result_left = {new ResultMap {}, new ResultMap {}};
    ResultMapPair temp_result_right = {new ResultMap {}, new ResultMap {}};

    QueueData temp_task {};

    int n;
    bool high_half_check, low_half_check;
    KmerData left_temp_k_mer, right_temp_k_mer;

    if (MAX_MER <= ABS_UINT64_MAX_MER) {
        CounterMap* k_mer_counter_map = nullptr;
        if (TABLE_MAX_MER < MAX_MER) {
            k_mer_counter_map = new CounterMap[MAX_MER - TABLE_MAX_MER];
        }

        while (1) {
            task_queue -> pop(temp_task);
            if (temp_task.loc_vector == nullptr) {
                break;
            }
            for (auto& [st, nd]: *(temp_task.loc_vector)) {
                std::pair<uint64_t, uint64_t> lef_seq, rht_seq;
                n = nd - st + 1;

                if (2 * MIN_MER <= n) {
                    left_temp_k_mer = {0, 0};
                    right_temp_k_mer = {0, 0};

                    if (4 * MIN_MER <= n) {
                        left_temp_k_mer = k_mer_check(temp_task.buffer, st, st + (n / 2) - 1, rot_table, extract_k_mer, k_mer_counter, k_mer_counter_map,
                                                      k_mer_data, k_mer_counter_list, repeat_check_table, temp_result_left, k_mer_total_cnt, MIN_MER, MIN(n / 4, MAX_MER), &lef_seq);

                        if (left_temp_k_mer.first > 0 || left_temp_k_mer.second > 0) {
                            right_temp_k_mer = k_mer_check(temp_task.buffer, nd - ((n + 1) / 2) + 1, nd, rot_table, extract_k_mer, k_mer_counter, k_mer_counter_map,
                                                           k_mer_data, k_mer_counter_list, repeat_check_table, ResultMapPair {left_temp_k_mer.first > 0 ? nullptr : temp_result_right.first, left_temp_k_mer.second > 0 ? nullptr : temp_result_right.second},
                                                           k_mer_total_cnt, MIN_MER, MIN(n / 4, MAX_MER), &rht_seq);

                            if ((left_temp_k_mer.first == right_temp_k_mer.first) && (left_temp_k_mer.first > 0)) {
                                k_mer_target(temp_task.buffer, st, nd, rot_table, extract_k_mer, k_mer_counter, k_mer_counter_map,
                                             k_mer_data, k_mer_counter_list, repeat_check_table, {result.both.first, nullptr}, k_mer_total_cnt, left_temp_k_mer.first);
                            } else {
                                for (auto& [seq, cnt] : *(temp_result_left.first)) {
                                    (*(result.forward.first))[seq] += cnt;
                                }

                                for (auto& [seq, cnt] : *(temp_result_right.first)) {
                                    (*(result.backward.first))[seq] += cnt;
                                }
                            }

                            if ((left_temp_k_mer.second == right_temp_k_mer.second) && (left_temp_k_mer.second > 0)) {
                                k_mer_target(temp_task.buffer, st, nd, rot_table, extract_k_mer, k_mer_counter, k_mer_counter_map,
                                             k_mer_data, k_mer_counter_list, repeat_check_table, {nullptr, result.both.second}, k_mer_total_cnt, left_temp_k_mer.second);
                            } else {
                                for (auto& [seq, cnt] : *(temp_result_left.second)) {
                                    (*(result.forward.second))[seq] += cnt;
                                }

                                for (auto& [seq, cnt] : *(temp_result_right.second)) {
                                    (*(result.backward.second))[seq] += cnt;
                                }
                            }

                            temp_result_right.first -> clear();
                            temp_result_right.second -> clear();
                        } else {
                            right_temp_k_mer = k_mer_check(temp_task.buffer, nd - ((n + 1) / 2) + 1, nd, rot_table, extract_k_mer, k_mer_counter, k_mer_counter_map,
                                                           k_mer_data, k_mer_counter_list, repeat_check_table, result.backward, k_mer_total_cnt, MIN_MER, MIN(n / 4, MAX_MER), &rht_seq);
                        }

                        temp_result_left.first -> clear();
                        temp_result_left.second -> clear();
                    }

                    high_half_check = left_temp_k_mer.first == 0 && right_temp_k_mer.first == 0;
                    low_half_check = left_temp_k_mer.second == 0 && right_temp_k_mer.second == 0;

                    if (4 * MAX_MER > n && (high_half_check || low_half_check)) {
                        k_mer_check(temp_task.buffer, st, nd, rot_table, extract_k_mer, k_mer_counter, k_mer_counter_map,
                                    k_mer_data, k_mer_counter_list, repeat_check_table, {high_half_check ? result.both.first : nullptr, low_half_check ? result.both.second : nullptr}, k_mer_total_cnt, MAX(n / 4 + 1, MIN_MER), MIN(n / 2, MAX_MER));
                    }
                }
            }

            delete temp_task.loc_vector;
            free(temp_task.buffer);
        }

        delete[] k_mer_counter_map;
    } else {
        CounterMap_128* k_mer_counter_map = new CounterMap_128[MAX_MER - TABLE_MAX_MER];

        while (1) {
            task_queue -> pop(temp_task);
            if (temp_task.loc_vector == nullptr) {
                break;
            }
            for (auto& [st, nd]: *(temp_task.loc_vector)) {
                std::pair<uint128_t, uint128_t> lef_seq, rht_seq;
                n = nd - st + 1;

                if (2 * MIN_MER <= n) {
                    left_temp_k_mer = {0, 0};
                    right_temp_k_mer = {0, 0};

                    if (4 * MIN_MER <= n) {
                        left_temp_k_mer = k_mer_check_128(temp_task.buffer, st, st + (n / 2) - 1, rot_table, extract_k_mer_128, k_mer_counter, k_mer_counter_map,
                                                          k_mer_data_128, k_mer_counter_list, repeat_check_table, temp_result_left, k_mer_total_cnt, MIN_MER, MIN(n / 4, MAX_MER), &lef_seq);

                        if (left_temp_k_mer.first > 0 || left_temp_k_mer.second > 0) {
                            right_temp_k_mer = k_mer_check_128(temp_task.buffer, nd - ((n + 1) / 2) + 1, nd, rot_table, extract_k_mer_128, k_mer_counter, k_mer_counter_map,
                                                               k_mer_data_128, k_mer_counter_list, repeat_check_table, ResultMapPair {left_temp_k_mer.first > 0 ? nullptr : temp_result_right.first, left_temp_k_mer.second > 0 ? nullptr : temp_result_right.second},
                                                               k_mer_total_cnt, MIN_MER, MIN(n / 4, MAX_MER), &rht_seq);

                            if ((left_temp_k_mer.first == right_temp_k_mer.first) && (left_temp_k_mer.first > 0)) {
                                k_mer_target_128(temp_task.buffer, st, nd, rot_table, extract_k_mer_128, k_mer_counter, k_mer_counter_map,
                                                 k_mer_data_128, k_mer_counter_list, repeat_check_table, {result.both.first, nullptr}, k_mer_total_cnt, left_temp_k_mer.first);
                            } else {
                                for (auto& [seq, cnt] : *(temp_result_left.first)) {
                                    (*(result.forward.first))[seq] += cnt;
                                }

                                for (auto& [seq, cnt] : *(temp_result_right.first)) {
                                    (*(result.backward.first))[seq] += cnt;
                                }
                            }

                            if ((left_temp_k_mer.second == right_temp_k_mer.second) && (left_temp_k_mer.second > 0)) {
                                k_mer_target_128(temp_task.buffer, st, nd, rot_table, extract_k_mer_128, k_mer_counter, k_mer_counter_map,
                                                 k_mer_data_128, k_mer_counter_list, repeat_check_table, {nullptr, result.both.second}, k_mer_total_cnt, left_temp_k_mer.second);
                            } else {
                                for (auto& [seq, cnt] : *(temp_result_left.second)) {
                                    (*(result.forward.second))[seq] += cnt;
                                }

                                for (auto& [seq, cnt] : *(temp_result_right.second)) {
                                    (*(result.backward.second))[seq] += cnt;
                                }
                            }

                            temp_result_right.first -> clear();
                            temp_result_right.second -> clear();
                        } else {
                            right_temp_k_mer = k_mer_check_128(temp_task.buffer, nd - ((n + 1) / 2) + 1, nd, rot_table, extract_k_mer_128, k_mer_counter, k_mer_counter_map,
                                                               k_mer_data_128, k_mer_counter_list, repeat_check_table, result.backward, k_mer_total_cnt, MIN_MER, MIN(n / 4, MAX_MER), &rht_seq);
                        }

                        temp_result_left.first -> clear();
                        temp_result_left.second -> clear();
                    }

                    high_half_check = left_temp_k_mer.first == 0 && right_temp_k_mer.first == 0;
                    low_half_check = left_temp_k_mer.second == 0 && right_temp_k_mer.second == 0;

                    if (4 * MAX_MER > n && (high_half_check || low_half_check)) {
                        k_mer_check_128(temp_task.buffer, st, nd, rot_table, extract_k_mer_128, k_mer_counter, k_mer_counter_map,
                                        k_mer_data_128, k_mer_counter_list, repeat_check_table, {high_half_check ? result.both.first : nullptr, low_half_check ? result.both.second : nullptr}, k_mer_total_cnt, MAX(n / 4 + 1, MIN_MER), MIN(n / 2, MAX_MER));
                    }
                }
            }

            delete temp_task.loc_vector;
            free(temp_task.buffer);
        }

        delete[] k_mer_counter_map;
    }

    delete temp_result_left.first;
    delete temp_result_left.second;
    delete temp_result_right.first;
    delete temp_result_right.second;

    free(k_mer_total_cnt);
    return result;
}

ResultMapData buffer_task_pair(TBBPairQueue* task_queue, ThreadData* thread_data, uint32_t** rot_table, const uint64_t *extract_k_mer, const uint128_t *extract_k_mer_128, uint8_t** repeat_check_table) {
    auto [k_mer_counter, k_mer_data, k_mer_data_128, k_mer_counter_list] = thread_data -> init_check();

    int16_t* k_mer_total_cnt = (int16_t*) malloc(sizeof(int16_t) * (MAX_MER - MIN_MER + 2));
    if (k_mer_total_cnt == nullptr) {
        fprintf(stderr, "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    ResultMapData result = ResultMapData {{new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}};

    ResultMapPair temp_result_left = {new ResultMap {}, new ResultMap {}};
    ResultMapPair temp_result_right = {new ResultMap {}, new ResultMap {}};

    PairQueueData temp_task {};

    int snum;

    int n, n1, n2;
    int st1, nd1, st2, nd2;

    KmerData left_temp_k_mer, right_temp_k_mer;
    int ti, tj;

    KmerData k_mer, temp_k_mer, lef_k_mer;
    KmerData si, sj;
    std::pair<bool, bool> repeat_end;

    uint128_t seq_rev;

    if (MAX_MER <= ABS_UINT64_MAX_MER) {
        CounterMap* k_mer_counter_map = nullptr;
        std::pair<uint64_t, uint64_t> k_mer_seq, temp_k_mer_seq;
        std::pair<uint64_t, uint64_t> lef_seq, rht_seq;

        if (TABLE_MAX_MER < MAX_MER) {
            k_mer_counter_map = new CounterMap[MAX_MER - TABLE_MAX_MER];
        }

        auto get_dir_seq = [](int i, int k, uint64_t seq, bool is_for) {
            if (i <= 2 == is_for) {
                return seq;
            } else {
                return get_rot_seq(reverse_complement_64(seq) >> (2 * (32 - k)), k);
            }
        };

        while (1) {
            task_queue -> pop(temp_task);
            if (temp_task.loc_vector1 == nullptr) {
                break;
            }

            size_t min_size = std::min(temp_task.loc_vector1 -> size(), temp_task.loc_vector2 -> size());
            for (size_t i = 0; i < min_size; i++) {
                st1 = (*temp_task.loc_vector1)[i].first;
                nd1 = (*temp_task.loc_vector1)[i].second;
                n1 = nd1 - st1 + 1;

                st2 = (*temp_task.loc_vector2)[i].first;
                nd2 = (*temp_task.loc_vector2)[i].second;
                n2 = nd2 - st2 + 1;

                n = MIN(n1, n2);

                if (2 * MIN_MER <= n) {
                    lef_k_mer = {0, 0};
                    k_mer = {0, 0};

                    if (4 * MIN_MER <= n) {
                        std::vector<const char*> buffer_vector= {temp_task.buffer1, temp_task.buffer1, temp_task.buffer2, temp_task.buffer2};
                        std::vector<std::pair<int, int>> index_vector= {{st1, st1 + (n1 / 2) - 1}, {nd1 - ((n1 + 1) / 2) + 1, nd1},
                                                                        {nd2 - ((n2 + 1) / 2) + 1, nd2}, {st2, st2 + (n2 / 2) - 1}};

                        snum = 4;
                        si = {1, 1};
                        k_mer = {0, 0};

                        for (ti = 1; ti <= snum && (!repeat_end.first || !repeat_end.second) ; ti++) {
                            temp_k_mer = k_mer_check(buffer_vector[ti - 1], index_vector[ti - 1].first, index_vector[ti - 1].second, rot_table, extract_k_mer, k_mer_counter,
                                                         k_mer_counter_map, k_mer_data, k_mer_counter_list, repeat_check_table,
                                                         {repeat_end.first ? nullptr : (ti <= 2 ? temp_result_left.first : temp_result_right.first),
                                                                     repeat_end.second ? nullptr : (ti <= 2 ? temp_result_left.second : temp_result_right.second)},
                                                                     k_mer_total_cnt, MIN_MER, MIN(n / 4, MAX_MER), &temp_k_mer_seq);

                            if (!repeat_end.first && temp_k_mer.first > 0 && ((k_mer.first == temp_k_mer.first and k_mer_seq.first == get_dir_seq(ti, temp_k_mer.first, temp_k_mer_seq.first, true)) || ti == 1)) {
                                si.first += 1;
                                k_mer.first = temp_k_mer.first;
                                if (ti == 1) {
                                    k_mer_seq.first = temp_k_mer_seq.first;
                                }
                                repeat_end.first = false;
                            } else {
                                repeat_end.first = true;
                            }
                            if (!repeat_end.second && temp_k_mer.second > 0 && ((k_mer.second == temp_k_mer.second and k_mer_seq.second == get_dir_seq(ti, temp_k_mer.second, temp_k_mer_seq.second, true)) || ti == 1)) {
                                si.second += 1;
                                k_mer.second = temp_k_mer.second;
                                if (ti == 1) {
                                    k_mer_seq.second = temp_k_mer_seq.second;
                                }
                                repeat_end.second = false;
                            } else {
                                repeat_end.second = true;
                            }
                        }
                        lef_k_mer = k_mer;

                        if (si.first == snum + 1) {
                            for (auto& [seq, cnt] : *(temp_result_left.first)) {
                                seq_rev = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
                                (*(result.both.first))[KmerSeq {seq.first, MIN(seq.second, seq_rev)}] += cnt;
                            }

                            for (auto& [seq, cnt] : *(temp_result_right.first)) {
                                seq_rev = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
                                (*(result.both.first))[KmerSeq {seq.first, MIN(seq.second, seq_rev)}] += cnt;
                            }
                        }
                        if (si.second == snum + 1) {
                            for (auto& [seq, cnt] : *(temp_result_left.second)) {
                                seq_rev = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
                                (*(result.both.second))[KmerSeq {seq.first, MIN(seq.second, seq_rev)}] += cnt;
                            }

                            for (auto& [seq, cnt] : *(temp_result_right.second)) {
                                seq_rev = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
                                (*(result.both.second))[KmerSeq {seq.first, MIN(seq.second, seq_rev)}] += cnt;
                            }
                        }

                        if (si.first <= snum or si.second <= snum) {
                            sj = {snum, snum};
                            k_mer = {0, 0};
                            repeat_end = {false, false};
                            for (tj = snum; !repeat_end.first || !repeat_end.second; tj--) {
                                temp_k_mer = k_mer_check(buffer_vector[tj - 1], index_vector[tj - 1].first, index_vector[tj - 1].second, rot_table, extract_k_mer, k_mer_counter,
                                                             k_mer_counter_map, k_mer_data, k_mer_counter_list, repeat_check_table,
                                                {repeat_end.first ? nullptr : (tj <= 2 ? temp_result_right.first : temp_result_left.first),
                                                            repeat_end.second ? nullptr : (tj <= 2 ? temp_result_right.second : temp_result_left.second)},
                                                            k_mer_total_cnt, MIN_MER, MIN(n / 4, MAX_MER), &temp_k_mer_seq);

                                if (sj.first >= si.first && !repeat_end.first && temp_k_mer.first > 0 && ((k_mer.first == temp_k_mer.first and
                                    k_mer_seq.first == get_dir_seq(tj, temp_k_mer.first, temp_k_mer_seq.first, false)) || tj == snum)) {
                                    sj.first -= 1;
                                    k_mer.first = temp_k_mer.first;
                                    if (tj == snum) {
                                        k_mer_seq.first = temp_k_mer_seq.first;
                                    }
                                    repeat_end.first = false;
                                } else {
                                    repeat_end.first = true;
                                }
                                if (sj.second >= si.second && !repeat_end.second && temp_k_mer.second > 0 && ((k_mer.second == temp_k_mer.second and
                                    k_mer_seq.second == get_dir_seq(tj, temp_k_mer.second, temp_k_mer_seq.second, false)) || tj == snum)) {
                                    sj.second -= 1;
                                    k_mer.second = temp_k_mer.second;
                                    if (tj == snum) {
                                        k_mer_seq.second = temp_k_mer_seq.second;
                                    }
                                    repeat_end.second = false;
                                } else {
                                    repeat_end.second = true;
                                }
                            }
                        }

                        if (si.first <= snum) {
                            for (auto& [seq, cnt] : *(temp_result_left.first)) {
                                (*(result.forward.first))[seq] += cnt;
                            }

                            for (auto& [seq, cnt] : *(temp_result_right.first)) {
                                (*(result.backward.first))[seq] += cnt;
                            }
                        }
                        if (si.second <= snum) {
                            for (auto& [seq, cnt] : *(temp_result_left.second)) {
                                (*(result.forward.second))[seq] += cnt;
                            }

                            for (auto& [seq, cnt] : *(temp_result_right.second)) {
                                (*(result.backward.second))[seq] += cnt;
                            }
                        }

                        temp_result_left.first -> clear();
                        temp_result_left.second -> clear();

                        temp_result_right.first -> clear();
                        temp_result_right.second -> clear();
                    }

                    left_temp_k_mer = {0, 0};
                    right_temp_k_mer = {0, 0};

                    if (4 * MAX_MER > n) {
                        if (lef_k_mer.first == 0 or lef_k_mer.second == 0) {
                            left_temp_k_mer = k_mer_check(temp_task.buffer1, st1, nd1, rot_table, extract_k_mer, k_mer_counter, k_mer_counter_map,
                                                          k_mer_data, k_mer_counter_list, repeat_check_table,
                                                          {lef_k_mer.first == 0 ? temp_result_left.first : nullptr,
                                                          lef_k_mer.second == 0 ? temp_result_left.second : nullptr},
                                                          k_mer_total_cnt, MAX(n / 4 + 1, MIN_MER), MIN(n / 2, MAX_MER), &lef_seq);
                        }
                        if (k_mer.first == 0 or k_mer.second == 0) {
                            right_temp_k_mer = k_mer_check(temp_task.buffer2, st2, nd2, rot_table, extract_k_mer, k_mer_counter, k_mer_counter_map,
                                                        k_mer_data, k_mer_counter_list, repeat_check_table,
                                                        {k_mer.first == 0 ? temp_result_left.first : nullptr,
                                                        k_mer.second == 0 ? temp_result_left.second : nullptr},
                                                        k_mer_total_cnt, MAX(n / 4 + 1, MIN_MER), MIN(n / 2, MAX_MER), &rht_seq);

                        }

                        if (lef_k_mer.first == 0 and k_mer.first == 0 and left_temp_k_mer.first == right_temp_k_mer.first
                            and left_temp_k_mer.first > 0 and lef_seq.first == get_rot_seq(reverse_complement_64(rht_seq.first) >> (2 * (32 - right_temp_k_mer.first)), right_temp_k_mer.first)) {
                            for (auto& [seq, cnt] : *(temp_result_left.first)) {
                                seq_rev = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
                                (*(result.both.first))[KmerSeq {seq.first, MIN(seq.second, seq_rev)}] += cnt;
                            }
                        }
                        if (lef_k_mer.second == 0 and k_mer.second == 0 and left_temp_k_mer.second == right_temp_k_mer.second
                            and left_temp_k_mer.second > 0 and lef_seq.second == get_rot_seq(reverse_complement_64(rht_seq.second) >> (2 * (32 - right_temp_k_mer.second)), right_temp_k_mer.second)) {
                            for (auto& [seq, cnt] : *(temp_result_left.second)) {
                                seq_rev = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
                                (*(result.both.second))[KmerSeq {seq.first, MIN(seq.second, seq_rev)}] += cnt;
                            }
                        }

                        for (auto& [seq, cnt] : *(temp_result_left.first)) {
                            (*(result.forward.first))[seq] += cnt;
                        }
                        for (auto& [seq, cnt] : *(temp_result_left.second)) {
                            (*(result.forward.second))[seq] += cnt;
                        }
                    }
                }
            }

            delete temp_task.loc_vector1;
            free(temp_task.buffer1);

            delete temp_task.loc_vector2;
            free(temp_task.buffer2);
        }

        delete[] k_mer_counter_map;
    } else {
        CounterMap_128* k_mer_counter_map = new CounterMap_128[MAX_MER - TABLE_MAX_MER];
        std::pair<uint128_t, uint128_t> k_mer_seq, temp_k_mer_seq;
        std::pair<uint128_t, uint128_t> lef_seq, rht_seq;

        auto get_dir_seq = [](int i, int k, uint128_t seq, bool is_for) {
            if (i <= 2 == is_for) {
                return seq;
            } else {
                return get_rot_seq_128(reverse_complement_128(seq) >> (2 * (64 - k)), k);
            }
        };

        while (1) {
            task_queue -> pop(temp_task);
            if (temp_task.loc_vector1 == nullptr) {
                break;
            }

            size_t min_size = std::min(temp_task.loc_vector1 -> size(), temp_task.loc_vector2 -> size());
            for (size_t i = 0; i < min_size; i++) {
                st1 = (*temp_task.loc_vector1)[i].first;
                nd1 = (*temp_task.loc_vector1)[i].second;
                n1 = nd1 - st1 + 1;

                st2 = (*temp_task.loc_vector2)[i].first;
                nd2 = (*temp_task.loc_vector2)[i].second;
                n2 = nd2 - st2 + 1;

                n = MIN(n1, n2);

                if (2 * MIN_MER <= n) {
                    lef_k_mer = {0, 0};
                    k_mer = {0, 0};

                    if (4 * MIN_MER <= n) {
                        std::vector<const char*> buffer_vector= {temp_task.buffer1, temp_task.buffer1, temp_task.buffer2, temp_task.buffer2};
                        std::vector<std::pair<int, int>> index_vector= {{st1, st1 + (n1 / 2) - 1}, {nd1 - ((n1 + 1) / 2) + 1, nd1},
                                                                        {nd2 - ((n2 + 1) / 2) + 1, nd2}, {st2, st2 + (n2 / 2) - 1}};

                        snum = 4;
                        si = {1, 1};
                        k_mer = {0, 0};

                        for (ti = 1; ti <= snum && (!repeat_end.first || !repeat_end.second) ; ti++) {
                            temp_k_mer = k_mer_check_128(buffer_vector[ti - 1], index_vector[ti - 1].first, index_vector[ti - 1].second, rot_table, extract_k_mer_128, k_mer_counter,
                                                         k_mer_counter_map, k_mer_data_128, k_mer_counter_list, repeat_check_table,
                                                         {repeat_end.first ? nullptr : (ti <= 2 ? temp_result_left.first : temp_result_right.first),
                                                                     repeat_end.second ? nullptr : (ti <= 2 ? temp_result_left.second : temp_result_right.second)},
                                                                     k_mer_total_cnt, MIN_MER, MIN(n / 4, MAX_MER), &temp_k_mer_seq);

                            if (!repeat_end.first && temp_k_mer.first > 0 && ((k_mer.first == temp_k_mer.first and k_mer_seq.first == get_dir_seq(ti, temp_k_mer.first, temp_k_mer_seq.first, true)) || ti == 1)) {
                                si.first += 1;
                                k_mer.first = temp_k_mer.first;
                                if (ti == 1) {
                                    k_mer_seq.first = temp_k_mer_seq.first;
                                }
                                repeat_end.first = false;
                            } else {
                                repeat_end.first = true;
                            }
                            if (!repeat_end.second && temp_k_mer.second > 0 && ((k_mer.second == temp_k_mer.second and k_mer_seq.second == get_dir_seq(ti, temp_k_mer.second, temp_k_mer_seq.second, true)) || ti == 1)) {
                                si.second += 1;
                                k_mer.second = temp_k_mer.second;
                                if (ti == 1) {
                                    k_mer_seq.second = temp_k_mer_seq.second;
                                }
                                repeat_end.second = false;
                            } else {
                                repeat_end.second = true;
                            }
                        }
                        lef_k_mer = k_mer;

                        if (si.first == snum + 1) {
                            for (auto& [seq, cnt] : *(temp_result_left.first)) {
                                seq_rev = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
                                (*(result.both.first))[KmerSeq {seq.first, MIN(seq.second, seq_rev)}] += cnt;
                            }

                            for (auto& [seq, cnt] : *(temp_result_right.first)) {
                                seq_rev = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
                                (*(result.both.first))[KmerSeq {seq.first, MIN(seq.second, seq_rev)}] += cnt;
                            }
                        }
                        if (si.second == snum + 1) {
                            for (auto& [seq, cnt] : *(temp_result_left.second)) {
                                seq_rev = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
                                (*(result.both.second))[KmerSeq {seq.first, MIN(seq.second, seq_rev)}] += cnt;
                            }

                            for (auto& [seq, cnt] : *(temp_result_right.second)) {
                                seq_rev = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
                                (*(result.both.second))[KmerSeq {seq.first, MIN(seq.second, seq_rev)}] += cnt;
                            }
                        }

                        if (si.first <= snum or si.second <= snum) {
                            sj = {snum, snum};
                            k_mer = {0, 0};
                            repeat_end = {false, false};
                            for (tj = snum; !repeat_end.first || !repeat_end.second; tj--) {
                                temp_k_mer = k_mer_check_128(buffer_vector[tj - 1], index_vector[tj - 1].first, index_vector[tj - 1].second, rot_table, extract_k_mer_128, k_mer_counter,
                                                             k_mer_counter_map, k_mer_data_128, k_mer_counter_list, repeat_check_table,
                                                {repeat_end.first ? nullptr : (tj <= 2 ? temp_result_right.first : temp_result_left.first),
                                                            repeat_end.second ? nullptr : (tj <= 2 ? temp_result_right.second : temp_result_left.second)},
                                                            k_mer_total_cnt, MIN_MER, MIN(n / 4, MAX_MER), &temp_k_mer_seq);

                                if (sj.first >= si.first && !repeat_end.first && temp_k_mer.first > 0 && ((k_mer.first == temp_k_mer.first and
                                    k_mer_seq.first == get_dir_seq(tj, temp_k_mer.first, temp_k_mer_seq.first, false)) || tj == snum)) {
                                    sj.first -= 1;
                                    k_mer.first = temp_k_mer.first;
                                    if (tj == snum) {
                                        k_mer_seq.first = temp_k_mer_seq.first;
                                    }
                                    repeat_end.first = false;
                                } else {
                                    repeat_end.first = true;
                                }
                                if (sj.second >= si.second && !repeat_end.second && temp_k_mer.second > 0 && ((k_mer.second == temp_k_mer.second and
                                    k_mer_seq.second == get_dir_seq(tj, temp_k_mer.second, temp_k_mer_seq.second, false)) || tj == snum)) {
                                    sj.second -= 1;
                                    k_mer.second = temp_k_mer.second;
                                    if (tj == snum) {
                                        k_mer_seq.second = temp_k_mer_seq.second;
                                    }
                                    repeat_end.second = false;
                                } else {
                                    repeat_end.second = true;
                                }
                            }
                        }

                        if (si.first <= snum) {
                            for (auto& [seq, cnt] : *(temp_result_left.first)) {
                                (*(result.forward.first))[seq] += cnt;
                            }

                            for (auto& [seq, cnt] : *(temp_result_right.first)) {
                                (*(result.backward.first))[seq] += cnt;
                            }
                        }
                        if (si.second <= snum) {
                            for (auto& [seq, cnt] : *(temp_result_left.second)) {
                                (*(result.forward.second))[seq] += cnt;
                            }

                            for (auto& [seq, cnt] : *(temp_result_right.second)) {
                                (*(result.backward.second))[seq] += cnt;
                            }
                        }

                        temp_result_left.first -> clear();
                        temp_result_left.second -> clear();

                        temp_result_right.first -> clear();
                        temp_result_right.second -> clear();
                    }

                    left_temp_k_mer = {0, 0};
                    right_temp_k_mer = {0, 0};

                    if (4 * MAX_MER > n and (lef_k_mer.first == 0 or lef_k_mer.second == 0 or k_mer.first == 0 or k_mer.second == 0)) {
                        if (lef_k_mer.first == 0 or lef_k_mer.second == 0) {
                            left_temp_k_mer = k_mer_check_128(temp_task.buffer1, st1, nd1, rot_table, extract_k_mer_128, k_mer_counter, k_mer_counter_map,
                                                          k_mer_data_128, k_mer_counter_list, repeat_check_table,
                                                          {lef_k_mer.first == 0 ? temp_result_left.first : nullptr,
                                                          lef_k_mer.second == 0 ? temp_result_left.second : nullptr},
                                                          k_mer_total_cnt, MAX(n / 4 + 1, MIN_MER), MIN(n / 2, MAX_MER), &lef_seq);
                        }

                        if (k_mer.first == 0 or k_mer.second == 0) {
                            right_temp_k_mer = k_mer_check_128(temp_task.buffer2, st2, nd2, rot_table, extract_k_mer_128, k_mer_counter, k_mer_counter_map,
                                                        k_mer_data_128, k_mer_counter_list, repeat_check_table,
                                                        {k_mer.first == 0 ? temp_result_left.first : nullptr,
                                                        k_mer.second == 0 ? temp_result_left.second : nullptr},
                                                        k_mer_total_cnt, MAX(n / 4 + 1, MIN_MER), MIN(n / 2, MAX_MER), &rht_seq);

                        }

                        if (lef_k_mer.first == 0 and k_mer.first == 0 and left_temp_k_mer.first == right_temp_k_mer.first
                            and left_temp_k_mer.first > 0 and lef_seq.first == get_rot_seq_128(reverse_complement_128(rht_seq.first) >> (2 * (64 - right_temp_k_mer.first)), right_temp_k_mer.first)) {
                            for (auto& [seq, cnt] : *(temp_result_left.first)) {
                                seq_rev = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
                                (*(result.both.first))[KmerSeq {seq.first, MIN(seq.second, seq_rev)}] += cnt;
                            }
                        }

                        if (lef_k_mer.second == 0 and k_mer.second == 0 and left_temp_k_mer.second == right_temp_k_mer.second
                            and left_temp_k_mer.second > 0 and lef_seq.second == get_rot_seq_128(reverse_complement_128(rht_seq.second) >> (2 * (64 - right_temp_k_mer.second)), right_temp_k_mer.second)) {
                            for (auto& [seq, cnt] : *(temp_result_left.second)) {
                                seq_rev = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
                                (*(result.both.second))[KmerSeq {seq.first, MIN(seq.second, seq_rev)}] += cnt;
                            }
                        }

                        for (auto& [seq, cnt] : *(temp_result_left.first)) {
                            (*(result.forward.first))[seq] += cnt;
                        }
                        for (auto& [seq, cnt] : *(temp_result_left.second)) {
                            (*(result.forward.second))[seq] += cnt;
                        }

                        temp_result_left.first -> clear();
                        temp_result_left.second -> clear();
                    }
                }
            }

            delete temp_task.loc_vector1;
            free(temp_task.buffer1);

            delete temp_task.loc_vector2;
            free(temp_task.buffer2);
        }

        delete[] k_mer_counter_map;
    }

    delete temp_result_left.first;
    delete temp_result_left.second;
    delete temp_result_right.first;
    delete temp_result_right.second;

    free(k_mer_total_cnt);
    return result;
}

ResultMapData buffer_task_long(TBBQueue* task_queue, ThreadData* thread_data, uint32_t** rot_table, const uint64_t *extract_k_mer, const uint128_t *extract_k_mer_128, uint8_t** repeat_check_table) {
    auto [k_mer_counter, k_mer_data, k_mer_data_128, k_mer_counter_list] = thread_data -> init_check();

    int16_t* k_mer_total_cnt = (int16_t*) malloc(sizeof(int16_t) * (MAX_MER - MIN_MER + 2));
    if (k_mer_total_cnt == nullptr) {
        fprintf(stderr, "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    ResultMapData result = ResultMapData {{new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}};

    ResultMapPair temp_result_left = {new ResultMap {}, new ResultMap {}};

    QueueData temp_task {};

    KmerData k_mer, temp_k_mer;
    int tst, tnd, ti, tj;
    KmerData si, sj;

    std::pair<bool, bool> repeat_end;

    int snum;
    int mid;
    int sl, mid_bonus_sl;

    uint128_t seq_rev;

    if (MAX_MER <= ABS_UINT64_MAX_MER) {
        CounterMap* k_mer_counter_map = nullptr;
        if (TABLE_MAX_MER < MAX_MER) {
            k_mer_counter_map = new CounterMap[MAX_MER - TABLE_MAX_MER];
        }

        while (1) {
            task_queue -> pop(temp_task);
            if (temp_task.loc_vector == nullptr) {
                break;
            }
            for (auto& [st, nd]: *(temp_task.loc_vector)) {
                std::pair<uint128_t, uint128_t> lef_seq, rht_seq;
                tst = st;
                tnd = nd;

                snum = (tnd - tst + 1) / SLICE_LENGTH;
                mid = (snum + 1) / 2;
                mid_bonus_sl = (tnd - tst + 1) % SLICE_LENGTH;

                si = {1, 1};
                k_mer = {0, 0};
                repeat_end = {false, false};
                for (ti = 1; ti <= snum && (!repeat_end.first || !repeat_end.second) ; ti++, tst += sl) {
                    sl = SLICE_LENGTH + (ti == mid ? mid_bonus_sl : 0);
                    temp_k_mer = k_mer_check(temp_task.buffer, tst, tst + sl - 1, rot_table, extract_k_mer, k_mer_counter,
                                                 k_mer_counter_map, k_mer_data, k_mer_counter_list, repeat_check_table,
                                                 {repeat_end.first ? nullptr : temp_result_left.first, repeat_end.second ? nullptr : temp_result_left.second}, k_mer_total_cnt, MIN_MER, MAX_MER);

                    if (!repeat_end.first && temp_k_mer.first > 0 && (k_mer.first == temp_k_mer.first|| ti == 1)) {
                        si.first += 1;
                        k_mer.first = temp_k_mer.first;
                        repeat_end.first = false;
                    } else {
                        repeat_end.first = true;
                    }
                    if (!repeat_end.second && temp_k_mer.second > 0 && (k_mer.second == temp_k_mer.second || ti == 1)) {
                        si.second += 1;
                        k_mer.second = temp_k_mer.second;
                        repeat_end.second = false;
                    } else {
                        repeat_end.second = true;
                    }
                }

                if (si.first == snum + 1) {
                    for (auto& [seq, cnt] : *(temp_result_left.first)) {
                        seq_rev = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
                        (*(result.both.first))[KmerSeq {seq.first, MIN(seq.second, seq_rev)}] += cnt;
                    }
                }
                if (si.second == snum + 1) {
                    for (auto& [seq, cnt] : *(temp_result_left.second)) {
                        seq_rev = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
                        (*(result.both.second))[KmerSeq {seq.first, MIN(seq.second, seq_rev)}] += cnt;
                    }
                }

                if (si.first <= snum or si.second <= snum) {
                    sj = {snum, snum};
                    k_mer = {0, 0};
                    repeat_end = {false, false};
                    for (tj = snum; !repeat_end.first || !repeat_end.second; tj--, tnd -= sl) {
                        sl = SLICE_LENGTH + (tj == mid ? mid_bonus_sl : 0);
                        temp_k_mer = k_mer_check(temp_task.buffer, tnd - sl + 1, tnd, rot_table, extract_k_mer, k_mer_counter,
                                                     k_mer_counter_map, k_mer_data, k_mer_counter_list, repeat_check_table,
                                                     {repeat_end.first ? nullptr : result.backward.first, repeat_end.second ? nullptr : result.backward.second}, k_mer_total_cnt, MIN_MER, MAX_MER);

                        if (sj.first >= si.first && !repeat_end.first && temp_k_mer.first > 0 && (k_mer.first == temp_k_mer.first || tj == snum)) {
                            sj.first -= 1;
                            k_mer.first = temp_k_mer.first;
                            repeat_end.first = false;
                        } else {
                            repeat_end.first = true;
                        }
                        if (sj.second >= si.second && !repeat_end.second && temp_k_mer.second > 0 && (k_mer.second == temp_k_mer.second || tj == snum)) {
                            sj.second -= 1;
                            k_mer.second = temp_k_mer.second;
                            repeat_end.second = false;
                        } else {
                            repeat_end.second = true;
                        }
                    }

                    if (si.first <= snum) {
                        for (auto& [seq, cnt] : *(temp_result_left.first)) {
                            (*(result.forward.first))[seq] += cnt;
                        }
                    }
                    if (si.second <= snum) {
                        for (auto& [seq, cnt] : *(temp_result_left.second)) {
                            (*(result.forward.second))[seq] += cnt;
                        }
                    }
                }
                temp_result_left.first -> clear();
                temp_result_left.second -> clear();
            }

            delete temp_task.loc_vector;
            free(temp_task.buffer);
        }

        delete[] k_mer_counter_map;
    } else {
        CounterMap_128* k_mer_counter_map = new CounterMap_128[MAX_MER - TABLE_MAX_MER];

        while (1) {
            task_queue -> pop(temp_task);
            if (temp_task.loc_vector == nullptr) {
                break;
            }
            for (auto& [st, nd]: *(temp_task.loc_vector)) {
                std::pair<uint128_t, uint128_t> lef_seq, rht_seq;
                tst = st;
                tnd = nd;

                snum = (tnd - tst + 1) / SLICE_LENGTH;
                mid = (snum + 1) / 2;
                mid_bonus_sl = (tnd - tst + 1) % SLICE_LENGTH;

                si = {1, 1};
                k_mer = {0, 0};
                repeat_end = {false, false};
                for (ti = 1; ti <= snum && (!repeat_end.first || !repeat_end.second) ; ti++, tst += sl) {
                    sl = SLICE_LENGTH + (ti == mid ? mid_bonus_sl : 0);
                    temp_k_mer = k_mer_check_128(temp_task.buffer, tst, tst + sl - 1, rot_table, extract_k_mer_128, k_mer_counter,
                                                 k_mer_counter_map, k_mer_data_128, k_mer_counter_list, repeat_check_table,
                                                 {repeat_end.first ? nullptr : temp_result_left.first, repeat_end.second ? nullptr : temp_result_left.second}, k_mer_total_cnt, MIN_MER, MAX_MER, &lef_seq);

                    if (!repeat_end.first && temp_k_mer.first > 0 && (k_mer.first == temp_k_mer.first|| ti == 1)) {
                        si.first += 1;
                        k_mer.first = temp_k_mer.first;
                        repeat_end.first = false;
                    } else {
                        repeat_end.first = true;
                    }
                    if (!repeat_end.second && temp_k_mer.second > 0 && (k_mer.second == temp_k_mer.second || ti == 1)) {
                        si.second += 1;
                        k_mer.second = temp_k_mer.second;
                        repeat_end.second = false;
                    } else {
                        repeat_end.second = true;
                    }
                }

                if (si.first == snum + 1) {
                    for (auto& [seq, cnt] : *(temp_result_left.first)) {
                        seq_rev = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
                        (*(result.both.first))[KmerSeq {seq.first, MIN(seq.second, seq_rev)}] += cnt;
                    }
                }
                if (si.second == snum + 1) {
                    for (auto& [seq, cnt] : *(temp_result_left.second)) {
                        seq_rev = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
                        (*(result.both.second))[KmerSeq {seq.first, MIN(seq.second, seq_rev)}] += cnt;
                    }
                }

                if (si.first <= snum || si.second <= snum) {
                    sj = {snum, snum};
                    k_mer = {0, 0};
                    repeat_end = {false, false};
                    for (tj = snum; !repeat_end.first || !repeat_end.second; tj--, tnd -= sl) {
                        sl = SLICE_LENGTH + (tj == mid ? mid_bonus_sl : 0);
                        temp_k_mer = k_mer_check_128(temp_task.buffer, tnd - sl + 1, tnd, rot_table, extract_k_mer_128, k_mer_counter,
                                                     k_mer_counter_map, k_mer_data_128, k_mer_counter_list, repeat_check_table,
                                                     {repeat_end.first ? nullptr : result.backward.first, repeat_end.second ? nullptr : result.backward.second}, k_mer_total_cnt, MIN_MER, MAX_MER, &rht_seq);

                        if (sj.first >= si.first && !repeat_end.first && temp_k_mer.first > 0 && (k_mer.first == temp_k_mer.first || tj == snum)) {
                            sj.first -= 1;
                            k_mer.first = temp_k_mer.first;
                            repeat_end.first = false;
                        } else {
                            repeat_end.first = true;
                        }
                        if (sj.second >= si.second && !repeat_end.second && temp_k_mer.second > 0 && (k_mer.second == temp_k_mer.second || tj == snum)) {
                            sj.second -= 1;
                            k_mer.second = temp_k_mer.second;
                            repeat_end.second = false;
                        } else {
                            repeat_end.second = true;
                        }
                    }

                    if (si.first <= snum) {
                        for (auto& [seq, cnt] : *(temp_result_left.first)) {
                            (*(result.forward.first))[seq] += cnt;
                        }
                    }
                    if (si.second <= snum) {
                        for (auto& [seq, cnt] : *(temp_result_left.second)) {
                            (*(result.forward.second))[seq] += cnt;
                        }
                    }
                }
                temp_result_left.first -> clear();
                temp_result_left.second -> clear();
            }

            delete temp_task.loc_vector;
            free(temp_task.buffer);
        }

        delete[] k_mer_counter_map;
    }

    delete temp_result_left.first;
    delete temp_result_left.second;
    free(k_mer_total_cnt);
    return result;
}

void read_fastq_thread(FileReader& file_reader, TBBQueue* buffer_task_queue) {
    int num = 0;
    int shift = 0;
    int bytes_read;
    int idx = -1;
    char* buffer;
    char* buffer_new;

    buffer = (char*) malloc(sizeof(char) * LENGTH);
    while (1) {
        LocationVector* loc_vector = new LocationVector{};

        bytes_read = file_reader.read(buffer + shift, LENGTH - 1 - shift);
        buffer[bytes_read + shift] = '\0';

        for (int i = 0; i < bytes_read + shift; i++) {
            if (buffer[i] == '\n') {
                num += 1;
                if ((num & 3) == 2) {
                    if ((i - 1) - (idx + 1) + 1 > MAX_SEQ) {
                        fprintf(stderr, "This mode is designed for short-read sequencing. Please use 'trew long'.\n");
                        exit(EXIT_FAILURE);
                    }
                    loc_vector->emplace_back(idx + 1, i - 1);
                }
                idx = i;
            }
        }

        if (bytes_read <= 0) {
            buffer_task_queue->push(QueueData{buffer, loc_vector});
            if (file_reader.eof()) {
                break;
            } else {
                fprintf(stderr, "File-IO Error: %s.\n", file_reader.error());
                exit(EXIT_FAILURE);
            }
        } else {
            buffer_new = (char*) malloc(sizeof(char) * LENGTH);
            if ((num & 3) == 1) {
                strcpy(buffer_new, buffer + idx + 1);
                shift = bytes_read + shift - idx - 1;
                idx = -1;
            } else {
                shift = 0;
            }

            buffer_task_queue->push(QueueData{buffer, loc_vector});
            buffer = buffer_new;
        }
    }
}

void read_pair_fastq_thread(FileReader& file_reader1, FileReader& file_reader2, TBBPairQueue* buffer_task_queue) {
    int num1 = 0, num2 = 0;
    int shift1 = 0, shift2 = 0;
    int bytes_read1, bytes_read2;
    int idx1 = -1, idx2 = -1;
    bool is_end1 = false, is_end2 = false;
    char* buffer1;
    char* buffer2;
    char* buffer_new1;
    char* buffer_new2;

    buffer1 = (char*) malloc(sizeof(char) * LENGTH);
    buffer2 = (char*) malloc(sizeof(char) * LENGTH);

    while (true) {
        LocationVector* loc_vector1 = new LocationVector{};
        LocationVector* loc_vector2 = new LocationVector{};

        if (!is_end1) {
            bytes_read1 = file_reader1.read(buffer1 + shift1, LENGTH - 1 - shift1);
            if (bytes_read1 <= 0) {
                if (file_reader1.eof()) {
                    is_end1 = true;
                } else {
                    fprintf(stderr, "File 1 IO Error: %s.\n", file_reader1.error());
                    exit(EXIT_FAILURE);
                }
            }
        } else {
            bytes_read1 = 0;
        }

        if (!is_end2) {
            bytes_read2 = file_reader2.read(buffer2 + shift2, LENGTH - 1 - shift2);
            if (bytes_read2 <= 0) {
                if (file_reader2.eof()) {
                    is_end2 = true;
                } else {
                    fprintf(stderr, "File 2 IO Error: %s.\n", file_reader2.error());
                    exit(EXIT_FAILURE);
                }
            }
        } else {
            bytes_read2 = 0;
        }

        buffer1[bytes_read1 + shift1] = '\0';
        buffer2[bytes_read2 + shift2] = '\0';

        for (int i = 0; i < bytes_read1 + shift1; i++) {
            if (buffer1[i] == '\n') {
                num1 += 1;
                if ((num1 & 3) == 2) {
                    loc_vector1->emplace_back(idx1 + 1, i - 1);
                }
                idx1 = i;
            }
        }

        for (int i = 0; i < bytes_read2 + shift2; i++) {
            if (buffer2[i] == '\n') {
                num2 += 1;
                if ((num2 & 3) == 2) {
                    loc_vector2->emplace_back(idx2 + 1, i - 1);
                }
                idx2 = i;
            }
        }

        if (is_end1 and is_end2) {
            if (num1 != num2) {
                fprintf(stderr, "Error: Mismatched record counts between files (num1: %d, num2: %d).\n", num1, num2);
                exit(EXIT_FAILURE);
            }
            buffer_task_queue->push(PairQueueData{buffer1, buffer2, loc_vector1, loc_vector2});
            break;
        }

        if (loc_vector1 -> size() == 0 and loc_vector2 -> size() > 0 or
            loc_vector1 -> size() > 0 and loc_vector2 -> size() == 0) {
            fprintf(stderr, "Paired-end error\n");
            exit(EXIT_FAILURE);
        }

        buffer_new1 = (char*) malloc(sizeof(char) * LENGTH);
        buffer_new2 = (char*) malloc(sizeof(char) * LENGTH);

        size_t min_size = std::min(loc_vector1->size(), loc_vector2->size());

        if (loc_vector1->size() > min_size) {
            idx1 = (*loc_vector1)[min_size].first - 1;
            strcpy(buffer_new1, buffer1 + idx1 + 1);
            shift1 = bytes_read1 + shift1 - idx1 - 1;
            num1 = ((num1 - 2) / 4) * 4 + 2 - 4 * (loc_vector1->size() - min_size);
            idx1 = -1;
        } else if ((num1 & 3) == 1) {
            strcpy(buffer_new1, buffer1 + idx1 + 1);
            shift1 = bytes_read1 + shift1 - idx1 - 1;
            idx1 = -1;
        } else {
            shift1 = 0;
        }

        if (loc_vector2->size() > min_size) {
            idx2 = (*loc_vector2)[min_size].first - 1;
            strcpy(buffer_new2, buffer2 + idx2 + 1);
            shift2 = bytes_read2 + shift2 - idx2 - 1;
            num2 = ((num2 - 2) / 4) * 4 + 2 - 4 * (loc_vector2->size() - min_size);
            idx2 = -1;
        } else if ((num2 & 3) == 1) {
            strcpy(buffer_new2, buffer2 + idx2 + 1);
            shift2 = bytes_read2 + shift2 - idx2 - 1;
            idx2 = -1;
        } else {
            shift2 = 0;
        }

        buffer_task_queue->push(PairQueueData{buffer1, buffer2, loc_vector1, loc_vector2});
        buffer1 = buffer_new1;
        buffer2 = buffer_new2;
    }

}

void read_fastq_long_thread(FileReader& file_reader, TBBQueue* buffer_task_queue) {
    int num = 0;
    int shift = 0;
    int bytes_read;
    int idx = -1;
    char* buffer;
    char* buffer_new;

    buffer = (char*) malloc(sizeof(char) * LENGTH);
    while (1) {
        LocationVector* loc_vector = new LocationVector{};

        bytes_read = file_reader.read(buffer + shift, LENGTH - 1 - shift);
        buffer[bytes_read + shift] = '\0';

        for (int i = 0; i < bytes_read + shift; i++) {
            if (buffer[i] == '\n') {
                num += 1;
                if ((num & 3) == 2 && (i - 1) - (idx + 1) + 1 >= SLICE_LENGTH) {
                    loc_vector -> emplace_back(idx + 1, i - 1);
                }
                idx = i;
            }
        }

        if (bytes_read <= 0) {
            buffer_task_queue->push(QueueData{buffer, loc_vector});
            if (file_reader.eof()) {
                break;
            } else {
                fprintf(stderr, "File-IO Error: %s.\n", file_reader.error());
                exit(EXIT_FAILURE);
            }
        } else {
            buffer_new = (char*) malloc(sizeof(char) * LENGTH);
            if ((num & 3) == 1) {
                strcpy(buffer_new, buffer + idx + 1);
                shift = bytes_read + shift - idx - 1;
                idx = -1;
            } else {
                shift = 0;
            }

            buffer_task_queue->push(QueueData{buffer, loc_vector});
            buffer = buffer_new;
        }
    }
}

void read_fastq_gz_thread_long(gzFile fp, TBBQueue* buffer_task_queue) {
    int num = 0;
    int shift = 0;
    int bytes_read;
    int idx;
    char* buffer;
    char* buffer_new;

    buffer = (char*) malloc(sizeof(char) * LENGTH);
    while (1) {
        LocationVector* loc_vector = new LocationVector{};

        bytes_read = gzread(fp, buffer + shift, LENGTH - 1 - shift);
        buffer[bytes_read + shift] = '\0';

        for (int i = 0; i < bytes_read + shift; i++) {
            if (buffer[i] == '\n') {
                num += 1;
                if ((num & 3) == 2 && (i - 1) - (idx + 1) + 1 >= SLICE_LENGTH) {
                    loc_vector -> emplace_back(idx + 1, i - 1);
                }
                idx = i;
            }
        }

        if (bytes_read < LENGTH - 1 - shift) {
            buffer_task_queue -> push(QueueData{buffer, loc_vector});
            if (gzeof(fp)) {
                break;
            }
            else {
                fprintf(stderr, "File-IO Error: %s.\n", strerror(errno));
                exit(EXIT_FAILURE);
            }
        }
        else {
            buffer_new = (char*) malloc(sizeof(char) * LENGTH);
            if ((num & 3) == 1) {
                strcpy(buffer_new, buffer + idx + 1);
                shift = bytes_read + shift - idx - 1;
                idx = -1;
            }
            else {
                shift = 0;
            }
            buffer_task_queue -> push(QueueData{buffer, loc_vector});
            buffer = buffer_new;
        }
    }
}

FinalFastqOutput process_kmer(const char* file_name, uint8_t **repeat_check_table, uint32_t **rot_table,
                                  uint64_t *extract_k_mer, uint128_t *extract_k_mer_128, uint128_t *extract_k_mer_ans,
                                  ThreadData* thread_data_list, bool is_gz) {
    ResultMapData* result_list = (ResultMapData*) malloc(sizeof(ResultMapData) * NUM_THREAD);

    TBBQueue buffer_task_queue {};
    TaskGroup tasks;

    if (QUEUE_SIZE >= 4) {
        buffer_task_queue.set_capacity(QUEUE_SIZE / 4);
    }

    for (int i = 0; i < NUM_THREAD - 1; i++) {
        tasks.run([&result_list, i, &buffer_task_queue, &thread_data_list, &rot_table, &extract_k_mer, &extract_k_mer_128, &repeat_check_table]{
            result_list[i] = buffer_task(&buffer_task_queue, thread_data_list + i, rot_table, extract_k_mer, extract_k_mer_128, repeat_check_table);
        });
    }

    FileReader reader;
    if (is_gz) {
        gzFile fp = gzopen(file_name, "r");
        if (fp == nullptr) {
            fprintf(stderr, "File open failed\n");
            exit (EXIT_FAILURE);
        }
        reader = fp;
    } else {
        FILE* fp = fopen(file_name, "r");
        if (fp == nullptr) {
            fprintf(stderr, "File open failed\n");
            exit (EXIT_FAILURE);
        }
        reader = fp;
    }

    read_fastq_thread(reader, &buffer_task_queue);
    reader.close();

    for (int i = 1; i < NUM_THREAD; i ++) {
        buffer_task_queue.push(QueueData{nullptr, nullptr});
    }

    if (!buffer_task_queue.empty()) {
        buffer_task_queue.push(QueueData{nullptr, nullptr});

        int i = NUM_THREAD - 1;
        if (NUM_THREAD > 1) {
            tasks.run([&result_list, i, &buffer_task_queue, &thread_data_list, &rot_table, &extract_k_mer, &extract_k_mer_128, &repeat_check_table]{
                result_list[i] = buffer_task(&buffer_task_queue, thread_data_list + i, rot_table, extract_k_mer, extract_k_mer_128, repeat_check_table);
            });
        } else {
            result_list[i] = buffer_task(&buffer_task_queue, thread_data_list + i, rot_table, extract_k_mer, extract_k_mer_128, repeat_check_table);
        }
    } else {
        result_list[NUM_THREAD - 1] = ResultMapData{{new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}};
    }

    if (NUM_THREAD > 1) {
        tasks.wait();
    }

    return process_output(file_name, result_list, rot_table, extract_k_mer_ans);
}

FinalFastqOutput process_kmer_pair(const char* file_name1, const char* file_name2, uint8_t **repeat_check_table, uint32_t **rot_table,
                                   uint64_t *extract_k_mer, uint128_t *extract_k_mer_128, uint128_t *extract_k_mer_ans,
                                   ThreadData* thread_data_list, bool is_gz1, bool is_gz2) {
    ResultMapData* result_list = (ResultMapData*) malloc(sizeof(ResultMapData) * NUM_THREAD);

    TBBPairQueue buffer_task_queue {};
    TaskGroup tasks;

    if (QUEUE_SIZE >= 4) {
        buffer_task_queue.set_capacity(QUEUE_SIZE / 4);
    }

    for (int i = 0; i < NUM_THREAD - 1; i++) {
        tasks.run([&result_list, i, &buffer_task_queue, &thread_data_list, &rot_table, &extract_k_mer, &extract_k_mer_128, &repeat_check_table]{
            result_list[i] = buffer_task_pair(&buffer_task_queue, thread_data_list + i, rot_table, extract_k_mer, extract_k_mer_128, repeat_check_table);
        });
    }

    FileReader reader1;
    if (is_gz1) {
        gzFile fp = gzopen(file_name1, "r");
        if (fp == nullptr) {
            fprintf(stderr, "File open failed\n");
            exit (EXIT_FAILURE);
        }
        reader1 = fp;
    } else {
        FILE* fp = fopen(file_name1, "r");
        if (fp == nullptr) {
            fprintf(stderr, "File open failed\n");
            exit (EXIT_FAILURE);
        }
        reader1 = fp;
    }

    FileReader reader2;
    if (is_gz2) {
        gzFile fp = gzopen(file_name2, "r");
        if (fp == nullptr) {
            fprintf(stderr, "File open failed\n");
            exit (EXIT_FAILURE);
        }
        reader2 = fp;
    } else {
        FILE* fp = fopen(file_name2, "r");
        if (fp == nullptr) {
            fprintf(stderr, "File open failed\n");
            exit (EXIT_FAILURE);
        }
        reader2 = fp;
    }

    read_pair_fastq_thread(reader1, reader2, &buffer_task_queue);
    reader1.close();
    reader2.close();

    for (int i = 1; i < NUM_THREAD; i ++) {
        buffer_task_queue.push(PairQueueData{nullptr, nullptr});
    }

    if (!buffer_task_queue.empty()) {
        buffer_task_queue.push(PairQueueData{nullptr, nullptr});

        int i = NUM_THREAD - 1;
        if (NUM_THREAD > 1) {
            tasks.run([&result_list, i, &buffer_task_queue, &thread_data_list, &rot_table, &extract_k_mer, &extract_k_mer_128, &repeat_check_table]{
                result_list[i] = buffer_task_pair(&buffer_task_queue, thread_data_list + i, rot_table, extract_k_mer, extract_k_mer_128, repeat_check_table);
            });
        } else {
            result_list[i] = buffer_task_pair(&buffer_task_queue, thread_data_list + i, rot_table, extract_k_mer, extract_k_mer_128, repeat_check_table);
        }
    } else {
        result_list[NUM_THREAD - 1] = ResultMapData{{new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}};
    }

    if (NUM_THREAD > 1) {
        tasks.wait();
    }

    return process_output(file_name1, result_list, rot_table, extract_k_mer_ans);
}

FinalFastqOutput process_kmer_long(const char* file_name, uint8_t **repeat_check_table, uint32_t **rot_table,
                                       uint64_t *extract_k_mer, uint128_t *extract_k_mer_128, uint128_t *extract_k_mer_ans,
                                       ThreadData* thread_data_list, bool is_gz) {
    ResultMapData* result_list = (ResultMapData*) malloc(sizeof(ResultMapData) * NUM_THREAD);

    TBBQueue buffer_task_queue {};
    TaskGroup tasks;

    if (QUEUE_SIZE >= 4) {
        buffer_task_queue.set_capacity(QUEUE_SIZE / 4);
    }

    for (int i = 0; i < NUM_THREAD - 1; i++) {
        tasks.run([&result_list, i, &buffer_task_queue, &thread_data_list, &rot_table, &extract_k_mer, &extract_k_mer_128, &repeat_check_table]{
            result_list[i] = buffer_task_long(&buffer_task_queue, thread_data_list + i, rot_table, extract_k_mer, extract_k_mer_128, repeat_check_table);
        });
    }

    FileReader reader;
    if (is_gz) {
        FILE* fp = fopen(file_name, "rb");
        if (fp == nullptr) {
            fprintf(stderr, "File open failed\n");
            exit (EXIT_FAILURE);
        }

        reader = fp;
    } else {
        FILE* fp = fopen(file_name, "r");
        if (fp == nullptr) {
            fprintf(stderr, "File open failed\n");
            exit (EXIT_FAILURE);
        }

        reader = fp;
    }

    read_fastq_long_thread(reader, &buffer_task_queue);
    reader.close();

    for (int i = 1; i < NUM_THREAD; i ++) {
        buffer_task_queue.push(QueueData{nullptr, nullptr});
    }

    if (!buffer_task_queue.empty()) {
        buffer_task_queue.push(QueueData{nullptr, nullptr});

        int i = NUM_THREAD - 1;
        if (NUM_THREAD > 1) {
            tasks.run([&result_list, i, &buffer_task_queue, &thread_data_list, &rot_table, &extract_k_mer, &extract_k_mer_128, &repeat_check_table]{
                result_list[i] = buffer_task_long(&buffer_task_queue, thread_data_list + i, rot_table, extract_k_mer, extract_k_mer_128, repeat_check_table);
            });
        } else {
            result_list[i] = buffer_task_long(&buffer_task_queue, thread_data_list + i, rot_table, extract_k_mer, extract_k_mer_128, repeat_check_table);
        }
    } else {
        result_list[NUM_THREAD - 1] = ResultMapData{{new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}};
    }

    if (NUM_THREAD > 1) {
        tasks.wait();
    }

    return process_output(file_name, result_list, rot_table, extract_k_mer_ans);
}

FinalFastqOutput process_output(const char* file_name, ResultMapData* result_list, uint32_t **rot_table, uint128_t *extract_k_mer_ans) {
    uint128_t _t;
    uint128_t kseq;
    int64_t _tcnt;


    char* buffer = (char*) malloc(sizeof(char) * (ABS_MAX_MER + 1));

    ResultMapData result_data = result_list[0];
    for (int i = 1; i < NUM_THREAD; i++) {
        for (auto &[k, v]: *(result_list[i].forward.first)) {
            (*(result_data.forward.first))[k] += v;
        }
        for (auto &[k, v]: *(result_list[i].backward.first)) {
            (*(result_data.backward.first))[k] += v;
        }
        for (auto &[k, v]: *(result_list[i].both.first)) {
            (*(result_data.both.first))[k] += v;
        }

        for (auto &[k, v]: *(result_list[i].forward.second)) {
            (*(result_data.forward.second))[k] += v;
        }
        for (auto &[k, v]: *(result_list[i].backward.second)) {
            (*(result_data.backward.second))[k] += v;
        }
        for (auto &[k, v]: *(result_list[i].both.second)) {
            (*(result_data.both.second))[k] += v;
        }

        delete result_list[i].forward.first;
        delete result_list[i].backward.first;
        delete result_list[i].both.first;

        delete result_list[i].forward.second;
        delete result_list[i].backward.second;
        delete result_list[i].both.second;
    }
    // free(result_list);

    for (auto& [seq, cnt] : *(result_data.backward.first)) {
        (*(result_data.forward.first))[KmerSeq {seq.first, get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first)}] += cnt;
    }
    for (auto& [seq, cnt] : *(result_data.backward.second)) {
        (*(result_data.forward.second))[KmerSeq {seq.first, get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first)}] += cnt;
    }

    FinalFastqData final_result_low = FinalFastqData {};
    for (auto& [seq, cnt] : *(result_data.forward.second)) {
        _t = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
        kseq = MIN(_t, seq.second);

        if (not final_result_low.contains(KmerSeq {seq.first, kseq})) {
            final_result_low[KmerSeq {seq.first, kseq}] = FinalData<int64_t> {0, _t == seq.second ? -1 : 0, 0};
        }

        if (kseq == seq.second) {
            final_result_low[KmerSeq {seq.first, kseq}].forward = cnt;
        }
        else {
            final_result_low[KmerSeq {seq.first, kseq}].backward = cnt;
        }
    }
    for (auto& [seq, cnt] : *(result_data.both.second)) {
        _t = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
        if (final_result_low.contains(seq)) {
            final_result_low[seq].both = cnt;
        }
        else {
            final_result_low[seq] = FinalData<int64_t> {0, _t == seq.second ? -1 : 0, cnt};
        }
    }

    delete result_data.forward.second;
    delete result_data.backward.second;
    delete result_data.both.second;

    FinalFastqData final_result_high = FinalFastqData {};
    for (auto& [seq, cnt] : *(result_data.forward.first)) {
        _t = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
        kseq = MIN(_t, seq.second);

        if (not final_result_high.contains(KmerSeq {seq.first, kseq})) {
            final_result_high[KmerSeq {seq.first, kseq}] = FinalData<int64_t> {0, _t == seq.second ? -1 : 0, 0};
        }

        if (kseq == seq.second) {
            final_result_high[KmerSeq {seq.first, kseq}].forward = cnt;
        }
        else {
            final_result_high[KmerSeq {seq.first, kseq}].backward = cnt;
        }
    }
    for (auto& [seq, cnt] : *(result_data.both.first)) {
        _t = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
        if (final_result_high.contains(seq)) {
            final_result_high[seq].both = cnt;
        }
        else {
            final_result_high[seq] = FinalData<int64_t> {0, _t == seq.second ? -1 : 0, cnt};
        }
    }

    delete result_data.forward.first;
    delete result_data.backward.first;
    delete result_data.both.first;

    FinalFastqVector* final_result_low_vector = new FinalFastqVector {};
    for (auto& [k, v] : final_result_low) {
        if (check_ans_seq(k, extract_k_mer_ans, rot_table)) {
            final_result_low_vector -> emplace_back(k, v);
        }
    }

    std::sort(final_result_low_vector -> begin(), final_result_low_vector -> end(), [](auto &a, auto &b) {
        if (a.second.forward == b.second.forward) {
            return a.second.both > b.second.both;
        } else {
            return a.second.forward > b.second.forward;
        }
    });

    FinalFastqVector* final_result_high_vector = new FinalFastqVector {};
    for (auto& [k, v] : final_result_high) {
        if (check_ans_seq(k, extract_k_mer_ans, rot_table)) {
            final_result_high_vector -> emplace_back(k, v);
        }
    }

    std::sort(final_result_high_vector -> begin(), final_result_high_vector -> end(), [](auto &a, auto &b) {
        if (a.second.forward == b.second.forward) {
            return a.second.both > b.second.both;
        } else {
            return a.second.forward > b.second.forward;
        }
    });

    fprintf(stdout, ">H:%s\n", file_name);
    for (auto& [k, v] : *final_result_high_vector) {
        if (v.forward + v.backward + v.both >= ABS_MIN_PRINT_COUNT) {
            int_to_four(buffer, k.second, k.first);
            fprintf(stdout, "%d,%s,%" PRId64",%" PRId64",%" PRId64",%c\n", k.first, buffer, MAX(v.forward, v.backward), MIN(v.forward, v.backward), v.both,
                    v.forward > v.backward ? '+' : (v.forward < v.backward ? '-' : '?'));
        }
    }

    fprintf(stdout, ">L:%s\n", file_name);
    for (auto& [k, v] : *final_result_low_vector) {
        if (v.forward + v.backward + v.both >= ABS_MIN_PRINT_COUNT) {
            int_to_four(buffer, k.second, k.first);
            fprintf(stdout, "%d,%s,%" PRId64",%" PRId64",%" PRId64",%c\n", k.first, buffer, MAX(v.forward, v.backward), MIN(v.forward, v.backward), v.both,
                    v.forward > v.backward ? '+' : (v.forward < v.backward ? '-' : '?'));
        }
    }

    return FinalFastqOutput {final_result_high_vector, final_result_low_vector};
}

uint8_t** set_repeat_check_table() {
    uint8_t** repeat_check_table = (uint8_t**) malloc(sizeof(uint8_t*) * (MIN(MAX_MER, TABLE_MAX_MER) - ABS_MIN_MER + 1));
    if (repeat_check_table == nullptr) {
        fprintf(stderr, "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    uint32_t max_ind;
    for (int i = 0; i <= MIN(MAX_MER, TABLE_MAX_MER) - ABS_MIN_MER; i++) {
        max_ind = 1 << (2 * (i + ABS_MIN_MER));

        repeat_check_table[i] = (uint8_t*) calloc(max_ind, sizeof(uint16_t));
        if (repeat_check_table[i] == nullptr) {
            fprintf(stderr, "memory allocation failure\n");
            exit(EXIT_FAILURE);
        }
    }

    return repeat_check_table;
}


uint64_t* set_extract_k_mer() {
    uint64_t* extract_k_mer = (uint64_t*) malloc(sizeof(uint64_t) * (MAX_MER - MIN_MER + 1));
    for (int i = 0; i <= MAX_MER - MIN_MER; i ++) {
        extract_k_mer[i] = MIN_MER + i < 32 ? (1ULL << (2 * (i + MIN_MER))) - 1 : -1ULL;
    }
    return extract_k_mer;
}

uint128_t* set_extract_k_mer_128() {
    uint128_t* extract_k_mer = (uint128_t*) malloc(sizeof(uint128_t) * (MAX_MER - MIN_MER + 1));
    for (int i = 0; i <= MAX_MER - MIN_MER; i ++) {
        extract_k_mer[i] = MIN_MER + i < 64 ? (uint128_t(1) << (2 * (i + MIN_MER))) - 1 : uint128_t(-1);
    }
    return extract_k_mer;
}

uint128_t* set_extract_k_mer_ans() {
    uint128_t* extract_k_mer = (uint128_t*) malloc(sizeof(uint128_t) * (MIN_MER - ABS_MIN_MER));
    for (int i = 0; i < MIN_MER - ABS_MIN_MER; i ++) {
        extract_k_mer[i] = ABS_MIN_MER + i < 64 ? (uint128_t(1) << (2 * (i + ABS_MIN_MER))) - 1 : uint128_t(-1);
    }
    return extract_k_mer;
}

uint32_t** set_rotation_table(uint8_t** repeat_check_table) {
    uint32_t** rot_table = (uint32_t**) malloc(sizeof(uint32_t*) * (MIN(MAX_MER, TABLE_MAX_MER) - ABS_MIN_MER + 1));
    if (rot_table == nullptr) {
        fprintf(stderr, "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    uint32_t max_ind;
    for (int i = 0; i <= MIN(MAX_MER, TABLE_MAX_MER) - ABS_MIN_MER; i++) {
        max_ind = 1 << (2 * (i + ABS_MIN_MER));

        rot_table[i] = (uint32_t*) calloc(max_ind, sizeof(uint32_t));
        if (rot_table[i] == nullptr)  {
            fprintf(stderr, "memory allocation failure\n");
            exit(EXIT_FAILURE);
        }

        for (uint32_t j = 0; j < max_ind; j++) {
            if (rot_table[i][j] == 0) {
                fill_rotation_table(rot_table[i], i + ABS_MIN_MER, j, repeat_check_table);
            }
        }
    }

    return rot_table;
}

uint16_t** set_k_mer_counter() {
    if (MIN_MER <= TABLE_MAX_MER) {
        uint16_t** cnt = (uint16_t**) malloc(sizeof(uint16_t*) * (MIN(MAX_MER, TABLE_MAX_MER) - MIN_MER + 1));
        if (cnt == nullptr) {
            fprintf(stderr, "memory allocation failure\n");
            exit(EXIT_FAILURE);
        }

        uint32_t max_ind;
        for (int i = 0; i <= MIN(MAX_MER, TABLE_MAX_MER) - MIN_MER; i++) {
            max_ind =  1 << (2 * (i + MIN_MER));

            cnt[i] = (uint16_t*) calloc(max_ind, sizeof(uint16_t));
            if (cnt[i] == nullptr)  {
                fprintf(stderr, "memory allocation failure\n");
                exit(EXIT_FAILURE);
            }
        }

        return cnt;
    } else {
        return nullptr;
    }
}

uint64_t** set_k_mer_data() {
    uint64_t** k_mer_data = (uint64_t**) malloc(sizeof(uint64_t*) * (MAX_MER - MIN_MER + 1));
    if (k_mer_data == nullptr) {
        fprintf(stderr, "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < MAX_MER - MIN_MER + 1; i++) {
        k_mer_data[i] = (uint64_t*) calloc(3, sizeof(uint64_t));
        if (k_mer_data[i] == nullptr)  {
            fprintf(stderr, "memory allocation failure\n");
            exit(EXIT_FAILURE);
        }
    }

    return k_mer_data;
}

uint128_t** set_k_mer_data_128() {
    uint128_t** k_mer_data = (uint128_t**) malloc(sizeof(uint128_t*) * (MAX_MER - MIN_MER + 1));
    if (k_mer_data == nullptr) {
        fprintf(stderr, "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < MAX_MER - MIN_MER + 1; i++) {
        k_mer_data[i] = (uint128_t*) calloc(3, sizeof(uint128_t));
        if (k_mer_data[i] == nullptr)  {
            fprintf(stderr, "memory allocation failure\n");
            exit(EXIT_FAILURE);
        }
    }

    return k_mer_data;
}

uint32_t** set_k_mer_counter_list() {
    if (MIN_MER <= TABLE_MAX_MER) {
        uint32_t** k_mer_counter_list = (uint32_t**) malloc(sizeof(uint32_t*) * (MIN(MAX_MER, TABLE_MAX_MER) - MIN_MER + 1));
        if (k_mer_counter_list == nullptr) {
            fprintf(stderr, "memory allocation failure\n");
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < MIN(MAX_MER, TABLE_MAX_MER) - MIN_MER + 1; i++) {
            k_mer_counter_list[i] = (uint32_t*) calloc(MAX_SEQ, sizeof(uint32_t));
            if (k_mer_counter_list[i] == nullptr)  {
                fprintf(stderr, "memory allocation failure\n");
                exit(EXIT_FAILURE);
            }
        }
        return k_mer_counter_list;
    } else {
        return nullptr;
    }
}

void fill_rotation_table(uint32_t* rot, int k, uint32_t seq, uint8_t** repeat_check_table) {
    uint32_t code;
    uint32_t tmp = seq;
    uint32_t cnt = 0;
    uint16_t nut_cnt[4]={0, 0, 0, 0};
    rot[tmp] = seq;
    for (int i = 0; i < k - 1; i++) {
        code = tmp & 0x3;
        nut_cnt[code]++;

        tmp = ((tmp & 0x3) << (2 * (k - 1))) + (tmp >> 2);
        rot[tmp] = seq;
    }
    nut_cnt[tmp & 0x3]++;

    for (uint16_t i : nut_cnt) {
        if (i > 0) {
            cnt++;
        }
    }

    repeat_check_table[k - ABS_MIN_MER][seq] = cnt <= ABS_MIN_DNA_COUNT ? 1 : 0;
}

uint64_t get_rot_seq(uint64_t seq, int k) {
    uint64_t tmp = seq;
    uint64_t ans = seq;
    for (int i = 0; i < k - 1; i++) {
        tmp = ((tmp & 0x3) << (2 * (k - 1))) + (tmp >> 2);
        ans = MIN(tmp, ans);
    }
    return ans;
}

uint128_t get_rot_seq_128(const uint128_t& seq, int k) {
    uint128_t tmp = seq;
    uint128_t ans = seq;
    for (int i = 0; i < k - 1; i++) {
        tmp = ((tmp & 0x3) << (2 * (k - 1))) + (tmp >> 2);
        ans = MIN(tmp, ans);
    }
    return ans;
}

uint8_t get_repeat_check(uint64_t seq, int k) {
    uint32_t cnt = 0;
    uint16_t nut_cnt[4]={0, 0, 0, 0};

    for (int i = 0; i < k; i++, seq >>= 2) {
        nut_cnt[seq & 0x3]++;
    }

    for (uint16_t i : nut_cnt) {
        if (i > 0) {
            cnt++;
        }
    }

    return cnt <= ABS_MIN_DNA_COUNT ? 1 : 0;
}

uint8_t get_repeat_check(uint128_t seq, int k) {
    uint32_t cnt = 0;
    uint16_t nut_cnt[4]={0, 0, 0, 0};

    for (int i = 0; i < k; i++, seq >>= 2) {
        nut_cnt[(uint8_t) (seq & 0x3)]++;
    }

    for (uint16_t i : nut_cnt) {
        if (i > 0) {
            cnt++;
        }
    }

    return cnt <= ABS_MIN_DNA_COUNT ? 1 : 0;
}

int get_dna_count(uint128_t seq, int k) {
    int cnt = 0;
    uint16_t nut_cnt[4]={0, 0, 0, 0};

    for (int i = 0; i < k; i++, seq >>= 2) {
        nut_cnt[(uint8_t) (seq & 0x3)]++;
    }

    for (uint16_t i : nut_cnt) {
        if (i > 0) {
            cnt++;
        }
    }

    return cnt;
}

void int_to_four(char* buffer, uint128_t seq, int n) {
    for (int i = 0; i < n; i++) {
        buffer[n - 1 - i] = trans_arr[(uint8_t) (seq & 0x3)];
        seq >>= 2;
    }
    buffer[n] = '\0';
}

void k_mer_target(const char* seq, int st, int nd, uint32_t** rot_table, const uint64_t *extract_k_mer,
                  uint16_t** k_mer_counter, CounterMap* k_mer_counter_map,
                  uint64_t** k_mer_data, uint32_t** k_mer_counter_list, uint8_t** repeat_check_table,
                  ResultMapPair result_pair,
                  int16_t* k_mer_total_cnt,
                  int k) {

    memset(k_mer_total_cnt, 0, sizeof(int16_t) * (MAX_MER - MIN_MER + 2));

    int len;
    int chr_code;
    int err_loc;
    uint16_t cur_cnt;
    uint32_t estimated_count;
    uint64_t tmp;
    uint64_t val;
    uint32_t _a, _b;
    uint64_t _c;

    double baseline = result_pair.first == nullptr ? LOW_BASELINE : HIGH_BASELINE;
    double target_frequency = 0;
    int target_k;
    double frequency;

    err_loc = st - 1;
    for (int i = st ; i <= nd; i++) {
        chr_code = codes[seq[i]];
        if (chr_code >= 0) {
            len = MIN(i - err_loc, MAX_MER);
            if (len >= MIN_MER) {
                k_mer_total_cnt[0]++;
                k_mer_total_cnt[len - MIN_MER + 1]--;
            }
        } else {
            err_loc = i;
        }
    }

    for (int i = 1; i <= MAX_MER - MIN_MER; i++) {
        k_mer_total_cnt[i] = (int16_t)(k_mer_total_cnt[i] + k_mer_total_cnt[i - 1]);
    }

    tmp = 0; err_loc = st - 1;
    for (int i = st ; i <= nd; i++) {
        tmp <<= 2;
        chr_code = codes[seq[i]];
        if (chr_code >= 0) {
            tmp += chr_code;
            if (i - err_loc >= k) {
                if (k <= TABLE_MAX_MER) {
                    val = rot_table[k - ABS_MIN_MER][tmp & extract_k_mer[k - MIN_MER]];
                    cur_cnt = ++k_mer_counter[k - MIN_MER][val];

                    k_mer_counter_list[k - MIN_MER][k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]++] = val;
                } else {
                    val = get_rot_seq(tmp & extract_k_mer[k - MIN_MER], k);
                    cur_cnt = ++k_mer_counter_map[k - TABLE_MAX_MER - 1][val];
                    k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]++;
                }

                if (k_mer_data[k - MIN_MER][K_MER_DATA_MAX] < cur_cnt) {
                    k_mer_data[k - MIN_MER][K_MER_DATA_MAX] = cur_cnt;
                    k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ] = val;
                }

                estimated_count = k_mer_data[k - MIN_MER][K_MER_DATA_MAX] + k_mer_total_cnt[k - MIN_MER] - k_mer_data[k - MIN_MER][K_MER_DATA_COUNT];
                if (estimated_count < k_mer_total_cnt[k - MIN_MER] * baseline) {
                    break;
                }
            }
        } else {
            err_loc = i;
        }
    }

    target_k = k;
    target_frequency = (k <= TABLE_MAX_MER ? repeat_check_table[k - ABS_MIN_MER][k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ]] : get_repeat_check(k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ], k)) ? 0 : (double) k_mer_data[k - MIN_MER][K_MER_DATA_MAX] / (double) k_mer_data[k - MIN_MER][K_MER_DATA_COUNT];

    if (target_frequency >= baseline) {
        ResultMap* result = result_pair.first == nullptr ? result_pair.second : result_pair.first;
        if (target_k == k) {
            if (k <= TABLE_MAX_MER) {
                for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                    val = k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]];
                    if (val > 0) {
                        _a = rot_table[k - ABS_MIN_MER][reverse_complement_32(k_mer_counter_list[k - MIN_MER][i]) >> (2 * (16 - k))];
                        _b = k_mer_counter_list[k - MIN_MER][i];
                        (*result)[KmerSeq {k, MIN(_a, _b)}] += val;
                        k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                    }
                }
            } else {
                for (auto& [key, value] : k_mer_counter_map[k - TABLE_MAX_MER - 1]) {
                    _c = get_rot_seq(reverse_complement_64(key) >> (2 * (32 - k)), k);
                    (*result)[KmerSeq {k, MIN(key, _c)}] += value;
                }
                k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
            }
        } else {
            if (k <= TABLE_MAX_MER) {
                for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                    k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                }
            }
            else {
                k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
            }
        }
        k_mer_data[k - MIN_MER][K_MER_DATA_COUNT] = 0;
        k_mer_data[k - MIN_MER][K_MER_DATA_MAX] = 0;
    }
    else {
        if (k <= TABLE_MAX_MER) {
            for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
            }
        }
        else {
            k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
        }
        k_mer_data[k - MIN_MER][K_MER_DATA_COUNT] = 0;
        k_mer_data[k - MIN_MER][K_MER_DATA_MAX] = 0;
    }
}

void k_mer_target_128(const char* seq, int st, int nd, uint32_t** rot_table, const uint128_t *extract_k_mer,
                      uint16_t** k_mer_counter, CounterMap_128* k_mer_counter_map,
                      uint128_t** k_mer_data, uint32_t** k_mer_counter_list, uint8_t** repeat_check_table,
                      ResultMapPair result_pair,
                      int16_t* k_mer_total_cnt,
                      int k) {

    memset(k_mer_total_cnt, 0, sizeof(int16_t) * (MAX_MER - MIN_MER + 2));

    int len;
    int chr_code;
    int err_loc;
    uint16_t cur_cnt;
    uint32_t estimated_count;
    uint128_t tmp;
    uint128_t val;
    uint32_t _a, _b;
    uint128_t _c;

    double baseline = result_pair.first == nullptr ? LOW_BASELINE : HIGH_BASELINE;
    double target_frequency = 0;
    int target_k;
    double frequency;

    err_loc = st - 1;
    for (int i = st ; i <= nd; i++) {
        chr_code = codes[seq[i]];
        if (chr_code >= 0) {
            len = MIN(i - err_loc, MAX_MER);
            if (len >= MIN_MER) {
                k_mer_total_cnt[0]++;
                k_mer_total_cnt[len - MIN_MER + 1]--;
            }
        } else {
            err_loc = i;
        }
    }

    for (int i = 1; i <= MAX_MER - MIN_MER; i++) {
        k_mer_total_cnt[i] = (int16_t)(k_mer_total_cnt[i] + k_mer_total_cnt[i - 1]);
    }

    tmp = 0; err_loc = st - 1;
    for (int i = st ; i <= nd; i++) {
        tmp <<= 2;
        chr_code = codes[seq[i]];
        if (chr_code >= 0) {
            tmp += chr_code;
            if (i - err_loc >= k) {
                if (k <= TABLE_MAX_MER) {
                    val = (uint64_t) rot_table[k - ABS_MIN_MER][(uint32_t) (tmp & extract_k_mer[k - MIN_MER])];
                    cur_cnt = ++k_mer_counter[k - MIN_MER][(uint32_t) val];

                    k_mer_counter_list[k - MIN_MER][(uint32_t) (k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]++)] = (uint32_t) val;
                } else {
                    val = get_rot_seq_128(tmp & extract_k_mer[k - MIN_MER], k);
                    cur_cnt = ++k_mer_counter_map[k - TABLE_MAX_MER - 1][val];
                    k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]++;
                }

                if (k_mer_data[k - MIN_MER][K_MER_DATA_MAX] < cur_cnt) {
                    k_mer_data[k - MIN_MER][K_MER_DATA_MAX] = cur_cnt;
                    k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ] = val;
                }

                estimated_count = (uint32_t) k_mer_data[k - MIN_MER][K_MER_DATA_MAX] + k_mer_total_cnt[k - MIN_MER] - (uint32_t) k_mer_data[k - MIN_MER][K_MER_DATA_COUNT];
                if (estimated_count < k_mer_total_cnt[k - MIN_MER] * LOW_BASELINE) {
                    break;
                }
            }
        } else {
            err_loc = i;
        }
    }

    target_k = k;
    target_frequency = (k <= TABLE_MAX_MER ? repeat_check_table[k - ABS_MIN_MER][(uint32_t) k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ]] : get_repeat_check(k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ], k)) ? 0 : (double) k_mer_data[k - MIN_MER][K_MER_DATA_MAX] / (double) k_mer_data[k - MIN_MER][K_MER_DATA_COUNT];

    if (target_frequency >= baseline) {
        ResultMap* result = result_pair.first == nullptr ? result_pair.second : result_pair.first;
        if (target_k == k) {
            if (k <= TABLE_MAX_MER) {
                for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                    val = k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]];
                    if (val > 0) {
                        _a = rot_table[k - ABS_MIN_MER][reverse_complement_32(k_mer_counter_list[k - MIN_MER][i]) >> (2 * (16 - k))];
                        _b = k_mer_counter_list[k - MIN_MER][i];
                        (*result)[KmerSeq {k, MIN(_a, _b)}] += (uint32_t) val;
                        k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                    }
                }
            } else {
                for (auto& [key, value] : k_mer_counter_map[k - TABLE_MAX_MER - 1]) {
                    _c = get_rot_seq_128(reverse_complement_128(key) >> (2 * (64 - k)), k);
                    (*result)[KmerSeq {k, MIN(key, _c)}] += value;
                }
                k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
            }
        } else {
            if (k <= TABLE_MAX_MER) {
                for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                    k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                }
            }
            else {
                k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
            }
        }
        k_mer_data[k - MIN_MER][K_MER_DATA_COUNT] = 0;
        k_mer_data[k - MIN_MER][K_MER_DATA_MAX] = 0;
    }
    else {
        if (k <= TABLE_MAX_MER) {
            for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
            }
        }
        else {
            k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
        }
        k_mer_data[k - MIN_MER][K_MER_DATA_COUNT] = 0;
        k_mer_data[k - MIN_MER][K_MER_DATA_MAX] = 0;
    }
}

KmerData k_mer_check(const char* seq, int st, int nd, uint32_t** rot_table, const uint64_t *extract_k_mer,
                     uint16_t** k_mer_counter, CounterMap* k_mer_counter_map,
                     uint64_t** k_mer_data, uint32_t** k_mer_counter_list, uint8_t** repeat_check_table,
                     ResultMapPair result_pair,
                     int16_t* k_mer_total_cnt, int min_mer, int max_mer, std::pair<uint64_t, uint64_t>* repeat_seq) {

    memset(k_mer_total_cnt, 0, sizeof(int16_t) * (max_mer - min_mer + 2));

    int len;
    int chr_code;
    int err_loc;
    uint16_t cur_cnt;
    uint32_t estimated_count;
    uint64_t tmp;
    uint64_t val;

    double target_frequency_low = 0;
    double target_frequency_high = 0;
    int target_k_low = 0, target_k_high = 0;
    double frequency;

    err_loc = st - 1;
    for (int i = st ; i <= nd; i++) {
        chr_code = codes[seq[i]];
        if (chr_code >= 0) {
            len = MIN(i - err_loc, max_mer);
            if (len >= min_mer) {
                k_mer_total_cnt[0]++;
                k_mer_total_cnt[len - min_mer + 1]--;
            }
        } else {
            err_loc = i;
        }
    }

    for (int i = 1; i <= max_mer - min_mer; i++) {
        k_mer_total_cnt[i] = (int16_t)(k_mer_total_cnt[i] + k_mer_total_cnt[i - 1]);
    }

    for (int k = min_mer; k <= max_mer; k++) {
        tmp = 0; err_loc = st - 1;
        for (int i = st ; i <= nd; i++) {
            tmp <<= 2;
            chr_code = codes[seq[i]];
            if (chr_code >= 0) {
                tmp += chr_code;
                if (i - err_loc >= k) {
                    if (k <= TABLE_MAX_MER) {
                        val = rot_table[k - ABS_MIN_MER][tmp & extract_k_mer[k - MIN_MER]];
                        cur_cnt = ++k_mer_counter[k - MIN_MER][val];

                        k_mer_counter_list[k - MIN_MER][k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]++] = val;
                    } else {
                        val = get_rot_seq(tmp & extract_k_mer[k - MIN_MER], k);
                        cur_cnt = ++k_mer_counter_map[k - TABLE_MAX_MER - 1][val];
                        k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]++;
                    }

                    if (k_mer_data[k - MIN_MER][K_MER_DATA_MAX] < cur_cnt) {
                        k_mer_data[k - MIN_MER][K_MER_DATA_MAX] = cur_cnt;
                        k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ] = val;
                    }

                    estimated_count = k_mer_data[k - MIN_MER][K_MER_DATA_MAX] + k_mer_total_cnt[k - min_mer] - k_mer_data[k - MIN_MER][K_MER_DATA_COUNT];
                    if (estimated_count < k_mer_total_cnt[k - min_mer] * LOW_BASELINE) {
                        break;
                    }
                }
            } else {
                err_loc = i;
            }
        }
    }

    bool repeat_k;
    std::vector<int> target_k_vector;

    for (int k = min_mer; k <= max_mer; k++) {
        repeat_k = false;
        frequency = (double) k_mer_data[k - MIN_MER][K_MER_DATA_MAX] / (double) k_mer_data[k - MIN_MER][K_MER_DATA_COUNT];
        if (frequency >= MAX(LOW_BASELINE, target_frequency_low) && !(k <= TABLE_MAX_MER ? repeat_check_table[k - ABS_MIN_MER][k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ]] : get_repeat_check(k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ], k))) {
            for (auto& tk : target_k_vector) {
                if (k % tk == 0) {
                    repeat_k = true;
                    break;
                }
            }

            if (!repeat_k) {
                target_k_low = k;
                target_frequency_low = frequency;
                target_k_vector.push_back(k);
            }
        }
    }

    target_k_vector.clear();
    for (int k = min_mer; k <= max_mer; k++) {
        repeat_k = false;
        frequency = (double) k_mer_data[k - MIN_MER][K_MER_DATA_MAX] / (double) k_mer_data[k - MIN_MER][K_MER_DATA_COUNT];
        if (frequency >= MAX(HIGH_BASELINE, target_frequency_high) && !(k <= TABLE_MAX_MER ? repeat_check_table[k - ABS_MIN_MER][k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ]] : get_repeat_check(k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ], k))) {
            for (auto& tk : target_k_vector) {
                if (k % tk == 0) {
                    repeat_k = true;
                    break;
                }
            }

            if (!repeat_k) {
                target_k_high = k;
                target_frequency_high = frequency;
                target_k_vector.push_back(k);
            }
        }
    }

    if (repeat_seq != nullptr) {
        *repeat_seq = {target_k_high == 0 ? 0 : k_mer_data[target_k_high - MIN_MER][K_MER_DATA_MAX_SEQ], target_k_low == 0 ? 0 : k_mer_data[target_k_low - MIN_MER][K_MER_DATA_MAX_SEQ]};
    }

    if (target_k_low > 0 || target_k_high > 0) {
        ResultMap* result;
        for (int k = min_mer; k <= max_mer; k++) {
            if ((target_k_low == k && result_pair.second != nullptr) && (target_k_high == k && result_pair.first != nullptr)) {
                if (k <= TABLE_MAX_MER) {
                    for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                        val = k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]];
                        if (val > 0) {
                            (*result_pair.first)[KmerSeq {k, k_mer_counter_list[k - MIN_MER][i]}] += (uint32_t) val;
                            (*result_pair.second)[KmerSeq {k, k_mer_counter_list[k - MIN_MER][i]}] += (uint32_t) val;
                            k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                        }
                    }
                } else {
                    for (auto& [key, value] : k_mer_counter_map[k - TABLE_MAX_MER - 1]) {
                        (*result_pair.first)[KmerSeq {k, key}] += value;
                        (*result_pair.second)[KmerSeq {k, key}] += value;
                    }
                    k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
                }
            } else if (target_k_high == k && result_pair.first != nullptr) {
                if (k <= TABLE_MAX_MER) {
                    for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                        val = k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]];
                        if (val > 0) {
                            (*result_pair.first)[KmerSeq {k, k_mer_counter_list[k - MIN_MER][i]}] += (uint32_t) val;
                            k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                        }
                    }
                } else {
                    for (auto& [key, value] : k_mer_counter_map[k - TABLE_MAX_MER - 1]) {
                        (*result_pair.first)[KmerSeq {k, key}] += value;
                    }
                    k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
                }
            } else if (target_k_low == k && result_pair.second != nullptr) {
                if (k <= TABLE_MAX_MER) {
                    for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                        val = k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]];
                        if (val > 0) {
                            (*result_pair.second)[KmerSeq {k, k_mer_counter_list[k - MIN_MER][i]}] += (uint32_t) val;
                            k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                        }
                    }
                } else {
                    for (auto& [key, value] : k_mer_counter_map[k - TABLE_MAX_MER - 1]) {
                        (*result_pair.second)[KmerSeq {k, key}] += value;
                    }
                    k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
                }
            } else {
                if (k <= TABLE_MAX_MER) {
                    for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                        k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                    }
                }
                else {
                    k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
                }
            }
            k_mer_data[k - MIN_MER][K_MER_DATA_COUNT] = 0;
            k_mer_data[k - MIN_MER][K_MER_DATA_MAX] = 0;
        }
        return {target_k_high, target_k_low};
    }
    else {
        for (int k = min_mer; k <= max_mer; k++) {
            if (k <= TABLE_MAX_MER) {
                for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                    k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                }
            }
            else {
                k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
            }
            k_mer_data[k - MIN_MER][K_MER_DATA_COUNT] = 0;
            k_mer_data[k - MIN_MER][K_MER_DATA_MAX] = 0;
        }
        return {0, 0};
    }
}

KmerData k_mer_check_128(const char* seq, int st, int nd, uint32_t** rot_table, const uint128_t *extract_k_mer,
                         uint16_t** k_mer_counter, CounterMap_128* k_mer_counter_map,
                         uint128_t** k_mer_data, uint32_t** k_mer_counter_list, uint8_t** repeat_check_table,
                         ResultMapPair result_pair,
                         int16_t* k_mer_total_cnt, int min_mer, int max_mer, std::pair<uint128_t, uint128_t>* repeat_seq) {

    memset(k_mer_total_cnt, 0, sizeof(int16_t) * (max_mer - min_mer + 2));

    int len;
    int chr_code;
    int err_loc;
    uint16_t cur_cnt;
    uint32_t estimated_count;
    uint128_t tmp;
    uint128_t val;

    double target_frequency_low = 0;
    double target_frequency_high = 0;
    int target_k_low = 0, target_k_high = 0;
    double frequency;

    err_loc = st - 1;
    for (int i = st ; i <= nd; i++) {
        chr_code = codes[seq[i]];
        if (chr_code >= 0) {
            len = MIN(i - err_loc, max_mer);
            if (len >= min_mer) {
                k_mer_total_cnt[0]++;
                k_mer_total_cnt[len - min_mer + 1]--;
            }
        } else {
            err_loc = i;
        }
    }

    for (int i = 1; i <= max_mer - min_mer; i++) {
        k_mer_total_cnt[i] = (int16_t)(k_mer_total_cnt[i] + k_mer_total_cnt[i - 1]);
    }

    for (int k = min_mer; k <= max_mer; k++) {
        tmp = 0; err_loc = st - 1;
        for (int i = st ; i <= nd; i++) {
            tmp <<= 2;
            chr_code = codes[seq[i]];
            if (chr_code >= 0) {
                tmp += chr_code;
                if (i - err_loc >= k) {
                    if (k <= TABLE_MAX_MER) {
                        val = (uint64_t) rot_table[k - ABS_MIN_MER][(uint32_t) (tmp & extract_k_mer[k - MIN_MER])];
                        cur_cnt = ++k_mer_counter[k - MIN_MER][(uint32_t) val];

                        k_mer_counter_list[k - MIN_MER][(uint32_t) (k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]++)] = (uint32_t) val;
                    } else {
                        val = get_rot_seq_128(tmp & extract_k_mer[k - MIN_MER], k);
                        cur_cnt = ++k_mer_counter_map[k - TABLE_MAX_MER - 1][val];
                        k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]++;
                        // std::cout << k << "!!" << k_mer_data[k - MIN_MER][K_MER_DATA_COUNT] << std::endl;
                    }

                    if (k_mer_data[k - MIN_MER][K_MER_DATA_MAX] < cur_cnt) {
                        k_mer_data[k - MIN_MER][K_MER_DATA_MAX] = cur_cnt;
                        k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ] = val;
                    }

                    estimated_count = (uint32_t) k_mer_data[k - MIN_MER][K_MER_DATA_MAX] + k_mer_total_cnt[k - min_mer] - (uint32_t) k_mer_data[k - MIN_MER][K_MER_DATA_COUNT];
                    if (estimated_count < k_mer_total_cnt[k - min_mer] * LOW_BASELINE) {
                        break;
                    }
                }
            } else {
                err_loc = i;
            }
        }
    }

    bool repeat_k;
    std::vector<int> target_k_vector;

    for (int k = min_mer; k <= max_mer; k++) {
        repeat_k = false;
        frequency = (double) k_mer_data[k - MIN_MER][K_MER_DATA_MAX] / (double) k_mer_data[k - MIN_MER][K_MER_DATA_COUNT];
        if (frequency >= MAX(LOW_BASELINE, target_frequency_low) && !(k <= TABLE_MAX_MER ? repeat_check_table[k - ABS_MIN_MER][(uint32_t) k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ]] : get_repeat_check(k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ], k))) {
            for (auto& tk : target_k_vector) {
                if (k % tk == 0) {
                    repeat_k = true;
                    break;
                }
            }

            if (!repeat_k) {
                target_k_low = k;
                target_frequency_low = frequency;
                target_k_vector.push_back(k);
            }
        }
    }

    target_k_vector.clear();
    for (int k = min_mer; k <= max_mer; k++) {
        repeat_k = false;
        frequency = (double) k_mer_data[k - MIN_MER][K_MER_DATA_MAX] / (double) k_mer_data[k - MIN_MER][K_MER_DATA_COUNT];
        if (frequency >= MAX(HIGH_BASELINE, target_frequency_high) && !(k <= TABLE_MAX_MER ? repeat_check_table[k - ABS_MIN_MER][(uint32_t) k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ]] : get_repeat_check(k_mer_data[k - MIN_MER][K_MER_DATA_MAX_SEQ], k))) {
            for (auto& tk : target_k_vector) {
                if (k % tk == 0) {
                    repeat_k = true;
                    break;
                }
            }

            if (!repeat_k) {
                target_k_high = k;
                target_frequency_high = frequency;
                target_k_vector.push_back(k);
            }
        }
    }

    if (repeat_seq != nullptr) {
        *repeat_seq = {target_k_high == 0 ? 0 : k_mer_data[target_k_high - MIN_MER][K_MER_DATA_MAX_SEQ], target_k_low == 0 ? 0 : k_mer_data[target_k_low - MIN_MER][K_MER_DATA_MAX_SEQ]};
    }

    if (target_k_low > 0 || target_k_high > 0) {
        ResultMap* result;
        for (int k = min_mer; k <= max_mer; k++) {
            if ((target_k_low == k && result_pair.second != nullptr) && (target_k_high == k && result_pair.first != nullptr)) {
                if (k <= TABLE_MAX_MER) {
                    for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                        val = k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]];
                        if (val > 0) {
                            (*result_pair.first)[KmerSeq {k, k_mer_counter_list[k - MIN_MER][i]}] += (uint32_t) val;
                            (*result_pair.second)[KmerSeq {k, k_mer_counter_list[k - MIN_MER][i]}] += (uint32_t) val;
                            k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                        }
                    }
                } else {
                    for (auto& [key, value] : k_mer_counter_map[k - TABLE_MAX_MER - 1]) {
                        (*result_pair.first)[KmerSeq {k, key}] += value;
                        (*result_pair.second)[KmerSeq {k, key}] += value;
                    }
                    k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
                }
            } else if (target_k_high == k && result_pair.first != nullptr) {
                if (k <= TABLE_MAX_MER) {
                    for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                        val = k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]];
                        if (val > 0) {
                            (*result_pair.first)[KmerSeq {k, k_mer_counter_list[k - MIN_MER][i]}] += (uint32_t) val;
                            k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                        }
                    }
                } else {
                    for (auto& [key, value] : k_mer_counter_map[k - TABLE_MAX_MER - 1]) {
                        (*result_pair.first)[KmerSeq {k, key}] += value;
                    }
                    k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
                }
            } else if (target_k_low == k && result_pair.second != nullptr) {
                if (k <= TABLE_MAX_MER) {
                    for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                        val = k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]];
                        if (val > 0) {
                            (*result_pair.second)[KmerSeq {k, k_mer_counter_list[k - MIN_MER][i]}] += (uint32_t) val;
                            k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                        }
                    }
                } else {
                    for (auto& [key, value] : k_mer_counter_map[k - TABLE_MAX_MER - 1]) {
                        (*result_pair.second)[KmerSeq {k, key}] += value;
                    }
                    k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
                }
            } else {
                if (k <= TABLE_MAX_MER) {
                    for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                        k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                    }
                }
                else {
                    k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
                }
            }
            k_mer_data[k - MIN_MER][K_MER_DATA_COUNT] = 0;
            k_mer_data[k - MIN_MER][K_MER_DATA_MAX] = 0;
        }
        return {target_k_high, target_k_low};
    }
    else {
        for (int k = min_mer; k <= max_mer; k++) {
            if (k <= TABLE_MAX_MER) {
                for (uint32_t i = 0; i < k_mer_data[k - MIN_MER][K_MER_DATA_COUNT]; i++) {
                    k_mer_counter[k - MIN_MER][k_mer_counter_list[k - MIN_MER][i]] = 0;
                }
            }
            else {
                k_mer_counter_map[k - TABLE_MAX_MER - 1].clear();
            }
            k_mer_data[k - MIN_MER][K_MER_DATA_COUNT] = 0;
            k_mer_data[k - MIN_MER][K_MER_DATA_MAX] = 0;
        }
        return {0, 0};
    }
}

bool check_ans_seq(const KmerSeq& seq, uint128_t* extract_k_mer_ans, uint32_t** rot_table) {
    int i;
    uint128_t num_seq;
    uint128_t temp_seq, temp_bef_seq;
    for (int k = ABS_MIN_MER; k < MIN_MER; k++) {
        num_seq = seq.second;
        for (i = 0; i < seq.first - k + 1; i++) {
            temp_seq = k <= TABLE_MAX_MER ? rot_table[k - ABS_MIN_MER][(uint32_t) (num_seq & extract_k_mer_ans[k - ABS_MIN_MER])] : get_rot_seq_128(num_seq & extract_k_mer_ans[k - ABS_MIN_MER], k);
            if (i > 0 && temp_seq != temp_bef_seq) {
                break;
            }
            temp_bef_seq = temp_seq;
            num_seq >>= 2;
        }

        if (i == seq.first - k + 1) {
            return false;
        }
    }
    return true;
}

void final_process_output(FinalFastqData* total_result_high, FinalFastqData* total_result_low) {
    TRMDirMap trm_map = {};

    bool max_cnt_check = false;
    for (auto& [k, v] : *total_result_high) {
        if (v.forward + v.backward + v.both >= ABS_MIN_ANS_COUNT) {
            max_cnt_check = true;
            break;
        }
    }
    for (auto& [k, v] : *total_result_low) {
        if (v.forward + v.backward + v.both >= ABS_MIN_ANS_COUNT) {
            max_cnt_check = true;
            break;
        }
    }

    fprintf(stdout, ">Putative_TRM\n");
    if (max_cnt_check) {
        ResultMap* score_result_map_high = get_score_map(total_result_high);
        ResultMap* score_result_map = get_score_map(total_result_low);

        for (auto& [k, v] : *score_result_map_high) {
            (*score_result_map)[k] += v;
        }
        delete score_result_map_high;

        // Bonus point
        int dna_cnt;
        std::vector<std::pair<KmerSeq, std::pair<uint32_t, int>>> score_result_vector {};
        for (auto&[k, v] : *(score_result_map)) {
            auto low_result = (*total_result_low)[k];
            auto high_result = (*total_result_high)[k];

            int bonus = 0;
            int high_dir, low_dir;

            if (high_result.forward > high_result.backward) {
                high_dir = 1;
            } else if (high_result.forward < high_result.backward) {
                high_dir = -1;
            } else {
                high_dir = 0;
            }

            if (low_result.forward > low_result.backward) {
                low_dir = 1;
            } else if (low_result.forward < low_result.backward) {
                low_dir = -1;
            } else {
                low_dir = 0;
            }

            int final_dir;
            if (low_dir != 0 and (low_dir == high_dir)) {
                bonus += 1;
                final_dir = low_dir;
            } else if (low_dir == 0 and high_dir != 0) {
                final_dir = high_dir;
            } else if (low_dir != 0 and high_dir == 0) {
                final_dir = low_dir;
            } else if (low_dir != high_dir and (low_result.forward > 0 or low_result.backward > 0 or high_result.forward > 0 or high_result.backward > 0)) {
                if (low_result.forward < low_result.backward) {
                    std::swap(low_result.forward, low_result.backward);
                }

                if (high_result.forward < high_result.backward) {
                    std::swap(high_result.forward, high_result.backward);
                }

                if (low_result.backward * high_result.forward == high_result.backward * low_result.forward) {
                    if (low_result.forward + low_result.backward > high_result.forward + high_result.backward) {
                        final_dir = low_dir;
                    } else {
                        final_dir = high_dir;
                    }
                } else if (low_result.backward * high_result.forward < high_result.backward * low_result.forward) {
                    final_dir = low_dir;
                } else {
                    final_dir = high_dir;
                }
            } else {
                final_dir = 0;
            }

            dna_cnt = get_dna_count(k.second, k.first);
            if (dna_cnt > 2) {
                bonus += 1;
            }

            trm_map[k] = final_dir;
            score_result_vector.push_back({k, {v + bonus, dna_cnt}});
        }

        std::sort(score_result_vector.begin(), score_result_vector.end(), [](auto &a, auto &b) {
            if (a.second.first != b.second.first) {
                return a.second.first > b.second.first;
            } else if (a.second.second != b.second.second) {
                return a.second.second > b.second.second;
            } else {
                return a.first.first < b.first.first;
            }
        });

        char buffer[ABS_MAX_MER + 1];
        for (int i = 0; i < MIN(ABS_MAX_ANS_NUM, score_result_vector.size()); i++) {
            int_to_four(buffer, score_result_vector[i].first.second, score_result_vector[i].first.first);

            int dir = trm_map[score_result_vector[i].first];
            char sign = (dir == 1) ? '+' :
                        (dir == -1) ? '-' : '?';
            fprintf(stdout, "%d,%s,%" PRIu32",%c\n", score_result_vector[i].first.first, buffer, score_result_vector[i].second.first, sign);
        }
        delete score_result_map;
    } else {
        fprintf(stdout, "NO_PUTATIVE_TRM,-1\n");
    }

    delete total_result_low;
    delete total_result_high;
}

ResultMap* get_score_map(FinalFastqData* total_result) {
    int cnt;
    FinalFastqVector total_result_vector {};

    for (auto& [k, v] : *total_result) {
        if (v.forward + v.backward + v.both >= ABS_MIN_PRINT_COUNT) {
            if (v.backward > v.forward) {
                total_result_vector.emplace_back(k , FinalData<int64_t>{v.backward, v.forward, v.both});
            } else {
                total_result_vector.emplace_back(k ,v);
            }
        }
    }

    FinalFastqData ratio_result {};
    ResultMap* score_result_map = new ResultMap {};

    std::sort(total_result_vector.begin(), total_result_vector.end(), [](auto &a, auto &b) {
        return a.second.forward > b.second.forward;
    });

    cnt = 0;
    for (auto& [k, v] : total_result_vector) {
        if (v.forward == 0 || cnt >= NUM_RAT_CAND) {
            break;
        }
        if (v.backward >= 0) {
            cnt += 1;
            ratio_result[k] = v;
        }
    }

    for (int i = 0; i < MIN(NUM_FOR_MAX_COUNT, total_result_vector.size()); i++) {
        if (total_result_vector[i].second.forward == 0) {
            break;
        }
        (*score_result_map)[total_result_vector[i].first] += 1;
    }

    std::sort(total_result_vector.begin(), total_result_vector.end(), [](auto &a, auto &b) {
        return a.second.forward + a.second.backward + a.second.both > b.second.forward + b.second.backward + b.second.both;
    });

    cnt = 0;
    for (auto& [k, v] : total_result_vector) {
        if (cnt >= NUM_RAT_CAND) {
            break;
        }
        if (v.forward > 0 && v.backward >= 0) {
            cnt += 1;
            ratio_result[k] = v;
        }
    }

    for (int i = 0; i < MIN(NUM_TOT_MAX_COUNT, total_result_vector.size()); i++) {
        (*score_result_map)[total_result_vector[i].first] += 1;
    }

    FinalFastqVector ratio_result_vector(ratio_result.begin(), ratio_result.end());
    std::sort(ratio_result_vector.begin(), ratio_result_vector.end(), [](auto &a, auto &b) {
        return (double) a.second.backward / a.second.forward < (double) b.second.backward / b.second.forward;
    });

    for (int i = 0; i < MIN(NUM_RAT_MAX_COUNT, ratio_result_vector.size()); i++) {
        (*score_result_map)[ratio_result_vector[i].first] += 1;
    }

    return score_result_map;
}
