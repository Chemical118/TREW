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

ptrdiff_t file_extract(FILE *in, off_t offset, unsigned char *buf, size_t len) {
#ifdef _WIN32
    int ret = _fseeki64(in, offset, SEEK_SET);
#else
    int ret = fseeko(in, offset, SEEK_SET);
#endif
    if (ret != 0) {
        fprintf(stderr, "fseek failed.\n");
        return -1;
    }

    // Read 'len' bytes from the file into the buffer 'buf'
    size_t bytes_read = fread(buf, 1, len, in);
    if (bytes_read != len) {
        if (feof(in)) {
            // End of file reached
            return bytes_read;  // Return the number of bytes read before EOF
        } else {
            // File read error
            fprintf(stderr, "fread failed.\n");
            return -1;
        }
    }

    // Return the number of bytes successfully read
    return bytes_read;
}

// from zran.c
ptrdiff_t deflate_index_extract(FILE *in, gz_index *index,
                                off_t offset, unsigned char *buf, size_t len) {
    if (index == NULL)
        return file_extract(in, offset, buf, len);

    // Do a quick sanity check on the index.
    if (index->have < 1 || index->list[0].out != 0 ||
        index->strm.state == Z_NULL)
        return Z_STREAM_ERROR;

    // If nothing to extract, return zero bytes extracted.
    if (len == 0 || offset < 0 || offset >= index->length)
        return 0;

    // Find the access point closest to but not after offset.
    int lo = -1, hi = index->have;
    point_t *point = index->list;
    while (hi - lo > 1) {
        int mid = (lo + hi) >> 1;
        if (offset < point[mid].out)
            hi = mid;
        else
            lo = mid;
    }
    point += lo;

    // Initialize the input file and prime the inflate engine to start there.
#ifdef _WIN32
  	int ret = _fseeki64(in, point->in - (point->bits ? 1 : 0), SEEK_SET);
#else
  	int ret = fseeko(in, point->in - (point->bits ? 1 : 0), SEEK_SET);
#endif
    if (ret == -1)
        return Z_ERRNO;
    int ch = 0;
    if (point->bits && (ch = getc(in)) == EOF)
        return ferror(in) ? Z_ERRNO : Z_BUF_ERROR;
    index->strm.avail_in = 0;
    ret = inflateReset2(&index->strm, RAW);
    if (ret != Z_OK)
        return ret;
    if (point->bits)
        inflatePrime(&index->strm, point->bits, ch >> (8 - point->bits));
    inflateSetDictionary(&index->strm, point->window, point->dict);

    // Skip uncompressed bytes until offset reached, then satisfy request.
    unsigned char input[CHUNK];
    unsigned char discard[WINSIZE];
    offset -= point->out;       // number of bytes to skip to get to offset
    size_t left = len;          // number of bytes left to read after offset
    do {
        if (offset) {
            // Discard up to offset uncompressed bytes.
            index->strm.avail_out = offset < WINSIZE ? (unsigned)offset :
                                                       WINSIZE;
            index->strm.next_out = discard;
        }
        else {
            // Uncompress up to left bytes into buf.
            index->strm.avail_out = left < (unsigned)-1 ? (unsigned)left :
                                                          (unsigned)-1;
            index->strm.next_out = buf + len - left;
        }

        // Uncompress, setting got to the number of bytes uncompressed.
        if (index->strm.avail_in == 0) {
            // Assure available input.
            index->strm.avail_in = fread(input, 1, CHUNK, in);
            if (index->strm.avail_in < CHUNK && ferror(in)) {
                ret = Z_ERRNO;
                break;
            }
            index->strm.next_in = input;
        }
        unsigned got = index->strm.avail_out;
        ret = inflate(&index->strm, Z_NO_FLUSH);
        got -= index->strm.avail_out;

        // Update the appropriate count.
        if (offset)
            offset -= got;
        else {
            left -= got;
            if (left == 0)
                // Request satisfied.
                break;
        }

        // If we're at the end of a gzip member and there's more to read,
        // continue to the next gzip member.
        if (ret == Z_STREAM_END && index->mode == GZIP) {
            // Discard the gzip trailer.
            unsigned drop = 8;              // length of gzip trailer
            if (index->strm.avail_in >= drop) {
                index->strm.avail_in -= drop;
                index->strm.next_in += drop;
            }
            else {
                // Read and discard the remainder of the gzip trailer.
                drop -= index->strm.avail_in;
                index->strm.avail_in = 0;
                do {
                    if (getc(in) == EOF)
                        // The input does not have a complete trailer.
                        return ferror(in) ? Z_ERRNO : Z_BUF_ERROR;
                } while (--drop);
            }

            if (index->strm.avail_in || ungetc(getc(in), in) != EOF) {
                // There's more after the gzip trailer. Use inflate to skip the
                // gzip header and resume the raw inflate there.
                inflateReset2(&index->strm, GZIP);
                do {
                    if (index->strm.avail_in == 0) {
                        index->strm.avail_in = fread(input, 1, CHUNK, in);
                        if (index->strm.avail_in < CHUNK && ferror(in)) {
                            ret = Z_ERRNO;
                            break;
                        }
                        index->strm.next_in = input;
                    }
                    index->strm.avail_out = WINSIZE;
                    index->strm.next_out = discard;
                    ret = inflate(&index->strm, Z_BLOCK);  // stop after header
                } while (ret == Z_OK && (index->strm.data_type & 0x80) == 0);
                if (ret != Z_OK)
                    break;
                inflateReset2(&index->strm, RAW);
            }
        }

        // Continue until we have the requested data, the deflate data has
        // ended, or an error is encountered.
    } while (ret == Z_OK);

    // Return the number of uncompressed bytes read into buf, or the error.
    return ret == Z_OK || ret == Z_STREAM_END ? len - left : ret;
}

// Add an access point to the list. If out of memory, deallocate the existing
// list and return NULL. index->mode is temporarily the allocated number of
// access points, until it is time for deflate_index_build() to return. Then
// index->mode is set to the mode of inflation.
static gz_index *add_point(gz_index *index, off_t in,
                                       off_t out, off_t beg,
                                       unsigned char *buffer,
                                       unsigned char *bef_buffer,
                                       int buffer_size,
                                       int bef_buffer_size) {
    if (index->have == index->mode) {
        // The list is full. Make it bigger.
        index->mode = index->mode ? index->mode << 1 : 8;
        point_t *next = (point_t *) realloc(index->list, sizeof(point_t) * index->mode);
        if (next == NULL) {
            deflate_index_free(index);
            return NULL;
        }
        index->list = next;
    }

    // Fill in the access point and increment how many we have.
    point_t *next = (point_t *)(index->list) + index->have++;
    if (index->have < 0) {
        // Overflowed the int!
        deflate_index_free(index);
        return NULL;
    }
    next->out = out;
    next->in = in;
    next->bits = index->strm.data_type & 7;
    // next->dict always small then WINSIZE
    next->dict = out - beg > WINSIZE ? WINSIZE : (unsigned)(out - beg);
    next->window = (unsigned char*) malloc(next->dict);
    if (next->window == NULL) {
        deflate_index_free(index);
        return NULL;
    }
    unsigned recent = buffer_size - index->strm.avail_out;

    // enough recent check
    unsigned copy = recent > next->dict ? next->dict : recent;

    // copy nearby winsize
    memcpy(next->window + next->dict - copy, buffer + recent - copy, copy);
    copy = next->dict - copy;
    memcpy(next->window, bef_buffer + bef_buffer_size - copy, copy);

    // Return the index, which may have been newly allocated or destroyed.
    return index;
}

void deflate_index_free(gz_index *index) {
    if (index != NULL) {
        size_t i = index->have;
        while (i)
            free(index->list[--i].window);
        free(index->list);
        inflateEnd(&index->strm);
        free(index);
    }
}

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


FinalData<int64_t> add_data(FinalData<int64_t> a, FinalData<int64_t> b) {
    return FinalData<int64_t> {a.forward + b.forward, a.backward + b.backward, a.both + b.both};
}

ResultMapPairData buffer_task(TBBQueue* task_queue, ThreadData* thread_data, uint32_t** rot_table, const uint64_t *extract_k_mer, const uint128_t *extract_k_mer_128, uint8_t** repeat_check_table) {
    auto [k_mer_counter, k_mer_data, k_mer_data_128, k_mer_counter_list] = thread_data -> init_check();

    int16_t* k_mer_total_cnt = (int16_t*) malloc(sizeof(int16_t) * (MAX_MER - MIN_MER + 2));
    if (k_mer_total_cnt == nullptr) {
        fprintf(stderr, "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    ResultMapData result = ResultMapData {{new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}};
    std::vector<FastqLocData>* loc_vector = new std::vector<FastqLocData> {};

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

                    if (INDEX and left_temp_k_mer.first != right_temp_k_mer.first or left_temp_k_mer.second != right_temp_k_mer.second) {
                        loc_vector -> emplace_back(std::pair<int64_t, int64_t> {temp_task.buf_off + st, n}, right_temp_k_mer, left_temp_k_mer, rht_seq, lef_seq);
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

                    if (INDEX and left_temp_k_mer.first != right_temp_k_mer.first or left_temp_k_mer.second != right_temp_k_mer.second) {
                        loc_vector -> emplace_back(std::pair<int64_t, int64_t> {temp_task.buf_off + st, n}, right_temp_k_mer, left_temp_k_mer, rht_seq, lef_seq);
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
    return {result, loc_vector};
}

ResultMapPairData buffer_task_long(TBBQueue* task_queue, ThreadData* thread_data, uint32_t** rot_table, const uint64_t *extract_k_mer, const uint128_t *extract_k_mer_128, uint8_t** repeat_check_table) {
    auto [k_mer_counter, k_mer_data, k_mer_data_128, k_mer_counter_list] = thread_data -> init_check();

    int16_t* k_mer_total_cnt = (int16_t*) malloc(sizeof(int16_t) * (MAX_MER - MIN_MER + 2));
    if (k_mer_total_cnt == nullptr) {
        fprintf(stderr, "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    ResultMapData result = ResultMapData {{new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}};
    std::vector<FastqLocData>* loc_vector = new std::vector<FastqLocData> {};

    ResultMapPair temp_result_left = {new ResultMap {}, new ResultMap {}};
    ResultMapPair temp_result_right = {new ResultMap {}, new ResultMap {}};

    QueueData temp_task {};

    KmerData k_mer, temp_k_mer, lef_k_mer;
    int tst, tnd, ti, tj;
    KmerData si, sj;

    std::pair<bool, bool> repeat_end;

    int snum;
    int mid;
    int sl, mid_bonus_sl;
    bool is_total;

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
                temp_result_left.first -> clear();
                temp_result_left.second -> clear();

                lef_k_mer = k_mer;
                if (si.first <= snum or si.second <= snum) {
                    sj = {snum, snum};
                    k_mer = {0, 0};
                    repeat_end = {false, false};
                    for (tj = snum; !repeat_end.first || !repeat_end.second; tj--, tnd -= sl) {
                        sl = SLICE_LENGTH + (tj == mid ? mid_bonus_sl : 0);
                        temp_k_mer = k_mer_check(temp_task.buffer, tnd - sl + 1, tnd, rot_table, extract_k_mer, k_mer_counter,
                                                     k_mer_counter_map, k_mer_data, k_mer_counter_list, repeat_check_table,
                                                     {repeat_end.first ? nullptr : temp_result_right.first, repeat_end.second ? nullptr : temp_result_right.second}, k_mer_total_cnt, MIN_MER, MAX_MER);

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

                    uint32_t max_cnt = 0;
                    if (si.first <= snum) {
                        max_cnt = 0;
                        lef_seq.first = 0;
                        for (auto& [seq, cnt] : *(temp_result_left.first)) {
                            if (cnt > max_cnt) {
                                lef_seq.first = seq.second;
                                max_cnt = cnt;
                            }
                            (*(result.forward.first))[seq] += cnt;
                        }
                    }
                    if (si.second <= snum) {
                        max_cnt = 0;
                        lef_seq.second = 0;
                        for (auto& [seq, cnt] : *(temp_result_left.second)) {
                            if (cnt > max_cnt) {
                                lef_seq.second = seq.second;
                                max_cnt = cnt;
                            }
                            (*(result.forward.second))[seq] += cnt;
                        }
                    }

                    max_cnt = 0;
                    rht_seq.first = 0;
                    for (auto& [seq, cnt] : *(temp_result_right.first)) {
                        if (cnt > max_cnt) {
                            rht_seq.first = seq.second;
                            max_cnt = cnt;
                        }
                        (*(result.backward.first))[seq] += cnt;
                    }

                    max_cnt = 0;
                    rht_seq.second = 0;
                    for (auto& [seq, cnt] : *(temp_result_right.second)) {
                        if (cnt > max_cnt) {
                            rht_seq.second = seq.second;
                            max_cnt = cnt;
                        }
                        (*(result.backward.second))[seq] += cnt;
                    }

                    temp_result_right.first -> clear();
                    temp_result_right.second -> clear();
                }
                temp_result_left.first -> clear();
                temp_result_left.second -> clear();

                if (INDEX and (si.first <= snum or si.second <= snum) and (k_mer.first != 0 or k_mer.second != 0 or lef_k_mer.first != 0 or lef_k_mer.second != 0)) {
                    loc_vector -> emplace_back(std::pair<int64_t, int64_t> {temp_task.buf_off + (int64_t) st, nd - st + 1}, k_mer, lef_k_mer, rht_seq, lef_seq);
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
                                                     {repeat_end.first ? nullptr : temp_result_right.first, repeat_end.second ? nullptr : temp_result_right.second}, k_mer_total_cnt, MIN_MER, MAX_MER, &rht_seq);

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

                    uint32_t max_cnt = 0;
                    if (si.first <= snum) {
                        max_cnt = 0;
                        lef_seq.first = 0;
                        for (auto& [seq, cnt] : *(temp_result_left.first)) {
                            if (cnt > max_cnt) {
                                lef_seq.first = seq.second;
                                max_cnt = cnt;
                            }
                            (*(result.forward.first))[seq] += cnt;
                        }
                    }
                    if (si.second <= snum) {
                        max_cnt = 0;
                        lef_seq.second = 0;
                        for (auto& [seq, cnt] : *(temp_result_left.second)) {
                            if (cnt > max_cnt) {
                                lef_seq.second = seq.second;
                                max_cnt = cnt;
                            }
                            (*(result.forward.second))[seq] += cnt;
                        }
                    }

                    max_cnt = 0;
                    rht_seq.first = 0;
                    for (auto& [seq, cnt] : *(temp_result_right.first)) {
                        if (cnt > max_cnt) {
                            rht_seq.first = seq.second;
                            max_cnt = cnt;
                        }
                        (*(result.backward.first))[seq] += cnt;
                    }

                    max_cnt = 0;
                    rht_seq.second = 0;
                    for (auto& [seq, cnt] : *(temp_result_right.second)) {
                        if (cnt > max_cnt) {
                            rht_seq.second = seq.second;
                            max_cnt = cnt;
                        }
                        (*(result.backward.second))[seq] += cnt;
                    }

                    temp_result_right.first -> clear();
                    temp_result_right.second -> clear();
                }
                temp_result_left.first -> clear();
                temp_result_left.second -> clear();

                if (INDEX and (si.first != snum + 1 or si.second == snum + 1) and (k_mer.first != 0 or k_mer.second != 0 or lef_k_mer.first != 0 or lef_k_mer.second != 0)) {
                    loc_vector -> emplace_back(std::pair<int64_t, int64_t> {temp_task.buf_off + (int64_t) st, nd - st + 1}, k_mer, lef_k_mer, rht_seq, lef_seq);
                }
            }

            delete temp_task.loc_vector;
            free(temp_task.buffer);
        }

        delete[] k_mer_counter_map;
    }

    delete temp_result_left.first;
    delete temp_result_left.second;
    free(k_mer_total_cnt);
    return {result, loc_vector};
}

void read_fastq_thread(FILE* fp, TBBQueue* buffer_task_queue) {
    int num = 0;
    int shift = 0;
    int bef_shift;
    int idx;

    off_t totout = 0;
    off_t bytes_read;
    off_t avail_out = 0;

    unsigned char* buffer = (unsigned char*) malloc(sizeof(unsigned char) * (LENGTH + 1));
    unsigned char* bef_buffer;
    unsigned char* next_out = buffer;

    int buffer_size = LENGTH;
    QueueData bef_queue_data = {nullptr, nullptr, -1};

    while (!feof(fp)) {
        if (avail_out == 0) {
            avail_out = LENGTH;
            next_out = buffer + shift + buffer_size - LENGTH;
        }

        bytes_read = fread(next_out, 1, avail_out, fp);
        totout += bytes_read;
        avail_out -= bytes_read;

        if ((avail_out == 0 and totout > 0) or feof(fp)) {
            int seq_cnt = 0;
            int bef_num = num;

            LocationVector* loc_vector = new LocationVector{};

            buffer[LENGTH - avail_out + shift] = '\0';
            for (int i = 0; i < LENGTH - avail_out + shift; i++) {
                if (buffer[i] == '\n') {
                    num += 1;
                    if ((num & 3) == 2) {
                        if ((i - 1) - (idx + 1) + 1 > MAX_SEQ) {
                            fprintf(stderr, "This mode is designed for short-read sequencing. Please use 'trew long'.\n");
                            exit (EXIT_FAILURE);
                        }
                        seq_cnt += 1;
                        loc_vector -> emplace_back(idx + 1, i - 1);
                    }
                    idx = i;
                }
            }

            if ((num & 3) == 1 and bef_num == num and seq_cnt == 0) {
                delete loc_vector;

                buffer_size += LENGTH;
                buffer = (unsigned char*) realloc(buffer, buffer_size + LENGTH + 1);
            } else {
                bef_shift = shift;
                if ((num & 3) == 1) {
                    shift = LENGTH - avail_out + shift - idx - 1;
                }
                else {
                    shift = 0;
                }

                bef_buffer = (unsigned char*) malloc(sizeof(unsigned char) * (LENGTH + shift + 1));
                if ((num & 3) == 1) {
                    strcpy(reinterpret_cast<char*>(bef_buffer), reinterpret_cast<char*>(buffer + idx + 1));
                    idx = -1;
                }

                if (bef_queue_data.buffer != nullptr) {
                    buffer_task_queue -> push(bef_queue_data);
                }

                bef_queue_data = QueueData{reinterpret_cast<char*>(buffer), loc_vector,
                                (int64_t) totout - (int64_t) bef_shift - (int64_t) buffer_size + (int64_t) avail_out};


                std::swap(bef_buffer, buffer);
                buffer_size = LENGTH;
            }
        }
    }

    if (bef_queue_data.buffer != NULL) {
        buffer_task_queue -> push(bef_queue_data);
    }
}

void read_fastq_thread_long(FILE* fp, TBBQueue* buffer_task_queue) {
    int num = 0;
    int shift = 0;
    int bef_shift;
    int idx;

    off_t totout = 0;
    off_t bytes_read;
    off_t avail_out = 0;

    unsigned char* buffer = (unsigned char*) malloc(sizeof(unsigned char) * (LENGTH + 1));
    unsigned char* bef_buffer;
    unsigned char* next_out = buffer;

    int buffer_size = LENGTH;
    QueueData bef_queue_data = {nullptr, nullptr, -1};

    while (!feof(fp)) {
        if (avail_out == 0) {
            avail_out = LENGTH;
            next_out = buffer + shift + buffer_size - LENGTH;
        }

        bytes_read = fread(next_out, 1, avail_out, fp);
        totout += bytes_read;
        avail_out -= bytes_read;

        if ((avail_out == 0 and totout > 0) or feof(fp)) {
            int seq_cnt = 0;
            int bef_num = num;

            LocationVector* loc_vector = new LocationVector{};

            buffer[LENGTH - avail_out + shift] = '\0';
            for (int i = 0; i < LENGTH - avail_out + shift; i++) {
                if (buffer[i] == '\n') {
                    num += 1;
                    if ((num & 3) == 2) {
                        seq_cnt += 1;
                        if ((i - 1) - (idx + 1) + 1 >= SLICE_LENGTH) {
                            loc_vector -> emplace_back(idx + 1, i - 1);
                        }
                    }
                    idx = i;
                }
            }

            if ((num & 3) == 1 and bef_num == num and seq_cnt == 0) {
                delete loc_vector;

                buffer_size += LENGTH;
                buffer = (unsigned char*) realloc(buffer, buffer_size + LENGTH + 1);
            } else {
                bef_shift = shift;
                if ((num & 3) == 1) {
                    shift = LENGTH - avail_out + shift - idx - 1;
                }
                else {
                    shift = 0;
                }

                bef_buffer = (unsigned char*) malloc(sizeof(unsigned char) * (LENGTH + shift + 1));
                if ((num & 3) == 1) {
                    strcpy(reinterpret_cast<char*>(bef_buffer), reinterpret_cast<char*>(buffer + idx + 1));
                    idx = -1;
                }

                if (bef_queue_data.buffer != nullptr) {
                    buffer_task_queue -> push(bef_queue_data);
                }

                bef_queue_data = QueueData{reinterpret_cast<char*>(buffer), loc_vector,
                                (int64_t) totout - (int64_t) bef_shift - (int64_t) buffer_size + (int64_t) avail_out};


                std::swap(bef_buffer, buffer);
                buffer_size = LENGTH;
            }
        }
    }

    if (bef_queue_data.buffer != NULL) {
        buffer_task_queue -> push(bef_queue_data);
    }
}

void read_fastq_gz_thread(FILE* fp, gz_index **built, TBBQueue* buffer_task_queue) {
    int num = 0;
    int shift = 0;
    int bef_shift;
    int idx;

    // ZLIB index
    *built = NULL;
    off_t span = SPAN;

    // Create and initialize the index list.
    gz_index *index = (gz_index*)malloc(sizeof(gz_index));
    if (index == NULL) {
        fprintf(stderr,  "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    index->have = 0;
    index->mode = 0;            // entries in index->list allocation
    index->list = NULL;
    index->strm.state = Z_NULL; // so inflateEnd() can work

    // Set up the inflation state.
    index->strm.avail_in = 0;
    index->strm.avail_out = 0;

    off_t totin = 0;            // total bytes read from input
    off_t totout = 0;           // total bytes uncompressed
    off_t beg = 0;              // starting offset of last history reset
    int mode = 0;               // mode: RAW, ZLIB, or GZIP (0 => not set yet)

    unsigned char* input_buf = (unsigned char*) malloc(sizeof(unsigned char) * CHUNK);

    unsigned char* buffer = (unsigned char*) malloc(sizeof(unsigned char) * (LENGTH + 1));
    unsigned char* bef_buffer;

    int buffer_size = LENGTH;
    int bef_buffer_size = LENGTH;

    QueueData bef_queue_data = {nullptr, nullptr, -1};

    int ret;                    // the return value from zlib, or Z_ERRNO
    off_t last;                 // last access point uncompressed offset

    do {
        // Assure available input, at least until reaching EOF.
        if (index->strm.avail_in == 0) {
            index->strm.avail_in = fread(input_buf, 1, sizeof(unsigned char) * CHUNK, fp);
            totin += index->strm.avail_in;
            index->strm.next_in = input_buf;
            if (index->strm.avail_in < sizeof(unsigned char) * CHUNK && ferror(fp)) {
                ret = Z_ERRNO;
                break;
            }

            if (mode == 0) {
                // At the start of the input -- determine the type. Assume raw
                // if it is neither zlib nor gzip. This could in theory result
                // in a false positive for zlib, but in practice the fill bits
                // after a stored block are always zeros, so a raw stream won't
                // start with an 8 in the low nybble.
                mode = index->strm.avail_in == 0 ? RAW :    // will fail
                       (index->strm.next_in[0] & 0xf) == 8 ? ZLIB :
                       index->strm.next_in[0] == 0x1f ? GZIP :
                       /* else */ RAW;
                index->strm.zalloc = Z_NULL;
                index->strm.zfree = Z_NULL;
                index->strm.opaque = Z_NULL;
                ret = inflateInit2(&index->strm, mode);
                if (ret != Z_OK)
                    break;
            }
        }

        // Assure available output. This rotates the output through, for use as
        // a sliding window on the uncompressed data.
        if (index->strm.avail_out == 0) {
            index->strm.avail_out = LENGTH;
            index->strm.next_out = buffer + shift + buffer_size - LENGTH;
        }

        if (mode == RAW && index->have == 0)
            // We skip the inflate() call at the start of raw deflate data in
            // order generate an access point there. Set data_type to imitate
            // the end of a header.
            index->strm.data_type = 0x80;
        else {
            // Inflate and update the number of uncompressed bytes.
            unsigned before = index->strm.avail_out;
            ret = inflate(&index->strm, Z_BLOCK);
            totout += before - index->strm.avail_out;
        }

        if (INDEX && (index->strm.data_type & 0xc0) == 0x80 &&
            (index->have == 0 || totout - last >= span)) {
            // We are at the end of a header or a non-last deflate block, so we
            // can add an access point here. Furthermore, we are either at the
            // very start for the first access point, or there has been span or
            // more uncompressed bytes since the last access point, so we want
            // to add an access point here.
            index = add_point(index, totin - index->strm.avail_in, totout, beg,
                              buffer + shift, bef_buffer + bef_shift, buffer_size, bef_buffer_size);
            if (index == NULL) {
                ret = Z_MEM_ERROR;
                break;
            }
            last = totout;
        }

         if ((index->strm.avail_out == 0 and totout > 0) or ret == Z_STREAM_END) {
            int seq_cnt = 0;
            int bef_num = num;

            LocationVector* loc_vector = new LocationVector{};

            buffer[LENGTH - index->strm.avail_out + shift] = '\0';
            for (int i = 0; i < LENGTH - index->strm.avail_out + shift; i++) {
                if (buffer[i] == '\n') {
                    num += 1;
                    if ((num & 3) == 2) {
                        if ((i - 1) - (idx + 1) + 1 > MAX_SEQ) {
                            fprintf(stderr, "This mode is designed for short-read sequencing. Please use 'trew long'.\n");
                            exit (EXIT_FAILURE);
                        }
                        seq_cnt += 1;
                        loc_vector -> emplace_back(idx + 1, i - 1);
                    }
                    idx = i;
                }
            }

             if ((num & 3) == 1 and bef_num == num and seq_cnt == 0) {
                 delete loc_vector;

                 buffer_size += LENGTH;
                 buffer = (unsigned char*) realloc(buffer, buffer_size + LENGTH + 1);
             } else {
                 bef_shift = shift;
                 if ((num & 3) == 1) {
                     shift = LENGTH - index->strm.avail_out + shift - idx - 1;
                 }
                 else {
                     shift = 0;
                 }

                 bef_buffer = (unsigned char*) malloc(sizeof(unsigned char) * (LENGTH + shift + 1));
                 if ((num & 3) == 1) {
                     strcpy(reinterpret_cast<char*>(bef_buffer), reinterpret_cast<char*>(buffer + idx + 1));
                     idx = -1;
                 }

                 if (bef_queue_data.buffer != nullptr) {
                     buffer_task_queue -> push(bef_queue_data);
                 }

                 bef_queue_data = QueueData{reinterpret_cast<char*>(buffer), loc_vector,
                                 (int64_t) totout - (int64_t) bef_shift - (int64_t) buffer_size + (int64_t) index->strm.avail_out};


                 std::swap(bef_buffer, buffer);
                 bef_buffer_size = buffer_size;
                 buffer_size = LENGTH;
             }
        }

        if (ret == Z_STREAM_END && mode == GZIP &&
            (index->strm.avail_in || ungetc(getc(fp), fp) != EOF)) {
            // There is more input after the end of a gzip member. Reset the
            // inflate state to read another gzip member. On success, this will
            // set ret to Z_OK to continue decompressing.
            ret = inflateReset2(&index->strm, GZIP);
            beg = totout;           // reset history
        }

        // Keep going until Z_STREAM_END or error. If the compressed data ends
        // prematurely without a file read error, Z_BUF_ERROR is returned.
    } while (ret == Z_OK);

    if (ret != Z_STREAM_END) {
        // An error was encountered. Discard the index and return a negative
        // error code.

        deflate_index_free(index);
        fprintf(stderr, "ZLIB error");
        exit(ret);
    }

    if (bef_queue_data.buffer != NULL) {
        buffer_task_queue -> push(bef_queue_data);
    }

    index->mode = mode;
    index->length = totout;
    if (INDEX) {
        *built = index;
    } else {
        *built = NULL;
        deflate_index_free(index);
    }
}

void read_fastq_gz_thread_long(FILE* fp, gz_index **built, TBBQueue* buffer_task_queue) {
    int num = 0;
    int shift = 0;
    int bef_shift = 0;
    int idx;

    // ZLIB index
    *built = NULL;
    off_t span = SPAN;

    // Create and initialize the index list.
    gz_index *index = (gz_index*)malloc(sizeof(gz_index));
    if (index == NULL) {
        fprintf(stderr,  "memory allocation failure\n");
        exit(EXIT_FAILURE);
    }

    index->have = 0;
    index->mode = 0;            // entries in index->list allocation
    index->list = NULL;
    index->strm.state = Z_NULL; // so inflateEnd() can work

    // Set up the inflation state.
    index->strm.avail_in = 0;
    index->strm.avail_out = 0;

    off_t totin = 0;            // total bytes read from input
    off_t totout = 0;           // total bytes uncompressed
    off_t beg = 0;              // starting offset of last history reset
    int mode = 0;               // mode: RAW, ZLIB, or GZIP (0 => not set yet)

    unsigned char* input_buf = (unsigned char*) malloc(sizeof(unsigned char) * CHUNK);

    unsigned char* buffer = (unsigned char*) malloc(sizeof(unsigned char) * (LENGTH + 1));
    unsigned char* bef_buffer;

    int buffer_size = LENGTH;
    int bef_buffer_size = LENGTH;

    QueueData bef_queue_data = {nullptr, nullptr, -1};

    int ret;                    // the return value from zlib, or Z_ERRNO
    off_t last;                 // last access point uncompressed offset

    do {
        // Assure available input, at least until reaching EOF.
        if (index->strm.avail_in == 0) {
            index->strm.avail_in = fread(input_buf, 1, sizeof(unsigned char) * CHUNK, fp);
            totin += index->strm.avail_in;
            index->strm.next_in = input_buf;
            if (index->strm.avail_in < sizeof(unsigned char) * CHUNK && ferror(fp)) {
                ret = Z_ERRNO;
                break;
            }

            if (mode == 0) {
                // At the start of the input -- determine the type. Assume raw
                // if it is neither zlib nor gzip. This could in theory result
                // in a false positive for zlib, but in practice the fill bits
                // after a stored block are always zeros, so a raw stream won't
                // start with an 8 in the low nybble.
                mode = index->strm.avail_in == 0 ? RAW :    // will fail
                       (index->strm.next_in[0] & 0xf) == 8 ? ZLIB :
                       index->strm.next_in[0] == 0x1f ? GZIP :
                       /* else */ RAW;
                index->strm.zalloc = Z_NULL;
                index->strm.zfree = Z_NULL;
                index->strm.opaque = Z_NULL;
                ret = inflateInit2(&index->strm, mode);
                if (ret != Z_OK)
                    break;
            }
        }

        // Assure available output. This rotates the output through, for use as
        // a sliding window on the uncompressed data.
        if (index->strm.avail_out == 0) {
            index->strm.avail_out = LENGTH;
            index->strm.next_out = buffer + shift + buffer_size - LENGTH;
        }

        if (mode == RAW && index->have == 0)
            // We skip the inflate() call at the start of raw deflate data in
            // order generate an access point there. Set data_type to imitate
            // the end of a header.
            index->strm.data_type = 0x80;
        else {
            // Inflate and update the number of uncompressed bytes.
            unsigned before = index->strm.avail_out;
            ret = inflate(&index->strm, Z_BLOCK);
            totout += before - index->strm.avail_out;
        }

        if (INDEX && (index->strm.data_type & 0xc0) == 0x80 &&
            (index->have == 0 || totout - last >= span)) {
            // We are at the end of a header or a non-last deflate block, so we
            // can add an access point here. Furthermore, we are either at the
            // very start for the first access point, or there has been span or
            // more uncompressed bytes since the last access point, so we want
            // to add an access point here.
            index = add_point(index, totin - index->strm.avail_in, totout, beg,
                              buffer + shift, bef_buffer + bef_shift, buffer_size, bef_buffer_size);
            if (index == NULL) {
                ret = Z_MEM_ERROR;
                break;
            }
            last = totout;
        }

        if ((index->strm.avail_out == 0 and totout > 0) or ret == Z_STREAM_END) {
            int seq_cnt = 0;
            int bef_num = num;

            LocationVector* loc_vector = new LocationVector{};

            buffer[buffer_size - index->strm.avail_out + shift] = '\0';
            for (int i = 0; i < buffer_size - index->strm.avail_out + shift; i++) {
                if (buffer[i] == '\n') {
                    num += 1;
                    if ((num & 3) == 2) {
                        seq_cnt += 1;
                        if ((i - 1) - (idx + 1) + 1 >= SLICE_LENGTH) {
                            loc_vector -> emplace_back(idx + 1, i - 1);
                        } else {
                            int t = 1;
                        }
                    }
                    idx = i;
                }
            }

            if ((num & 3) == 1 and bef_num == num and seq_cnt == 0) {
                delete loc_vector;

                buffer_size += LENGTH;
                buffer = (unsigned char*) realloc(buffer, buffer_size + LENGTH + 1);

                printf("HELP!\n");
            } else {
                bef_shift = shift;
                if ((num & 3) == 1) {
                    shift = buffer_size - index->strm.avail_out + shift - idx - 1;
                }
                else {
                    shift = 0;
                }

                bef_buffer = (unsigned char*) malloc(sizeof(unsigned char) * (LENGTH + shift + 1));
                if ((num & 3) == 1) {
                    strcpy(reinterpret_cast<char*>(bef_buffer), reinterpret_cast<char*>(buffer + idx + 1));
                    idx = -1;
                }

                if (bef_queue_data.buffer != nullptr) {
                    buffer_task_queue -> push(bef_queue_data);
                }

                // printf("%lld\n", totout);

                bef_queue_data = QueueData{reinterpret_cast<char*>(buffer), loc_vector,
                                    (int64_t) totout - (int64_t) bef_shift - (int64_t) buffer_size + (int64_t) index->strm.avail_out};


                std::swap(bef_buffer, buffer);
                bef_buffer_size = buffer_size;
                buffer_size = LENGTH;
            }
        }

        if (ret == Z_STREAM_END && mode == GZIP &&
            (index->strm.avail_in || ungetc(getc(fp), fp) != EOF)) {
            // There is more input after the end of a gzip member. Reset the
            // inflate state to read another gzip member. On success, this will
            // set ret to Z_OK to continue decompressing.
            ret = inflateReset2(&index->strm, GZIP);
            beg = totout;           // reset history
        }

        // Keep going until Z_STREAM_END or error. If the compressed data ends
        // prematurely without a file read error, Z_BUF_ERROR is returned.
    } while (ret == Z_OK);

    if (ret != Z_STREAM_END) {
        // An error was encountered. Discard the index and return a negative
        // error code.

        deflate_index_free(index);
        fprintf(stderr, "ZLIB error");
        exit(ret);
    }

    if (bef_queue_data.buffer != NULL) {
        buffer_task_queue -> push(bef_queue_data);
    }

    index->mode = mode;
    index->length = totout;
    if (INDEX) {
        *built = index;
    } else {
        *built = NULL;
        deflate_index_free(index);
    }
}


FinalFastqOutput process_kmer(const char* file_name, uint8_t **repeat_check_table, uint32_t **rot_table,
                                  uint64_t *extract_k_mer, uint128_t *extract_k_mer_128, uint128_t *extract_k_mer_ans,
                                  ThreadData* thread_data_list, bool is_gz, gz_index **built) {
    ResultMapPairData* result_list = (ResultMapPairData*) malloc(sizeof(ResultMapPairData) * NUM_THREAD);

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

    if (is_gz) {
        FILE* fp = fopen(file_name, "rb");
        if (fp == nullptr) {
            fprintf(stderr, "File open failed\n");
            exit (EXIT_FAILURE);
        }

        read_fastq_gz_thread(fp, built, &buffer_task_queue);
        fclose(fp);
    } else {
        FILE* fp = fopen(file_name, "r");
        if (fp == nullptr) {
            fprintf(stderr, "File open failed\n");
            exit (EXIT_FAILURE);
        }

        read_fastq_thread(fp, &buffer_task_queue);
        fclose(fp);
    }

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
        result_list[NUM_THREAD - 1] = ResultMapPairData {{{new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}}, new std::vector<FastqLocData> {}};
    }

    if (NUM_THREAD > 1) {
        tasks.wait();
    }

    return process_output(file_name, result_list, rot_table, extract_k_mer_ans);
}

FinalFastqOutput process_kmer_long(const char* file_name, uint8_t **repeat_check_table, uint32_t **rot_table,
                                       uint64_t *extract_k_mer, uint128_t *extract_k_mer_128, uint128_t *extract_k_mer_ans,
                                       ThreadData* thread_data_list, bool is_gz, gz_index **built) {
    ResultMapPairData* result_list = (ResultMapPairData*) malloc(sizeof(ResultMapPairData) * NUM_THREAD);

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

    if (is_gz) {
        FILE* fp = fopen(file_name, "rb");
        if (fp == nullptr) {
            fprintf(stderr, "File open failed\n");
            exit (EXIT_FAILURE);
        }

        read_fastq_gz_thread_long(fp, built, &buffer_task_queue);
        fclose(fp);
    } else {
        FILE* fp = fopen(file_name, "r");
        if (fp == nullptr) {
            fprintf(stderr, "File open failed\n");
            exit (EXIT_FAILURE);
        }

        read_fastq_thread_long(fp, &buffer_task_queue);
        fclose(fp);
    }

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
        result_list[NUM_THREAD - 1] = ResultMapPairData {{{new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}, {new ResultMap {}, new ResultMap {}}}, new std::vector<FastqLocData> {}};
    }

    if (NUM_THREAD > 1) {
        tasks.wait();
    }

    return process_output(file_name, result_list, rot_table, extract_k_mer_ans);
}

FinalFastqOutput process_output(const char* file_name, ResultMapPairData* result_list, uint32_t **rot_table, uint128_t *extract_k_mer_ans) {
    uint128_t _t;
    uint128_t kseq;
    int64_t _tcnt;


    std::vector<FastqLocData>* fastq_loc_data = result_list[0].second;
    for (int i = 1; i < NUM_THREAD; i++) {
        for (auto &t : *(result_list[i].second)) {
            fastq_loc_data -> emplace_back(t);
        }

        delete result_list[i].second;
    }

    char* buffer = (char*) malloc(sizeof(char) * (ABS_MAX_MER + 1));

    ResultMapData result_data = result_list[0].first;
    for (int i = 1; i < NUM_THREAD; i++) {
        for (auto &[k, v]: *(result_list[i].first.forward.first)) {
            (*(result_data.forward.first))[k] += v;
        }
        for (auto &[k, v]: *(result_list[i].first.backward.first)) {
            (*(result_data.backward.first))[k] += v;
        }
        for (auto &[k, v]: *(result_list[i].first.both.first)) {
            (*(result_data.both.first))[k] += v;
        }

        for (auto &[k, v]: *(result_list[i].first.forward.second)) {
            (*(result_data.forward.second))[k] += v;
        }
        for (auto &[k, v]: *(result_list[i].first.backward.second)) {
            (*(result_data.backward.second))[k] += v;
        }
        for (auto &[k, v]: *(result_list[i].first.both.second)) {
            (*(result_data.both.second))[k] += v;
        }

        delete result_list[i].first.forward.first;
        delete result_list[i].first.backward.first;
        delete result_list[i].first.both.first;

        delete result_list[i].first.forward.second;
        delete result_list[i].first.backward.second;
        delete result_list[i].first.both.second;
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

    return FinalFastqOutput {final_result_high_vector, final_result_low_vector, fastq_loc_data};
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

TRMDirVector* final_process_output(FinalFastqData* total_result_high, FinalFastqData* total_result_low) {
    const auto final_output = new TRMDirVector {};

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

            final_output -> emplace_back(k, final_dir);
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

            int dir = 0;
            for (auto& [k, v] : *final_output) {
                if (k == score_result_vector[i].first) {
                    dir = v;
                    break;
                }
            }
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

    return final_output;
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

void get_trm_read(const std::filesystem::path &fastq_path, TRMDirVector* put_trm, FinalFastqOutput fastq_file_data, gz_index* index,
                   size_t st, size_t nd, char* temp_path) {
    TRMDirMap trm_dir = {};
    for (auto&[k, a] : *put_trm) {
        trm_dir[k] = a;
    }

    FILE *in = fopen(std::filesystem::canonical(fastq_path).string().c_str(), "rb");
    int buf_size = LENGTH;
    unsigned char* buf = (unsigned char*) malloc(sizeof(unsigned char) * (buf_size + 1));

    FILE *fp = fopen(temp_path, "w");
    for (size_t i = st; i < nd; i++) {
        auto &[pos, rht, lef, rht_seq, lef_seq] = (*fastq_file_data.k_mer_loc_vector)[i];

        if (buf_size < pos.second) {
            free(buf);
            buf_size = pos.second;
            buf = (unsigned char*) malloc(sizeof(unsigned char) * (buf_size + 1));
        }

        std::vector<std::pair<std::pair<bool, bool>, KmerSeq>> seq_vector = {};

        if (rht.first != rht.second) {
            if (rht.first != 0) {
                seq_vector.push_back({{false, true}, {rht.first, rht_seq.first}});
            }
            if (rht.second != 0) {
                seq_vector.push_back({{false, false}, {rht.second, rht_seq.second}});
            }
        }

        if (lef.first != lef.second) {
            if (lef.first != 0) {
                seq_vector.push_back({{true, true}, {lef.first, lef_seq.first}});
            }
            if (lef.second != 0) {
                seq_vector.push_back({{true, false}, {lef.second, lef_seq.second}});
            }
        }

        uint128_t _t;
        absl::flat_hash_set<KmerSeq> output_set {};
        bool seq_dir, aln_dir, loc_dir, is_high;

        for (auto &[data, seq] : seq_vector) {
            _t = get_rot_seq_128(reverse_complement_128(seq.second) >> (2 * (64 - seq.first)), seq.first);
            KmerSeq kseq = KmerSeq {seq.first, MIN(_t, seq.second)};
            if (trm_dir.contains(KmerSeq {seq.first, MIN(_t, seq.second)}) and (not output_set.contains(kseq))) {
                loc_dir = data.first;
                is_high = data.second;

                seq_dir = MIN(_t, seq.second) == seq.second;

                if (trm_dir[kseq] != 0) {
                    aln_dir = trm_dir[kseq] == 1;

                    char buffer[ABS_MAX_MER + 1];
                    int_to_four(buffer, MIN(_t, seq.second), seq.first);
                    fprintf(fp, ">%s.%c.%c\n", buffer, loc_dir ? 'F' : 'B', (seq_dir == loc_dir) == aln_dir ? 'T' : 'F');
                    ptrdiff_t got = deflate_index_extract(in, index, pos.first, buf, pos.second);
                    if (got <= 0)
                        fprintf(stderr, "zran: extraction failed: %s error\n",
                                got == Z_MEM_ERROR ? "out of memory" : "input corrupted");
                    else {
                        fwrite(buf, 1, got, fp);
                        fprintf(fp, "\n");
                    }
                    output_set.insert(kseq);
                } else {
                    fprintf(stderr, "Score direction error!\n");
                }
            }
        }
    }
    fclose(in);
    fclose(fp);
}

gz_index *get_thread_safe_index(gz_index* index) {
    if (index == NULL) {
        return NULL;
    }

    gz_index *safe_index = (gz_index*) malloc(sizeof(gz_index));
    safe_index->have = index->have;
    safe_index->mode = index->mode;
    safe_index->length = index->length;
    safe_index->list = index->list;
    safe_index->strm = index->strm;

    inflateInit2(&safe_index->strm, index->mode);
    return safe_index;
}