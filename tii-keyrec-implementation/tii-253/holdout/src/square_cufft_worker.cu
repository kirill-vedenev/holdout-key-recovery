// Optimized square-compressed Krylov supplier for the TII-249 holdout instance.
//
// This embeds Saarinen's generic holdout CUDA worker for request parsing,
// colex arithmetic, projection, and the support-complement M kernels.  The
// recurrence here is A = R M, where R is a deterministic systematic Toeplitz
// row compressor from the literal row space to the degree-5 monomial space.

#define main tii249_embedded_holdout_worker_main
#include "../third_party/saarinen/holdout_cuda_worker.cu"
#undef main

#include <cufft.h>

namespace tii249_square {

constexpr const char *STATE_SCHEMA = "tii249-square-cufft-krylov-state-v3";
constexpr const char *RESULT_SCHEMA = "tii249-square-cufft-krylov-result-v2";
constexpr const char *HEARTBEAT_SCHEMA = "tii249-square-cufft-heartbeat-v1";
constexpr const char *RECOVERY_SCHEMA = "tii249-square-cufft-recovery-result-v1";
constexpr u64 DEFAULT_COMPRESSOR_SEED = UINT64_C(0x544949323439524d);

#define CUFFT_CHECK(call) do {                                                   \
    const cufftResult status_ = (call);                                           \
    if (status_ != CUFFT_SUCCESS) {                                               \
        throw std::runtime_error("cuFFT failure code "                            \
            + std::to_string(static_cast<int>(status_)));                         \
    }                                                                             \
} while (false)

__device__ __forceinline__ u64 splitmix_device(u64 seed, u64 index) {
    u64 value = seed + (index + 1) * UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

__device__ __forceinline__ unsigned diagonal_bit(u64 seed, u64 index) {
    return static_cast<unsigned>((splitmix_device(seed, index >> 6) >> (index & 63)) & 1U);
}

__global__ void fill_toeplitz_circulant_kernel(
    double *circulant,
    u64 seed,
    u64 output_rows,
    u64 input_rows,
    u64 fft_length,
    u64 real_stride) {
    const u64 stride = static_cast<u64>(blockDim.x) * gridDim.x;
    for (u64 index = static_cast<u64>(blockIdx.x) * blockDim.x + threadIdx.x;
         index < real_stride;
         index += stride) {
        double value = 0.0;
        if (index < fft_length) {
            u64 diagonal_index = 0;
            bool present = false;
            if (index < output_rows) {
                diagonal_index = input_rows - 1 + index;
                present = true;
            } else if (index > fft_length - input_rows) {
                const u64 distance = fft_length - index;
                diagonal_index = input_rows - 1 - distance;
                present = true;
            }
            if (present) {
                value = static_cast<double>(diagonal_bit(seed, diagonal_index));
            }
        }
        circulant[index] = value;
    }
}

__global__ void unpack_tail_words_kernel(
    double *buffer,
    const u64 *tail,
    u64 tail_rows,
    int panel_words,
    int word_base,
    int local_words,
    u64 real_stride) {
    const u64 items = tail_rows * static_cast<u64>(local_words);
    const u64 stride = static_cast<u64>(blockDim.x) * gridDim.x;
    for (u64 item = static_cast<u64>(blockIdx.x) * blockDim.x + threadIdx.x;
         item < items;
         item += stride) {
        const int local_word = static_cast<int>(item / tail_rows);
        const u64 row = item - static_cast<u64>(local_word) * tail_rows;
        const u64 word = tail[row * panel_words + word_base + local_word];
        for (int bit = 0; bit < 64; ++bit) {
            buffer[(static_cast<u64>(local_word) * 64 + bit) * real_stride + row]
                = static_cast<double>((word >> bit) & 1U);
        }
    }
}

__global__ void multiply_spectra_kernel(
    cufftDoubleComplex *panels,
    const cufftDoubleComplex *circulant,
    u64 frequencies,
    u64 lanes) {
    const u64 items = frequencies * lanes;
    const u64 stride = static_cast<u64>(blockDim.x) * gridDim.x;
    for (u64 item = static_cast<u64>(blockIdx.x) * blockDim.x + threadIdx.x;
         item < items;
         item += stride) {
        const u64 frequency = item % frequencies;
        const cufftDoubleComplex left = panels[item];
        const cufftDoubleComplex right = circulant[frequency];
        panels[item] = {
            left.x * right.x - left.y * right.y,
            left.x * right.y + left.y * right.x,
        };
    }
}

__global__ void pack_xor_words_kernel(
    u64 *output,
    const double *buffer,
    u64 output_rows,
    int panel_words,
    int word_base,
    int local_words,
    u64 real_stride,
    double inverse_length) {
    const u64 items = output_rows * static_cast<u64>(local_words);
    const u64 stride = static_cast<u64>(blockDim.x) * gridDim.x;
    for (u64 item = static_cast<u64>(blockIdx.x) * blockDim.x + threadIdx.x;
         item < items;
         item += stride) {
        const int local_word = static_cast<int>(item / output_rows);
        const u64 row = item - static_cast<u64>(local_word) * output_rows;
        u64 word = 0;
        for (int bit = 0; bit < 64; ++bit) {
            const u64 lane = static_cast<u64>(local_word) * 64 + bit;
            const double scaled = buffer[lane * real_stride + row] * inverse_length;
            const long long rounded = llrint(scaled);
            word |= static_cast<u64>(rounded & 1LL) << bit;
        }
        output[row * panel_words + word_base + local_word] ^= word;
    }
}

__global__ void fill_csr_slice_class_kernel(
    std::uint32_t *columns,
    const u64 *row_offsets,
    u64 total_rows,
    u64 global_row_offset,
    u64 rows_per_slice,
    int complement_size,
    int support_column_size,
    int complement_weight,
    int support_weight,
    const unsigned char *__restrict__ complement_positions,
    const unsigned char *__restrict__ support_positions,
    int row_size_count,
    int size0,
    int size1,
    int size2,
    int size3,
    u64 offset0,
    u64 offset1,
    u64 offset2,
    u64 offset3) {
    const int sizes[4] = {size0, size1, size2, size3};
    const u64 offsets[4] = {offset0, offset1, offset2, offset3};
    const u64 stride = static_cast<u64>(blockDim.x) * gridDim.x;
    for (u64 row = static_cast<u64>(blockIdx.x) * blockDim.x + threadIdx.x;
         row < total_rows;
         row += stride) {
        const u64 slice = row / rows_per_slice;
        const u64 local_row = row - slice * rows_per_slice;
        int row_size_index = -1;
        for (int index = 0; index < row_size_count; ++index) {
            const u64 count = device_choose(support_weight, sizes[index]);
            if (local_row >= offsets[index] && local_row - offsets[index] < count) {
                row_size_index = index;
                break;
            }
        }
        if (row_size_index < 0) continue;
        const int row_size = sizes[row_size_index];
        const u64 local_derivative = colex_unrank_device(
            local_row - offsets[row_size_index], row_size, support_weight);
        const u64 local_complement = colex_unrank_device(
            slice, complement_size, complement_weight);
        const u64 support_mask = support_weight == 64
            ? ~u64{0}
            : ((u64{1} << support_weight) - 1);
        const u64 available = support_mask & ~local_derivative;
        const int extension_size = support_column_size - row_size;
        unsigned char available_positions[MAX_K]{};
        const int available_count = mask_positions(available, available_positions);
        u64 cursor = row_offsets[global_row_offset + row];
        if (extension_size <= available_count) {
            int indices[MAX_D]{};
            for (int index = 0; index < extension_size; ++index) indices[index] = index;
            bool more = true;
            while (more) {
                u64 extension = 0;
                for (int index = 0; index < extension_size; ++index) {
                    extension |= u64{1} << available_positions[indices[index]];
                }
                const u64 global_column = colex_rank_merged_local_masks_device(
                    local_complement,
                    complement_positions,
                    local_derivative | extension,
                    support_positions);
                columns[cursor++] = static_cast<std::uint32_t>(global_column);
                more = next_indices(indices, extension_size, available_count);
            }
        }
    }
}

template<int WORD_CAP>
__global__ void csr_spmv_panel_kernel(
    const u64 *__restrict__ row_offsets,
    const std::uint32_t *__restrict__ columns,
    const u64 *__restrict__ input,
    u64 *__restrict__ output,
    u64 rows,
    int words) {
    const unsigned lane = threadIdx.x & 31U;
    const u64 warp = (static_cast<u64>(blockIdx.x) * blockDim.x + threadIdx.x) >> 5;
    const u64 warp_stride = (static_cast<u64>(gridDim.x) * blockDim.x) >> 5;
    for (u64 row = warp; row < rows; row += warp_stride) {
        u64 accum[WORD_CAP]{};
        const u64 begin = row_offsets[row];
        const u64 end = row_offsets[row + 1];
        for (u64 edge = begin + lane; edge < end; edge += 32) {
            const u64 *source = input + static_cast<u64>(columns[edge]) * words;
#pragma unroll
            for (int word = 0; word < WORD_CAP; ++word) {
                if (word < words) accum[word] ^= source[word];
            }
        }
#pragma unroll
        for (int word = 0; word < WORD_CAP; ++word) {
            if (word < words) {
                u64 value = accum[word];
                for (int offset = 16; offset; offset >>= 1) {
                    value ^= static_cast<u64>(__shfl_xor_sync(
                        0xffffffffU,
                        static_cast<unsigned long long>(value),
                        offset));
                }
                if (lane == 0) output[row * words + word] = value;
            }
        }
    }
}


template<int WORD_CAP>
__global__ void accumulate_relation16_kernel(
    const u64 *__restrict__ current,
    const u64 *__restrict__ relation16,
    u64 *__restrict__ candidates,
    u64 rows,
    int words) {
    const u64 stride = static_cast<u64>(blockDim.x) * gridDim.x;
    for (u64 row = static_cast<u64>(blockIdx.x) * blockDim.x + threadIdx.x;
         row < rows;
         row += stride) {
        u64 accum[WORD_CAP]{};
        const u64 *source = current + row * static_cast<u64>(words);
        for (int in_word = 0; in_word < WORD_CAP; ++in_word) {
            if (in_word >= words) break;
            const u64 value = source[in_word];
            for (int chunk = 0; chunk < 4; ++chunk) {
                const unsigned index = static_cast<unsigned>((value >> (16 * chunk)) & 0xffffU);
                const u64 *entry = relation16
                    + (((static_cast<u64>(in_word) * 4 + chunk) << 16) + index)
                        * static_cast<u64>(words);
#pragma unroll
                for (int out_word = 0; out_word < WORD_CAP; ++out_word) {
                    if (out_word < words) accum[out_word] ^= entry[out_word];
                }
            }
        }
        u64 *target = candidates + row * static_cast<u64>(words);
#pragma unroll
        for (int out_word = 0; out_word < WORD_CAP; ++out_word) {
            if (out_word < words) target[out_word] ^= accum[out_word];
        }
    }
}

__global__ void reduce_or_kernel(
    const u64 *__restrict__ values,
    u64 *output,
    u64 word_count) {
    u64 accum = 0;
    const u64 stride = static_cast<u64>(blockDim.x) * gridDim.x;
    for (u64 index = static_cast<u64>(blockIdx.x) * blockDim.x + threadIdx.x;
         index < word_count;
         index += stride) {
        accum |= values[index];
    }
    if (accum != 0) {
        atomicOr(reinterpret_cast<unsigned long long *>(output),
                 static_cast<unsigned long long>(accum));
    }
}

u64 next_power_of_two(u64 value) {
    u64 output = 1;
    while (output < value) {
        if (output > (UINT64_MAX >> 1)) {
            throw std::runtime_error("FFT length overflow");
        }
        output <<= 1;
    }
    return output;
}

struct ToeplitzCompressor {
    u64 output_rows;
    u64 input_rows;
    u64 fft_length;
    u64 real_stride;
    u64 frequencies;
    int panel_words;
    int batch_words;
    u64 seed;
    double *buffer = nullptr;
    double *circulant_real = nullptr;
    cufftDoubleComplex *circulant = nullptr;
    cufftHandle forward{};
    cufftHandle inverse{};

    ToeplitzCompressor(
        u64 output_rows_value,
        u64 input_rows_value,
        int panel_words_value,
        int batch_words_value,
        u64 seed_value)
        : output_rows(output_rows_value),
          input_rows(input_rows_value),
          fft_length(next_power_of_two(output_rows_value + input_rows_value - 1)),
          real_stride(fft_length + 2),
          frequencies(fft_length / 2 + 1),
          panel_words(panel_words_value),
          batch_words(batch_words_value),
          seed(seed_value) {
        if (input_rows == 0 || output_rows == 0 || panel_words <= 0) {
            throw std::runtime_error("invalid Toeplitz compressor dimensions");
        }
        if (batch_words <= 0 || batch_words > panel_words) {
            throw std::runtime_error("invalid Toeplitz batch word count");
        }
        CUDA_CHECK(cudaMalloc(&circulant_real, real_stride * sizeof(double)));
        fill_toeplitz_circulant_kernel<<<launch_blocks(real_stride), THREADS>>>(
            circulant_real, seed, output_rows, input_rows, fft_length, real_stride);
        CUDA_CHECK(cudaGetLastError());
        cufftHandle single{};
        CUFFT_CHECK(cufftPlan1d(&single, static_cast<int>(fft_length), CUFFT_D2Z, 1));
        CUFFT_CHECK(cufftExecD2Z(
            single,
            reinterpret_cast<cufftDoubleReal *>(circulant_real),
            reinterpret_cast<cufftDoubleComplex *>(circulant_real)));
        CUFFT_CHECK(cufftDestroy(single));
        circulant = reinterpret_cast<cufftDoubleComplex *>(circulant_real);

        const int lanes = batch_words * 64;
        CUDA_CHECK(cudaMalloc(
            &buffer,
            static_cast<std::size_t>(lanes) * real_stride * sizeof(double)));
        int shape[1] = {static_cast<int>(fft_length)};
        int real_embed[1] = {static_cast<int>(real_stride)};
        int complex_embed[1] = {static_cast<int>(frequencies)};
        CUFFT_CHECK(cufftPlanMany(
            &forward, 1, shape,
            real_embed, 1, static_cast<int>(real_stride),
            complex_embed, 1, static_cast<int>(frequencies),
            CUFFT_D2Z, lanes));
        CUFFT_CHECK(cufftPlanMany(
            &inverse, 1, shape,
            complex_embed, 1, static_cast<int>(frequencies),
            real_embed, 1, static_cast<int>(real_stride),
            CUFFT_Z2D, lanes));
    }

    ~ToeplitzCompressor() {
        if (inverse) cufftDestroy(inverse);
        if (forward) cufftDestroy(forward);
        cudaFree(buffer);
        cudaFree(circulant_real);
    }

    void apply(const u64 *range, u64 *output) {
        const u64 head_bytes = output_rows * static_cast<u64>(panel_words) * sizeof(u64);
        CUDA_CHECK(cudaMemcpy(output, range, head_bytes, cudaMemcpyDeviceToDevice));
        const u64 *tail = range + output_rows * static_cast<u64>(panel_words);
        const double inverse_length = 1.0 / static_cast<double>(fft_length);
        for (int word_base = 0; word_base < panel_words; word_base += batch_words) {
            const int local_words = std::min(batch_words, panel_words - word_base);
            const int lanes = local_words * 64;
            CUDA_CHECK(cudaMemset(
                buffer, 0, static_cast<std::size_t>(lanes) * real_stride * sizeof(double)));
            unpack_tail_words_kernel<<<
                launch_blocks(input_rows * static_cast<u64>(local_words)), THREADS>>>(
                buffer, tail, input_rows, panel_words, word_base, local_words, real_stride);
            CUDA_CHECK(cudaGetLastError());
            CUFFT_CHECK(cufftExecD2Z(
                forward,
                reinterpret_cast<cufftDoubleReal *>(buffer),
                reinterpret_cast<cufftDoubleComplex *>(buffer)));
            multiply_spectra_kernel<<<
                launch_blocks(frequencies * static_cast<u64>(lanes)), THREADS>>>(
                reinterpret_cast<cufftDoubleComplex *>(buffer),
                circulant,
                frequencies,
                lanes);
            CUDA_CHECK(cudaGetLastError());
            CUFFT_CHECK(cufftExecZ2D(
                inverse,
                reinterpret_cast<cufftDoubleComplex *>(buffer),
                reinterpret_cast<cufftDoubleReal *>(buffer)));
            pack_xor_words_kernel<<<
                launch_blocks(output_rows * static_cast<u64>(local_words)), THREADS>>>(
                output,
                buffer,
                output_rows,
                panel_words,
                word_base,
                local_words,
                real_stride,
                inverse_length);
            CUDA_CHECK(cudaGetLastError());
        }
    }
};

struct State {
    u64 completed_steps = 0;
    std::vector<u64> panel;
};


u64 splitmix_host_word(u64 seed, u64 index) {
    u64 value = seed + (index + 1) * UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

std::vector<u64> independent_random_words(std::size_t count, u64 seed) {
    std::vector<u64> output(count);
    for (std::size_t index = 0; index < count; ++index) {
        output[index] = splitmix_host_word(seed, static_cast<u64>(index));
    }
    return output;
}

std::string operator_fingerprint(const Request &request) {
    std::vector<u64> binding{static_cast<u64>(request.k), static_cast<u64>(request.degree),
        static_cast<u64>(request.words), request.expected_columns, request.expected_rows};
    for (const Point &point : request.points) {
        binding.push_back(point.support);
        binding.push_back(point.support_high);
        binding.push_back(point.levels.size());
        for (int level : point.levels) binding.push_back(static_cast<u64>(level));
    }
    return hex_u64(checksum_words(binding));
}

std::vector<u64> build_row_offsets_host(
    const Request &request,
    const SlicePlan &plan) {
    std::vector<u64> offsets(static_cast<std::size_t>(request.expected_rows) + 1);
    u64 cursor = 0;
    u64 expected_row = 0;
    for (std::size_t point_index = 0; point_index < request.points.size(); ++point_index) {
        const Point &point = request.points[point_index];
        const int support_weight = point_weight(point);
        for (const SliceClass &slice_class : plan.points[point_index].classes) {
            const u64 total_rows = slice_class.slice_count * slice_class.rows_per_slice;
            if (slice_class.range_offset != expected_row) {
                throw std::runtime_error("slice ranges are not contiguous in CSR builder");
            }
            for (u64 row = 0; row < total_rows; ++row) {
                const u64 local_row = row % slice_class.rows_per_slice;
                int row_size = -1;
                for (int index = 0; index < slice_class.row_size_count; ++index) {
                    const u64 count = HOST_CHOOSE.get(support_weight, slice_class.row_sizes[index]);
                    if (local_row >= slice_class.row_offsets[index]
                        && local_row - slice_class.row_offsets[index] < count) {
                        row_size = slice_class.row_sizes[index];
                        break;
                    }
                }
                if (row_size < 0) {
                    throw std::runtime_error("CSR row size lookup failed");
                }
                offsets[static_cast<std::size_t>(expected_row + row)] = cursor;
                cursor += HOST_CHOOSE.get(
                    support_weight - row_size,
                    slice_class.support_column_size - row_size);
            }
            expected_row += total_rows;
            offsets[static_cast<std::size_t>(expected_row)] = cursor;
        }
    }
    if (expected_row != request.expected_rows || cursor != request.expected_nonzeros) {
        throw std::runtime_error("CSR host offsets do not match request dimensions");
    }
    return offsets;
}

struct DeviceCsr {
    u64 *row_offsets = nullptr;
    std::uint32_t *columns = nullptr;
    u64 nonzeros = 0;
    std::size_t bytes = 0;

    DeviceCsr() = default;
    DeviceCsr(const DeviceCsr &) = delete;
    DeviceCsr &operator=(const DeviceCsr &) = delete;
    DeviceCsr(DeviceCsr &&other) noexcept
        : row_offsets(other.row_offsets),
          columns(other.columns),
          nonzeros(other.nonzeros),
          bytes(other.bytes) {
        other.row_offsets = nullptr;
        other.columns = nullptr;
        other.nonzeros = 0;
        other.bytes = 0;
    }
    DeviceCsr &operator=(DeviceCsr &&other) noexcept {
        if (this != &other) {
            cudaFree(columns);
            cudaFree(row_offsets);
            row_offsets = other.row_offsets;
            columns = other.columns;
            nonzeros = other.nonzeros;
            bytes = other.bytes;
            other.row_offsets = nullptr;
            other.columns = nullptr;
            other.nonzeros = 0;
            other.bytes = 0;
        }
        return *this;
    }

    ~DeviceCsr() {
        cudaFree(columns);
        cudaFree(row_offsets);
    }
};

DeviceCsr build_device_csr(
    const Request &request,
    const SlicePlan &plan,
    const unsigned char *device_positions,
    const std::vector<u64> &host_offsets) {
    DeviceCsr csr;
    csr.nonzeros = host_offsets.back();
    const std::size_t offset_bytes = host_offsets.size() * sizeof(u64);
    const std::size_t column_bytes =
        static_cast<std::size_t>(csr.nonzeros) * sizeof(std::uint32_t);
    csr.bytes = offset_bytes + column_bytes;
    CUDA_CHECK(cudaMalloc(&csr.row_offsets, offset_bytes));
    CUDA_CHECK(cudaMalloc(&csr.columns, column_bytes));
    CUDA_CHECK(cudaMemcpy(
        csr.row_offsets, host_offsets.data(), offset_bytes, cudaMemcpyHostToDevice));
    for (std::size_t point_index = 0; point_index < request.points.size(); ++point_index) {
        const Point &point = request.points[point_index];
        const int support_weight = point_weight(point);
        const int complement_weight = request.k - support_weight;
        const unsigned char *support_positions = device_positions + point_index * 2 * MAX_K;
        const unsigned char *complement_positions = support_positions + MAX_K;
        for (const SliceClass &slice_class : plan.points[point_index].classes) {
            const u64 total_rows = slice_class.slice_count * slice_class.rows_per_slice;
            fill_csr_slice_class_kernel<<<launch_blocks(total_rows), THREADS>>>(
                csr.columns,
                csr.row_offsets,
                total_rows,
                slice_class.range_offset,
                slice_class.rows_per_slice,
                slice_class.complement_size,
                slice_class.support_column_size,
                complement_weight,
                support_weight,
                complement_positions,
                support_positions,
                slice_class.row_size_count,
                slice_class.row_sizes[0],
                slice_class.row_sizes[1],
                slice_class.row_sizes[2],
                slice_class.row_sizes[3],
                slice_class.row_offsets[0],
                slice_class.row_offsets[1],
                slice_class.row_offsets[2],
                slice_class.row_offsets[3]);
            CUDA_CHECK(cudaGetLastError());
        }
    }
    CUDA_CHECK(cudaDeviceSynchronize());
    return csr;
}

void write_state(
    const std::filesystem::path &path,
    const Request &request,
    u64 seed_y,
    u64 compressor_seed,
    const State &state) {
    std::ofstream stream(path, std::ios::binary | std::ios::trunc);
    if (!stream) throw std::runtime_error("cannot create square state");
    stream << STATE_SCHEMA << "\n"
           << "public_sha256 " << request.public_sha256 << "\n"
           << "operator_fnv1a64_le " << operator_fingerprint(request) << "\n"
           << "panel_prng splitmix64-indexed-v1\n"
           << "k " << request.k << "\n"
           << "degree " << request.degree << "\n"
           << "words " << request.words << "\n"
           << "columns " << request.expected_columns << "\n"
           << "rows " << request.expected_rows << "\n"
           << "seed_y " << seed_y << "\n"
           << "compressor_seed " << compressor_seed << "\n"
           << "completed_steps " << state.completed_steps << "\n"
           << "panel_words " << state.panel.size() << "\n"
           << "panel_fnv1a64_le " << hex_u64(checksum_words(state.panel)) << "\n"
           << "data_le\n";
    stream.write(
        reinterpret_cast<const char *>(state.panel.data()),
        static_cast<std::streamsize>(state.panel.size() * sizeof(u64)));
    if (!stream) throw std::runtime_error("cannot write square state");
}

State read_or_create_state(
    const std::filesystem::path &path,
    const Request &request,
    u64 seed_y,
    u64 compressor_seed) {
    const std::size_t panel_words =
        static_cast<std::size_t>(request.expected_columns) * request.words;
    if (path == "-") {
        return State{0, independent_random_words(panel_words, seed_y)};
    }
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("cannot open square state");
    std::string schema;
    std::getline(stream, schema);
    if (schema != STATE_SCHEMA) throw std::runtime_error("square state schema differs");
    std::string key;
    auto read_named = [&](const char *expected, auto &value) {
        if (!(stream >> key) || key != expected || !(stream >> value)) {
            throw std::runtime_error(std::string("expected square state key ") + expected);
        }
    };
    std::string public_sha256, operator_hash, panel_prng;
    int k = 0, degree = 0, words = 0;
    u64 columns = 0, rows = 0, stored_seed = 0, stored_compressor = 0;
    u64 completed = 0, declared_words = 0;
    std::string expected_fnv;
    read_named("public_sha256", public_sha256);
    read_named("operator_fnv1a64_le", operator_hash);
    read_named("panel_prng", panel_prng);
    read_named("k", k);
    read_named("degree", degree);
    read_named("words", words);
    read_named("columns", columns);
    read_named("rows", rows);
    read_named("seed_y", stored_seed);
    read_named("compressor_seed", stored_compressor);
    read_named("completed_steps", completed);
    read_named("panel_words", declared_words);
    read_named("panel_fnv1a64_le", expected_fnv);
    if (!(stream >> key) || key != "data_le") {
        throw std::runtime_error("square state data marker differs");
    }
    char newline = 0;
    stream.get(newline);
    if (public_sha256 != request.public_sha256 || operator_hash != operator_fingerprint(request)
        || panel_prng != "splitmix64-indexed-v1" || k != request.k
        || degree != request.degree || words != request.words
        || columns != request.expected_columns || rows != request.expected_rows
        || stored_seed != seed_y || stored_compressor != compressor_seed
        || declared_words != panel_words) {
        throw std::runtime_error("square state binding differs");
    }
    State state;
    state.completed_steps = completed;
    state.panel.resize(panel_words);
    stream.read(
        reinterpret_cast<char *>(state.panel.data()),
        static_cast<std::streamsize>(panel_words * sizeof(u64)));
    if (!stream) throw std::runtime_error("square state payload truncated");
    if (hex_u64(checksum_words(state.panel)) != expected_fnv) {
        throw std::runtime_error("square state checksum differs");
    }
    return state;
}

void apply_m_device(
    const Request &request,
    const SlicePlan &plan,
    const u64 *input,
    u64 *range,
    const unsigned char *positions) {
    for (std::size_t point_index = 0; point_index < request.points.size(); ++point_index) {
        const Point &point = request.points[point_index];
        for (const SliceClass &slice_class : plan.points[point_index].classes) {
            launch_m_slice_recompute_class_raw(
                request,
                point,
                slice_class,
                input,
                range + slice_class.range_offset * static_cast<u64>(request.words),
                positions,
                point_index);
        }
    }
}

void run_segment(
    const std::filesystem::path &request_path,
    const std::filesystem::path &state_path,
    u64 step_count,
    const std::filesystem::path &output_directory,
    int batch_words,
    u64 seed_y,
    u64 compressor_seed) {
    if (step_count == 0) throw std::runtime_error("step count must be positive");
    if (std::filesystem::exists(output_directory)) {
        throw std::runtime_error("refusing existing output directory: " + output_directory.string());
    }
    const auto wall_started = Clock::now();
    Request request = read_request(request_path);
    if (request.words <= 0 || request.words > 16) {
        throw std::runtime_error("unsupported block width");
    }
    if (request.expected_rows <= request.expected_columns) {
        throw std::runtime_error("square compression requires more rows than columns");
    }
    CUDA_CHECK(cudaMemcpyToSymbol(
        DEVICE_CHOOSE, HOST_CHOOSE.values.data(), sizeof(HOST_CHOOSE.values)));
    State state = read_or_create_state(state_path, request, seed_y, compressor_seed);
    const u64 start_step = state.completed_steps;
    const u64 end_step = start_step + step_count;

    const auto preprocessing_started = Clock::now();
    const SlicePlan plan = build_slice_plan(request, false);
    const std::vector<u64> host_row_offsets = build_row_offsets_host(request, plan);
    const std::vector<unsigned char> positions = point_position_tables(request);
    const std::size_t panel_words =
        static_cast<std::size_t>(request.expected_columns) * request.words;
    const std::size_t panel_bytes = panel_words * sizeof(u64);
    const std::size_t range_words =
        static_cast<std::size_t>(request.expected_rows) * request.words;
    const std::size_t range_bytes = range_words * sizeof(u64);
    const std::size_t position_bytes = request.points.size() * 2 * MAX_K;
    const int projection_chunks = projection_chunk_count(request.expected_columns);
    const std::size_t term_words = static_cast<std::size_t>(64) * request.words * request.words;
    const std::size_t term_bytes = term_words * sizeof(u64);
    const std::size_t partial_bytes = term_bytes * projection_chunks;

    std::size_t free_bytes = 0, total_bytes = 0;
    CUDA_CHECK(cudaMemGetInfo(&free_bytes, &total_bytes));
    const std::size_t csr_bytes =
        host_row_offsets.size() * sizeof(u64)
        + static_cast<std::size_t>(host_row_offsets.back()) * sizeof(std::uint32_t);
    const std::size_t allocated_estimate =
        3 * panel_bytes + range_bytes + position_bytes + partial_bytes + term_bytes
        + csr_bytes;
    if (allocated_estimate > free_bytes || free_bytes - allocated_estimate < (std::size_t{8} << 30)) {
        throw std::runtime_error("square CUDA allocation would leave less than 8 GiB free");
    }

    u64 *device_state = nullptr;
    u64 *device_next = nullptr;
    u64 *device_z = nullptr;
    u64 *device_range = nullptr;
    u64 *device_partials = nullptr;
    u64 *device_term = nullptr;
    unsigned char *device_positions = nullptr;
    struct Cleanup {
        u64 *&a; u64 *&b; u64 *&z; u64 *&r; u64 *&p; u64 *&t; unsigned char *&pos;
        ~Cleanup() {
            cudaFree(pos); cudaFree(t); cudaFree(p); cudaFree(r);
            cudaFree(z); cudaFree(b); cudaFree(a);
        }
    } cleanup{device_state, device_next, device_z, device_range, device_partials, device_term, device_positions};
    CUDA_CHECK(cudaMalloc(&device_state, panel_bytes));
    CUDA_CHECK(cudaMalloc(&device_next, panel_bytes));
    CUDA_CHECK(cudaMalloc(&device_z, panel_bytes));
    CUDA_CHECK(cudaMalloc(&device_range, range_bytes));
    CUDA_CHECK(cudaMalloc(&device_partials, partial_bytes));
    CUDA_CHECK(cudaMalloc(&device_term, term_bytes));
    CUDA_CHECK(cudaMalloc(&device_positions, position_bytes));

    std::vector<u64> z = independent_random_words(panel_words, seed_y ^ PROJECTION_SEED_XOR);
    CUDA_CHECK(cudaMemcpy(device_state, state.panel.data(), panel_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(device_z, z.data(), panel_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(device_positions, positions.data(), position_bytes, cudaMemcpyHostToDevice));
    z.clear();
    z.shrink_to_fit();

    DeviceCsr csr = build_device_csr(request, plan, device_positions, host_row_offsets);
    ToeplitzCompressor compressor(
        request.expected_columns,
        request.expected_rows - request.expected_columns,
        request.words,
        batch_words,
        compressor_seed);
    CUDA_CHECK(cudaDeviceSynchronize());
    const double preprocessing_seconds = std::chrono::duration<double>(
        Clock::now() - preprocessing_started).count();

    if (!std::filesystem::create_directories(output_directory)) {
        throw std::runtime_error("cannot create output directory");
    }
    const auto sequence_tmp = output_directory / "sequence.raw.bin.tmp";
    const auto sequence_path = output_directory / "sequence.raw.bin";
    const auto state_tmp = output_directory / "state.bin.tmp";
    const auto state_out = output_directory / "state.bin";
    std::ofstream sequence(sequence_tmp, std::ios::binary | std::ios::trunc);
    if (!sequence) throw std::runtime_error("cannot create square sequence");

    cudaEvent_t before_m = nullptr, after_m = nullptr, after_r = nullptr, before_p = nullptr, after_p = nullptr;
    CUDA_CHECK(cudaEventCreate(&before_m));
    CUDA_CHECK(cudaEventCreate(&after_m));
    CUDA_CHECK(cudaEventCreate(&after_r));
    CUDA_CHECK(cudaEventCreate(&before_p));
    CUDA_CHECK(cudaEventCreate(&after_p));
    struct EventCleanup {
        cudaEvent_t &a; cudaEvent_t &b; cudaEvent_t &c; cudaEvent_t &d; cudaEvent_t &e;
        ~EventCleanup() { cudaEventDestroy(e); cudaEventDestroy(d); cudaEventDestroy(c); cudaEventDestroy(b); cudaEventDestroy(a); }
    } event_cleanup{before_m, after_m, after_r, before_p, after_p};

    std::vector<u64> term(term_words);
    u64 sequence_checksum = UINT64_C(1469598103934665603);
    double m_seconds = 0.0, r_seconds = 0.0, projection_seconds = 0.0;
    for (u64 local_step = 0; local_step < step_count; ++local_step) {
        // CADO's homogeneous lingen path discards coefficient zero. Supply
        // S_0=Z^T Y so that its coefficient one is Z^T A Y, not Z^T A^2 Y.
        // Shifting this by one produces ker(A^2) candidates and can accept
        // an extra Jordan-chain direction if only compressed checks are used.
        CUDA_CHECK(cudaEventRecord(before_p));
        launch_projection_device(
            device_z, device_state, device_partials, device_term,
            request.expected_columns, request.words, projection_chunks);
        CUDA_CHECK(cudaEventRecord(after_p));
        CUDA_CHECK(cudaEventRecord(before_m));
        if (request.words <= 8) {
            csr_spmv_panel_kernel<8><<<launch_blocks(request.expected_rows * 32), 256>>>(
                csr.row_offsets,
                csr.columns,
                device_state,
                device_range,
                request.expected_rows,
                request.words);
        } else {
            csr_spmv_panel_kernel<16><<<launch_blocks(request.expected_rows * 32), 256>>>(
                csr.row_offsets,
                csr.columns,
                device_state,
                device_range,
                request.expected_rows,
                request.words);
        }
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaEventRecord(after_m));
        compressor.apply(device_range, device_next);
        CUDA_CHECK(cudaEventRecord(after_r));
        CUDA_CHECK(cudaEventSynchronize(after_r));
        float m_ms = 0, r_ms = 0, p_ms = 0;
        CUDA_CHECK(cudaEventElapsedTime(&m_ms, before_m, after_m));
        CUDA_CHECK(cudaEventElapsedTime(&r_ms, after_m, after_r));
        CUDA_CHECK(cudaEventElapsedTime(&p_ms, before_p, after_p));
        m_seconds += m_ms / 1000.0;
        r_seconds += r_ms / 1000.0;
        projection_seconds += p_ms / 1000.0;
        CUDA_CHECK(cudaMemcpy(term.data(), device_term, term_bytes, cudaMemcpyDeviceToHost));
        sequence.write(
            reinterpret_cast<const char *>(term.data()),
            static_cast<std::streamsize>(term_bytes));
        if (!sequence) throw std::runtime_error("cannot append square sequence term");
        sequence_checksum = update_checksum_words(sequence_checksum, term);
        std::swap(device_state, device_next);
        const u64 completed = start_step + local_step + 1;
        std::cout << "{\"schema\":\"" << HEARTBEAT_SCHEMA
                  << "\",\"completed_steps\":" << completed
                  << ",\"segment_end_step\":" << end_step
                  << ",\"elapsed_seconds\":"
                  << std::chrono::duration<double>(Clock::now() - wall_started).count()
                  << ",\"m_seconds_total\":" << m_seconds
                  << ",\"compress_seconds_total\":" << r_seconds
                  << "}" << std::endl;
    }
    sequence.flush();
    sequence.close();
    if (!sequence) throw std::runtime_error("cannot finalize square sequence");

    state.completed_steps = end_step;
    state.panel.resize(panel_words);
    CUDA_CHECK(cudaMemcpy(state.panel.data(), device_state, panel_bytes, cudaMemcpyDeviceToHost));
    write_state(state_tmp, request, seed_y, compressor_seed, state);
    fsync_path(sequence_tmp, false);
    fsync_path(state_tmp, false);
    publish_temporary_new(sequence_tmp, sequence_path);
    publish_temporary_new(state_tmp, state_out);
    fsync_path(output_directory, true);

    const auto result_tmp = output_directory / "result.json.tmp";
    const auto result_path = output_directory / "result.json";
    std::ostringstream result;
    result << std::setprecision(12)
           << "{\n"
           << "  \"schema\": \"" << RESULT_SCHEMA << "\",\n"
           << "  \"terminal\": \"tii249_square_cufft_krylov_segment_complete\",\n"
           << "  \"public_sha256\": \"" << request.public_sha256 << "\",\n"
           << "  \"operator_fnv1a64_le\": \"" << operator_fingerprint(request) << "\",\n"
           << "  \"operator\": \"A=R*M with deterministic systematic Toeplitz R\",\n"
           << "  \"sequence_convention\": \"S_i=Z^T*A^i*Y; includes S_0\",\n"
           << "  \"first_power\": " << start_step << ",\n"
           << "  \"columns\": " << request.expected_columns << ",\n"
           << "  \"rows\": " << request.expected_rows << ",\n"
           << "  \"tail_rows\": " << (request.expected_rows - request.expected_columns) << ",\n"
           << "  \"words\": " << request.words << ",\n"
           << "  \"block_bits\": " << request.words * 64 << ",\n"
           << "  \"fft_length\": " << compressor.fft_length << ",\n"
           << "  \"batch_words\": " << batch_words << ",\n"
           << "  \"start_step\": " << start_step << ",\n"
           << "  \"term_count\": " << step_count << ",\n"
           << "  \"end_step\": " << end_step << ",\n"
           << "  \"seed_y\": " << seed_y << ",\n"
           << "  \"seed_z\": " << (seed_y ^ PROJECTION_SEED_XOR) << ",\n"
           << "  \"panel_prng\": \"splitmix64-indexed-v1\",\n"
           << "  \"compressor_seed\": " << compressor_seed << ",\n"
           << "  \"free_device_bytes_before\": " << free_bytes << ",\n"
           << "  \"total_device_bytes\": " << total_bytes << ",\n"
           << "  \"allocated_bytes_estimate\": " << allocated_estimate << ",\n"
           << "  \"csr_bytes\": " << csr.bytes << ",\n"
           << "  \"csr_nonzeros\": " << csr.nonzeros << ",\n"
           << "  \"preprocessing_seconds\": " << preprocessing_seconds << ",\n"
           << "  \"m_seconds_total\": " << m_seconds << ",\n"
           << "  \"compress_seconds_total\": " << r_seconds << ",\n"
           << "  \"projection_seconds_total\": " << projection_seconds << ",\n"
           << "  \"wall_seconds\": "
           << std::chrono::duration<double>(Clock::now() - wall_started).count() << ",\n"
           << "  \"sequence_fnv1a64_le\": \"" << hex_u64(sequence_checksum) << "\",\n"
           << "  \"state_fnv1a64_le\": \"" << hex_u64(checksum_words(state.panel)) << "\",\n"
           << "  \"sequence_file\": \"sequence.raw.bin\",\n"
           << "  \"sequence_bytes\": " << std::filesystem::file_size(sequence_path) << ",\n"
           << "  \"state_file\": \"state.bin\",\n"
           << "  \"state_bytes\": " << std::filesystem::file_size(state_out) << "\n"
           << "}\n";
    std::ofstream out(result_tmp);
    if (!out) throw std::runtime_error("cannot create result json");
    out << result.str();
    out.close();
    publish_temporary_new(result_tmp, result_path);
    fsync_path(result_path, false);
    fsync_path(output_directory, true);
    std::cout << result.str();
}



void apply_square_csr(
    const Request &request,
    const DeviceCsr &csr,
    ToeplitzCompressor &compressor,
    const u64 *device_state,
    u64 *device_next,
    u64 *device_range) {
    if (request.words <= 8) {
        csr_spmv_panel_kernel<8><<<launch_blocks(request.expected_rows * 32), 256>>>(
            csr.row_offsets,
            csr.columns,
            device_state,
            device_range,
            request.expected_rows,
            request.words);
    } else {
        csr_spmv_panel_kernel<16><<<launch_blocks(request.expected_rows * 32), 256>>>(
            csr.row_offsets,
            csr.columns,
            device_state,
            device_range,
            request.expected_rows,
            request.words);
    }
    CUDA_CHECK(cudaGetLastError());
    compressor.apply(device_range, device_next);
}

std::vector<u64> maybe_transpose_coefficient(
    const std::vector<u64> &coefficient,
    int words,
    bool transpose) {
    const int bits = words * 64;
    if (!transpose) return coefficient;
    std::vector<u64> output(static_cast<std::size_t>(bits) * words, 0);
    for (int row = 0; row < bits; ++row) {
        const int row_word = row >> 6;
        const int row_bit = row & 63;
        for (int out_word = 0; out_word < words; ++out_word) {
            u64 word = 0;
            for (int bit = 0; bit < 64; ++bit) {
                const int source_row = out_word * 64 + bit;
                const u64 source = coefficient[static_cast<std::size_t>(source_row) * words + row_word];
                word |= ((source >> row_bit) & 1ULL) << bit;
            }
            output[static_cast<std::size_t>(row) * words + out_word] = word;
        }
    }
    return output;
}



u64 checksum_bytes_vector(const std::vector<unsigned char> &values) {
    u64 hash = UINT64_C(1469598103934665603);
    for (const unsigned char value : values) {
        hash ^= static_cast<u64>(value);
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

std::vector<u64> build_relation16_table(
    const std::vector<u64> &coefficient,
    int words) {
    const int chunks = words * 4;
    std::vector<u64> table(static_cast<std::size_t>(chunks) * 65536ULL * words, 0);
    for (int chunk = 0; chunk < chunks; ++chunk) {
        u64 *chunk_table = table.data() + static_cast<std::size_t>(chunk) * 65536ULL * words;
        const int base_bit = chunk * 16;
        for (unsigned index = 1; index < 65536U; ++index) {
            const unsigned low = index & (0U - index);
            const unsigned previous = index ^ low;
            const int bit = __builtin_ctz(low);
            const u64 *previous_entry = chunk_table + static_cast<std::size_t>(previous) * words;
            const u64 *row = coefficient.data() + static_cast<std::size_t>(base_bit + bit) * words;
            u64 *entry = chunk_table + static_cast<std::size_t>(index) * words;
            for (int word = 0; word < words; ++word) {
                entry[word] = previous_entry[word] ^ row[word];
            }
        }
    }
    return table;
}

struct PolynomialStats {
    std::vector<u64> weights;
    std::vector<u64> checksums;
    std::string thk1_checksum;
    std::string panel_checksum;
    std::uintmax_t thk1_bytes = 0;
    std::uintmax_t panel_bytes = 0;
};

// Pivot columns of the row space are independent columns of the original
// panel. Keep the original polynomials, not row-reduced combinations. When a
// target is given, stop as soon as it is reached (production only needs 280).
std::vector<int> independent_panel_columns(
    const std::vector<u64> &panel, u64 rows, int words, u64 target) {
    const int bits = words * 64;
    std::vector<u64> basis(static_cast<std::size_t>(bits) * words, 0);
    std::vector<bool> present(bits, false);
    std::vector<int> selected;
    std::vector<u64> row(words);
    for (u64 r = 0; r < rows; ++r) {
        std::copy_n(panel.data() + static_cast<std::size_t>(r) * words, words, row.data());
        for (int word = 0; word < words; ++word) {
            while (row[word]) {
                const int pivot = word * 64 + __builtin_ctzll(row[word]);
                u64 *stored = basis.data() + static_cast<std::size_t>(pivot) * words;
                if (!present[pivot]) {
                    std::copy(row.begin(), row.end(), stored);
                    present[pivot] = true;
                    selected.push_back(pivot);
                    break;
                }
                for (int w = word; w < words; ++w) row[w] ^= stored[w];
            }
            if (row[word]) break; // new pivot installed
        }
        if (selected.size() == (target ? target : static_cast<u64>(bits))) break;
    }
    std::sort(selected.begin(), selected.end());
    return selected;
}

void select_panel_columns(std::vector<u64> &panel, int words, const std::vector<int> &selected) {
    std::vector<u64> row(words);
    for (std::size_t begin = 0; begin < panel.size(); begin += words) {
        std::copy_n(panel.data() + begin, words, row.data());
        std::fill_n(panel.data() + begin, words, u64{0});
        for (std::size_t v = 0; v < selected.size(); ++v) {
            const int source = selected[v];
            panel[begin + v / 64] |= ((row[source / 64] >> (source % 64)) & 1ULL) << (v % 64);
        }
    }
}

void write_u64_le(std::ostream &output, u64 value) {
    for (int byte = 0; byte < 8; ++byte) {
        output.put(static_cast<char>((value >> (8 * byte)) & 0xff));
    }
}

PolynomialStats write_polynomial_artifacts(
    const std::filesystem::path &output_directory,
    const Request &request,
    u64 count,
    const std::vector<u64> &panel) {
    if (count == 0 || count > static_cast<u64>(request.words) * 64ULL) {
        throw std::runtime_error("invalid polynomial vector count");
    }
    const auto polynomial_directory = output_directory / "polynomials";
    if (!std::filesystem::create_directories(polynomial_directory)) {
        throw std::runtime_error("cannot create polynomial artifact directory");
    }
    const auto thk1_path = output_directory / "kernel.thk1";
    const auto panel_path = polynomial_directory / "kernel-panel.u64le";
    const auto index_path = polynomial_directory / "index.json";
    const auto readme_path = polynomial_directory / "README.txt";
    const auto monomial_path = polynomial_directory / "monomials-colex.u64";
    if (std::filesystem::exists(thk1_path) || std::filesystem::exists(panel_path)
        || std::filesystem::exists(index_path) || std::filesystem::exists(readme_path)
        || std::filesystem::exists(monomial_path)) {
        throw std::runtime_error("refusing existing polynomial artifacts");
    }

    const auto panel_tmp = panel_path.string() + ".tmp";
    {
        std::ofstream panel_out(panel_tmp, std::ios::binary | std::ios::trunc);
        if (!panel_out) throw std::runtime_error("cannot create packed panel artifact");
        panel_out.write(
            reinterpret_cast<const char *>(panel.data()),
            static_cast<std::streamsize>(panel.size() * sizeof(u64)));
        panel_out.flush();
        panel_out.close();
        if (!panel_out) throw std::runtime_error("packed panel write failed");
    }
    publish_temporary_new(panel_tmp, panel_path);
    fsync_path(panel_path, false);

    PolynomialStats stats;
    stats.weights.assign(static_cast<std::size_t>(count), 0);
    stats.checksums.assign(static_cast<std::size_t>(count), UINT64_C(1469598103934665603));
    const std::size_t stride = static_cast<std::size_t>((request.expected_columns + 7) / 8);
    std::vector<unsigned char> vector_bytes(stride);
    const auto thk1_tmp = thk1_path.string() + ".tmp";
    std::ofstream output(thk1_tmp, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("cannot create THK1 output");
    auto put32 = [&](std::uint32_t value) {
        for (int byte = 0; byte < 4; ++byte) output.put(static_cast<char>(value >> (8 * byte)));
    };
    output.write("THK1", 4);
    put32(static_cast<std::uint32_t>(request.k));
    put32(static_cast<std::uint32_t>(request.degree));
    put32(static_cast<std::uint32_t>(request.expected_columns));
    put32(static_cast<std::uint32_t>(count));
    for (u64 vector = 0; vector < count; ++vector) {
        std::fill(vector_bytes.begin(), vector_bytes.end(), 0);
        const int source_word = static_cast<int>(vector >> 6);
        const int source_bit = static_cast<int>(vector & 63);
        u64 weight = 0;
        for (u64 row = 0; row < request.expected_columns; ++row) {
            const u64 bit = (panel[static_cast<std::size_t>(row) * request.words + source_word]
                >> source_bit) & 1ULL;
            vector_bytes[static_cast<std::size_t>(row >> 3)] |=
                static_cast<unsigned char>(bit << (row & 7));
            weight += bit;
        }
        stats.weights[static_cast<std::size_t>(vector)] = weight;
        stats.checksums[static_cast<std::size_t>(vector)] = checksum_bytes_vector(vector_bytes);
        output.write(
            reinterpret_cast<const char *>(vector_bytes.data()),
            static_cast<std::streamsize>(vector_bytes.size()));
        if (!output) throw std::runtime_error("THK1 write failed");
    }
    output.flush();
    output.close();
    if (!output) throw std::runtime_error("THK1 finalize failed");
    publish_temporary_new(thk1_tmp, thk1_path);
    fsync_path(thk1_path, false);
    stats.thk1_bytes = std::filesystem::file_size(thk1_path);
    stats.panel_bytes = std::filesystem::file_size(panel_path);
    stats.panel_checksum = hex_u64(checksum_words(panel));

    // Compute a file-level checksum over THK1 through the vector byte checksums and header fields.
    std::vector<u64> checksum_material;
    checksum_material.reserve(4 + stats.checksums.size());
    checksum_material.push_back(static_cast<u64>(request.k));
    checksum_material.push_back(static_cast<u64>(request.degree));
    checksum_material.push_back(request.expected_columns);
    checksum_material.push_back(count);
    checksum_material.insert(checksum_material.end(), stats.checksums.begin(), stats.checksums.end());
    stats.thk1_checksum = hex_u64(checksum_words(checksum_material));

    const auto monomial_tmp = monomial_path.string() + ".tmp";
    {
        std::ofstream monomials(monomial_tmp, std::ios::binary | std::ios::trunc);
        if (!monomials) throw std::runtime_error("cannot create polynomial monomial map");
        for (u64 row = 0; row < request.expected_columns; ++row) {
            write_u64_le(monomials, colex_unrank_host(row, request.degree, request.k));
        }
        monomials.flush(); monomials.close();
        if (!monomials) throw std::runtime_error("monomial map write failed");
    }
    publish_temporary_new(monomial_tmp, monomial_path);
    fsync_path(monomial_path, false);

    const auto readme_tmp = readme_path.string() + ".tmp";
    {
        std::ofstream readme(readme_tmp, std::ios::trunc);
        if (!readme) throw std::runtime_error("cannot create polynomial README");
        readme
            << "Verified independent single-position holdout-kernel polynomials.\n"
            << "kernel.thk1 is the vector-major packed polynomial set accepted by build/holdout-cpu verify.\n"
            << "kernel-panel.u64le is row-major/colex packed as " << request.expected_columns
            << " rows times " << request.words << " u64 words; bit v in row r is coefficient of monomial r in polynomial v.\n"
            << "monomials-colex.u64 maps row indices to degree-" << request.degree << " colex monomial masks.\n"
            << "Use scripts/decode_kernel_polynomial.py to expand one polynomial if a text support list is needed.\n";
        readme.flush();
        readme.close();
        if (!readme) throw std::runtime_error("polynomial README write failed");
    }
    publish_temporary_new(readme_tmp, readme_path);
    fsync_path(readme_path, false);

    const auto index_tmp = index_path.string() + ".tmp";
    {
        std::ofstream index(index_tmp, std::ios::trunc);
        if (!index) throw std::runtime_error("cannot create polynomial index");
        index << "{\n"
              << "  \"schema\": \"tii249-holdout-polynomial-index-v1\",\n"
              << "  \"public_sha256\": \"" << request.public_sha256 << "\",\n"
              << "  \"operator_fnv1a64_le\": \"" << operator_fingerprint(request) << "\",\n"
              << "  \"k\": " << request.k << ",\n"
              << "  \"degree\": " << request.degree << ",\n"
              << "  \"monomial_order\": \"colex\",\n"
              << "  \"monomials\": " << request.expected_columns << ",\n"
              << "  \"count\": " << count << ",\n"
              << "  \"independent_rank\": " << count << ",\n"
              << "  \"panel_words\": " << request.words << ",\n"
              << "  \"kernel_file\": \"../kernel.thk1\",\n"
              << "  \"panel_file\": \"kernel-panel.u64le\",\n"
              << "  \"monomial_map\": \"monomials-colex.u64\",\n"
              << "  \"panel_fnv1a64_le\": \"" << stats.panel_checksum << "\",\n"
              << "  \"thk1_derived_fnv1a64_le\": \"" << stats.thk1_checksum << "\",\n"
              << "  \"polynomials\": [\n";
        for (u64 vector = 0; vector < count; ++vector) {
            index << "    {\"index\": " << vector
                  << ", \"weight\": " << stats.weights[static_cast<std::size_t>(vector)]
                  << ", \"packed_bytes_fnv1a64_le\": \""
                  << hex_u64(stats.checksums[static_cast<std::size_t>(vector)]) << "\"}";
            if (vector + 1 != count) index << ',';
            index << "\n";
        }
        index << "  ]\n}\n";
        index.flush();
        index.close();
        if (!index) throw std::runtime_error("polynomial index write failed");
    }
    publish_temporary_new(index_tmp, index_path);
    fsync_path(index_path, false);
    fsync_path(polynomial_directory, true);
    return stats;
}


void run_recovery(
    const std::filesystem::path &request_path,
    const std::filesystem::path &ffile_path,
    const std::filesystem::path &output_directory,
    int batch_words,
    u64 seed_y,
    u64 compressor_seed,
    bool transpose_coefficients,
    u64 output_count) {
    if (std::filesystem::exists(output_directory)) {
        throw std::runtime_error("refusing existing output directory: " + output_directory.string());
    }
    const auto wall_started = Clock::now();
    Request request = read_request(request_path);
    if (request.words <= 0 || request.words > 16) {
        throw std::runtime_error("unsupported block width");
    }
    const u64 block_bits = static_cast<u64>(request.words) * 64ULL;
    if (output_count > block_bits) {
        throw std::runtime_error("requested more vectors than block width");
    }
    const std::uintmax_t coefficient_bytes =
        static_cast<std::uintmax_t>(block_bits) * request.words * sizeof(u64);
    const std::uintmax_t ffile_bytes = std::filesystem::file_size(ffile_path);
    if (coefficient_bytes == 0 || ffile_bytes == 0 || ffile_bytes % coefficient_bytes) {
        throw std::runtime_error("generator file has partial block-matrix coefficients");
    }
    const u64 coefficient_count = static_cast<u64>(ffile_bytes / coefficient_bytes);

    CUDA_CHECK(cudaMemcpyToSymbol(
        DEVICE_CHOOSE, HOST_CHOOSE.values.data(), sizeof(HOST_CHOOSE.values)));
    const auto preprocessing_started = Clock::now();
    const SlicePlan plan = build_slice_plan(request, false);
    const std::vector<u64> host_row_offsets = build_row_offsets_host(request, plan);
    const std::vector<unsigned char> positions = point_position_tables(request);
    const std::size_t panel_words =
        static_cast<std::size_t>(request.expected_columns) * request.words;
    const std::size_t panel_bytes = panel_words * sizeof(u64);
    const std::size_t range_words =
        static_cast<std::size_t>(request.expected_rows) * request.words;
    const std::size_t range_bytes = range_words * sizeof(u64);
    const std::size_t position_bytes = request.points.size() * 2 * MAX_K;
    const std::size_t table_words =
        static_cast<std::size_t>(request.words) * 4ULL * 65536ULL * request.words;
    const std::size_t table_bytes = table_words * sizeof(u64);

    std::size_t free_bytes = 0, total_bytes = 0;
    CUDA_CHECK(cudaMemGetInfo(&free_bytes, &total_bytes));
    const std::size_t csr_bytes =
        host_row_offsets.size() * sizeof(u64)
        + static_cast<std::size_t>(host_row_offsets.back()) * sizeof(std::uint32_t);
    const std::size_t allocated_estimate =
        3 * panel_bytes + range_bytes + position_bytes + table_bytes + csr_bytes;
    if (allocated_estimate > free_bytes || free_bytes - allocated_estimate < (std::size_t{8} << 30)) {
        throw std::runtime_error("recovery CUDA allocation would leave less than 8 GiB free");
    }

    u64 *device_current = nullptr;
    u64 *device_next = nullptr;
    u64 *device_candidates = nullptr;
    u64 *device_range = nullptr;
    u64 *device_relation16 = nullptr;
    u64 *device_residual_or = nullptr;
    unsigned char *device_positions = nullptr;
    struct Cleanup {
        u64 *&current; u64 *&next; u64 *&candidates; u64 *&range; u64 *&table;
        u64 *&residual; unsigned char *&positions;
        ~Cleanup() {
            cudaFree(positions); cudaFree(residual); cudaFree(table); cudaFree(range);
            cudaFree(candidates); cudaFree(next); cudaFree(current);
        }
    } cleanup{device_current, device_next, device_candidates, device_range,
        device_relation16, device_residual_or, device_positions};
    CUDA_CHECK(cudaMalloc(&device_current, panel_bytes));
    CUDA_CHECK(cudaMalloc(&device_next, panel_bytes));
    CUDA_CHECK(cudaMalloc(&device_candidates, panel_bytes));
    CUDA_CHECK(cudaMalloc(&device_range, range_bytes));
    CUDA_CHECK(cudaMalloc(&device_relation16, table_bytes));
    CUDA_CHECK(cudaMalloc(&device_residual_or, sizeof(u64)));
    CUDA_CHECK(cudaMalloc(&device_positions, position_bytes));
    CUDA_CHECK(cudaMemset(device_candidates, 0, panel_bytes));
    std::vector<u64> initial = independent_random_words(panel_words, seed_y);
    CUDA_CHECK(cudaMemcpy(device_current, initial.data(), panel_bytes, cudaMemcpyHostToDevice));
    initial.clear();
    initial.shrink_to_fit();
    CUDA_CHECK(cudaMemcpy(device_positions, positions.data(), position_bytes, cudaMemcpyHostToDevice));

    DeviceCsr csr = build_device_csr(request, plan, device_positions, host_row_offsets);
    ToeplitzCompressor compressor(
        request.expected_columns,
        request.expected_rows - request.expected_columns,
        request.words,
        batch_words,
        compressor_seed);
    CUDA_CHECK(cudaDeviceSynchronize());
    const double preprocessing_seconds = std::chrono::duration<double>(
        Clock::now() - preprocessing_started).count();

    if (!std::filesystem::create_directories(output_directory)) {
        throw std::runtime_error("cannot create output directory");
    }
    std::ifstream ffile(ffile_path, std::ios::binary);
    if (!ffile) throw std::runtime_error("cannot open generator file");
    std::vector<u64> coefficient(static_cast<std::size_t>(block_bits) * request.words);
    double table_seconds = 0.0;
    double accumulate_seconds = 0.0;
    double apply_seconds = 0.0;
    for (u64 degree = 0; degree < coefficient_count; ++degree) {
        ffile.read(
            reinterpret_cast<char *>(coefficient.data()),
            static_cast<std::streamsize>(coefficient.size() * sizeof(u64)));
        if (!ffile) throw std::runtime_error("short generator read");
        const auto table_started = Clock::now();
        std::vector<u64> effective = maybe_transpose_coefficient(
            coefficient, request.words, transpose_coefficients);
        std::vector<u64> relation16 = build_relation16_table(effective, request.words);
        table_seconds += std::chrono::duration<double>(Clock::now() - table_started).count();
        CUDA_CHECK(cudaMemcpy(
            device_relation16,
            relation16.data(),
            table_bytes,
            cudaMemcpyHostToDevice));
        relation16.clear();
        relation16.shrink_to_fit();
        effective.clear();
        effective.shrink_to_fit();
        const auto accumulate_started = Clock::now();
        if (request.words <= 8) {
            accumulate_relation16_kernel<8><<<launch_blocks(request.expected_columns), THREADS>>>(
                device_current,
                device_relation16,
                device_candidates,
                request.expected_columns,
                request.words);
        } else {
            accumulate_relation16_kernel<16><<<launch_blocks(request.expected_columns), THREADS>>>(
                device_current,
                device_relation16,
                device_candidates,
                request.expected_columns,
                request.words);
        }
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
        accumulate_seconds += std::chrono::duration<double>(Clock::now() - accumulate_started).count();
        if (degree + 1 < coefficient_count) {
            const auto apply_started = Clock::now();
            apply_square_csr(request, csr, compressor, device_current, device_next, device_range);
            CUDA_CHECK(cudaDeviceSynchronize());
            std::swap(device_current, device_next);
            apply_seconds += std::chrono::duration<double>(Clock::now() - apply_started).count();
        }
        std::cout << "{\"schema\":\"" << RECOVERY_SCHEMA
                  << "\",\"phase\":\"coefficient\",\"completed_coefficients\":"
                  << (degree + 1) << ",\"coefficient_count\":" << coefficient_count
                  << ",\"elapsed_seconds\":"
                  << std::chrono::duration<double>(Clock::now() - wall_started).count()
                  << "}" << std::endl;
    }

    const auto verify_started = Clock::now();
    if (request.words <= 8) {
        csr_spmv_panel_kernel<8><<<launch_blocks(request.expected_rows * 32), 256>>>(
            csr.row_offsets,
            csr.columns,
            device_candidates,
            device_range,
            request.expected_rows,
            request.words);
    } else {
        csr_spmv_panel_kernel<16><<<launch_blocks(request.expected_rows * 32), 256>>>(
            csr.row_offsets,
            csr.columns,
            device_candidates,
            device_range,
            request.expected_rows,
            request.words);
    }
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemset(device_residual_or, 0, sizeof(u64)));
    reduce_or_kernel<<<launch_blocks(range_words), THREADS>>>(
        device_range, device_residual_or, range_words);
    CUDA_CHECK(cudaGetLastError());
    u64 residual_or = 0;
    CUDA_CHECK(cudaMemcpy(&residual_or, device_residual_or, sizeof(u64), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaDeviceSynchronize());
    const double verify_seconds = std::chrono::duration<double>(
        Clock::now() - verify_started).count();

    if (residual_or != 0) {
        const auto path = output_directory / "recovery-result.json";
        write_new(path, "{\"schema\":\"" + std::string(RECOVERY_SCHEMA)
            + "\",\"terminal\":\"rejected\",\"exact_kernel\":false,\"residual_or\":\""
            + hex_u64(residual_or) + "\",\"output_count\":0}\n");
        throw std::runtime_error("literal Hasse residual is nonzero; no kernel artifacts published");
    }

    std::vector<u64> panel(panel_words);
    const auto copy_started = Clock::now();
    CUDA_CHECK(cudaMemcpy(panel.data(), device_candidates, panel_bytes, cudaMemcpyDeviceToHost));
    const double copy_seconds = std::chrono::duration<double>(
        Clock::now() - copy_started).count();
    const auto rank_started = Clock::now();
    const u64 requested_count = output_count;
    const std::vector<int> selected = independent_panel_columns(
        panel, request.expected_columns, request.words, requested_count);
    if (selected.empty() || (requested_count && selected.size() < requested_count)) {
        write_new(output_directory / "recovery-result.json",
            "{\"schema\":\"" + std::string(RECOVERY_SCHEMA)
            + "\",\"terminal\":\"rejected\",\"exact_kernel\":true,\"independent_rank\":"
            + std::to_string(selected.size()) + ",\"output_count\":0}\n");
        throw std::runtime_error("insufficient independent nonzero kernel polynomials; no kernel artifacts published");
    }
    output_count = selected.size();
    select_panel_columns(panel, request.words, selected);
    const double rank_seconds = std::chrono::duration<double>(Clock::now() - rank_started).count();
    const auto export_started = Clock::now();
    PolynomialStats polynomial_stats = write_polynomial_artifacts(
        output_directory, request, output_count, panel);
    panel.clear();
    panel.shrink_to_fit();
    const double export_seconds = std::chrono::duration<double>(
        Clock::now() - export_started).count();

    const auto result_tmp = output_directory / "recovery-result.json.tmp";
    const auto result_path = output_directory / "recovery-result.json";
    std::ostringstream result;
    result << std::setprecision(12)
           << "{\n"
           << "  \"schema\": \"" << RECOVERY_SCHEMA << "\",\n"
           << "  \"terminal\": \"tii249_square_cufft_recovery_complete\",\n"
           << "  \"public_sha256\": \"" << request.public_sha256 << "\",\n"
           << "  \"operator_fnv1a64_le\": \"" << operator_fingerprint(request) << "\",\n"
           << "  \"generator_file\": \"" << ffile_path.string() << "\",\n"
           << "  \"generator_bytes\": " << ffile_bytes << ",\n"
           << "  \"coefficient_count\": " << coefficient_count << ",\n"
           << "  \"transpose_coefficients\": " << (transpose_coefficients ? "true" : "false") << ",\n"
           << "  \"columns\": " << request.expected_columns << ",\n"
           << "  \"rows\": " << request.expected_rows << ",\n"
           << "  \"words\": " << request.words << ",\n"
           << "  \"output_count\": " << output_count << ",\n"
           << "  \"requested_count\": " << requested_count << ",\n"
           << "  \"independent_rank\": " << output_count << ",\n"
           << "  \"rank_seconds\": " << rank_seconds << ",\n"
           << "  \"seed_y\": " << seed_y << ",\n"
           << "  \"panel_prng\": \"splitmix64-indexed-v1\",\n"
           << "  \"compressor_seed\": " << compressor_seed << ",\n"
           << "  \"csr_bytes\": " << csr.bytes << ",\n"
           << "  \"table_bytes\": " << table_bytes << ",\n"
           << "  \"residual_or\": \"" << hex_u64(residual_or) << "\",\n"
           << "  \"exact_kernel\": " << (residual_or == 0 ? "true" : "false") << ",\n"
           << "  \"preprocessing_seconds\": " << preprocessing_seconds << ",\n"
           << "  \"table_seconds_total\": " << table_seconds << ",\n"
           << "  \"accumulate_seconds_total\": " << accumulate_seconds << ",\n"
           << "  \"apply_seconds_total\": " << apply_seconds << ",\n"
           << "  \"verify_seconds\": " << verify_seconds << ",\n"
           << "  \"device_to_host_seconds\": " << copy_seconds << ",\n"
           << "  \"export_seconds\": " << export_seconds << ",\n"
           << "  \"wall_seconds\": "
           << std::chrono::duration<double>(Clock::now() - wall_started).count() << ",\n"
           << "  \"candidate_panel_fnv1a64_le\": \"" << polynomial_stats.panel_checksum << "\",\n"
           << "  \"kernel_file\": \"kernel.thk1\",\n"
           << "  \"polynomial_index_file\": \"polynomials/index.json\",\n"
           << "  \"polynomial_panel_file\": \"polynomials/kernel-panel.u64le\",\n"
           << "  \"kernel_thk1_bytes\": " << polynomial_stats.thk1_bytes << ",\n"
           << "  \"polynomial_panel_bytes\": " << polynomial_stats.panel_bytes << "\n"
           << "}\n";
    std::ofstream out(result_tmp);
    if (!out) throw std::runtime_error("cannot create recovery result json");
    out << result.str();
    out.close();
    publish_temporary_new(result_tmp, result_path);
    fsync_path(result_path, false);
    fsync_path(output_directory, true);
    std::cout << result.str();
}

u64 parse_u64_arg(const char *text, const char *label) {
    std::string value(text);
    if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos) {
        throw std::runtime_error(std::string("invalid ") + label);
    }
    return std::stoull(value);
}

} // namespace tii249_square

int main(int argc, char **argv) {
    try {
        if (argc >= 2 && std::string(argv[1]) == "segment") {
            if (argc < 6) {
                std::cerr
                    << "usage: " << argv[0]
                    << " segment REQUEST.txt STATE.bin|- STEPS OUTPUT_DIR"
                    << " [BATCH_WORDS=2] [SEED_Y=24920260930]"
                    << " [COMPRESSOR_SEED]\n";
                return 2;
            }
            const u64 steps = tii249_square::parse_u64_arg(argv[4], "step count");
            int batch_words = 2;
            if (argc >= 7) {
                batch_words = static_cast<int>(tii249_square::parse_u64_arg(argv[6], "batch words"));
            }
            u64 seed_y = 24920260930ULL;
            if (argc >= 8) {
                seed_y = tii249_square::parse_u64_arg(argv[7], "seed_y");
            }
            u64 compressor_seed = tii249_square::DEFAULT_COMPRESSOR_SEED;
            if (argc >= 9) {
                compressor_seed = tii249_square::parse_u64_arg(argv[8], "compressor seed");
            }
            tii249_square::run_segment(
                argv[2], argv[3], steps, argv[5], batch_words, seed_y, compressor_seed);
            return 0;
        }
        if (argc >= 2 && std::string(argv[1]) == "recover") {
            if (argc < 5) {
                std::cerr
                    << "usage: " << argv[0]
                    << " recover REQUEST.txt CADO_FFILE OUTPUT_DIR"
                    << " [BATCH_WORDS=2] [SEED_Y=24920260930]"
                    << " [COMPRESSOR_SEED] [TRANSPOSE=1] [OUTPUT_COUNT=0(max independent)]\n";
                return 2;
            }
            int batch_words = 2;
            if (argc >= 6) {
                batch_words = static_cast<int>(tii249_square::parse_u64_arg(argv[5], "batch words"));
            }
            u64 seed_y = 24920260930ULL;
            if (argc >= 7) {
                seed_y = tii249_square::parse_u64_arg(argv[6], "seed_y");
            }
            u64 compressor_seed = tii249_square::DEFAULT_COMPRESSOR_SEED;
            if (argc >= 8) {
                compressor_seed = tii249_square::parse_u64_arg(argv[7], "compressor seed");
            }
            bool transpose = true; // CADO lingen_output_to_singlefile writes transposed F.
            if (argc >= 9) {
                const u64 value = tii249_square::parse_u64_arg(argv[8], "transpose flag");
                if (value > 1) throw std::runtime_error("transpose flag must be 0 or 1");
                transpose = value != 0;
            }
            u64 output_count = 0;
            if (argc >= 10) {
                output_count = tii249_square::parse_u64_arg(argv[9], "output count");
            }
            tii249_square::run_recovery(
                argv[2], argv[3], argv[4], batch_words, seed_y,
                compressor_seed, transpose, output_count);
            return 0;
        }
        std::cerr
            << "usage: " << argv[0]
            << " segment REQUEST.txt STATE.bin|- STEPS OUTPUT_DIR [BATCH_WORDS=2]"
            << " [SEED_Y=24920260930] [COMPRESSOR_SEED]\n"
            << "   or: " << argv[0]
            << " recover REQUEST.txt CADO_FFILE OUTPUT_DIR [BATCH_WORDS=2]"
            << " [SEED_Y=24920260930] [COMPRESSOR_SEED] [TRANSPOSE=1] [OUTPUT_COUNT=0(max independent)]\n";
        return 2;
    } catch (const std::exception &error) {
        std::cerr << "tii249 square worker refused: " << error.what() << "\n";
        return 1;
    }
}
