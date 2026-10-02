CXX ?= c++
NVCC ?= nvcc
CUDA_ARCH ?= 90
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic -pthread
NVCCFLAGS ?= -O3 -std=c++17 -lineinfo -Xcompiler=-pthread

.PHONY: all cpu gpu
all: cpu
cpu: build/holdout-cpu
gpu: build/square-worker

build:
	mkdir -p build

build/holdout-cpu: src/cpu.cpp src/common.hpp src/engine.hpp third_party/saarinen/colex.hpp | build
	$(CXX) $(CXXFLAGS) src/cpu.cpp -o $@

build/square-worker: src/square_cufft_worker.cu third_party/saarinen/holdout_cuda_worker.cu | build
	$(NVCC) $(NVCCFLAGS) -gencode arch=compute_$(CUDA_ARCH),code=sm_$(CUDA_ARCH) src/square_cufft_worker.cu -lcufft -o $@
