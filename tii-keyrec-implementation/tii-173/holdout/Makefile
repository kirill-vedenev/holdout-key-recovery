CXX ?= c++
NVCC ?= nvcc
CUDA_ARCH ?= 90
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic -pthread
NVCCFLAGS ?= -O3 -std=c++17 -lineinfo -Xcompiler=-pthread

.PHONY: all cpu cuda
all: cpu
cpu: build/holdout-cpu
cuda: build/holdout-cuda

build:
	mkdir -p build

build/holdout-cpu: src/cpu.cpp src/common.hpp src/engine.hpp third_party/saarinen/colex.hpp | build
	$(CXX) $(CXXFLAGS) src/cpu.cpp -o $@

build/holdout-cuda: src/cuda.cu src/common.hpp src/engine.hpp third_party/saarinen/colex.hpp | build
	$(NVCC) $(NVCCFLAGS) -gencode arch=compute_$(CUDA_ARCH),code=sm_$(CUDA_ARCH) src/cuda.cu -o $@
