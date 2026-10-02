CXX ?= c++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra
ifdef M4RI_OMP
M4RI_CFLAGS := -isystem $(M4RI_OMP)/src/m4ri-master
M4RI_LIBS := $(M4RI_OMP)/libm4ri_omp.a -fopenmp
CXXFLAGS += -fopenmp
else
ifdef M4RI_PREFIX
M4RI_CFLAGS := -isystem $(M4RI_PREFIX)/include
M4RI_LIBS := -L$(M4RI_PREFIX)/lib -Wl,-rpath,$(M4RI_PREFIX)/lib -lm4ri
else
M4RI_CFLAGS := $(subst @SIMD_CFLAGS@,,$(shell pkg-config --cflags m4ri 2>/dev/null))
M4RI_LIBS := $(shell pkg-config --libs m4ri 2>/dev/null)
endif
endif
.PHONY: all check
all: build/holdout build/wild-goppa
build:
	mkdir -p build
build/holdout: src/holdout.cpp include/algebra.hpp include/public.hpp | build
	$(CXX) $(CXXFLAGS) $(M4RI_CFLAGS) -Iinclude $< -o $@ $(M4RI_LIBS)
build/wild-goppa: src/wild_goppa.cpp src/series.cpp include/algebra.hpp include/public.hpp include/series_api.hpp | build
	$(CXX) $(CXXFLAGS) $(M4RI_CFLAGS) -Iinclude src/wild_goppa.cpp src/series.cpp -o $@ $(M4RI_LIBS)
check: build/wild-goppa
	./build/wild-goppa self-test
