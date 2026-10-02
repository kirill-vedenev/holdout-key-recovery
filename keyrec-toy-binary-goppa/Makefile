CXX ?= c++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra

ifdef M4RI_PREFIX
M4RI_CFLAGS := -isystem $(M4RI_PREFIX)/include
M4RI_LIBS := -L$(M4RI_PREFIX)/lib -Wl,-rpath,$(M4RI_PREFIX)/lib -lm4ri
else
M4RI_CFLAGS := $(shell pkg-config --cflags m4ri 2>/dev/null)
M4RI_LIBS := $(shell pkg-config --libs m4ri 2>/dev/null)
endif

.PHONY: all check
all: build/holdout build/toy-goppa

build:
	mkdir -p build

build/holdout: src/holdout.cpp | build
	$(CXX) $(CXXFLAGS) $(M4RI_CFLAGS) $< -o $@ $(M4RI_LIBS)

build/toy-goppa: src/toy_goppa.cpp src/holdout.cpp include/algebra.hpp include/series_api.hpp | build
	$(CXX) $(CXXFLAGS) $(M4RI_CFLAGS) -Iinclude -DSHARED_ONLY src/toy_goppa.cpp src/holdout.cpp -o $@ $(M4RI_LIBS)

check: build/toy-goppa
	./build/toy-goppa self-test
