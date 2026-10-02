CXX ?= c++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -pthread

.PHONY: all check clean
all: build/product_rank

build:
	mkdir -p build

build/product_rank: src/product_rank.cpp | build
	$(CXX) $(CXXFLAGS) $< -o $@

check: build/product_rank
	./build/product_rank self-test

clean:
	rm -rf build
