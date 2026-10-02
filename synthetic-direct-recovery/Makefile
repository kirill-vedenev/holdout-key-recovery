CXX ?= c++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra

.PHONY: all check clean
all: build/synthetic_direct

build:
	mkdir -p build

build/synthetic_direct: src/synthetic_direct.cpp | build
	$(CXX) $(CXXFLAGS) $< -o $@

check: build/synthetic_direct
	./build/synthetic_direct self-test

clean:
	rm -rf build
