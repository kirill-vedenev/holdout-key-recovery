SAGE ?= sage
SAGE_LOCAL ?= $(shell $(SAGE) -python -c 'from sage.env import SAGE_LOCAL; print(SAGE_LOCAL)')
CXXFLAGS = -O3 -std=c++17 -Wno-expansion-to-defined -I$(SAGE_LOCAL)/include
LDFLAGS = -L$(SAGE_LOCAL)/lib -Wl,-rpath,$(SAGE_LOCAL)/lib -lm4ri
all: holdout libseries.so libincremental.so
holdout: native.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)
libseries.so: native.cpp
	$(CXX) $(CXXFLAGS) -DSHARED_ONLY -shared -fPIC $< -o $@ $(LDFLAGS)
libincremental.so: native.cpp
	$(CXX) $(CXXFLAGS) -DSHARED_ONLY -shared -fPIC $< -o $@ $(LDFLAGS)
