#pragma once
#include <cstdint>
extern "C" {
int local_rows(int k,int d,std::uint64_t point,std::uint8_t* output);
int series_products(int k,int d,int m,int modulus,int depth,const std::uint8_t* arcs,std::uint8_t* output);
void* series_cache_create(int k,int d,int m,int modulus,int depth);
void series_cache_free(void* ptr);
int series_cache_set(void* ptr,int r,const std::uint8_t* values);
int series_cache_eval(void* ptr,int r,std::uint8_t* output);
}
