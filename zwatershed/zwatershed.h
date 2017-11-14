#ifndef TESTLIB_H
#define TESTLIB_H

#include <iostream>
#include <list>
#include <string>
#include <map>
#include <tuple>
#include <vector>
#include <utility>
#include "zwatershed_util/types.hpp"


struct ZWShedResult {
    region_graph_ptr<uint64_t, uint8_t> rg;
    volume_ptr<uint64_t> seg_ref;
    std::vector<std::size_t> counts_ref;
    uint8_t getWeight(size_t idx) { return std::get<0>((*rg)[idx]); }
    uint64_t getEdge1(size_t idx) { return std::get<1>((*rg)[idx]); }
    uint64_t getEdge2(size_t idx) { return std::get<2>((*rg)[idx]); }
    size_t rg_size() { return rg->size(); }
};

ZWShedResult zwshed_initial_c(const int dx, const int dy, const int dz, uint8_t* affs, uint8_t LOW, uint8_t HIGH);

void merge_no_stats(
    int dimX, int dimY, int dimZ, 
    ZWShedResult &result, int thresh, uint64_t *seg_ptr);


#endif