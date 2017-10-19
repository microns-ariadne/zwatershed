#ifndef TESTLIB_H
#define TESTLIB_H

#include <iostream>
#include <list>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include "zwatershed_util/types.hpp"

struct ZWShedResult {
    std::vector<uint64_t> edge_1;
    std::vector<uint64_t> edge_2;
    std::vector<float> weight;
    volume_ptr<uint64_t> seg_ref;
    std::vector<std::size_t> counts_ref;
};

ZWShedResult zwshed_initial_c(const int dx, const int dy, const int dz, float* affs);

void merge_no_stats(
    int dimX, int dimY, int dimZ, 
    ZWShedResult &result, int thresh, uint64_t *seg_ptr);


#endif