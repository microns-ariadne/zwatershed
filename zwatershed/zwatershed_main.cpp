/* Connected components
 * developed and maintained by Srinivas C. Turaga <sturaga@mit.edu>
 * do not distribute without permission.
 */
#include "zwatershed.h"
//#pragma once
#include "zwatershed_util/agglomeration.hpp"
#include "zwatershed_util/region_graph.hpp"
#include "zwatershed_util/basic_watershed.hpp"
#include "zwatershed_util/limit_functions.hpp"
#include "zwatershed_util/types.hpp"
#include "zwatershed_util/main_helper.hpp"

#include <memory>
#include <type_traits>

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <vector>
#include <algorithm>
#include <tuple>
#include <map>
#include <list>
#include <set>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>
#include <boost/make_shared.hpp>
using namespace std;
bool RECREATE_RG = true;
using seg_t = uint64_t;

ZWShedResult zwshed_initial_c(const int dimX, const int dimY, const int dimZ, uint8_t* affs, uint8_t LOW, uint8_t HIGH)
{
    std::cout << "calculating basic watershed..." << std::endl;

    // read data
    ZWShedResult result;
    affinity_graph_ptr<uint8_t> aff(new affinity_graph<uint8_t>
                              (affs,
			       boost::extents[dimX][dimY][dimZ][3],
                               boost::c_storage_order()));
    std::tie(result.seg_ref , result.counts_ref) = 
        watershed<seg_t>(aff, LOW, HIGH);


    // calculate region graph
    std::cout << "calculating rgn graph..." << std::endl;
    auto rg = get_region_graph(aff, result.seg_ref , result.counts_ref.size()-1);

    // save and return
    for ( const auto& e: *rg ){
        result.edge_1.push_back(std::get<1>(e));
        result.edge_2.push_back(std::get<2>(e));
        result.weight.push_back(std::get<0>(e));
    }
    return result;
 }


void merge_no_stats(int dimX, int dimY, int dimZ, 
                    ZWShedResult &result, int thresh,
		    seg_t *seg_ptr){
    std::cout << "evaluating..." << std::endl;

    // read data
    std::vector<std::size_t> &counts = result.counts_ref;
    region_graph_ptr<seg_t, uint8_t> rg( new region_graph<seg_t, uint8_t> );
    for(int i=0;i<result.weight.size();i++)
	(*rg).emplace_back(result.weight[i],result.edge_1[i],result.edge_2[i]);

    // merge
    std::cout << "thresh: " << thresh << "\n";
    double t = (double) thresh;
	merge_segments_with_function(
	    result.seg_ref, rg, counts, square(t), 10,RECREATE_RG);

    //
    // The segmentation is in format X, Y, Z and the numpy memory
    // is in Z, Y, X
    //
    for (size_t i = 0; i < dimX * dimY * dimZ; i++)
	seg_ptr[i] = result.seg_ref->data()[i];
}
