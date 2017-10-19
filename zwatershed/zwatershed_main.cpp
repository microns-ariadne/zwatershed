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
// arb funcs
#include "zwatershed_util/region_graph_arb.hpp"
#include "zwatershed_util/basic_watershed_arb.hpp"
#include "zwatershed_util/main_helper_arb.hpp"


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
// these values based on 5% at iter = 10000
double LOW=  .0001;
double HIGH= .9999;
bool RECREATE_RG = true;

ZWShedResult zwshed_initial_c(const int dimX, const int dimY, const int dimZ, float* affs)
{
    std::cout << "calculating basic watershed..." << std::endl;

    // read data
    ZWShedResult result;
    affinity_graph_ptr<float> aff(new affinity_graph<float>
                              (boost::extents[dimX][dimY][dimZ][3],
                               boost::fortran_storage_order()));
    for(int i=0;i<dimX*dimY*dimZ*3;i++)
        aff->data()[i] = affs[i];
    std::tie(result.seg_ref , result.counts_ref) = 
        watershed<uint64_t>(aff, LOW, HIGH);


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
		    uint64_t *seg_ptr){
    std::cout << "evaluating..." << std::endl;

    // read data
    std::vector<std::size_t> &counts = result.counts_ref;
    region_graph_ptr<uint64_t,float> rg( new region_graph<uint64_t,float> );
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
    for (size_t z = 0; z < dimZ; z++)
	for (size_t y = 0; y < dimY; y++)
	    for (size_t x = 0; x < dimX; x++)
		seg_ptr[(z * dimY + y) * dimX + x] =
		    result.seg_ref->data()[(x * dimY + y) * dimZ + z];
}
