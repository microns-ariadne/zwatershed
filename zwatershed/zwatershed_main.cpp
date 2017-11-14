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

/*********************************
 *
 * zwshed_initial_c - run the watershed and make the region graph
 *
 * affs - an array laid out as Z varying first, then Y, then X and on the
 *        outside, the channel # (0 for X, 1 for Y, 2 for Z).
 * dimX - the dimensions along the X axis
 * dimY - the dimensions along the Y axis
 * dimZ - the dimensions along the Z axis
 * LOW - the value indicating "definitely not connected"
 * HIGH - the value indicating "definitely connected"
 *********************************/
ZWShedResult zwshed_initial_c(const int dimX, const int dimY, const int dimZ, uint8_t* affs, uint8_t LOW, uint8_t HIGH)
{
    /*
     * Some ideas taken from:
     * https://stackoverflow.com/questions/2168082/how-to-rewrite-array-from-row-order-to-column-order
     */
    typedef boost::general_storage_order<4> Storage;
    std::cout << "calculating basic watershed..." << std::endl;

    ZWShedResult result;
    /*
     * This transposes the array so that the channel is the last dimension
     * instead of the first.
     */
    int ordering[] = {2, 1, 0, 3};
    bool ascending[] = {true, true, true, true};
    affinity_graph_ptr<uint8_t> aff(new affinity_graph<uint8_t>
                              (affs,
			       boost::extents[dimX][dimY][dimZ][3],
                               Storage(ordering, ascending)));
    std::tie(result.seg_ref , result.counts_ref) = 
        watershed<seg_t>(aff, LOW, HIGH);


    // calculate region graph
    std::cout << "calculating rgn graph..." << std::endl;
    result.rg = get_region_graph(aff, result.seg_ref , result.counts_ref.size()-1);

    return result;
 }


void merge_no_stats(int dimX, int dimY, int dimZ, 
                    ZWShedResult &result, int thresh,
		    seg_t *seg_ptr){
    std::cout << "evaluating..." << std::endl;

    // read data
    std::vector<std::size_t> &counts = result.counts_ref;

    // merge
    std::cout << "thresh: " << thresh << "\n";
    double t = (double) thresh;
	merge_segments_with_function(
	    result.seg_ref, result.rg, counts, square(t), 10,RECREATE_RG);

    //
    // The segmentation is in format X, Y, Z and the numpy memory
    // is in Z, Y, X
    //
    for (size_t i = 0; i < dimX * dimY * dimZ; i++)
	seg_ptr[i] = result.seg_ref->data()[i];
}
