from libcpp.algorithm cimport copy
from libcpp.list cimport list
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.pair cimport pair
from libc.stdint cimport uint64_t, uint32_t, uint8_t
import numpy as np
cimport numpy as np

#-------------- interface methods --------------------------------------------------------------
def zwatershed(np.ndarray[np.uint8_t, ndim=4] affs, threshes, 
               uint8_t LOW=1, uint8_t HIGH=254, return_rg=False):
    '''Produce z-watersheds at one or more thresholds
    :param affs: a numpy array of UINT8 adjacent-voxel affinities (255 = merge,
    0 = split) arranged as a c, z, y, x array (c = 0 for z, 1 for y, 2 for x)
    :param threshes: a list of thresholds at which to evaluate the watershed.
    :param LOW: always exclude edges with this predicted value or lower
    :param HIGH: always join edges with this predicted value or higher
    :returns: a list of segmentations as UINT32 numpy arrays
    '''
    threshes.sort()
    affs = np.ascontiguousarray(affs)
    dims = affs.shape
    cdef ZWShedResult result
    cdef np.ndarray[np.uint64_t, ndim=1] seg
    cdef:
        np.ndarray[np.uint64_t] edge_1
        np.ndarray[np.uint64_t] edge_2
        np.ndarray[np.uint8_t] weight
    cdef int thresh
    segs = []
    with nogil:
        result = zwshed_initial_c(
            dims[1], dims[2], dims[3], <uint8_t *>affs.data, LOW, HIGH)

    # get segs, stats
    for i in range(len(threshes)):
        seg = \
            np.zeros((dims[1] * dims[2] * dims[3]), np.uint64)
        if result.edge_1.size() > 0:
            thresh = threshes[i]
            with nogil:
                merge_no_stats(
                    dims[1], dims[2], dims[3], result, thresh, 
                    <uint64_t *>(seg.data))
        segs.append(seg.reshape(dims[3], dims[2], dims[1]).transpose(2, 1, 0))
        del seg
    if return_rg:
        edge_1 = np.zeros(result.edge_1.size(), np.uint64)
        edge_2 = np.zeros(result.edge_2.size(), np.uint64)
        weight = np.zeros(result.weight.size(), np.uint8)
        copy(result.edge_1.begin(), result.edge_1.end(),
             <uint64_t *>edge_1.data)
        copy(result.edge_2.begin(), result.edge_2.end(),
             <uint64_t *>edge_2.data)
        copy(result.weight.begin(), result.weight.end(), 
             <uint8_t *>weight.data)
        return segs, (edge_1, edge_2, weight)
    return segs

#-------------- c++ methods --------------------------------------------------------------
cdef extern from "zwatershed.h" nogil:
    struct ZWShedResult:
        vector[uint64_t] edge_1
        vector[uint64_t] edge_2
        vector[uint8_t] weight
        
    ZWShedResult zwshed_initial_c(
        int dimX, int dimY, int dimZ, 
        np.uint8_t*affs, np.uint8_t LOW, np.uint8_t HIGH)
    void merge_no_stats(int dx, int dy, int dz,
                        ZWShedResult &result, int thresh, uint64_t *seg_ptr)
