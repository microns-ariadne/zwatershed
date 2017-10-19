from libcpp.list cimport list
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.pair cimport pair
from libc.stdint cimport uint64_t
import numpy as np
import os
cimport numpy as np
import h5py

#-------------- interface methods --------------------------------------------------------------
def zwatershed(affs, threshes):
    threshes.sort()
    return zwshed_no_stats(affs, threshes)

def zwshed_no_stats(np.ndarray[np.float32_t, ndim=4] affs, threshes):
    # get initial seg,rg
    affs = np.asfortranarray(np.transpose(affs, (1, 2, 3, 0)))
    dims = affs.shape
    cdef ZWShedResult result = zwshed_initial_c(
        dims[0], dims[1], dims[2], <float *>affs.data)
    cdef np.ndarray[np.uint64_t, ndim=3] seg
    segs = []

    # get segs, stats
    for i in range(len(threshes)):
        seg = \
            np.zeros((dims[0], dims[1], dims[2]), np.uint64)
        if(result.edge_1.size() > 0):
            map = merge_no_stats(
            dims[0], dims[1], dims[2], result, threshes[i], 
            <uint64_t *>(seg.data))
        segs.append(seg)
    return segs

#-------------- c++ methods --------------------------------------------------------------
cdef extern from "zwatershed.h":
    struct ZWShedResult:
        vector[uint64_t] edge_1
        vector[uint64_t] edge_2
        vector[float] weights
        
    ZWShedResult zwshed_initial_c(int dimX, int dimY, int dimZ, np.float32_t*affs)
    void merge_no_stats(int dx, int dy, int dz,
                        ZWShedResult &result, int thresh, uint64_t *seg_ptr)
