import numpy as np
cimport numpy as np
ctypedef np.float64_t dtype_t 
ctypedef np.uint16_t DTYPE_t
from genda.pysam_callbacks.allele_counter cimport AlleleCounter
from cpython cimport bool
from pysam.calignmentfile cimport AlignmentFile, IteratorRowRegion, AlignedSegment
import cython
from cython.parallel import prange
from . import *

STUFF = 'hi'

@cython.boundscheck(False)
cpdef aei_count_samples(np.ndarray[dtype_t, ndim=2] geno,
        np.ndarray geno_columns,
        np.ndarray[DTYPE_t, ndim=2] c_m, 
        np.ndarray chrom, 
        np.ndarray[DTYPE_t, ndim=1] pos,
        bool reduced = False):
    """ A parallel cython version, should make this part of a class
    """
    cdef int n_samples
    cdef int counter
    cdef AlignmentFile bamfile
    cdef int n_markers
    cdef int j 
    cdef int i
    cdef char * chrm
    cdef int p
    cdef IteratorRowRegion temp
    #cdef char *j
    #cdef int k 
    cdef AlleleCounter variant 
    n_samples = len(geno_columns)
    n_markers = len(chrom)
    # Assume positions are sorted
    #index = sample_genotype.index
    #for i in prange(n_samples, nogil=False):
    for i in xrange(n_samples):
    #for i in prange(n_samples, nogil=True):
        bamfile = Samfile(geno_columns[i], 'rb')
        counter = 0
        #for i, j, k in zip(chrom, pos, index):
        for j in xrange(n_markers):
            print(chrom[j], pos[j])
            variant = AlleleCounter(chrom[j], pos[j], phredThreshold=0)
            bamfile.fetch(variant.region, variant.position, variant.position +
                    1, callback=variant)
            #print(variant.A_n, variant.C_n, variant.G_n, variant.T_n)
            c_m[counter, i * 4: i * 4 + 4] = np.asarray([variant.A_n,
                variant.C_n, variant.G_n, variant.T_n]).T
            #print(c_m[counter, i * 4: i * 4 + 4])
            counter += 1
    return(c_m)
