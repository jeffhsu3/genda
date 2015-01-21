import numpy as np
cimport numpy as np

from cpython cimport bool
#from multiprocessing import Process
#import cython
cimport pysam.csamtools as csam
from pysam.calignmentfile cimport AlignedSegment
from libc.stdlib cimport malloc, free
from cpython cimport PyObject
#cimport csamtools
DTYPE = np.uint16
ctypedef np.uint16_t DTYPE_t
ctypedef np.float64_t dtype_t 

cdef inline int int_max(int i , int j): return i if i >= j else j


cdef class AlleleCounter:
    """ Counts the various alleles at a given position.
    usage:
    import pysam
    t = AlleleCounter(c_m, "chr1", snp-position - 1)
    bamfile.fetch("chr1", 124415155, 124415156, callback=t)
    print(t.A_n/float(t.G_n))  

    :TODO refactor out callback
    """
    #cdef DTYPE_t [:] counts
    def __cinit__(self, char *region, 
            int position, int phredThreshold=0, bool isIndel = False):
        #cdef dict BASE_INDEX
        #:TODO get rid of this dictionary make it in C?
        #self.BASE_INDEX = {'A':0, 'a':0, 'C':1, 'c':1, 'G':2, 'g':2, 'T':3, 't':3}
        #cdef np.ndarray counts = np.zeros(4, dtype=DTYPE_t)
        cdef int A_n = 0
        cdef int G_n = 0
        cdef int C_n = 0
        cdef int T_n = 0
        self.A_n = A_n
        self.G_n = G_n
        self.C_n = C_n
        self.T_n = T_n
        self.region = region
        self.position = int(position)
        self.phredThreshold = phredThreshold
        #self.counts = np.zeros(4, dtype=np.uint16)
        


    def  __call__(self, AlignedSegment alignment, int position = 0,
                 int phredThreshold = 0,
                 bool isIndel = False):
        cdef int index
        cdef int inserts
        cdef char *b_qual
        cdef char *base
        #cdef int iindex
        cdef int i_s
        cdef int i_e
        cdef char *base_pair
        #cdef np.ndarray[np.uint16, ndim=1] counts = self.counts
        #cdef bool alignment.is_duplicaate
        # Settting up variables
        position = self.position
        #qualT = self.phredThreshold
        #print(alignment.seq)
        
        #if alignment.is_duplicate:
        if False:
            pass
        else:
            inserts = 0
            i_e = alignment.pos
            i_s = alignment.pos

            for i in alignment.cigar:
                # 0 is cigar string match
                inserts += i[0] * 1/3 * i[1]
                i_e += i[1]  
                if i_s <= position < i_e and i[0] == 3:
                    break
                elif i_s <= position < i_e:
                    index = position - inserts - alignment.pos - 1 
                    # :TODO INDELS!
                    base_pair = alignment.seq[index]
                    if base_pair == b'A' or base_pair == b'a':
                        self.A_n += 1
                    elif base_pair == b'C' or base_pair == b'c':
                        self.C_n += 1
                    elif base_pair == b'G' or base_pair == b'g':
                        self.G_n += 1
                    elif base_pair == b'T' or base_pair == b't':
                        self.T_n += 1
                    else:
                        pass
                    break
                    """
                    try:
                        
                        iindex = self.BASE_INDEX[alignment.seq[index]]
                        self.counts[iindex] += 1
                        break
                    except KeyError:
                        pass
                    except IndexError:
                        pass
                    """
                else:
                    pass
                i_s = i_e


