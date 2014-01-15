import numpy as np

cimport numpy as np
DTYPE = np.int

from cpython cimport bool

cimport pysam.csamtools as csam
from pysam.csamtools cimport AlignedRead
#cimport csamtools

ctypedef np.int_t DTYPE_t

cdef inline int int_max(int i , int j): return i if i >= j else j



class AlleleCounter(object):
    """ Counts the various alleles at a given position.

    usage:
    import pysam
    t = AlleleCounter("chr1", snp-position - 1)

    bamfile.fetch("chr1", 124415155, 124415156)
    t.counts = [A counts, C counts, G counts, T counts]
    """
    BASE_INDEX = {'A':0, 'a':0, 'C':1, 'c':1, 'G':2, 'g':2, 'T':3, 't':3}
    def __init__(self, region, int position, int phredThreshold=0, bool isIndel = False):
        self.region = region
        self.position = int(position)
        self.phredThreshold = phredThreshold
        self.counts = np.zeros(4, dtype=np.uint8)

    def  __call__(self, alignment, int position = 0,
                 int phredThreshold = 0,
                 bool isIndel = False):

        cdef int index
        cdef int inserts
        cdef char *b_qual
        cdef char *base
        cdef int iindex
        cdef int i_s
        cdef int i_e

        # Settting up variables
        position = self.position
        qualT = self.phredThreshold
        """
        if alignment.is_duplicate or alignment.mapq <= 50:
            pass
        """
        inserts = 0
        i_e = alignment.pos
        i_s = alignment.pos
        for i in alignment.cigar:
            inserts += i[0] * 1/3 * i[1]
            i_e += i[1]  
            if i_s <= position < i_e and i[0] == 3:
                break
            elif i_s <= position < i_e:
                index = position - inserts - alignment.pos
                try:
                    iindex = self.BASE_INDEX[alignment.seq[index]]
                    self.counts[iindex] += 1
                except IndexError:
                    pass
            else:
                pass
            i_s = i_e
