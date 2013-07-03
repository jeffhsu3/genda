import numpy as np

cimport numpy as np
DTYPE = np.int

from cpython cimport bool

#cimport pysam.csamtools as csam
#from pysam.csamtools cimport AlignedRead
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

        # Settting up variables
        position = self.position
        qualT = self.phredThreshold
        """
        if alignment.is_duplicate or alignment.mapq <= 50:
            pass
        """
        print(alignment.cigar)
        try:
            inserts = 0
            for i in alignment.cigar:
                inserts += i[0] * 1/3 * i[1]
            index_b = position - alignment.pos
            index_a = position - inserts - alignment.pos - 1
            print(index_b)
            print(index_a)
            index = int_max(index_b, index_a)
        except IndexError:
            index = position - alignment.pos - 1
        try:
            base = alignment.seq[index]
            iindex = self.BASE_INDEX[base]
            self.counts[iindex] += 1
        except IndexError:
            pass


        """
        if 3 in [i[0] for i in alignment.cigar]:
            t = [i[1] for i in alignment.cigar if i[0] == 3]
            inserts = sum(t)
            index = position - inserts - alignment.pos - 1
        else:
            index = position - alignment.pos - 1
        if index >= 0:
            base = alignment.seq[index]
            b_qual = alignment.qual[index]
            if base != "N" and ord(b_qual)-33 > qualT:
                base =  self.BASE_INDEX[base]
                self.counts[base] += 1
            else: pass
        else: pass
        """
