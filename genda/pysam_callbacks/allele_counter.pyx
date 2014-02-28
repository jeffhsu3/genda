import numpy as np

cimport numpy as np

from cpython cimport bool

cimport pysam.csamtools as csam
from pysam.csamtools cimport AlignedRead
from cython.parallel cimport prange
#cimport csamtools
DTYPE = np.uint16

ctypedef np.uint16_t DTYPE_t

cdef inline int int_max(int i , int j): return i if i >= j else j


class AlleleCounter(object):
    """ Counts the various alleles at a given position.

    usage:
    import pysam
    t = AlleleCounter("chr1", snp-position - 1)

    bamfile.fetch("chr1", 124415155, 124415156)
    t.counts = [A counts, C counts, G counts, T counts]
    """
    def __init__(self, region, int position, int phredThreshold=0, bool isIndel = False):
        cdef dict BASE_INDEX
        self.BASE_INDEX = {'A':0, 'a':0, 'C':1, 'c':1, 'G':2, 'g':2, 'T':3, 't':3}
        self.region = region
        self.position = int(position)
        self.phredThreshold = phredThreshold
        self.counts = np.zeros(4, dtype=np.uint16)

    def  __call__(self, AlignedRead alignment, int position = 0,
                 int phredThreshold = 0,
                 bool isIndel = False):
        cdef int index
        cdef int inserts
        cdef char *b_qual
        cdef char *base
        cdef int iindex
        cdef int i_s
        cdef int i_e
        #cdef np.ndarray[np.uint16, ndim=1] counts = self.counts
        #cdef bool alignment.is_duplicaate

        # Settting up variables
        position = self.position
        #qualT = self.phredThreshold
        
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
                    try:
                        iindex = self.BASE_INDEX[alignment.seq[index]]
                        self.counts[iindex] += 1
                        break
                    except IndexError:
                        pass
                else:
                    pass
                i_s = i_e


'''
cpdef np.ndarray[DTYPE_t, DTYPE_t] aei_counts_samples(np.ndarray[DTYPE_t, ndim=2] c_m, bool reduced = False, 
        chrom = None, pos=None, path= None, bamfile_func = column_names_to_bams):
    """ A parallel cython version, should make this part of a class
    """
    cdef int n_samples
    n_samples = c_m.shape[0]
    sample_genotype.index
    for i in prange(n_samples, nogil=True):
        bamfile = pysam.Samfile(sample_genotype.name, 'rb')
        print(sample_genotype.name)
        for i, j, k in zip(chrom, pos, index):
            variant = AlleleCounter(str(i), j, phredThreshold=0)
            bamfile.fetch(variant.region, variant.position, variant.position +
                    1, callback=variant)
        with gil:
            c_m.ix[k, i * 4: i * 4 + 3] = variant.counts.T

    return(c_m)
'''



