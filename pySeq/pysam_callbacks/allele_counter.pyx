import numpy as np

cimport numpy as np
DTYPE = np.int

ctypedef np.int_t DTYPE_t

class AlleleCounter():
    """ Counts the various alleles at a given position.

    usage:
    import pysam
    t = AlleleCounter("chr1", snp-position - 1)

    bamfile.fetch("chr1", 124415155, 124415156)
    t.counts = [A counts, C counts, G counts, T counts]
    """
    BASE_INDEX = {'A':0, 'a':0, 'C':1, 'c':1, 'G':2, 'g':2, 'T':3, 't':3}
    def __init__(self, region, position, phredThreshold=0, isIndel = False):
        self.region = region
        self.position = int(position)
        self.phredThreshold = phredThreshold
        self.counts = np.zeros(4, dtype=np.uint32)

    def __call__(self, alignment, int position = None, 
                 int phredThreshold = None,
                 isIndel = None):

        cdef int index
        cdef int inserts

        # Settting up variables
        if position == None:
            position = self.position
        else: pass
        if phredThreshold == None:
            qualT = self.phredThreshold
        else: pass
        if alignment.is_duplicate or alignment.mapq <= 50:
            pass

        else:
            # :TODO Too many nestings, need to refactor
            if isIndel:
                pass
            else:
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
