import unittest

import numpy as np

from genda.pysam_callbacks.allele_counter import AlleleCounter
from pysam.calignmentfile import AlignedSegment

def make_read(aligned_read, pos, cigar):
    aligned_read.pos = pos
    aligned_read.cigar = cigar 
    return aligned_read


class FakeAlignment(AlignedSegment):
    """ Simulating sequencing reads
    """

    def __init__(self, position, seq, cigar):
        """ Only relevent information for this test
        """
        self.pos = position
        self.seq = seq
        self.cigar = cigar
        self.is_duplicate = False
        self.mapq = 255
        self.qual = len(seq) * 'I'


class testCounter(unittest.TestCase):
    #:TODO read actual bam file 
    def setUp(self):
        # Read positions in BAM are zero-indexed while SNP positions are almost
        # always 1-indexed
        self.read = AlignedSegment()
        self.read.seq = "ATTAGGATAG" 
        self.mapq = 255
        self.qual = len(self.read.seq) * 'I'


    def testRegular(self):
        # Position should be an G
        regular_read = make_read(self.read, 24, [(0,10)])
        test = AlleleCounter("chr1", 29)
        test(regular_read)
        np.testing.assert_equal(
                np.asarray([test.A_n, test.G_n, 
            test.G_n, test.T_n]), 
                np.asarray([0,0,1,0]))


    def testSNPinIntron(self):
        snp_in_intron = make_read(self.read, 24, [(0,5), (3,5), (0,5)])
        test = AlleleCounter("chr1", 31)
        test(snp_in_intron)
        t = np.asarray([test.A_n, test.G_n, 
            test.G_n, test.T_n])
        np.testing.assert_equal(t, np.asarray([0,0,0,0]))


    def testSecondExon(self):
        snp_in_second_exon = make_read(self.read, 24, [(0,5),(3,5),(0,5)])
        test = AlleleCounter("chr1", 37)
        test(snp_in_second_exon)
        t = np.asarray([test.A_n, test.G_n, 
            test.G_n, test.T_n])
        np.testing.assert_equal(t, np.asarray([0,0,0,1]))

    """
    Not yet implemented
    def testIndels(self):
        indel = FakeAlignment(25, self.seq, [(0,20), (3, 45), (0,10)])
        pass
    """



if __name__ == '__main__':
    unittest.main()
