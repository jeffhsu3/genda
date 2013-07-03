import unittest

import numpy as np

from genda.pysam_callbacks.allele_counter import AlleleCounter

class FakeAlignment(object):
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
    def setUp(self):
        self.seq = "ATTAGGATAGATGATGATGA"

    def testRegular(self):
        regular_read = FakeAlignment(25, self.seq, [(0,20)])
        # Position should be an A
        test = AlleleCounter("chr1", 35)
        test(regular_read)
        np.testing.assert_equal(test.counts, np.asarray([0,0,1,0]))

    def testSNPinIntron(self):
        snp_in_intron = FakeAlignment(25, self.seq, [(0,10), (3,100), (0,10)])
        test = AlleleCounter("chr1", 35)
        test(snp_in_intron)
        np.testing.assert_equal(test.counts, np.asarray([0,0,0,0]))

    def testSecondExon(self):
        snp_in_second_exon = FakeAlignment(25, self.seq,[(0,10),(3,5),(0,10)])
        test = AlleleCounter("chr1", 47)
        test(snp_in_second_exon)
        np.testing.assert_equal(test.counts, np.asarray([1,0,0,0]))

    """
    Not yet implemented
    def testIndels(self):
        indel = FakeAlignment(25, self.seq, [(0,20), (3, 45), (0,10)])
        pass
    """



if __name__ == '__main__':
    unittest.main()
