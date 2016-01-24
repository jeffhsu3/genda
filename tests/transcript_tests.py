import unittest
import pysam
from genda.transcripts import (Exon, Gene, Transcript, unique_sets,
        compare_two_transcripts)

class TestGeneBreaks(unittest.TestCase):

    def setUp(self):
        self.simple_gene = Gene([Transcript([Exon('chr1', 10, 20),
                                             Exon('chr1', 50, 60),
                                             Exon('chr1', 80, 120)]),
                                 Transcript([Exon('chr1', 10, 20),
                                             Exon('chr1', 80, 120)]),
                                 Transcript([Exon('chr1', 5, 20),
                                             Exon('chr1', 50, 60),
                                             Exon('chr1', 80,120)])])
        self.seperate_genes = 0

        self.read_length = 75

    def uniqueSets(self):
        """ Test to calculate the unique regions such that, given an
        annotation, regions are unique iff a set of s.  Note,
        regions defined in such a manner can be disjoint.

        """
        pass


class TestGene(unittest.TestCase):

    def setUp(self):
        self.simple_gene = Gene('ENSG01', chrom=2, start=0, end=20,
        symbol='HSU')

    def test_get_transcript(self):    
        gtf = 'bleh' # Need to generate a test gtf file




class TestGenerateSets(unittest.TestCase):

    def setUp(self):
        self.set_a = set(['a', 'b', 'c'])
        self.set_b = set(['a', 'b'])
        self.set_c = set(['a', 'b', 'd'])


    def find_sets(self):
        self.assertEqual(unique_sets([self.set_a, self.set_b, self.set_c]),
                                     [set(['a', 'b']), set(['c']), set(['d'])])



class TestCompareTwoTranscripts(unittest.TestCase):
    """:TODO add more test cases
    """
    def setUp(self):
        self.transcript_dict = {
                't1' : [(10, 20, 1), 
                    (50, 60, 2), 
                    (70, 90, 3),
                    (120, 180, 4),],
                't2' : [(12, 20, 1), 
                    (70, 90, 2)]
                }

    def test_compare_transcripts_skipped_exon_simple(self):
        exclusive_juncs, torder, matching_exons, skipped_exons =\
                compare_two_transcripts('t1', 't2', self.transcript_dict)
        se = skipped_exons[0]
        self.assertEqual((50, 60), (se.start, se.end))
        self.assertEqual('skipped_exon', se.event_type)
        self.assertEqual((2, None), se.exon_num)
        


if __name__ == '__main__':
    unittest.main()
