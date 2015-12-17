from genda.transcripts import (Exon, Gene, Transcript, unique_sets,
        compare_two_transcripts)
import unittest

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


class TestGenerateSets(unittest.TestCase):

    def setUp(self):
        self.set_a = set(['a', 'b', 'c'])
        self.set_b = set(['a', 'b'])
        self.set_c = set(['a', 'b', 'd'])


    def find_sets(self):
        self.assertEqual(unique_sets([self.set_a, self.set_b, self.set_c]),
                                     [set(['a', 'b']), set(['c']), set(['d'])])



class TestCompareTwoTranscripts(unittest.TestCase):
    """
    """
    def setUp(self):
        self.transcript_dict = {
                't1' : [(10, 20, 1), 
                    (50, 60, 2), 
                    (70, 90, 3)],
                't2' : [(10, 20, 1), (70, 90, 2)]
                }

    def test_compare_transcripts_skipped_exon(self):
        exclusive_juncs, torder, matching_exons, skipped_exons =\
                compare_two_transcripts('t1', 't2', self.transcript_dict)
        print(torder)
        print(exclusive_juncs)
        self.assertEqual([(50, 60, (None, 2), (0, 10)),], skipped_exons)
        


if __name__ == '__main__':
    unittest.main()
