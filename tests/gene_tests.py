
import unittest
import pysam
from genda.transcripts import Gene



class TestGene(unittest.TestCase):

    def setUp(self):
        self.simple_gene = Gene('ENSG01', chrom=2, start=0, end=20,
        symbol='HSU')

    def test_get_transcript(self):    
        gtf = pysam.Tabixfile('') # Need to generate a test gtf file
