import unittest, os
import numpy as np
import pandas as pd

from pySeq.formats.PED import PED

class TestLoadingSingleColumnData(unittest.TestCase):
    """ Testing loading VCF files into pandas
    """

    def setUp(self):
        self.pedfile = './data/test.ped'
        self.mapfile = './data/test.map'
        self.encoder = {'snp1':'A/C','snp2':'A/C','snp3':'C/A','snp4':'T/G','snp5':'C/A'}
        self.PED = PED(self.pedfile, self.mapfile, encoder = self.encoder)

    def testSamples(self):
        self.assertEqual(self.PED.geno.shape, (5,4))
        self.assertEqual(self.PED.geno.ix[0,0],0)
        self.assertEqual(np.isnan(self.PED.geno.ix[0,3]), True)
        self.assertEqual(self.PED.hardyweinberg('snp1'), True)

    def testInfo(self):
        pass

def main():
    unittest.main()

if __name__ == '__main__':
    print(__file__)
    main()


