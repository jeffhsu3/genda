import unittest, os
import pandas as pd

from genda.formats import Genotype

class TestGenotype(unittest.TestCase):
    """ Testing loading VCF files into pandas
    """

    def setUp(self):
        pass

    def testSamples(self):
        self.assertEqual(max(Genotype.chi2_association(pd.DataFrame([[0,0,2],[2,2,1],[1,2,0]]),\
                pd.DataFrame([[2,1,1],[2,1,0],[1,1,1]])))[1], 0.90560382312624199)

    def testInfo(self):
        pass



def main():
    unittest.main()

if __name__ == '__main__':
    print(__file__)
    main()


