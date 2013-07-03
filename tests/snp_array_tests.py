import unittest, os

from genda.formats.Snp_array import SNP_array

class TestLoadingSingleColumnData(unittest.TestCase):
    """ Testing loading VCF files into pandas
    """

    def setUp(self):
        self.filename = './data/one-column-test-data'
        self.array = SNP_array(self.filename, fileformat = 'one column', delim = '\t')

    def testSamples(self):
        self.assertEqual(self.array.df.shape, (9,4))
        self.assertEqual(self.array.df.ix[0,3],'AA')
        self.assertEqual(self.array.geno.shape,(9,1))
        self.assertEqual(self.array.geno.ix[2,0], 2.0)
        self.assertEqual(self.array.geno.index[3], 'rs12124819')

    def testInfo(self):
        pass


class TestLoadingDoubleColumnData(unittest.TestCase):
    """ Testing loading VCF files into pandas
    """

    def setUp(self):
        self.filename = './data/two-column-test-data'
        self.array = SNP_array(self.filename, fileformat = 'two column', delim = '\t')

    def testSamples(self):
        self.assertEqual(self.array.df.shape, (9,4))
        self.assertEqual(self.array.df.ix[0,3],'AA')
        self.assertEqual(self.array.geno.shape,(9,1))
        self.assertEqual(self.array.geno.ix[2,0], 2.0)
        self.assertEqual(self.array.geno.index[3], 'rs12124819')

    def testInfo(self):
        pass

def main():
    unittest.main()

if __name__ == '__main__':
    print(__file__)
    main()


