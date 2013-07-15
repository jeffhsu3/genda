import unittest, os

from genda.formats.Snp_array import SNP_array

class TestLoadingSingleColumnData(unittest.TestCase):
    """ Testing loading one column SNP array files into pandas
    """

    def setUp(self):
        self.filename = './data/one-column-test-data'
        self.array = SNP_array(self.filename, fileformat = 'one column', delim = '\t', encoding = {'rs4477212':'A/G','rs3094315':'A/G',\
                'rs3131972':'G/A','rs12124819':'A/G','rs11240777':'A/G','rs6681049':'C/T','rs4970383':'T/C','rs4475691':'T/C',\
                'rs7537756':'A/T'})

    def testSamples(self):
        self.assertEqual(self.array.df.shape, (9,4))
        self.assertEqual(self.array.df.ix[0,3],'AA')
        self.assertEqual(self.array.geno.shape,(9,1))
        self.assertEqual(self.array.geno.ix[2,0], 0)
        self.assertEqual(self.array.geno.index[3], 'rs12124819')

    def testInfo(self):
        pass


class TestLoadingDoubleColumnData(unittest.TestCase):
    """ Testing loading two column files into pandas
    """

    def setUp(self):
        self.filename = './data/two-column-test-data'
        self.array = SNP_array(self.filename, fileformat = 'two column', delim = '\t',\
                encoding = {'rs4477212':'A/G','rs3094315':'A/G','rs3131972':'G/A','rs12124819':'A/G','rs11240777':'A/G',\
                'rs6681049':'C/T','rs4970383':'T/C','rs4475691':'T/C','rs7537756':'A/T'})

    def testSamples(self):
        self.assertEqual(self.array.df.shape, (9,8))
        self.assertEqual(self.array.df.ix[0,1],'rs4477212')
        self.assertEqual(self.array.geno.shape,(9,2))
        self.assertEqual(self.array.geno.ix[2,0], 0)
        self.assertEqual(self.array.geno.index[3], 'rs12124819')

    def testInfo(self):
        pass

def main():
    unittest.main()

if __name__ == '__main__':
    print(__file__)
    main()


