import unittest, os

from pySeq.formats.panVCF import VCF

class TestLoadingVCF(unittest.TestCase):
    """ Testing loading VCF files into pandas
    """

    def setUp(self):
        self.filename = './data/chrY.test.vcf'
        self.VCF = VCF(self.filename)

    def testSamples(self):
        self.assertEqual(len(self.VCF.samples), 526)

    def testInfo(self):
        pass



def main():
    unittest.main()

if __name__ == '__main__':
    print(__file__)
    main()


