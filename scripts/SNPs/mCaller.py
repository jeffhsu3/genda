#!/usr/bin/python

""" Calls genotypes from target enriched regions.  Outputs a VCF with
genotyping calls.

usage:
    python mCaller.py [OPTIONS] in1.bam [in2.bam ..]
"""

import pysam

import pandas

class mPileup(object):
    def __init__(self, Samfiles):
        self.files = Samfile

def callGenotypes():
    pass

def main():
    # Right now requires intervals bed, genearlize to whole-genome later
    import optparse, sys
    p = optparse.OptionParser(__doc__)
    p.add_option("-L", "--intervals", dest="intervals", help="Target \
            Enriched Region or Region of Interest.  A Bedfile")
    p.add_option("-D", "--dbsnp", dest="G", help=\
            "Use imputed/genotypes if available, should be in VCF file format")
    p.add_option("-v",  "--vcf_file", action="store_true", dest="inputisvcfile",
                 help="the input is a VCF file")
    p.add_option("-q", "--quality_threshold", type="int", dest="qual",
                 help="base quality threshold to take allele counts from")
    p.add_option("-p", "--pileup", action="store_true", dest="p",
                 help= "Input files are pileup files")
    p.add_option("-D", "--debug", action="store_true", dest="D", help="debug")
    p.add_option("-c", "--count-threshold", action="store", type="int",
                 dest="c", help="Set the count threshold for making AEI calls")
    p.add_option("-V", "--output_vcf", action="store_true", dest="outputVCF",
                 help="Output the results to a VCF file")

    options, args = p.parse_args()

    try:
        bam_files = [pysam.Samfile(i, 'rb' for i in args)]
    except IOError:
        sys.exit()

    # Initialize Count Matrix

    if options.intervals:
        with open(options.intervals, 'rU') as intervals:
            for line in interavls:
                if line[0] == "#":
                    pass
                else:
                    call_genotypes()


if __name__ == '__main__':
    main()
