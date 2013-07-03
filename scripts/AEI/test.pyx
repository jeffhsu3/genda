#!/usr/bin/env python

""" Gets the read counts at each position specified by the first argument, which
is a VCF/BED/GFF/GTF file, from an indexed BAMfile/s or pileup/s and loads into a python pickle file:

Usage:
    python getReads.py bedfile/vcffile bamfile1 [bamfile2 ...] -o OUTPUT

:TODO Need to handle indels
:TODO add an option that includes a gene annotation file
:TODO add support for pileup files
:TODO also how does PCR duplicates affect?
:TODO makes this more modular and testable, especially test for changes due to pysam API changes

:TODO Output the results into a VCF file


Written by Jeffrey Hsu
2010
"""
# Standard Libraries
import optparse, pickle

import pysam
import numpy as np

cimport numpy as np
cimport csamtools as csam

DTYPE = np.int
ctypedef np.int_t DTYPE_t

class AlleleCounter():
    """ Gets the allele counts
    """
    BASE_INDEX = {'A':0, 'a':0, 'C':1, 'c':1, 'G':2, 'g':2, 'T':3, 't':3}
    def __init__(self, region, position, phredThreshold=0):
        self.region = region
        self.position = int(position)
        self.phredThreshold = phredThreshold
        self.counts = np.zeros(4, dtype=np.uint32)
        self.debugcounts = 0

    def __call__(self, alignment, position = None, phredThreshold=None):
        if position == None:
            position = self.position
        else: pass
        if phredThreshold == None:
            qualT = self.phredThreshold
        else: pass
        #if alignment.is_duplicate or alignment.mapq <= 50:
        if alignment.mapq <= 0:
            pass
        else:
            # This needs testing
            if 3 in [i[0] for i in alignment.cigar]:
                t = [i[1] for i in alignment.cigar if i[0] == 3]
                inserts = sum(t)
                # -1 to account for python 0-based index
                index = position - inserts - alignment.pos - 1
            else:
                index = position - alignment.pos - 1
            if index >= 0:
                base = alignment.seq[index]
                b_qual = alignment.qual[index]
                if base != "N" and ord(b_qual)-33 > qualT:
                    base =  self.BASE_INDEX[base]
                    self.counts[base] += 1
                else: pass
            else: pass

def main():


    #################################################################
    # Argument and Options Parsing
    #################################################################

    p = optparse.OptionParser(__doc__)
    p.add_option("-o", "--output", dest="filename", help="write \
            report to FILE")
    p.add_option("-I", "--input", dest="input",
            help="BAM files are listed in a file instead of command line\
                    arguments")
    p.add_option("-G", "--genotype", dest="G", help=\
            "Use imputed/genotypes if available, should be in VCF file format")
    p.add_option("-v",  "--vcf_file", action="store_true", dest="inputisvcfile",
                 help="the input is a VCF file")
    p.add_option("-q", "--quality_threshold", type="int", dest="qual",
                 help="base quality threshold to take allele counts from")
    p.add_option("-D", "--debug", action="store_true", dest="D", help="debug")
    p.add_option("-V", "--output_vcf", action="store_true", dest="outputVCF",
                 help="Output the results to a VCF file")

    options, args = p.parse_args()
    if options.qual: pass
    else: options.qual = 20


    # For testing purposes
    debug = 1
    output = open(options.filename, 'wb')
    file_a = open(args[0],"rb")

    if not options.input:
        bam_Names = args[1:]
    else:
        bam_Names = []
        inputs = open(options.input, 'rU')
        for line in inputs:
            bam_Names.append(line.strip("\n"))

    bam_files = []
    for filename in bam_Names:
        bam_files.append(pysam.Samfile(filename,"rb"))


    INDEX_BASE = ['A', 'C', 'G', 'T']
    counts_matrix = []
    c_m = []
    rsID = []
    t = 0
    for line in file_a:
        counts = []
        c = []

        line = line.strip('\n').split('\t')
        # Counts is a list of numpy arrays

        if options.inputisvcfile:
            region = str(line[0])
            position = int(line[1]) # Start for a VCF file
            if len(line[3]) > 1 or len(line[4]) > 1:
                isIndel = True
            else:
                rsID.append(line[2])
                isIndel = False
        else:
            region = str(line[0])
            position = int(line[2]) # End position
            start = int(line[1])
            if position - start != 1:
                isIndel = True
            else:
                # Unfortunately the BED file does not specifiy this
                rsID.append(line[-1])
                t += 1
                isIndel = False

        cA = np.zeros(len(bam_files), dtype=np.uint32)
        cC = np.zeros(len(bam_files), dtype=np.uint32)
        cG = np.zeros(len(bam_files), dtype=np.uint32)
        cT = np.zeros(len(bam_files), dtype=np.uint32)

        if not isIndel:
            for i, bamfile in enumerate(bam_files):
                p_v = AlleleCounter(region, position,
                        phredThreshold=options.qual)
                # -1 to convert 0-based which pysam uses
                bamfile.fetch(p_v.region, p_v.position-1, p_v.position-1,
                        callback=p_v)

                # c.append(p_v.counts)
                cA[i] = p_v.counts[0]
                cC[i] = p_v.counts[1]
                cG[i] = p_v.counts[2]
                cT[i] = p_v.counts[3]

            c_m.append(cA)
            c_m.append(cC)
            c_m.append(cG)
            c_m.append(cT)

        else:
            # For right now
            pass

        # For testing purposes
        if options.D:
            if debug > 10: break
            else: debug += 1
        else:pass
    #print(c_m) 
    INDEX_BASE = ['A', 'C', 'G', 'T']
    out = {}
    out['counts'] = c_m
    out['index'] = INDEX_BASE
    out['rsIDs'] = rsID
    pickle.dump(out, output)
    #counts_df = pd.DataFrame(counts_matrix)
    #counts_df.save(opttion.out)

if __name__ == '__main__':
    main()
