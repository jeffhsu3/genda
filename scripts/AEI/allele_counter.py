#!/usr/bin/env python

""" Gets the read bases at each position specified by the first argument, which
is a VCF/BED/GFF/GTF file, from an indexed BAMfile/s or pileup/s:
    1. Whether the loci is heterzygous or homozygous if the genotyping is not
    known
    2. The p-value of the loci exhibiting allelic imbalance if the loci is
    heterzygous.

Usage:
    python getReads.py bedfile/vcf bamfile1 [bamfile2 ...] -o OUTPUT

Limtations: Doesn't accurately account paralogous gene expression, which can
manifest itself in allelic imbalance.  One way would be to have an input file
of all paralogous genes and calculated the probability of a read mismatching
and then the liklihood that it is something not valulable.

:TODO Do comparison with reference allele
:TODO Account for mis-alignments of the minor allele
:TODO Need to handle indel
:TODO add an option that includes a gene annotation file
:TODO add support for pileup files
:TODO also how does PCR duplicates affect?
:TODO makes this more modular and testable
:TODO also make an R version?

:TODO Output the results into a VCF file

Written by Jeffrey Hsu
2010
"""
# Standard Libraries
import optparse, pickle

import pysam
import numpy as np

class AlleleCounter():
    """ Gets the allele counts 
    """
    BASE_INDEX = {'A':0, 'a':0, 'C':1, 'c':1, 'G':2, 'g':2, 'T':3, 't':3}
    def __init__(self, region, position, phredThreshold=0):
        self.region = region
        self.position = int(position)
        self.phredThreshold = phredThreshold
        self.counts = np.zeros(4, dtype=np.uint32)

    def __call__(self, alignment, position = None, phredThreshold=None):
        if position == None: 
            position = self.position
        else: pass
        if phredThreshold == None: 
            qualT = self.phredThreshold
        else: pass
        if alignment.is_duplicate or alignment.mapq <= 50:
            pass
        else:
            if 3 in [i[0] for i in alignment.cigar]:
                t = [i[1] for i in alignment.cigar if i[0] == 3]
                # print(len(t))
                inserts = sum(t)
                #print("Alignment Start: %i" % alignment.pos)
                #print(alignment.seq)
                index = position - inserts - alignment.pos - 1   
                #print(index)
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
    """ Main loop that iterates over the desired loci to calculate Aellic
    Imbalance

    """


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

    # Handling of multiple BAM/SAM inputs
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
        else:
            region = str(line[0])
            position = int(line[2]) # End position
            start = int(line[1])
            if position - start != 1:
                isIndel = True
            else: 
                rsID.append(line[-1])  
                t += 1
                isIndel = False
        
        if not isIndel:
            for bamfile, bamNames in map(None, bam_files, bam_Names):
                # :TODO in the VCF and bed files make sure to type the attributes
                p_v = AlleleCounter(region, position,
                        phredThreshold=options.qual)

                bamfile.fetch(p_v.region, p_v.position-1, p_v.position-1,
                        callback=p_v)

                c.append(p_v.counts)

            c_m.append(c)
        else: 
            # For right now
            pass

        # For testing purposes
        if options.D:
            if debug > 500: break
            else: debug += 1
        else:pass
    
    INDEX_BASE = ['A', 'C', 'G', 'T']
    out = {}
    out['counts'] = c_m
    out['rsIDs'] = rsID
    pickle.dump(out, output)
    #counts_df = pd.DataFrame(counts_matrix)
    #counts_df.save(opttion.out)

if __name__ == '__main__':
    main()
