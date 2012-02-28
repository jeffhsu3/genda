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
import cProfile

import pysam
import numpy as np

from pySeq.formats.VCF import VCFfile
import pySeq.stats.likelihood_funcs as lf


"""
class AlleleCounter():
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
                base =  self.BASE_INDEX[base]
                b_qual = alignment.qual[index]
                if base != "N" and ord(b_qual)-33 > qualT:
                    self.counts[base] += 1
                else: pass
            else: pass
            """



class lociInformation():
    """ Contains information about a given position in a BAM file
    """

    BASE_INDEX = {'A':0, 'a':0, 'C':1, 'c':1, 'G':2, 'g':2, 'T':3, 't':3}
    def __init__(self, region, position, phredThreshold=0):
        """ Location is one-based.  Internally this will get changed to zero
        base by pysam.
        """
        self.region = region
        self.position = int(position)
        # So allele_counts['sample'] = [0,14,0,15]
        # BASE_INDEX converts index to base, ie C, in 14
        self.phredThreshold = phredThreshold
        self.counts = np.zeros(4, dtype=np.uint32)


    def __call__(self, pileups, position=None,
        phredThreshold=None, genotype=None):
        """ Genotype file contains the infromation the actualy genotyping
        calls.

        """
        if position == None: position = self.position
        else:  pass
        if phredThreshold == None: qualT = self.phredThreshold
        if pileups.pos != position-1: pass
        else:
            for read in pileups.pileups:
                if read.alignment.is_duplicate or read.alignment.mapq <= 50:
                    pass
                else:
                    qpos = read.qpos
                    base_quality = read.alignment.qual[qpos]
                    if ord(base_quality)-33 > qualT:
                        base = read.alignment.seq[qpos]
                        base = self.BASE_INDEX[base]
                        #strand = read.alignment.is_reverse
                        self.counts[base] += 1
                    else: pass



def threshold_counts(counts, threshold=30, number=1):
    """ Returns True

    Makes sure that at least one of the samples meet a read count threshold.

    """
    counter = 0
    for i in counts:
        if sum(i) > threshold: counter += 1
    if counter >= number:return True
    else: return False


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
    p.add_option("-G", "--genotype", dest="G", help=\
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
    p.add_option("-A", "--auto_parse", action="store_true", dest="auto",
                 help="Autoparse readgroups, if set to false will assume a\
                 single sample in each file")

    options, args = p.parse_args()
    if options.qual: pass
    else: options.qual = 20


    # For testing purposes
    debug = 1
    output = open(options.filename, 'wb')

    file_a = open(args[0],"rb")

    # Handling of multiple BAM/SAM inputs
    bam_Names = args[1:]
    bam_files = []
    for filename in bam_Names:
        bam_files.append(pysam.Samfile(filename,"rb"))


    INDEX_BASE = ['A', 'C', 'G', 'T']
    if options.c:
        count_threshold = options.c
    else:
        count_threshold = 20

    counts_matrix = []
    c_m = []
    rsID = {}
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
                rsID[line[-1]] = t 
                t += 1
                isIndel = False

        if not isIndel:
            for bamfile, bamNames in map(None, bam_files, bam_Names):
                # :TODO in the VCF and bed files make sure to type the attributes
                variant = lociInformation(region, position,
                                          phredThreshold=options.qual)
                bamfile.pileup(variant.region, variant.position,
                               variant.position+1, callback=variant)


                counts.append(variant.counts)

            counts_matrix.append(counts)
        else: 
            # For right now
            pass

        # For testing purposes
        if options.D:
            if debug > 500: break
            else: debug += 1
        else:pass

    out = {}
    out['counts'] = counts_matrix
    out['rsIDs'] = rsID
    pickle.dump(out, output)
    #counts_df = pd.DataFrame(counts_matrix)
    #counts_df.save(opttion.out)

if __name__ == '__main__':
    cProfile.run('main()')
