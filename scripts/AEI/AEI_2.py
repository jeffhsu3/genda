#!/usr/bin/env python

""" Gets the read bases at each position specified by the first argument, which
is a VCF/BED/GFF/GTF file, from an indexed BAMfile/s or pileup/s:
    1. Whether the loci is heterzygous or homozygous if the genotyping is not
    known
    2. The p-value of the loci exhibiting allelic imbalance if the loci is
    heterzygous.

Usage:
    python getReads.py bedfile/vcf tab-delimited sample_mapping -o OUTPUT

Limtations: Doesn't accurately account paralogous gene expression, which can
manifest itself in allelic imbalance.  One way would be to have an input file
of all paralogous genes and calculated the probability of a read mismatching
and then the liklihood that it is something not valulable.

Written by Jeffrey Hsu
2010
"""

# Standard Libraries
import optparse, pickle
import cProfile
import functools
from collections import defaultdict

import pysam
import numpy as np
import pandas as pd

from pySeq.formats.panVCF import VCF
import pySeq.stats.likelihood_funcs as lf
from pySeq.pysam_callbacks.allele_counter import AlleleCounter


class lociInformation():
    """ Contains information about a given position in a BAM file

    Change into slots also make this into cython.
    """

    BASE_INDEX = {'A':0, 'a':0, 'C':1, 'c':1, 'G':2, 'g':2, 'T':3, 't':3}
    def __init__(self, region, position, phredThreshold=0):
        """ Location is one-based.  Internally this will get changed to zero
        base by pysam.
        """
        self.region = region
        self.position = int(position)
        self.phredThreshold = phredThreshold
        self.counts = np.zeros(4, dtype=np.uint8)


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
                qpos = read.qpos
                base_quality = read.alignment.qual[qpos]
                if ord(base_quality)-33 > qualT:
                    base = read.alignment.seq[qpos]
                    base = self.BASE_INDEX[base]
                    #strand = read.alignment.is_reverse
                    self.counts[base] += 1



def threshold_counts(counts, threshold=30, number=1):
    """ Returns True

    Makes sure that at least one of the samples meet a read count threshold.
    """
    counter = 0
    for i in counts:
        if sum(i) > threshold: counter += 1
    if counter >= number:return True
    else: return False


def column_names_to_bams(colname):
    """ Conversion of colnames to actual bam files
    """

    return colname + ".bam"


def counts_for_individuals(sample_genotype, c_m=None, chrm=None, pos=None, path=None,
                           bamfile_func = column_names_to_bams):
    """
    """
    # Restrict to heterozygotes
    reduced_chrm = chrm[sample_genotype == 1]
    reduced_pos = pos[sample_genotype == 1]
    reduced_index = sample_genotype.index[sample_genotype == 1]
    bamfile = pysam.Samfile(sample_genotype.name, 'rb')
    for i, j, k in zip(reduced_chrm, reduced_pos, reduced_index):
        """ Vectorize this?
        """
        print(i, j)
        variant = lociInformation(str(i), j,
                                  phredThreshold=0)
        bamfile.pileup(variant.region, variant.position,
                        variant.position+1, callback=variant)
        """
        variant_2 = AlleleCounter('chr' + str(i), j,
                                   phredThreshold=0)
        bamfile.fetch(variant.region, variant.position,
                        variant.position+1, callback=variant_2)
        print(variant.counts, variant_2.counts)
        """
        c_m.ix[k, [(sample_genotype.name, 0), (sample_genotype.name, 1), (sample_genotype.name, 2), (sample_genotype.name, 3),]] = variant.counts.T




def inner_apply(x, chrom, position):
    """ Try inner apply for the
    """
    pass


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

    vcf = VCF(args[0])

    # A tab delimited file mapping sample names to bams ############# 
    bam_inputs = open(args[1], 'rU')
    sample_to_file = {}
    for line in bam_inputs:
        line = line.split("\t")
        sample_to_file[line[0]] = line[1].rstrip("\n")


    # Handling of multiple BAM/SAM inputs


    INDEX_BASE = ['A', 'C', 'G', 'T']
    if options.c:
        count_threshold = options.c
    else:
        count_threshold = 20

    multi_tuples = []
    # Temproary change the column names
    vcf.vcf.rename(columns=sample_to_file, inplace=True)
    for i in vcf.vcf.columns[9:]:
        multi_tuples.append(i)
        multi_tuples.append(i)
        multi_tuples.append(i)
        multi_tuples.append(i)
    multi = zip(multi_tuples, [0,1,2,3]*len(vcf.vcf.columns[9:]))
    multi_index = pd.MultiIndex.from_tuples(multi, names=['sample', 'alleles'])

    counts_matrix = pd.DataFrame(np.zeros((vcf.vcf.shape[0], len(vcf.vcf.columns[9:])*4), dtype=np.int16),
                                 index=vcf.vcf.index, columns=multi_index)

    counts_fixed = functools.partial(counts_for_individuals, c_m=counts_matrix,
                                     chrm=vcf.vcf['CHROM'], pos=vcf.vcf['POS'])
    vcf.vcf.ix[:, 9:].apply(counts_fixed, axis=0)
    counts_matrix.to_csv(options.filename)
    pickle.p    

if __name__ == '__main__':
    main()
    #cProfile.run('main()')
