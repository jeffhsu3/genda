#!/usr/bin/env python

""" Gets the read bases at each position specified by the first argument, which
is a VCF/BED/GFF/GTF file, from an indexed BAMfile/s or pileup/s:
    1. Whether the loci is heterzygous or homozygous if the genotyping is not
    known
    2. The p-value of the loci exhibiting allelic imbalance if the loci is
    heterzygous.

Usage:
    python AEI_2.py bedfile/vcf tab-delimited sample_mapping -o OUTPUT

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

import requests
import pysam
import numpy as np
import pandas as pd

from genda.formats.panVCF import VCF
import genda.stats.likelihood_funcs as lf
from genda.pysam_callbacks.allele_counter import AlleleCounter


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
    """ Given a genotype, only run AEI on the known heterozygotes
    """
    # Restrict to heterozygotes
    reduced_chrm = chrm[sample_genotype == 1]
    reduced_pos = pos[sample_genotype == 1]
    reduced_index = sample_genotype.index[sample_genotype == 1]
    bamfile = pysam.Samfile(sample_genotype.name, 'rb')
    for i, j, k in zip(reduced_chrm, reduced_pos, reduced_index):
        """ Vectorize this?
        """
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
        print(variant.counts.T)
        print(k)
        print(sample_genotype.name)
        print(c_m.columns[0:4])
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
    p.add_option("-a", "--annot", dest="annot", help="Annotation file\
            for SNPs")
    p.add_option("-G", "--genotype", dest="G", help=\
            "Use imputed/genotypes if available, should be in VCF file format")
    p.add_option("-v",  "--vcf_file", action="store_true", dest="inputisvcfile",
                 help="the input is a VCF file")
    p.add_option("-q", "--quality_threshold", type="int", dest="qual",
                 help="base quality threshold to take allele counts from")
    p.add_option("-P", "--variant_positions", action="store", help="")
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

    if options.v:
        vcf = VCF(args[0])
        chrm = vcf.vcf['#CHROM']
        pos = vcf.vcf['POS']
    elif options.G:
        pass
        #file_a = open(options.G,  'rU') 
    else:
        vcf = pd.read_csv(args[0], sep=" ")
        # Read in annotation file
        try:
            vcf_annot = pd.read_csv(options.annot)
        except (IOError, AttributeError):
            ensembl_server = "http://beta.rest.ensembl.org"
            ensembl_extension = "/vep/human/id/%s/consequences?"
            chrm = []
            pos = []
            for i in vcf.ix[: ,0]:
                r = requests.get(ensembl_server + ensembl_extension % i)
                if r.status_code == 200:
                    r = r.json()
                    try:
                        chrm.append(r['data']['location']['name'])
                        pos.append(r['data']['location']['start'])
                    except IndexError:
                        print(r)

            chrm = np.array(chrm)
            pos = np.array(pos)


    # A tab delimited file mapping sample names to bams ############# 
    bam_inputs = open(args[1], 'rU')
    sample_to_file = {}
    for line in bam_inputs:
        line = line.split("\t")
        sample_to_file[line[0]] = line[1].rstrip("\n").rstrip("/n")


    INDEX_BASE = ['A', 'C', 'G', 'T']
    if options.c:
        count_threshold = options.c
    else:
        count_threshold = 20

    multi_tuples = []
    # Temproary change the column names
    #vcf.vcf.rename(columns=sample_to_file, inplace=True)
    subset_vcf = vcf.vcf.ix[:, sample_to_file.keys()].copy()
    subset_vcf.rename(columns=sample_to_file, inplace=True)
    subset_geno = vcf.geno.ix[:, sample_to_file.keys()].copy()
    subset_geno.rename(columns=sample_to_file, inplace=True)
    for i in subset_geno.columns:
        multi_tuples.append(i)
        multi_tuples.append(i)
        multi_tuples.append(i)
        multi_tuples.append(i)
    multi = zip(multi_tuples, [0,1,2,3]*len(subset_geno))
    multi_index = pd.MultiIndex.from_tuples(multi, names=['sample', 'alleles'])
    counts_matrix = pd.DataFrame(np.zeros((subset_geno.shape[0], 
                                           len(subset_geno.columns)*4), 
                                           dtype=np.int16),
                                 index=subset_vcf.index, columns=multi_index)

    counts_fixed = functools.partial(counts_for_individuals, c_m=counts_matrix,
                                     chrm=chrm, pos=pos)
    subset_geno.apply(counts_fixed, axis=0)
    #counts_matrix.to_csv(options.filename)
    counts_matrix.save('options.filename')

if __name__ == '__main__':
    main()
    #cProfile.run('main()')
