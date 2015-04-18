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


from IPython import embed

from genda.formats.panVCF import VCF
import genda.stats.likelihood_funcs as lf
from genda.pysam_callbacks.allele_counter import AlleleCounter
from genda.stats.aei_count_samples import aei_count_samples

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


def reget_sample_names(x, mapping=None):
    """
    """
    try:
        return mapping[x]
    except Nonetype:
        print('BAAAAAAAAAAAAD')


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
    p.add_option("-G", "--genotype", dest="G", action='store', help=\
            "Use imputed/genotypes if available, should be in VCF file format", 
            default=None)
    p.add_option("-v",  "--vcf_file", action="store_true", dest="inputisvcffile",
                 help="the input is a VCF file", default=False)
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
    p.add_option("-R", "--reduced", action="store_true", dest="reduced",
                 default=False)
    print('Hello')

    options, args = p.parse_args()
    if options.qual: pass
    else: options.qual = 20

    bam_inputs = open(args[1], 'rU')
    sample_to_file = {}
    for line in bam_inputs:
        line = line.split("\t")
        sample_to_file[line[0]] = line[1].rstrip("\n").rstrip("/n")
        file_to_sample = {v:k for k, v in sample_to_file.items()}

    # For testing purposes
    debug = 1
    output = open(options.filename, 'wb')
    if options.inputisvcffile:
        vcf = VCF(args[0])
        chrm = vcf.vcf['#CHROM']
        pos = vcf.vcf['POS']
        #geno = vcf.geno
        geno = pd.DataFrame(np.zeros((len(pos), len(sample_to_file))), 
                columns = pd.Index(sample_to_file.keys()))

    # :TODO fix different annotation options
    elif options.G:
        pass
    elif options.annot:
        chrom = 18
        print('Right annotation file')
        rsIDs = []
        pos = []
        file_a = pysam.Tabixfile(options.annot)
        # :TODO fix this
        # Currently counts at all SNPs
        a_iter = file_a.fetch('18')
        chrm = []
        '''
        base_path = '/proj/GenomicsHD/Atrial_RNASeq/'
        s_ann = pd.read_pickle(base_path + 'ref/snp_annot/' + str(chrom) +\
                '.pkl')
        '''
        debug = 0
        for i in a_iter:
            i = i.split("\t")
            if not int(i[-5]) == 0:
                rsIDs.append(i[3])
                chrm.append('chr' + str(i[0]))
                pos.append(int(i[1]))
                '''
                debug += 1
                if debug > 200:
                    break
                '''
        geno = pd.DataFrame(np.zeros((len(rsIDs), len(sample_to_file))), 
            index=pd.Index(rsIDs), columns = pd.Index(sample_to_file.keys()))
        #print(s_ann.head())
        pos = np.asarray(pos, dtype=np.uint32)
        chrm = np.asarray(chrm)
        print('Chrom length')
        print(len(chrm))
        print(len(pos))
        print(pos[-1])
    else:
        # grabbin at only a small subset of the genes
        vcf = pd.read_csv(args[0], sep=" ")
        # Read in annotation file
        try:
            vcf_annot = pd.read_csv(options.annot)
        except (IOError, AttributeError):
            ensembl_server = "http://beta.rest.ensembl.org"
            ensembl_extension = "/vep/human/id/%s/consequences?"
            headers = {'Content-Type' : 'application/json'}
            chrm = []
            pos = []
            for i in vcf.index:
                r = requests.get(ensembl_server + ensembl_extension % i,
                        headers=headers)
                if r.status_code == 200:
                    r = r.json()
                    try:
                        chrm.append(r['data'][0]['location']['name'])
                        pos.append(r['data'][0]['location']['start'])
                    except IndexError:
                        print(r)
                else: 
                    print('whelp')

            chrm = np.array(['chr' + j for j in chrm])
            pos = np.array(pos)
            geno = vcf

    # A tab delimited file mapping sample names to bams ############# 
    INDEX_BASE = ['A', 'C', 'G', 'T']
    if options.c:
        count_threshold = options.c
    else:
        count_threshold = 20
    print(len(pos))

    multi_tuples = []
    subset_geno = geno.ix[:, sample_to_file.keys()].copy()
    subset_geno.rename(columns=sample_to_file, inplace=True)
    for i in subset_geno.columns:
        multi_tuples.append(i)
        multi_tuples.append(i)
        multi_tuples.append(i)
        multi_tuples.append(i)
    multi = zip(multi_tuples, [0,1,2,3]*subset_geno.shape[1])
    multi_index = pd.MultiIndex.from_tuples(multi, names=['sample', 'alleles'])
    counts_matrix = pd.DataFrame(np.zeros((subset_geno.shape[0], 
                                           len(subset_geno.columns)*4), 
                                           dtype=np.uint32),
                                 index=subset_geno.index, columns=multi_index)

    c_m = aei_count_samples(subset_geno.values, subset_geno.columns.values,
            counts_matrix.values, np.asarray(chrm), pos)
    c_m = pd.DataFrame(c_m, index=subset_geno.index, columns=multi_index)
    c_m.to_pickle(options.filename)

if __name__ == '__main__':
    main()
    #cProfile.run('main()')
