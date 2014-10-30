#!/usr/local/bin

import pickle
import itertools
import heapq
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import binom
import matplotlib.pyplot as plt


def generate_mi(rsIDs):
    """ Usefuls for reindexing for a new subset of rsIDS
    """
    BASES = ['A', 'C', 'G', 'T']
    temp_index = itertools.chain(*[[i]*4 for i in rsIDs])
    tuples = zip(temp_index, itertools.cycle(BASES))
    return pd.MultiIndex.from_tuples(tuples, names=['ID', 'BASE'])


class AEIData(object):
    """ Multi-index implementation of CountMatrix
    """

    def __init__(self, filename):
        data = self._load(filename)
        rsIDs = data['rsIDs']
        # Might implement this in the future as a multi_index
        self.df = pd.DataFrame(data['counts'], index=generate_mi(rsIDs))
        self.df = self.df/float(1)
        self.rsIDs = rsIDs

    def _load(self, filename):
        """ Loads a count matrix that is outputed from allele_counter.py
        """
        pkl_file = open(filename, 'rb')
        return(pickle.load(pkl_file))

    def set_genotypes(self, genotypes):
        """ Genotypes is a pandas Dataframe containing genotype information
        from the samples
        """
        #:TODO try some column (sample) automatching.  Not always applicable however.  
        self.genotypes = genotypes.reindex(index = self.rsIDs)

    def set_annotations(self, annotationfile):
        """ Loads a VCF SNP annotation file from dbSNP and automatically aligns
        data from the new set with the old set.
        """
        # Check if it is already a DataFrame
        # VCFs and BEDFiles
        # VCF files usually have really long comment sections
        # Read in file first and find the header line
        # :TODO Need to parse VCF header
        annot = pd.read_csv(annotationfile, sep="\t")
        grouped = annot.groupby("ID")
        index = [gp_keys[0] for gp_keys in grouped.groups.values()]
        temp = annot.reindex(index)
        temp.index = temp["ID"]

        # Drop these extraneous columns
        temp = temp.drop(["QUAL", "FILTER", "INFO", "ID"], axis=1)

        # Reindex according to the count data
        self.annot = temp.reindex(self.rsIDs)

        # Deal with mismatches in the annotation file used to generate
        # Need to use CHROM since it's the only float value, and np.NaN is
        # used for missing values.
        self.annot = self.annot[np.logical_not(np.isnan(self.annot["CHROM"]))]
        self.df = self.df.reindex(generate_mi(self.annot.index))
        self.genotypes = self.genotypes.reindex(self.annot.index)
        self.rsIDs = list(self.annot.index)

    def mask(self, count_threshold = 20, impute_threshold = 0.5):
        """ Mask locations that aren't heterozygotes and the loctions that
        don't meet a read count threshold.
        """
        try:
            # Reducing the dataframe based on annotations
            ref_tup = zip(self.annot.index, self.annot["REF"])
            alt_tup = zip(self.annot.index, self.annot["ALT"])
            ref = self.df.ix[ref_tup]
            alt = self.df.ix[alt_tup]
            # Need to collapse the multi index
            # :TODO find a more rigorous way to do this
            ref.index = self.genotypes.index
            alt.index = self.genotypes.index
            hets = np.logical_or(self.genotypes < 0.5, self.genotypes > 1.5)
            sums = (ref + alt) < count_threshold
            ref[np.logical_or(hets, sums)] = np.NaN
            alt[np.logical_or(hets, sums)] = np.NaN
            self.binom_p = pd.DataFrame(binom.sf(alt - 1, ref + alt, 0.5), 
                    columns=self.df.columns,index=ref.index)
            self.ratio = ref/((ref+alt)/float(1))
        except AttributeError:
            print("Need to run set_annotations and set_genotyeps first")

    def binom_test(self):
        pass

    def to_R_dataframe(self):
        """ Converts
        """
        pass

    def to_UCSC_track(self):
        """ Creates a bedfile with the combined p-values at each of the SNPs
        """
        pass

    def to_SQL(self):
        """ Write the data frame to SQL.
        """
        pass


class CountMatrix(object):
    """ Count Matrix holds counts from sequencing data
    """


    def __init__(self, filename):
        data = self._load(filename)
        # Might implement this in the future as a multi_index
        self.df = pd.Panel(data['counts'], items=data['rsIDs'],
                           minor_axis=['A','C','G','T'])
        self.index = data['rsIDs']


    def _load(self, filename):
        """ Loads a count matrix that is outputed from allele_counter.py
        """
        pkl_file = open(filename, 'rb')
        return(pickle.load(pkl_file))

    def ratio_df(self, threshold = 20):
        """ Converts the allele count pd.Panel to a dataframe of the ratios
        """
        # This could be easily optimized
        def _second_largest(x):
            x = heapq.nlargest(2, x)
            return float(x[1])
        major = self.df.apply(max, axis=2)
        minor = self.df.apply(_second_largest, axis=2)
        # Make  sure that the indexes match up.  Remove the ones that don't
        # Check the other way as well

        total = major + minor
        ratio = major/total
        sum_matrix = self.df.apply(sum, axis=2)
        self.sums = sum_matrix.transpose()
        self.major = major.transpose()
        self.minor = minor.transpose()
        self.ratio = ratio.transpose()
        self.ratio[self.sums < threshold] = np.NaN

    def binom_test(self):
        # Probably slower, but certainly more flexible
        self.biallelic = self.df.apply(heapq.nlargest, n=2)

    def set_genotypes(self, genotypes):
        """ Genotypes is a pandas Dataframe containing genotype information
        from the samples
        """
        self.genotypes = genotypes.reindex(index = self.df.items)

        # Automatically check if the genotypes are the same

    def set_annotations(self, annotationfile):
        """ Sets what the major allele frequencies and the minor allele
        frequencies are using an VCF annotation file.
        """
        # Check if it is already a DataFrame
        # VCFs and BEDFiles
        # VCF files usually have really long comment sections
        # Read in file first and find the header line
        annot = pd.read_csv(annotationfile, sep="\t")
        grouped = anot.groupby("ID")
        index = [gp_keys[0] for gp_keys in grouped.groups.values()]
        temp = annot.reindex(index)
        temp.index = temp["ID"]
        # Drop these extraneous columns
        temp.drop(["QUAL", "FILTER", "INFO", "ID"])
        # Reindex according to
        self.annot = temp.reindex(self.df.items)
        # Deal with mismatches in the annotation file used to generate
        self.annot = self.annot[np.logical_not(np.isnan(self.annot["REF"]))]
        self.df.items



def collide_snps(snp_annotation, gff, genes_of_interest = None):
    """ Grabs all SNPs that reside within a gene.  

    Arguments
    ---------

    Returns
    -------
    A dictionary of genes with keys being the gene names and a list of SNPs
    being the values.
    """
    out_dict = defaultdict(list)
    gene_ids = gff[8].apply(lambda x: (x.split("; ")[-1]
                                       .lstrip('gene_id "')
                                       .rstrip('"')))
    if genes_of_interest:
        gene_ids = gene_ids.apply(lambda x: x in genes_of_interest)
    m_matrix = gff.ix[np.logical_and(gene_ids, gff.ix[:, 2] == 'exonic_part'),
            [2,3,4,8]]
    #gene_ids = gene_ids[np.logical_and(





