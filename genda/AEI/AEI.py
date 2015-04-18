from math import ceil
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from textwrap import wrap

from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from genda import calculate_minor_allele_frequency, calculate_ld
from genda.plotting import (should_not_plot, add_gene_bounderies,
        add_snp_arrow, plot_transcript)
from genda.eQTL import plot_eQTL


#:TODO move this to plotting utils
def remove_tr_spines(ax):
    invisible_spines = ['top', 'right']
    for i in invisible_spines:
        ax.spines[i].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
    
    return(ax)


def single_snp_aei_test(geno, outliers, allelic_ratio, num_threshold=2):
    """
    """
    geno = geno[np.logical_not(outliers)]
    het_combined = allelic_ratio[np.array(geno == 1)]
    homo_combined = allelic_ratio[np.array(np.logical_or(geno==0, geno==2))]
    if len(het_combined) < num_threshold or len(homo_combined) < num_threshold:
        return(1)
    else:
        return(ttest_ind(het_combined, homo_combined, equal_var=False)[1])


def single_snp_aei_test2(geno, allelic_ratio, num_threshold=2):
    het_combined = allelic_ratio[np.array(geno == 1)]
    homo_combined = allelic_ratio[np.array(np.logical_or(geno==0, geno==2))]
    # Het_combined must have a higher ratio than homo_combined
    if len(het_combined) < num_threshold or len(homo_combined) < num_threshold:
        return(1)
    elif (het_combined.mean() > homo_combined.mean()):
        return(1)
    else:
        return(ttest_ind(het_combined, homo_combined, equal_var=False)[1])


def dosage_round(geno, threshold = 0.5):
    """ Rounds dosage to threshold
    """
    geno[ geno < 1 - threshold] = 0
    geno[np.logical_and(geno >= 1-threshold, geno <= 1 + threshold)] = 1
    geno[ geno > 1 + threshold] = 2
    return geno


class AEI(object):
    def __init__(self, aei_dataframe, geno, snp_annot, gene_name,
            maf_threshold = 0.05):
        """ 
        Arguments
        ---------
        aei_dataframe - AEI count dataframe gotton from func :TODO list script
        to obtain aei from bam file in genda/scripts.  Must be subsetted for
        SNPs within gene.
        geno - genotype dataframe to do testing on
        surrounding SNP.
        snp_annot - dataframe of SNP annotations
        """
        self.aei = aei_dataframe
        # :TODO assert self.geno and self.snp_annot are same shape
        try:
            self.sample_ids = [i[0] for i in self.aei.columns][::4]
        except AttributeError:
            self.sample_ids = [i[0] for i in self.aei.index][::4]


        ids = (pd.Index(self.sample_ids)
                .intersection(geno.columns))
        idx = pd.IndexSlice
        try:
            self.aei.sort_index(axis=1, inplace=True)
            self.aei = self.aei.loc[:, idx[ids, :]]
            self.aei.sort_index(axis=1, inplace=True)
            self.geno = geno.ix[:, 
                    self.aei.columns.get_level_values(0)[::4]]
        except TypeError:
            # Case where aei_dataframe is just a series
            self.aei = self.aei.sort_index()
            self.aei = self.aei.loc[idx[ids, :]]
            self.aei = self.aei.sort_index()
            self.geno = geno.ix[:,
                    self.aei.index.get_level_values(0)[::4]]
        self.ids = self.geno.columns
        self.maf = calculate_minor_allele_frequency(self.geno.ix[:, ids])
        # Restrict to  > 5%
        self.geno = self.geno.ix[(self.maf >= maf_threshold) &\
                (self.maf <= 1 - maf_threshold), :]
        self.snp_annot = snp_annot.ix[self.geno.index, :]
        self.maf = self.maf[self.geno.index]
        self.gene_name = gene_name
        self.aei = self.aei.ix[self.geno.index.intersection(self.aei.index),:]
        # SNPs of interest
        #self.ld = calculate_ld(self.geno.ix[:,0], snps_interest) 
        # :TODO Need to merge aei ids into here


    def calc_aei(self, num_threshold=5, count_threshold = 30, 
            single_snp=None, outlier_thresh=-4.32):
        """
        :TODO split this up into two functions 1 for series and one for
        dataframe
        Arguments
        ---------
        eQTL - eQTL pvalues
        num_threshold - threshold of heterozygous samples at the tag SNP
        required
        single_snp - calculate the pvalue at just a single SNP
        count_threshold - number counts at the tag SNP within an 
        individual before that individual will be used in the fit
        """
        base_pair_to_index = {'A':0, 'C': 1, 'G': 2, 'T': 3}
        ids = self.ids
        # Focusing only on single-nucleotide polymorphisms as indicator
        if type(self.aei) == pd.DataFrame:
            not_indel = [i for i in self.aei.index if\
                    (len(self.snp_annot.ix[i, 'a1']) == 1)\
                     & len(self.snp_annot.ix[i, 'a0']) == 1]
        else:
            not_indel = [self.aei.name]
        geno_t = dosage_round(self.geno)
        hets = geno_t.ix[not_indel, :]
        m_size = (len(ids), len(not_indel))
        aei_ratios = pd.DataFrame(
                np.zeros(m_size, dtype=np.uint32),
                index=self.geno.columns, 
                columns=pd.Index(not_indel))
        overall_counts = pd.DataFrame(
                np.zeros((len(not_indel), 2), dtype=np.uint32), 
                index=pd.Index(not_indel))
        # 1 = False.  True for outliers_m means not outliers
        outliers_m = pd.DataFrame(np.ones(m_size,
                dtype=bool),
                index=self.geno.columns, 
                columns=pd.Index(not_indel))
        # sufficient hets is really het counts
        sufficient_hets = pd.Series(
                data=np.repeat(0, len(not_indel)), 
                index=pd.Index(not_indel),
                dtype=np.int32)
        hets_dict = {}

        #REF = [base_pair_to_index[i] for i in self.snp_annot.ix[not_indel, 'a0']]
        #ALT = [base_pair_to_index[i] for i in self.snp_annot.ix[not_indel, 'a1']]
        def _run_df(hetr):
            #hi = hetr[(hetr == 1.0) &  np.logical_not(pd.isnull(hetr))] 
            return((hetr == 1.0) &  np.logical_not(pd.isnull(hetr)))

        #:TODO deal with aei being class series

        hets = hets.apply(_run_df, axis=1)
        #print(hets)
        sufficient_hets = hets.sum(axis=1)
        hets = hets.ix[sufficient_hets >= num_threshold, :]
        aei_t = self.aei.ix[hets.index,:]
        #self.aei.apply(calc_pvalues, axis=1, args=(hets))
        pvalues_out = np.ones((self.geno.shape[0], aei_t.shape[0]),
            dtype=np.double)
        pvalues_out = pd.DataFrame(pvalues_out, index=self.geno.index,
                columns=aei_t.index)
        aei_ratios = pd.DataFrame(
                np.zeros((len(ids), len(aei_t.index)), dtype=np.uint32),
                index=self.geno.columns, 
                columns=aei_t.index)

        c = 0

        for i, row in hets.iterrows():
            try:
                REF = base_pair_to_index[self.snp_annot.ix[i, 'a0']]
                ALT = base_pair_to_index[self.snp_annot.ix[i, 'a1']]
            except KeyError:
                continue
            # hi = hets index
            hi = list(row[row].index)
            idx = pd.IndexSlice
            ref = np.asarray(aei_t.loc[i, idx[hi, REF]].values,
                    dtype=np.float64)
            alt = np.asarray(aei_t.loc[i, idx[hi, ALT]].values,
                    dtype=np.float64)
            s = ref + alt
            s = s >= count_threshold 
            allelic_ratio = alt/ref
            allelic_ratio = pd.Series(allelic_ratio,
                    index=hi)
            aei_ratios.ix[hi, c] = allelic_ratio
            allelic_ratio[allelic_ratio > 1] =\
                    1/allelic_ratio[allelic_ratio > 1]

            outliers = (np.log2(allelic_ratio) < outlier_thresh) |\
                    np.isinf(np.log2(allelic_ratio))
            outliers = np.logical_or(outliers, np.logical_not(s))
            outliers_m.ix[hi, i] = outliers
            allelic_ratio = allelic_ratio[np.logical_not(outliers)]
            geno_i = geno_t.ix[:, allelic_ratio.index]
            pvalues_out.iloc[:, c] = (geno_i.
                               apply(single_snp_aei_test2, axis=1, 
                               args=(allelic_ratio,)))
            hets_dict[i] = hi
            c += 1 
        self.pvalues = pvalues_out
        self.hets_dict = hets_dict
        self.ratios = aei_ratios
        self.outliers = outliers_m
        

    def aei_barplot(self, csnp, tsnp, ax=None, gene_name=None, title=True):
        """ AEI barplot

        Arguments
        ---------
        csnp - cis snp to test association
        tsnp - tag snp to run the AEI test on
        ax - axis to plot on if None creates a new figure
        gene_name - replace the default ID ie replace with more 
        readable gene symbol
        """
        color = dosage_round(self.geno.ix[csnp, self.hets_dict[tsnp]])
        nplots = 1
        if ax:
            pass
        else:
            fig, ax = plt.subplots(nrows=nplots, ncols=1, 
                    figsize=(12, 4*nplots), sharey=False,
                    sharex=True, subplot_kw=dict(axisbg='#FFFFFF'))
        title_str = ('AEI at tag SNP {tag} for {gene}\n shaded by'
                'genotype at {csnp}')
        if gene_name:
            title_str = title_str.format(tag=tsnp, 
                    gene=gene_name, csnp=csnp) 
        else:
            title_str = title_str.format(tag=tsnp, 
                    gene=self.gene_name, csnp=csnp) 
        if title:
            ax.set_title(title_str)
        else:  
            ax.set_title('')
        ax.set_xlabel('Samples')
        ax.set_ylabel('$log_{2}$ Allelic Fraction')
        width = 1
        allelic_ratio = self.ratios.ix[self.hets_dict[tsnp], tsnp]
        allelic_ratio_i = np.argsort(allelic_ratio.values)
        allelic_ratio = np.log2(allelic_ratio.iloc[allelic_ratio_i])
        outliers = np.logical_not(self.outliers.ix[self.hets_dict[tsnp], tsnp])
        color_geno = []
        color = color[allelic_ratio_i][outliers]
        for i in color:
            if i == 0 or i ==2:
                color_geno.append('green')
            else:
                color_geno.append('#FFAE00')
        allelic_ratio = allelic_ratio[outliers]
        ind = np.arange(len(allelic_ratio))
        rects1 = ax.bar(ind, allelic_ratio, width, color = color_geno, 
                linewidth=1)
        ax.set_xlim((-1, len(allelic_ratio+1)))
        ax = remove_tr_spines(ax)
        if ax: 
            return(ax)
        else:
            fig.tight_layout()
            return(fig)


    def aei_within_individual(self, ax=None):
        ''' A hidden Markov attempt at determining AEI within
        a transcript
        '''
        self.aei


    def barplot_across_genes(self, samples, other_samples, 
            ax=None, bar_width = 0.35):
        if ax:
            pass
        else:
            fig, ax = plt.subplots(nrows=2 , ncols=1, 
                    sharey=False, sharex=True, 
                    subplot_kw=dict(axisbg='#FFFFFF'))
        np.zeros((2, self.ratios.shape[1]), dtype=np.float64)
        ar_samp = defaultdict(list)
        ar_osamp = defaultdict(list)
        for i in self.ratios.columns:
            ar_samp[i] = []
            ar_osamp[i] = []
        for i in samples:
            ar = []
            for j in self.ratios.columns:
                if i in self.hets_dict[j]:
                    c_ar = self.ratios.ix[i, j]
                    if c_ar <= 1:
                        c_ar = 1/c_ar
                    else: pass
                    if self.outliers.ix[i, j]:
                        pass
                    else:
                        ar_samp[j].append(c_ar)
                else: pass
        for i in other_samples:
            for j in self.ratios.columns:
                if i in self.hets_dict[j]:
                    c_ar = self.ratios.ix[i, j]
                    if c_ar <= 1:
                        c_ar = 1/c_ar
                    else: pass
                    if self.outliers.ix[i, j]:
                        pass
                    else:
                        ar_osamp[j].append(c_ar)
                else: pass
        smean = []
        sstd = []
        for snps, ratios in ar_samp.iteritems():
            mdat = np.ma.masked_array(ratios, np.isnan(ratios))
            smean.append(np.mean(mdat))
            sstd.append(np.std(mdat))
        osm = []
        sosm = []
        for snps, ratios in ar_osamp.iteritems():
            mdat = np.ma.masked_array(ratios,np.isnan(ratios))
            osm.append(np.mean(mdat))
            sosm.append(np.std(mdat))
        ind = np.arange(len(ar_samp.keys()))
        rects1 = ax.bar(ind, smean, bar_width, color='r', yerr=sstd)
        rects2 = ax.bar(ind + bar_width, osm, bar_width,  color='y', yerr=sosm)
        if ax:
            return(ax)
        else: return(fig)


    def plot_aei_within_individual(self, samples, 
            other_samples=None, ax=None, exon_scale=1,
            intron_scale=20):
        if ax:
            pass
        else:
            fig, ax = plt.subplots(nrows=1 , ncols=1, 
                    sharey=False, sharex=True, 
                    subplot_kw=dict(axisbg='#FFFFFF'))
        counter = 0
        for i in samples:
            pos = []
            ar = []
            for j in self.ratios.columns:
                if i in self.hets_dict[j]:
                    c_ar = self.ratios.ix[i, j]
                    if c_ar <= 1:
                        c_ar = 1/c_ar
                    else: pass
                    if self.outliers.ix[i, j]:
                        pass
                    else:
                        pos.append(self.snp_annot.ix[j, 'pos'])
                        ar.append(c_ar)
                else: pass
            scatter = ax.plot(pos, ar, '-o')
            counter += 1
        if ax:
            return(ax)
        else: return(fig)



    def aei_plot_single(self, cissnp, ax=None, focus_snp=None):
        """
        Arguments
        ---------
        tsnp - a particular tag snp to plot association with """
        x_scale = 1e6
        cm = plt.cm.get_cmap('Blues')
        size_maf = ((200 * self.maf) + 20)
        pos = self.snp_annot.loc[:, 'pos']
        pos = np.asarray(pos, dtype=np.uint64)/x_scale
        if ax:
            pass
        else:
            fig, ax = plt.subplots(nrows=1 , ncols=1, 
                    sharey=False, sharex=True, 
                    subplot_kw=dict(axisbg='#FFFFFF'))
        adj_pvalue = -1*np.log10(self.pvalues.loc[:, cissnp])
        if focus_snp:
            print('focus snp')
            snp = focus_snp
        else:
            # :TODO fix for both use cases
            #snp = subset.iloc[np.nanargmax(adj_pv), 0]
            snp = self.pvalues.index[np.nanargmax(adj_pvalue)]
        print(snp)
        snp_iloc = [i for i, j in enumerate(adj_pvalue.index)\
                if j == snp][0]
        color1 = calculate_ld(self.geno,
            snp)[adj_pvalue.index].values
        scatter = ax.scatter(pos, 
                        adj_pvalue, s=size_maf, c=color1)
        ax.set_ylabel(r'-$log_{10}$ AEI p-value')
        ylim = (max(adj_pvalue) + max(adj_pvalue/6.0))
        ax.set_ylim((-0.01, ylim))
        ax = add_snp_arrow(adj_pvalue[snp_iloc], pos[snp_iloc], snp, ax)
        if ax:
            return(ax)
        else:
            return(fig)


    def aei_plot(self, meQTL, snp_plot=None, n_sufficient_hets=50, 
            common_only=False, mpld3plot=False, 
            focus_snp=None, ax=None, **kwargs):
        """
        Arguments
        ---------
        snp_plot - :TODO what is this again?
        n_sufficient - 
        common_only - Only common SNPs (>0.05)
        ax - matplotlib axis in which to add plot (for multiploting) NOT
        IMPLMENTED
        """
        x_scale = 1e6
        size_maf = ((200 * self.maf) + 20)
        cm = plt.cm.get_cmap('Blues')
        '''
        color1 = calculate_ld(self.geno,
            focus_snp)[self.geno.index].values
        '''
        if type(snp_plot) == pd.Series or type(snp_plot) == list:
            suff_hets = pd.Series(snp_plot, index = pd.Index(snp_plot))
        else:
            suff_hets = self.sufficient_hets[
                        np.logical_and(self.sufficient_hets >= n_sufficient_hets,
                        self.overall_counts.sum(axis=1)>=500)]
        nplots = len(suff_hets) + 2
        pos = self.snp_annot.loc[:, 'pos']
        pos = np.asarray(pos, dtype=np.uint64)/x_scale
        if ax:
            pass
        else:
            if len(suff_hets) > 1:
                ncols = 2
            else: 
                ncols = 1
            fig, ax = plt.subplots(nrows=int(ceil(nplots/ncols)) , ncols=2,
                    figsize=(20, 4*nplots/2), 
                    sharey=False, sharex=True, 
                    subplot_kw=dict(axisbg='#FFFFFF'))
        ko = 0
        for j in range(ncols):
            io = 0
            for i in range(int(ceil(len(suff_hets)/2.0))):
                if ko < len(suff_hets):
                    curr = suff_hets.index[ko]
                    adj_pvalue = -1*np.log10(self.pvalues.loc[:, curr])
                    scatter = ax[io, j].scatter(pos, 
                                    adj_pvalue, s=30)
                    ax[io, j].set_ylabel(r'-$log_{10}$ AEI p-value')
                    ylim = (max(adj_pvalue) + max(adj_pvalue/6.0))
                    ax[io, j].set_ylim((-0.01, ylim))
                    # Need to make the text relative positioning
                    #labels = list(self.annot_table.index)
                    #tooltip = mpld3.plugins.PointLabelTooltip(scatter, labels=labels)
                    ko += 1
                    io += 1
                else: 
                    pass
            scatter = ax[-1, j].scatter(pos,
                    -1*np.log10(meQTL.ix[self.snp_annot.index , 'p-value']),
                                     s=size_maf, cmap=cm)
            ax[-1, j].set_ylabel('-$log_{10}$ p-value')
            ax[-1, j].set_xlabel('Genomic Position (mb)')
            ax[-1, j] = plot_eQTL(meQTL, gene_name = self.gene_name,
                  annotation=self.snp_annot, dosage=self.geno, ax = ax[-1,j],
                  **kwargs)
            ax[-1, j] = remove_tr_spines(ax[-1, j])
            #ax[-1, j].set_title('%s eQTL plot' % (self.gene_name,), fontsize=25)
            '''
            labels = list(self.snp_table.index)
            tooltip = mpld3.plugins.PointLabelTooltip(scatter, labels=labels)
            mpld3.plugins.connect(fig, tooltip,  plugins.LinkedBrush(scatter))
            '''
        fig.tight_layout()
        return(fig)
