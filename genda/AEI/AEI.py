from math import ceil
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from textwrap import wrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from genda import calculate_minor_allele_frequency, calculate_ld
from genda.plotting import (should_not_plot, add_gene_bounderies,
        add_snp_arrow)
from genda.eQTL import plot_eQTL


#:TODO move this to plotting utils
def remove_tr_spines(ax):
    invisible_spines = ['top', 'right']
    for i in invisible_spines:
        ax.spines[i].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
    
    return(ax)


def single_snp_aei_test(geno, outliers, allelic_ratio, num_threshold=5):
    """
    """
    geno = geno[np.logical_not(outliers)]
    het_combined = allelic_ratio[np.array(geno == 1)]
    homo_combined = allelic_ratio[np.array(np.logical_or(geno==0, geno==2))]
    if len(het_combined) < num_threshold or len(homo_combined) < num_threshold:
        return(1)
    else:
        return(ttest_ind(het_combined, homo_combined, equal_var=False)[1])


def single_snp_aei_test2(geno, allelic_ratio, num_threshold=5):
    het_combined = allelic_ratio[np.array(geno == 1)]
    homo_combined = allelic_ratio[np.array(np.logical_or(geno==0, geno==2))]
    # Het_combined must have a higher ratio than homo_combined
    if len(het_combined) < num_threshold or len(homo_combined) < num_threshold:
        return(1)
    elif het_combined.mean() < homo_combined.mean():
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
    def __init__(self, aei_dataframe, geno, snp_annot, gene_name):
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
        self.geno = geno
        # :TODO calculate MAF of the snps in self.geno
        # :TODO assert self.geno and self.snp_annot are same shape
        self.snp_annot = snp_annot
        try:
            self.sample_ids = [i[0] for i in self.aei.columns][::4]
        except AttributeError:
            self.sample_ids = [i[0] for i in self.aei.index][::4]
        ids = (pd.Index(self.sample_ids)
                .intersection(self.geno.columns))
        self.ids = ids
        self.maf = calculate_minor_allele_frequency(self.geno.ix[:, ids])
        self.gene_name = gene_name
        # SNPs of interest
        #self.ld = calculate_ld(self.geno.ix[:,0], snps_interest) 
        # :TODO Need to merge aei ids into here


    def calc_aei(self, num_threshold=5, single_snp=None):
        """
        Arguments
        ---------
        eQTL - eQTL pvalues
        num_threshold - threshold of heterozygous samples before
        single_snp - calculate the pvalue at just a single SNP
        """
        base_pair_to_index = {'A':0, 'C': 1, 'G': 2, 'T': 3}
        ids = self.ids
        # Focusing only on single-nucleotide polymorphisms
        if type(self.aei) == pd.DataFrame:
            not_indel = [i for i in self.aei.index if\
                    len(self.snp_annot.ix[i, 'a1']) == 1]
        else:
            not_indel = [self.aei.name]
        hets = dosage_round(self.geno.ix[not_indel, :])[ids]
        pvalues_out = np.ones((self.geno.shape[0], len(not_indel)),
            dtype=np.double)
        m_size = (len(ids), len(not_indel))
        aei_ratios = pd.DataFrame(
                np.zeros(m_size, dtype=np.uint32),
                index=pd.Index(ids), 
                columns=pd.Index(not_indel))
        overall_counts = pd.DataFrame(
                np.zeros((len(not_indel), 2), dtype=np.uint32), 
                index=pd.Index(not_indel))
        # sufficient hets is really het counts
        sufficient_hets = pd.Series(
                data=np.repeat(0, len(not_indel)), 
                index=pd.Index(not_indel),
                dtype=np.int32)
        hets_dict = {}
        pvalues_out = pd.DataFrame(pvalues_out, index=self.geno.index,
                columns=pd.Index(not_indel))
        # :TODO make this into a function
        for i, j in enumerate(not_indel):
            try:
                REF = base_pair_to_index[self.snp_annot.ix[j, 'a0']]
                ALT = base_pair_to_index[self.snp_annot.ix[j, 'a1']]
            except KeyError:
                continue
            #  hi = hets for the current snp
            hi = hets.ix[i,:]
            hi = hi[np.logical_and(hi == 1.0, 
                                    np.logical_not(pd.isnull(hi)))]
            sufficient_hets[j] = len(hi)
            if len(hi) >= num_threshold:
                refi = [(k, REF) for k in hi.index]
                alti = [(k, ALT) for k in hi.index]
                try:
                    ref = np.asarray(self.aei.ix[j, refi].values,
                            dtype=np.float64)
                    alt = np.asarray(self.aei.ix[j, alti].values, 
                            dtype=np.float64)
                except pd.core.indexing.IndexingError:
                    ref = np.asarray(self.aei[refi].values, dtype=np.float64)
                    alt = np.asarray(self.aei[alti].values, dtype=np.float64)
                # THIS IS THE WRONG WAY TO MULTIINDEX (commented out below)
                overall_counts.ix[i, :] = [np.nansum(ref), np.nansum(alt)]
                allelic_ratio = alt/ref
                aei_ratios.ix[hi.index, i] = pd.Series(allelic_ratio, index=hi.index)
                # Flip the ratio
                allelic_ratio[allelic_ratio > 1] =\
                        1/allelic_ratio[allelic_ratio > 1]
                outliers = (np.log2(allelic_ratio) < -3.5)
                # Need to do something better for outliers
                allelic_ratio = allelic_ratio[np.logical_not(outliers)]
                if not single_snp:
                    geno_t = self.geno.ix[:, hi.index]
                    geno_t = dosage_round(geno_t.ix[:, np.logical_not(outliers)])
                    pvalues_out.iloc[:,i] = (geno_t.
                                       apply(single_snp_aei_test2, axis=1, 
                                       args=(allelic_ratio,)))
                else:
                    geno_t = self.geno.ix[single_snp, hi.index]
                    geno_t = dosage_round(geno_t[np.logical_not(outliers)])
                    pvalues_out.ix[single_snp, i] =\
                            single_snp_aei_test2(geno_t, allelic_ratio)
                hets_dict[j] = hi.index
            else:
                outliers = np.repeat(False, len(hi))
                pvalues_out.iloc[:, i] = np.repeat(np.nan, self.geno.shape[0])
        self.pvalues = pvalues_out
        self.hets_dict = hets_dict
        self.sufficient_hets = sufficient_hets
        overall_counts['Nhets'] = sufficient_hets
        self.overall_counts = overall_counts
        # This never actually gets used, but there is currently an error if you
        # of the local variable being unbound
        #self.outliers = outliers
        self.ratios = aei_ratios
        # Outliers need to be fixed only set on the most recent one atm
        


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
        outliers =  np.logical_not(np.logical_or(
                          allelic_ratio < -3.0 ,
                          allelic_ratio > 3.0
                          )) 
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


    def aei_plot_single(self, tsnp, ax=None, focus_snp=None):
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
        adj_pvalue = -1*np.log10(self.pvalues.loc[:, tsnp])
        if focus_snp:
            snp = focus_snp
        else:
            # :TODO fix for both use cases
            #snp = subset.iloc[np.nanargmax(adj_pv), 0]
            snp = self.pvalues.index[np.nanargmax(adj_pvalue)]
        snp_iloc = [i for i, j in enumerate(adj_pvalue.index)\
                if j == snp][0]
        color1 = calculate_ld(self.geno,
            snp)[adj_pvalue.index].values
        scatter = ax.scatter(pos, 
                        adj_pvalue, s=size_maf, c=color1)
        ax.set_ylabel(r'-$log_{10}$ AEI p-value')
        ylim = (max(adj_pvalue) + max(adj_pvalue/6.0))
        ax.set_ylim((-0.01, ylim))
        ax = add_snp_arrow(adj_pvalue[snp], pos[snp_iloc], snp, ax)
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
