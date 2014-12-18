from math import ceil
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from textwrap import wrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


from genda import calculate_minor_allele_frequency, calculate_ld
from genda.plotting import should_not_plot, add_gene_bounderies
from genda.eQTL import plot_eQTL


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


def single_snp_aei_test2(geno, allelic_ratio, num_threshold=10):
    het_combined = allelic_ratio[np.array(geno == 1)]
    homo_combined = allelic_ratio[np.array(np.logical_or(geno==0, geno==2))]
    if len(het_combined) < num_threshold or len(homo_combined) < num_threshold:
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
        self.sample_ids = [i[0] for i in self.aei][::4]
        ids = pd.Index(self.sample_ids).intersection(self.geno.columns)
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
        not_indel = [i for i in self.aei.index if\
                len(self.snp_annot.ix[i, 'a1']) == 1]
        hets = dosage_round(self.geno.ix[not_indel, :])[ids]
        pvalues_out = np.zeros((self.geno.shape[0], len(not_indel)),
            dtype=np.double)
        m_size = (len(ids), len(not_indel))
        aei_ratios= pd.DataFrame(
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
                ref = np.asarray(self.aei.ix[j, refi].values,
                        dtype=np.float64)
                alt = np.asarray(self.aei.ix[j, alti].values, 
                        dtype=np.float64)
                # THIS IS THE WRONG WAY TO MULTIINDEX (commented out below)
                #print(self.aei.ix[ j, hi.index][0:20])
                overall_counts.ix[i, :] = [np.nansum(ref), np.nansum(alt)]
                allelic_ratio = alt/ref
                aei_ratios.ix[hi.index, i] = pd.Series(allelic_ratio, index=hi.index)
                # Flip 
                allelic_ratio[allelic_ratio > 1] = 1/allelic_ratio[allelic_ratio > 1]
                outliers = (np.log2(allelic_ratio) < -4.5)
                allelic_ratio = allelic_ratio[np.logical_not(outliers)]
                if not single_snp:
                    pvalues_out.iloc[:,i] = (self.geno.ix[:, hi.index].
                                       apply(single_snp_aei_test, axis=1, 
                                       args=(outliers, allelic_ratio)))
                else:
                    try:
                        geno_t = self.geno.ix[single_snp, hi.index]
                        geno_t = geno_t[np.logical_not(outliers)]
                        pvalues_out.ix[single_snp, i] =\
                                single_snp_aei_test2(geno_t, allelic_ratio)
                    except ValueError:
                        geno_t = self.geno.ix[single_snp, hi.index]
                        geno_t = geno_t.iloc[0,:]
                        geno_t = geno_t[np.logical_not(outliers)]
                        pvalues_out.ix[single_snp, i] =\
                                single_snp_aei_test2(geno_t, allelic_ratio)
                        # Why are there two identical indexes for geno?
                hets_dict[j] = hi.index
            else:
                outliers = np.repeat(False, len(hi))
                pvalues_out.iloc[:, i] = np.repeat(np.nan, self.geno.shape[0])
        self.pvalues = pvalues_out
        self.hets_dict = hets_dict
        self.sufficient_hets = sufficient_hets
        overall_counts['Nhets'] = sufficient_hets
        self.overall_counts = overall_counts
        self.outliers = outliers
        


    def aei_bar_plot(self, csnp, tsnp, ax=None):
        """ AEI barplot
        """
        pass


    def aei_plot(self, meQTL, snp_plot=None, n_sufficient_hets=50, 
            common_only=False, mpld3plot=False, focus_snp=None, ax=None, **kwargs):
        """
        Arguments
        ---------
        snp_plot - :TODO what is this again?
        n_sufficient - 
        common_only - Only common SNPs (>0.05)
        ax - matplotlib axis in which to add plot (for multiploting) NOT
        IMPLMENTED
        """
        if ax:
            raise NotImplementedError
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
        #pos = self.meQTL.loc[:, 'pos']
        pos = np.asarray(pos, dtype=np.uint64)/x_scale

        fig, ax = plt.subplots(nrows=int(ceil(nplots/2.0)) , ncols=2,
                figsize=(20, 4*nplots/2), 
                sharey=False, sharex=True, 
                subplot_kw=dict(axisbg='#FFFFFF'))
        ko = 0
        for j in range(2):
            io = 0
            for i in range(int(ceil(len(suff_hets)/2.0))):
                if ko < len(suff_hets):
                    curr = suff_hets.index[ko]
                    adj_pvalue = -1*np.log10(self.pvalues.loc[:, curr])
                    scatter = ax[io, j].scatter(pos, 
                                    adj_pvalue, s=30)
                    ax[io, j].set_ylabel(r'-$log_{10}$ AEI p-value',
                            fontsize=18)
                    ylim = (max(adj_pvalue) + max(adj_pvalue/6.0))
                    ax[io, j].set_ylim((-0.01, ylim))
                    """ 
                    ax[io, j].set_ylabel("\n".join(wrap(r'AEI -$log_{10}$ p-value ',
                        'for {tag}'.format(tag=curr)), 30), fontsize=15)
                    """
                    #ax[io, j].set_xlabel('Genomic Position (mb)', fontsize=15)
                    """
                    ax[io, j].set_title('AEI plot for %s (N=%i)' %
                            (curr,
                        self.overall_counts.ix[suff_hets.index[ko], 'Nhets']), fontsize=25)
                    """
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
            ax[-1, j].set_ylabel('-$log_{10}$ p-value', fontsize=15)
            ax[-1, j].set_xlabel('Genomic Position (mb)', fontsize=15)
            ax[-1, j] = plot_eQTL(meQTL, gene_name = self.gene_name,
                  annotation=self.snp_annot, dosage=self.geno, ax = ax[-1,j],
                  **kwargs)
            invisible_spines = ['top', 'right']
            # This needs to be a function
            for i in invisible_spines:
                ax[-1, j].spines[i].set_visible(False)
                ax[-1, j].xaxis.set_ticks_position('bottom')
                ax[-1, j].yaxis.set_ticks_position('left')
            #ax[-1, j].set_title('%s eQTL plot' % (self.gene_name,), fontsize=25)
            '''
            labels = list(self.snp_table.index)
            tooltip = mpld3.plugins.PointLabelTooltip(scatter, labels=labels)
            mpld3.plugins.connect(fig, tooltip,  plugins.LinkedBrush(scatter))
            '''
        fig.tight_layout()
        return(fig)
