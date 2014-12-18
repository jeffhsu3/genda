import copy
import numpy as np
import pandas as pd
import pysam
from scipy.stats import ttest_ind, pearsonr
import matplotlib.pyplot as plt
import mpld3
from mpld3 import plugins
from math import ceil
import matplotlib.patches as patches
from genda.formats.dosages import grab_gene_location
from mpld3.plugins import PluginBase
import jinja2
import json
import matplotlib
from matplotlib.path import Path

def make_rectangle(start, end, y1, y2):
    verts = [(start, y1),
            (start, y2),
            (end, y2),
            (end, y1),
            (start, y1),
            ]
    codes = [
            Path.MOVETO,
            Path.LINETO,
            Path.LINETO,
            Path.LINETO,
            Path.CLOSEPOLY,
            ]
    return (Path(verts, codes))


def should_not_plot(x):
    if x is None:
        return True
    elif isinstance(x, np.ndarray):
        return x.size==0
    elif isinstance(x, pd.DataFrame):
        return True
    else:
        return(bool(x))

def dosage_round(geno, threshold = 0.5):
    """ Rounds dosage to threshold
    """
    geno[ geno < 1 - threshold] = 0
    geno[np.logical_and(geno >= 1-threshold, geno <= 1 + threshold)] = 1
    geno[ geno > 1 + threshold] = 2
    return geno


def multiple_snp_aei_test(geno, outliers, allelic_ration, num_threshold=5):
    """
    """
    pass


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

    
class AEI_object(object):
    """ A class storing the aei data for a particular gene. 
    """
    def __init__(self, aei_pvalues, gene_name, annot_table, 
                 sufficient_hets, matrix_eQTL):
        self.aei_pvalues = aei_pvalues
        self.gene_name = gene_name
        self.annot_table = annot_table.ix[aei_pvalues.index, :]
        self.sufficient_hets = sufficient_hets
        self.meQTL = matrix_eQTL
        self.gene_names = 'Nothing'
        
        
    def aei_bar_plot(self, dosage, cis_snp, tag_snp, gene_name=None):
        """ AEI barplot
        """
        nplots = 1
        color = dosage_round(dosage.ix[cis_snp, self.hets_dict[tag_snp]])
        fig, ax = plt.subplots(nrows=nplots, ncols=1, figsize=(12, 4*nplots), 
                               sharey=False, sharex=True, 
                               subplot_kw=dict(axisbg='#FFFFFF'))
        if gene_name:
            title = 'AEI at tag %s for %s and\ncolored by genotype at %s' % (gene_name, tag_snp, cis_snp)
        else:
            title = "AEI at tag %s and\ncolored by genotype at %s" % (tag_snp,
                    cis_snp)

        ax.set_title(title, fontsize=20)
        ax.set_xlabel('Samples', fontsize=15)
        ax.set_ylabel('Allelic Fraction ($log_{2}$)', fontsize=15)
        width = 0.5
        allelic_ratio = self.ratios.ix[self.hets_dict[tag_snp],
            tag_snp]
        allelic_ratio_i = np.argsort(allelic_ratio.values)
        allelic_ratio = np.log2(allelic_ratio.iloc[allelic_ratio_i])
        outliers =  np.logical_not(np.logical_or(
                          allelic_ratio < -3.0 ,
                          allelic_ratio > 3.0
                          )) 
        color_geno = []
        color = color[allelic_ratio_i][outliers]
        for i in color:
            if i == 0 or i == 2:
                color_geno.append('green')
            else:
                color_geno.append('orange')
        allelic_ratio = allelic_ratio[outliers]
        ind = np.arange(len(allelic_ratio))
        rects1 = ax.bar(ind, allelic_ratio, width, color = color_geno)
        ax.set_xlim((-1, len(allelic_ratio+1)))
        return(fig)

        
    def aei_plot(self, snp_plot=None, n_sufficient_hets=50, common_only=False):
        """ AEI plots in mpld3
        """
        x_scale=1e6
        size_maf =((200 * self.maf) + 20) 
        cm = plt.cm.get_cmap('winter')
        if type(snp_plot) == pd.Series or type(snp_plot) == list:
            suff_hets = pd.Series(snp_plot, index = pd.Index(snp_plot))
        else:
            suff_hets = self.sufficient_hets[
                        np.logical_and(self.sufficient_hets >= n_sufficient_hets,
                                       self.overall_counts.sum(axis=1)>=500)]
        nplots = len(suff_hets) + 2
        pos = self.meQTL.loc[:, 'pos']
        pos = np.asarray(pos, dtype=np.uint64)/x_scale
        min_x = np.min(pos)
        max_x = np.max(pos)
        #range_x = max_x - min_x
        text_x = min_x 
        fig, ax = plt.subplots(nrows=int(ceil(nplots/2.0)) , ncols=2, figsize=(24, 4*nplots/2), 
                               sharey=False, sharex=True, subplot_kw=dict(axisbg='#EEEEEE'))
        ko = 0
        print(int(ceil(len(suff_hets)/2.0)))
        for j in range(2):
            io = 0
            for i in range(int(ceil(len(suff_hets)/2.0))):
                if ko < len(suff_hets):
                    curr = suff_hets.index[ko]
                    adj_pvalue = -1*np.log10(self.aei_pvalues.loc[:, curr])
                    scatter = ax[io, j].scatter(pos, 
                                    adj_pvalue, s=30)
                    ax[io, j].set_ylabel(r'-1*$log_{10}$ p-value', fontsize=15)
                    ax[io, j].set_xlabel('Genomic Position (mb)', fontsize=15)
                    ax[io, j].set_title('AEI plot for %s (N=%i)' %
                            (curr,
                        self.overall_counts.ix[suff_hets.index[ko], 'Nhets']), fontsize=25)
                    # Need to make the text relative positioning
                    #labels = list(self.annot_table.index)
                    #tooltip = mpld3.plugins.PointLabelTooltip(scatter, labels=labels)
                    ko += 1
                    io += 1
                else: 
                    pass
            scatter = ax[-1, j].scatter(pos,
                    -1*np.log10(self.meQTL.ix[: , 'p-value']),
                                     c=self.ld,
                                     s=size_maf, cmap=cm)
            ax[-1, j].set_ylabel('-1*$log_{10}$ p-value', fontsize=15)
            ax[-1, j].set_xlabel('Genomic Position (mb)', fontsize=15)
            ax[-1, j].set_title('%s eQTL plot' % (self.gene_name,), fontsize=25)
            labels = list(self.annot_table.index)
            tooltip = mpld3.plugins.PointLabelTooltip(scatter, labels=labels)
            mpld3.plugins.connect(fig, tooltip,  plugins.LinkedBrush(scatter))
        fig.tight_layout()
        return(fig)
  
        
def combind_aei_with_dosage(aei, dosage, snp):
    """ aei - aei dataframe
    dosage = dosage dataframe
    """
    new_columns = [i[0] for i in aei][::4]
    hets = dosage_round(dosage.ix[snp,:])[new_columns]
    hets = hets[hets == 1]
    return aei.ix[snp, hets.index]


def aei_test_2(full, aei_df, annot_table, gene_snps, gene, num_threshold=5):
    """ Calculates aei for all heterozygous SNPs within a gene across all cis-eQTL SNPs.
    
    Paremeters 
    ----------
    matrixeQTL - full containing multiple population groups
    dosage - dosage object
    aei_df - aellic expression dataframe
    gene_snps - a dictionary with keys being the gene name and 
    
    Returns
    -------
    AEI object

    Refrences
    ---------
    """
    base_pair_to_index = {'A':0, 'C': 1, 'G': 2, 'T': 3}
    gene_i = full.euro.ix[:, 1] == gene
    g_meQTL = full.euro.ix[gene_i, :]
    g_meQTL = pd.merge(annot_table, g_meQTL, right_on="SNP", left_index=True,
            sort=False, how='inner')
    snps_cis = g_meQTL.loc[:, 'SNP']
    new_columns = [i[0] for i in aei_df][::4]
    not_indel = [i for i in gene_snps[gene] if len(annot_table.ix[i, 'a1']) ==1]
    #comb_dosage = pd.concat([full.euro_dos, full)
    hets = dosage_round(full.euro_dos.ix[not_indel,:])[new_columns]
    pvalues_out = np.zeros((len(snps_cis), len(not_indel)), dtype=np.float64)
    sufficient_hets = pd.Series(data=np.repeat(0, len(not_indel)), 
                                index=pd.Index(not_indel),
                                dtype=np.int32)
    aei_ratios= pd.DataFrame(np.zeros((len(new_columns), len(not_indel)), dtype=np.uint32),
                              index=pd.Index(new_columns), columns=pd.Index(not_indel))
    overall_counts = pd.DataFrame(np.zeros((len(not_indel), 2), dtype=np.uint32), 
                                  index=pd.Index(not_indel))
    hets_dict = {}
    maf = calculate_minor_allele_frequency(full.euro_dos.ix[snps_cis, new_columns])
    snp_interest = np.nanargmin(g_meQTL['p-value'])
    snp_interest = g_meQTL['SNP'].iloc[snp_interest]
    g_meQTL.index = pd.Index(g_meQTL['SNP'])
    outliers = []
    for j, i in enumerate(not_indel):
        try:
            REF = base_pair_to_index[annot_table.ix[i, 'a0']]
            ALT = base_pair_to_index[annot_table.ix[i, 'a1']]
        except KeyError:
            print('Indel skipping')
        hets_j = hets.ix[j,:]
        hets_j = hets_j[np.logical_and(hets_j == 1.0, 
                                np.logical_not(pd.isnull(hets_j)))]
        sufficient_hets[i] = len(hets_j)
        if len(hets_j) >= num_threshold:
            # Maybe we should include all samples.  If not maybe make in a sparse index. 
            ref = np.asarray(aei_df.ix[i, hets_j.index].ix[REF::4])
            alt = np.asarray(aei_df.ix[i, hets_j.index].ix[ALT::4], dtype=np.float64)
            overall_counts.ix[j, :] = [np.nansum(ref), np.nansum(alt)]
            allelic_ratio = alt/ref
            aei_ratios.ix[hets_j.index, j] = pd.Series(allelic_ratio, index= hets_j.index)
            allelic_ratio[allelic_ratio > 1] = 1/allelic_ratio[allelic_ratio > 1]
            outliers = np.log2(allelic_ratio) < -3.5
            allelic_ratio = allelic_ratio[np.logical_not(outliers)]
            pvalues_out[:,j] = full.euro_dos.ix[snps_cis, hets_j.index].apply(single_snp_aei_test, axis=1, 
                    args=(outliers, allelic_ratio))
            hets_dict[i] = hets_j.index
        else:
            pvalues_out[:, j] = np.repeat(np.nan, len(snps_cis))
    pvalues_out = pd.DataFrame(pvalues_out, index=pd.Index(snps_cis), 
                               columns=pd.Index(not_indel)) 
    aei_object = AEI_object(pvalues_out, gene, annot_table, sufficient_hets,
            g_meQTL)
    aei_object.ld = calculate_ld(full.euro_dos.ix[snps_cis, new_columns],
            snp_interest)
    overall_counts['Nhets'] = sufficient_hets
    aei_object.overall_counts = overall_counts
    aei_object.ratios = aei_ratios
    aei_object.hets_dict = hets_dict
    aei_object.maf = maf
    aei_object.outliers = outliers
    # This needs to change
    return(aei_object)



def aei_test(matrixeQTL, dosage, aei_df, annot_table, gene_snps, gene, num_threshold=5):
    """ Calculates aei for all heterozygous SNPs within a gene across all cis-eQTL SNPs.
    
    Paremeters 
    ----------
    matrixeQTL - matrixeQTL object
    dosage - dosage object
    aei_df - aellic expression dataframe
    annot_table - snp_annotation table
    gene_snps - a dictionary with keys being the gene name and 

    gene - gene of interest
    :TODO make it test against all genes
    
    Returns
    -------
    AEI object

    Refrences
    ---------
    """
    base_pair_to_index = {'A':0, 'C': 1, 'G': 2, 'T': 3}
    gene_i = matrixeQTL.ix[:, 1] == gene
    g_meQTL = matrixeQTL.ix[gene_i, :]
    g_meQTL = pd.merge(annot_table, g_meQTL, right_on="SNP", left_index=True,
            sort=False, how='inner')
    snps_cis = g_meQTL.loc[:, 'SNP']
    new_columns = [i[0] for i in aei_df][::4]
    not_indel = [i for i in gene_snps[gene] if len(annot_table.ix[i, 'a1']) == 1]
    not_indel = [i for i in not_indel if i in dosage.index and i in
            aei_df.index]
    hets = dosage_round(dosage.ix[not_indel,:])[new_columns]
    pvalues_out = np.zeros((len(snps_cis), len(not_indel)), dtype=np.float64)
    sufficient_hets = pd.Series(data=np.repeat(0, len(not_indel)), 
                                index=pd.Index(not_indel),
                                dtype=np.int32)
    aei_ratios= pd.DataFrame(np.zeros((len(new_columns), len(not_indel)), dtype=np.uint32),
                              index=pd.Index(new_columns), columns=pd.Index(not_indel))
    overall_counts = pd.DataFrame(np.zeros((len(not_indel), 2), dtype=np.uint32), 
                                  index=pd.Index(not_indel))
    hets_dict = {}
    maf = calculate_minor_allele_frequency(dosage.ix[snps_cis, new_columns])
    snp_interest = np.nanargmin(g_meQTL['p-value'])
    snp_interest = g_meQTL['SNP'].iloc[snp_interest]
    g_meQTL.index = pd.Index(g_meQTL['SNP'])
    outliers = []
    for j, i in enumerate(not_indel):
        try:
            REF = base_pair_to_index[annot_table.ix[i, 'a0']]
            ALT = base_pair_to_index[annot_table.ix[i, 'a1']]
        except KeyError:
            print('Indel skipping')
        hets_j = hets.ix[j,:]
        hets_j = hets_j[np.logical_and(hets_j == 1.0, 
                                np.logical_not(pd.isnull(hets_j)))]
        sufficient_hets[i] = len(hets_j)
        if len(hets_j) >= num_threshold:
            # Maybe we should include all samples.  If not maybe make in a sparse index. 
            ref = np.asarray(aei_df.ix[i, hets_j.index].ix[REF::4])
            alt = np.asarray(aei_df.ix[i, hets_j.index].ix[ALT::4], dtype=np.float64)
            overall_counts.ix[j, :] = [np.nansum(ref), np.nansum(alt)]
            allelic_ratio = alt/ref
            aei_ratios.ix[hets_j.index, j] = pd.Series(allelic_ratio, index= hets_j.index)
            allelic_ratio[allelic_ratio > 1] = 1/allelic_ratio[allelic_ratio > 1]
            outliers = np.log2(allelic_ratio) < -3.5
            allelic_ratio = allelic_ratio[np.logical_not(outliers)]
            pvalues_out[:,j] = dosage.ix[snps_cis, hets_j.index].apply(single_snp_aei_test, axis=1, 
                    args=(outliers, allelic_ratio))
            hets_dict[i] = hets_j.index
        else:
            pvalues_out[:, j] = np.repeat(np.nan, len(snps_cis))
    pvalues_out = pd.DataFrame(pvalues_out, index=pd.Index(snps_cis), 
                               columns=pd.Index(not_indel)) 
    aei_object = AEI_object(pvalues_out, gene, annot_table, sufficient_hets,
            g_meQTL)
    aei_object.ld = calculate_ld(dosage.ix[snps_cis, new_columns],
            snp_interest)
    overall_counts['Nhets'] = sufficient_hets
    aei_object.overall_counts = overall_counts
    aei_object.ratios = aei_ratios
    aei_object.hets_dict = hets_dict
    aei_object.maf = maf
    aei_object.outliers = outliers
    return(aei_object)


def calculate_concordance(aei_pvalues, eqtl_pvalues, threshold=0.05):
    """ Returns
    """

    eqtl_pvalues_i = np.nanargmin(eqtl_pvalues)
    print(eqtl_pvalues_i)
    print(aei_pvalues.iloc[eqtl_pvalues_i])
    if aei_pvalues.iloc[eqtl_pvalues_i] <= threshold:
        return(True)
    else:
        return(False)


def calculate_ld(genotypes, snp):
    """ Calculates SNP linkage disquilibrium

    Parameters 
    ----------
    genotypes : a genotype dataframe with polymorphisms going rowise and samples
    columnwise. 
    snp: snp-identifier

    Returns
    -------
    A vector of linkage disequilibrium 

    """
    ld = genotypes.apply(lambda x: pearsonr(genotypes.ix[snp, :], x)[0], axis=1)
    ld[np.isnan(ld.values)] = 0
    ld = np.sqrt(ld ** 2)
    return(ld)


def calculate_minor_allele_frequency(genotypes):
    """
    Parameters
    ----------
    genotypes : a genotype dataframe with polymorphisms going row-wise
    """
    try:
        maf = genotypes.sum(axis=1)/(2* genotypes.shape[1])
    except IndexError:
        maf = genotypes.sum()/(2 * len(genotypes))
    return(maf)


def single_individual_aei_test():
    """
    Parameters
    ----------
    """
    pass


def single_individual_aei_plot(aei_object, dosage, sample):
    """
    Parameters
    ----------
    """
    fig, ax = plt.subplots(nrows=nplots, figsize=(12, 4*nplots), 
                           sharey=False, sharex=True, 
                           subplot_kw=dict(axisbg='#FFFFFF'))
    genotype = dosage.round(dosage.ix[aei_object.ratios.columns, sample])


def subset_meQTL(meQTL, gene_name):
    """
    """
    try:
        index = meQTL.gene == gene_name
        subset = meQTL.ix[index, :]
    except AttributeError:
        # Case where eQTL matrix is already 1 gene
        subset = meQTL
    return(subset)


def add_gene_bounderies(ax, gene_annot, gene_name, x_scale):
    """ Add gene bounderies to an axis
    """
    gene_bounds = grab_gene_location(gene_name, cis_buffer=0)
    path_a = make_rectangle(float(gene_bounds[1])/x_scale, 
            float(gene_bounds[2])/x_scale,
            0, 200)
    patch = patches.PathPatch(path_a, facecolor='grey', alpha=0.25)
    return(patch)


def plot_eQTL(meQTL, gene_name, annotation, dosage, ax=None,
        symbol=None, focus_snp=None, gene_annot=None):
    """ Plot eQTL from a full_matrix object
    """
    subset = subset_meQTL(meQTL, gene_name)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    cm = plt.cm.get_cmap('Blues')
    x_scale= 1e6
    try:
        adj_pv = -1*np.log10(subset.ix[:,4])
        pos = np.asarray(annotation.ix[subset.ix[:, 0], 'pos'], dtype=np.uint32)/x_scale
        dosage_sub = dosage.ix[subset.ix[:, "SNP"],:]
    except IndexError:
        adj_pv = -1*np.log10(subset.ix[:,0])
        pos = np.asarray(annotation.ix[subset.index, 'pos'], dtype=np.uint32)/x_scale
        dosage_sub = dosage.ix[subset.index,:]
    if should_not_plot(dosage):
        dosage_maf =\
                calculate_minor_allele_frequency(dosage_sub)
        dosage_maf = ((150 * dosage_maf) + 20)
    else:
        dosage_maf = np.repeat(1, subset.shape[0])
    if focus_snp:
        snp = focus_snp
    else:
        snp = subset.iloc[np.nanargmax(adj_pv), 0]
    try:
        iix = [i for i, j in enumerate(subset["SNP"]) if j == snp] 
    except KeyError:
        iix = [i for i, j in enumerate(subset.index) if j == snp]
    snpx = pos[iix[0]] 
    snp_pv = adj_pv.iloc[iix[0]]
    color1 = calculate_ld(dosage_sub,
            snp)[dosage_sub.index].values
    print(np.sum(np.isfinite(color1)))
    if ax is None:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18, 6), 
                               sharey=False, sharex=True,
                               subplot_kw=dict(axisbg='#FFFFFF'))
        ax.tick_params(axis='both', which='major', labelsize=24)
        im = ax.scatter(pos, adj_pv, s=dosage_maf, c = color1)
        ax.set_ylabel(r'$-log_{10}$ p-value', fontsize=24)
        ax.set_xlabel(r'Position (Mb)', fontsize=24)
        ax.set_title(r'cis-eQTL for %s' % gene_name, fontsize=30)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        bar = fig.colorbar(im, cax=cbar_ax)
        bar.ax.tick_params(labelsize=18)  
        bar.set_label('r$^{2}$', fontsize=24)
        bar.ax.set_position((0.8125, 0.14, 0.02, 0.725))
        ax.xaxis.set_major_formatter(x_formatter)
        ax.arrow(snpx, snp_pv + 1.3, 0, -1.2, 
                head_width=0.01, head_length=0.1, fc='k',
                ec='k')
        ax.text(snpx-0.05 , snp_pv + 1.5, snp, 
                style='italic', fontsize=16)
        if gene_annot:
            patch = add_gene_bounderies(ax, gene_annot, 
                    gene_name, x_scale)
            ax.add_patch(patch)
        fig.subplots_adjust(right=0.8)
        ax.set_xlim((min(pos) -0.01, max(pos) + 0.01))
        ax.set_ylim((-0.01, max(adj_pv) + 2))
    else:
        pass
    if symbol:
        gene_name = symbol
    ax.set_ylabel(r'$-log_{10}$ p-value', fontsize=24)
    ax.set_xlabel(r'Position (Mb)', fontsize=24)
    ax.set_title(r'eQTL for %s' % gene_name, fontsize=30)
    return(fig)




def test_replication(matrix_a, matrix_b, gene_name):
    index_a = matrix_a.matrix.gene == gene_name
    index_b = matrix_b.matrix.gene == gene_name
    subset_a = matrix_a.matrix.ix[index_a, :] 
    subset_b = matrix_b.matrix.ix[index_b, :]
    try:
        min_a = subset_a.iloc[np.nanargmin(subset_a.ix[:, 4]), :]
        if min_a[2] < 0:
            coef_neg = True
        else: 
            coef_neg = False
        #min_b = subset_b.iloc[np.nanargmin(subset_b.ix[:, 4]), :]
    except ValueError:
        return(False)
        #print(subset_b.ix[: 4])
    #print(subset_b.ix[subset_b.ix[:, 0] == min_a.ix[0], :])
    if subset_b.ix[subset_b.ix[:, 0] == min_a.ix[0], 4].values < 0.05:
        coef_b = (subset_b.ix[subset_b.ix[:, 0] == min_a[0], 2].values[0] < 0)
        #print(subset_b.ix[subset_b.ix[:, 0] == min_a[0], 2])
        if coef_b == coef_neg:
            return(True)
        else:
            return(False)
    else: 
        return(False)


def plot_coefs(matrix_a, matrix_b, genes):
    """
    """
    matrix_a
    pass


def is_cis_eQTL(matrixeQTL, gene_name, num_snps=1):
    """
    """
    pvalues_gene = matrixeQTL.matrix.ix[matrixeQTL.matrix.ix[:, 1] ==
            gene_name, 4]
    log_pvalue = -1 * np.log10(pvalues_gene)
    if np.sum(log_pvalue >= 7) >= num_snps:
        return True
    else:
        return False


def get_genes_with_eqtl(matrixeQTL):
    """
    """
    with_eqtl = {}
    temp = set(matrixeQTL.matrix.gene)
    for i in temp:
        if is_cis_eQTL(matrixeQTL, i):
            with_eqtl.append(i)
        else: pass
    return(with_eqtl)


