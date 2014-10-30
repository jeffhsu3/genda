""" Plotting functions for eQTLs
"""
import numpy as np
#import pandas as pd
import matplotlib
print(matplotlib.get_backend())

import matplotlib.pyplot as plt

from genda.plotting import should_not_plot, add_gene_bounderies
from genda import calculate_minor_allele_frequency, calculate_ld


def var_boxplot(ax, x, colors=None):
    """
    Arguments:
    ----------
    ax - axis object
    """
    for i in x:
        pass


def plot_dosage_by_rsID(chrom, gene, rsID, cov_mat, title = None, ax = None):
    """
    """
    chrom = chrom - 1
    if ax:
        pass
    

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


def plot_eQTL(meQTL, gene_name, annotation, dosage, ax=None,
        symbol=None, focus_snp=None, gene_annot=None):
    """ Plot eQTL from a full_matrix object

    Arguments
    ---------
    meQTL - a matrix eQTL dataframe or a series of pvalues
    gene_name - gene name
    annotation - snp annotation
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
        adj_pv = -1*np.log10(subset.iloc[:,0])
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
        # :TODO fix for both use cases
        #snp = subset.iloc[np.nanargmax(adj_pv), 0]
        snp = subset.index[np.nanargmax(adj_pv)]
    try:
        iix = [i for i, j in enumerate(subset["SNP"]) if j == snp] 
    except KeyError:
        iix = [i for i, j in enumerate(subset.index) if j == snp]
    snpx = pos[iix[0]] 
    snp_pv = adj_pv.iloc[iix[0]]
    color1 = calculate_ld(dosage_sub,
            snp)[dosage_sub.index].values
    plt.ioff()
    if ax is None:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18, 6), 
                               sharey=False, sharex=True,
                               subplot_kw=dict(axisbg='#FFFFFF'))
        ax.tick_params(axis='both', which='major', labelsize=24)
        im = ax.scatter(pos, adj_pv, s=dosage_maf, c = color1)
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
        if should_not_plot(gene_annot):
            patch = add_gene_bounderies(ax, gene_annot, 
                    gene_name, x_scale)
            ax.add_patch(patch)
        else: pass
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
