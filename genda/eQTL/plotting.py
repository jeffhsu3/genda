""" Plotting functions for eQTLs
"""

import numpy as np
import pandas as pd
import matplotlib
import statsmodels.api as sm
import matplotlib.pyplot as plt
from statsmodels.graphics.regressionplots import abline_plot

from genda.plotting import (snp_arrow)
from genda import calculate_minor_allele_frequency, calculate_ld


class gene_reference(object):
    """
    """
    def __init__(self, chrom, gene, rsID = None):
        """
        """
        self.chrom = chrom
        self.gene = gene
        self.rsID = rsID


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


def plot_dosage_by_rsID(gene_reference, dos, cov_mat, counts,
        title = None, ax = None,
        additional_covar=None, 
        adjx=True, adjy=True):
    """
    Arguments:
    ---------
    gene_reference - a gene reference object
    meqtl - a list of matrix-eQTL objects, one for each chromosome
    cov_mat - covariate matrix
    counts - counts
    additional_covar - a matrix of same rows as cov_mat to add to the model
    """
    gr = gene_reference
    cov_mat_t = cov_mat.copy(deep=True)
    try:
        geno = dos.ix[gr.rsID, cov_mat_t.index]
    except pd.core.indexing.IndexingError:
        geno = dos.ix[cov_mat_t.index]
    cov_mat_t[gr.rsID] = geno
    c = counts.ix[gr.gene, cov_mat_t.index]
    if adjx:
        results = sm.OLS(geno, cov_mat).fit()
        adj_dos_mat = geno -\
                np.dot(results.params, cov_mat.T)
    else: 
        adj_dos_mat = dos.ix[gr.rsID,:]

    if adjy:
        results = sm.OLS(c, cov_mat).fit()
        adj_counts = c - np.dot(results.params, cov_mat.T)
        const = results.params.const
    else:
        adj_counts = counts.ix[gr.gene, cov_mat_t.index]
        const = 0
    # Need to grab original genotypes
    colors = []
    # Make this into a function
    color_dict = np.linspace(0, 1, 3)
    for i in geno:
        if i <= 0.5:
            colors.append(color_dict[0])
        elif i > 0.5 and i <= 1.5:
            colors.append(color_dict[1])
        else:
            colors.append(color_dict[2])
    if ax:
        ax_orig = True
    else:
        ax_orig = None
        fig, ax = plt.subplots(nrows=1, ncols=1, sharey=False,
                sharex=False, subplot_kw=dict(axisbg='#FFFFFF'))
    ax.scatter(adj_dos_mat, adj_counts + const, s=50,
             c = colors)
    if title:
        ax.set_title('%s partial regression\non %s' % (title, gr.rsID))
    else: pass
    ax.set_ylabel('$log_{2}$ CPM')
    ax.set_xlabel('Fitted Dosages')
    xticks = ax.get_xticks()
    yticks = ax.get_yticks()
    ax.set_xticks(xticks[1::2])
    ax.set_yticks(yticks[1::2])
    fitted_line = sm.OLS(adj_counts, adj_dos_mat).fit()
    abline_plot(const, fitted_line.params[0], color='k', ax=ax)
    test = sm.OLS(c, cov_mat_t).fit()
    if test.params[gr.rsID] > 0:
        annot_y = - 1
    else:
        annot_y = 1
    yrange = yticks[-1] - yticks[0]
    ax.text(xticks[0] + 0.025 , yticks[annot_y] + annot_y/2*yrange/5, 
            '$R^{2}$=%s' % str(test.rsquared)[0:4], 
            style='italic')
    if ax_orig:
        return(ax, test)
    else:
        return(fig, test)





def plot_eQTL(meQTL, gene_name, annotation, dosage, ax=None,
        symbol=None, focus_snp=None, gene_annot=None, **kwargs):
    """ Plot eQTL from a full_matrix object
    Arguments
    ---------
    meQTL - a matrix eQTL dataframe or a series of pvalues
    gene_name - gene name
    annotation - snp annotation dataframe, index is rsID
    dosage - a dosage dataframe
    ax - axis to plot into
    """
    subset = subset_meQTL(meQTL, gene_name)

    if isinstance(subset.index, pd.core.index.MultiIndex):
        subset.index = subset.index.get_level_values('SNP')
    else: pass

    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    x_scale= 1e6
    try:
        adj_pv = -1*np.log10(subset.ix[:,'p-value'])
    except IndexError:
        adj_pv = -1*np.log10(subset.iloc[:,0])
    except pd.core.indexing.IndexingError:
        adj_pv = -1*np.log10(subset.iloc[:])
    try:
        pos = np.asarray(annotation.ix[subset.index, 'pos'], 
                dtype=np.double)/x_scale
    except KeyError:
        pos = np.asarray(annotation.ix[subset.index, 1], 
                dtype=np.double)/x_scale
    dosage_sub = dosage.ix[subset.index,:]
    print('subset shape')
    print(subset.shape)

    dosage_maf =\
            calculate_minor_allele_frequency(dosage_sub)
    dosage_maf = ((150 * dosage_maf) + 20)
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
    # Need this since pos is a numpy array not pandas series
    snpx = pos[iix[0]] 
    snp_pv = adj_pv.iloc[iix[0]]
    color1 = calculate_ld(dosage_sub,
            snp)[dosage_sub.index].values
    if ax is None:
        ax_orig = False 
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16, 6), 
                               sharey=False, sharex=True,
                               subplot_kw=dict(axisbg='#FFFFFF'))
        fig.tight_layout() 
        fig.subplots_adjust(right=0.8, bottom=0.2)
        #bar.ax.set_position((0.85, 0.14, 0.02, 0.725))
    else:
        ax_orig = True
    ax.set_xlim((min(pos) -0.01, max(pos) + 0.01))
    ylim = (max(adj_pv) + max(adj_pv/6.0))
    ax.set_ylim((-0.01,  ylim))
    ax.xaxis.set_major_formatter(x_formatter)
    ### Actual scatter #############################
    im = ax.scatter(pos, adj_pv, s=dosage_maf, c = color1)
    #:TODO make the arrow into a funciton
    if focus_snp:
        ax = snp_arrow(snpx, snp_pv, snp, ax)
    else: pass
    ax.set_ylabel(r'$-log_{10}$ eQTL p-value')
    ax.set_xlabel(r'Position (Mb)')
    if symbol:
        gene_name = symbol
    if ax_orig:
        return(ax)
    else: 
        cbar_ax = fig.add_axes([0.87, 0.15, 0.05, 0.7])
        bar = fig.colorbar(im, cax=cbar_ax)
        bar.ax.tick_params(labelsize=18)  
        bar.set_label('r$^{2}$')
        return(fig)
