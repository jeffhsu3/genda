from matplotlib import figure
import matplotlib.pyplot as plt
from matplotlib.patches import patches
from scipy.stats import linregress, pearsonr
import pandas as pd
import numpy as np
import pysam

def should_not_plot(x):
    """ Check if a potential numpy array exists or not
    """
    if x is None:
        return True
    elif isinstance(x, np.ndarray):
        return x.size == 0
    else:
        return(bool(x))

def plot_eQTL(dosage, highlight_snp=None, paths=None, title='eQTL', chrm=None,
        zoom_in=None, ld_dataframe=None, calculate_LD=True, aei_pvalues=None,
        transcripts=None, plot_exons=True):
    """
    """
    pvalues = dosage.pvalues
    gene_name = dosage.gene_name
    fig = plt.figure(num=None, figsize=(22, 10), dpi=300, facecolor='w',
            edgecolor='k')
    # Todo add another argument
    gtf = pysam.Tabix('/proj/genetics/Projects/shared/Subject\
             Sources/External/Ensembl/Data/Ensembl74/gtf/homo_sapiens/Homo_sapiens.GRCh37.74.gtf.bed.gz')
    commoni = [i for i in dosage.annot.index if not np.isnan(pvalues[i])]
    commonb = [True for i in dosage.annot.index if not np.isnan(pvalues[i])]
    ax = fig.add_subplot(111)
    if highlight_snp:
        scale = dosage.dosages.apply(lambda x:\
                pearsonr(dosage.dosages.ix[highlight_snp, :], x)[0], axis=1)
    else: 
        mindex = np.nanargmin(pvalues)
        dosage.dosages.ix[mindex, :]
        scale = dosage.dosages.apply(lambda x:
               pearsonr(dosage.dosages.ix[mindex, :], x)[0], axis=1)
    scale = scale[dosage.annot.index]
    scale[np.isnan(scale.values)] = 1
    scale = scale ** 2
    cm = plt.cm.get_cmap('winter')
    cax = ax.scatter(np.asarray(dosage.annot.ix[commoni, 'pos'], dtype= np.uint32),
            np.asarray(-1*np.log10(pvalues[commoni])), 
            c=np.sqrt(scale[commoni]), cmap=cm, s=90)
    if should_not_plot(aei_pvalues):
        bar = plt.colorbar(cax)
    else:
        ax = ax.twinx()
        aei_pvalues = pd.Series([1 if np.isnan(i) else i for i in aei_pvalues],
                index=aei_pvalues.index)
        aei_pvalues = np.absolute(-1 * np.log10(aei_pvalues))
        ax2.scatter(np.asarray(dosage.annot.ix[aei_pvalues.index, 'pos'],
            dtype=uint32), np.asarray(aei_pvalues[:]), 
            s=60, c='orange')
        ax2.set_ylim(0, max(aei_pvalues)+3)
    yrange = (0,max(-1*np.log10(pvalues[commoni]))+ 3)
    ax.set_ylim(*yrange) 
    if zoom_in:
        plt.xlim(zoom_in)
    else: 
        plt.xlim((dosage.annot.ix[0, 'pos'], dosage.annot.ix[-1, 'pos']))
    if highlight_snp:
        ax.scatter(np.asarray(dosage.annot.ix[highlight_snp,'pos'],dtype=uint32), 
                np.asarray(-1*np.log10(pvalues[highlight_snp])),
                c='red', s=200, alpha=0.5)
    plt.ylabel('-1 * log10(p-value)', fontsize=30)
    plt.xlabel('Genomic Position', fontsize=30)
    bar.set_label('Linkage Disequillibrium', fontsize=30)
    bar.ax.tick_params(labelsize=20)
    if title:
        plt.title(title, fontsize=40)
    else:
        plt.title(gene_name + ' eQTL', fontsize=40)
    plt.ylabel('-1 * log10(pvalue)', fontsize=30)
    if gene_name and chrm:
        pls = gtf.fetch(chrm, dosage.annot.ix[0, 'pos'], 
                dosage.annot.ix[-1, 'pos'])
        paths = create_path(pls, gene_name)
        enum_i = int(yrange[1]/5)
        if plot_exons:
            if not transcripts:
                transcripts = paths.keys()
            else: pass
            for i, j in paths.iteritems():
                if i in transcripts:
                    plot_transcript(i, paths, ax, y = enum_i,
                            height=yrange[1]/6)
                    enum_i += yrange[1]/5
        else:
            gene_bounds = grab_gene_location(dosage.gene_name, cis_buffer=0)
            path = make_rectangle(gene_bounds[1], gene_bounds[2], 
                    0, yrange[1] + 3)
            patch = patches.PathPatch(path, facecolor='lightblue')
            ax.add_patch(patch)
    plt.tick_params(axis='both', which='major', labelsize=25)
    plt.show()



def plot_matrixeQTL(matrixeQTL, annot, gene_name):
    """
    """
    fig = figure(num=None, figsize=(22, 10), dpi=300, facecolor='w',
            edgecolor='k')
    ax = fig.add_subplot(111)
    eQTL_ax = ax.scatter( annot['pos'], 
            (-1*np.log10(matrixeQTL.ix[matrixeQTL.ix[:, 1]\
                    == gene_name, ])).values)
    ax.set_ylim(0, 20)
    plt.show()
    
        


