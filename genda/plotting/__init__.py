from matplotlib import figure
form matplotlib import pyplot as plt
from scipy.stats import linregress, pearsonr
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
    fig = figure(num=None, figsize=(22, 10), dpi=300, facecolor='w',
            edgecolor='k')
    # Todo add another argument
    gtf = pysam.Tabix('/proj/genetics/Projects/shared/Subject\
             Sources/External/Ensembl/Data/Ensembl74/gtf/homo_sapiens/Homo_sapiens.GRCh37.74.gtf.bed.gz')
    commoni = [i for i in dosage.annot.index if not np.isnan(pvalues[i])]
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
