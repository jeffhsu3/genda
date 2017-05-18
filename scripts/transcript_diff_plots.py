# Plot CAV1 transcript/ratio eQTL
try:
    import configparser 
except ImportError:
    import ConfigParser as configparser
import sys
import timeit
from collections import defaultdict
#from interval import Interval, IntervalSet
import re

import pandas as pd
from sklearn.decomposition import PCA
import matplotlib
import statsmodels.api as sm
from statsmodels.graphics.regressionplots import plot_fit
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from matplotlib import cm
from matplotlib.pyplot import colorbar
from IPython import embed

plt.style.use('ggplot')
#from statsmodels.formula.api import ols


from genda.AEI import remove_tr_spines
from genda import calculate_minor_allele_frequency, calculate_ld
from genda.plotting import snp_arrow, add_size_legends
from genda.transcripts import (compare_two_transcripts, 
        pairwise_transcript_comparison, 
        get_transcript_ids, 
        Gene, DiffEvent)

from genda.transcripts.get_estimation import count_reads


from genda.eQTL import (plot_eQTL, plot_dosage_by_rsID,
        gene_reference)
from genda.plotting import (create_path,
        coverage_hist, plot_transcript, get_path_max_and_min,
        draw_arrows, draw_arc_label, 
        draw_junction_arcs,)


import pysam


from matplotlib import rcParams
rcParams['axes.labelsize'] = 14
rcParams['xtick.labelsize'] = 12
rcParams['ytick.labelsize'] = 12
rcParams['legend.fontsize'] = 14




def single_apply(snp, covariates, exog):
    '''
    '''
    cov = sm.add_constant(covariates)
    cov['snp'] = snp
    fit = sm.OLS(exog, cov, missing='drop').fit()
    return(fit.pvalues['snp'])





def get_genotype(chrom, rsid):
    """
    """
    geno_path = ('/home/hsuj/lustre/geno/'
            'CCF_1000G_Aug2013_Chr{0}.dose.double.ATB.RNASeq_MEQTL.txt')

    geno_gen = pd.read_csv(geno_path.format(str(chrom)), 
            sep=" ", chunksize = 10000)
    for i in geno_gen:
        if rsid in i.index:
            break
        else: pass
    return(i)


def combine_transcript():
    """
    """
    pass


        

def get_mmseq(redo = False):
    ''' # Move to utils
    '''
    if not redo:
        mmseq = pd.read_pickle('remapped/AF_mmseq.pkl')
        #mmseq_sd = pd.read_pickle('remapped/AF_mmseq_sd.pkl')
        mmseq_uh = pd.read_pickle('remapped/AF_mmseq_uh.pkl')
        mmseq_prop = pd.read_pickle('remapped/AF_mmseq_prop.pkl')
    else:
        raise IOError
    mmseq_sd = None
    #mmseq_uh = None
    #mmseq_prop = None
    return(mmseq, mmseq_sd, mmseq_uh, mmseq_prop)


def apply_normalization(y, uh, ufrac=0.2):
    """ DESEQ style apply normalization
    """
    good = ((uh > 0).sum(axis=1)/float(uh.shape[1]) >= ufrac)
    if good.sum() < 100:
        sys.exit('Not enough unique hit features!')
    y = y.ix[good.values, :]
    row_means = y.mean(axis=1)

    for i in y.columns:
        med_ydiff = np.median(y.ix[:, i] - row_means)
        y.loc[: ,i] = y.loc[:,i] - med_ydiff

    return(y)



def generate_count_frames(gene):
    """ Generate dataframe base on the event type
    """
    if de.event_type == 'AFE':
        pass
    else:
        pass




def main(gene, de, rsid, expr, cov=None):
    """rsid is simply a SNP within the region

    # Refactor transcript order shouldn't matter 
    # :TODO refactor so that it can handle multiple transcripts

    Arguments
    ---------
    gene : genda.transcripts.gene object
    de : diffevent
    covariates : add covriates to the fit
    """
    # :TODO add this to configure parser
    chrom = gene.chrom
    base_path = '/home/hsuj/Afib/'
    sann = pd.read_pickle(base_path + 'ref/snp_annot/{0!s}.pkl'.format(chrom))
    sample_mappings = pd.read_csv('/home/hsuj/lustre/sample_mapping.txt',
            sep="\t", index_col=0,)
    sample_mappings['new_name'] = [i.split("/")[-1].rstrip(".bam") \
            for i in sample_mappings.ix[:, 2]]
    coverages = pd.read_csv(('/home/hsuj/Afib/eQTL/gene_data.txt'),
            sep = ",", index_col=0)
    coverages = coverages.sum(axis=0)
    # Normalize
    ## PCA covariates
    cov_mat = pd.read_csv('/home/hsuj/Afib/eQTL/Exons/pca_final.csv', sep=",",
            index_col=0)
    pheno2 = cov_mat
    # Need at least 2 transcripts to compare.   
    #graf = gr.genotype_reader_h5py('/home/hsuj/lustre/AF_miso_AFE.hdf')
    path_dict = gene.transcripts
    ## Handle genotypes ############################################
    sann['pos'] = sann['pos'].astype(int)
    i = get_genotype(chrom, rsid)
    geno = i.ix[:, expr.columns]
    gaf = calculate_minor_allele_frequency(geno)
    gaf = ((gaf >= 0.10) & (gaf <= 0.90)).values
    geno = geno.ix[gaf,:]
    ################################################################

    plot_dict = {}
    fig2, ax2  = plt.subplots(figsize=(10, 10), nrows=4,
            sharex=True)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)  
    if de.exon_num[0]:
        sea = 0
    else:
        sea = 1
    eoi = [de.exon_num[sea]]
    cigar_skipped = getattr(de, 'cigar{0!s}'.format(sea*1 + 1))
    if len(cigar_skipped) > 3:
        eoi.append(eoi[0] + 1)
    else: pass
    cpath = path_dict[de.tid[sea]]
    eoi_intron_lengths = [i[1] for i in cigar_skipped if i[0] == 3]
    eoi2 = getattr(de,'cigar{0!s}'.format(2 - 1*sea))
    eoi2 = [abs(i[1]) for i in eoi2 if i[0] == 3 ]
    # Move generation of to plot into diffevents?
    to_plot1 = [i for i in cpath if (i[2] >= min(eoi) - 1) and (i[2] <= max(eoi)+ 1)]
    plot_dict[de.tid[sea]] = to_plot1

    try:
        to_plot2 = [i for i in path_dict[de.tid[1-sea]] if i[2] in
                [de.exon2[0], de.exon2[1]]]
    except TypeError:
        # For alternate first exons
        #to_plot2 = [i for i in path_dict[de.tid[1-sea]]]
        to_plot2 = [[i[1], i[2], 2] for i in getattr(de,
            'cigar{0!s}'.format(2-1*sea)) if i[0] == 0]
    plot_dict[de.tid[1-sea]] = to_plot2
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    xmin = min(to_plot1[0][0], to_plot2[0][0]) 
    try:
        xmax = max(to_plot1[2][1], to_plot2[1][1])
    except IndexError:
        # AFE exception 
        xmax = to_plot1[1][1]
    sample_mappings = sample_mappings.ix[
            np.logical_not(sample_mappings.index.duplicated()),:]
    nsamples =  sample_mappings.shape[0]
    buffer_bp = 0 
    t1 = de.tid[0]
    t2 = de.tid[1]
    all_juncs = [eoi_intron_lengths, eoi2]
    SCN5A = 'ENSG00000183873'
    transcript_min = min([i[0] for i in gi.transcripts[de.tid[0]]])
    transcript_max = max([i[1] for i in gi.transcripts[de.tid[1]]])
    out_frame = pd.DataFrame({de.tid[sea]:np.zeros(nsamples),
        de.tid[sea-1]: np.zeros(sample_mappings.shape[0]),
        }, 
        index = sample_mappings.index)
    ofc = pd.DataFrame(
            np.zeros((sample_mappings.shape[0],len(gi.transcripts))),
            index=sample_mappings.index, columns=gi.transcripts.keys())
    # Move this to function
    for bami in sample_mappings.index:
        # :TODO Make this more generalizable
        fname = sample_mappings.ix[bami,2].split("/")[-1]
        bi = '/home/hsuj/lustre/AFdata/{0!s}'.format(fname)
        bamf = pysam.Samfile(bi)
        bamiter = bamf.fetch('chr' + chrom,  xmin-buffer_bp,
                xmax + buffer_bp)
        #################### Getting intron junctions counts
        # :TODO convert to cython
        # Convert this to series
        intron_counts = np.zeros(2, dtype=np.int32)
        for i in bamiter:
            #exons = [j[1] for j in i.cigar if j[0] == 0]
            introns = [j[1] for j in i.cigar if j[0] == 3]
            # This depends on there not being other exact intron sizes
            for knum, ijunc in enumerate(all_juncs):
                try:
                    matches = [zi for zi in ijunc if zi in introns]
                except TypeError:
                    matches = [zi for zi in [ijunc] if zi in introns]
                if len(matches) > 0:
                    intron_counts[knum] += 1
        out_frame.ix[bami, 0:3] = intron_counts
        bleh = gi.transcripts[de.tid[0]][:-3]
        hmm = [bleh, gi.transcripts[de.tid[1]]]
        for i, j in zip(hmm, de.tid):
            bamiter = bamf.fetch('chr' + chrom,  transcript_min,
                    transcript_max)
            ofc.ix[bami ,j] = count_reads(i, bamiter)
    # Filter low counts
    count_threshold = 0 
    out_frame = out_frame.ix[out_frame.sum(axis=1) > count_threshold, :]
    read_length = 100
    intron_factor = 1
    # Refactor this to a function
    propi = (out_frame.ix[:, t1]+1)/intron_factor/(out_frame.ix[:,t1]/intron_factor +
                out_frame.ix[:, t2])
    X = cov_mat.T
    X = sm.add_constant(X)
    from lin_test import test_transcript
    # Let's get all the SNPs that fall within a certain region
    # Adding cav1_beta1 for plotting purposes
    de.tid = [de.tid[0], de.tid[1]]
    out_frame.columns = ['ENST00000502471', 'ENST00000033079',]
    for i, j in enumerate(de.tid):
        ax2[0] = plot_transcript(j, plot_dict, ax2[0], y=i*2.5, 
                height=2.)
        t_xmin = min([k[0] for k in plot_dict[j]])
        t_xmax = max([k[1] for k in plot_dict[j]])
        ax2[0].hlines((i*2.5 + 2) - 1, t_xmin, t_xmax, colors='darkgrey', lw=2)
        ax2[0].xaxis.set_major_formatter(x_formatter)
    ax2[0].get_yaxis().set_ticks([])
    ax2[0] = remove_tr_spines(ax2[0])
    goi = geno
    goi = goi.ix[:, out_frame.index]
    gfits = goi.apply(test_transcript, axis=1, args=(X, propi[X.index]))
    pvalues = [i.pvalues['geno'] for i in gfits]
    best_snp = 'rs17171731'
    pvalues = pd.Series(pvalues, index=geno.index)
    color = plt.rcParams['axes.color_cycle'][0]
    example_ylims = []
    for i in range(3):
        if i == 0: geno_string = sann.ix[best_snp, 'a0'] * 2
        elif i ==1:
            geno_string = sann.ix[best_snp, 'a0'] +\
                    sann.ix[best_snp, 'a1']
        elif i == 2: geno_string = sann.ix[best_snp, 'a1'] * 2
        hist = np.zeros(xmax - xmin, dtype=np.uint64)  
        c_geno = (goi.ix[best_snp, :] == i)
        # Random from out_Frame
        try:
            b_example = np.random.choice(goi.columns[c_geno.values], size=1)[0]
        except ValueError:
            continue
        het = pysam.Samfile(sample_mappings.ix[b_example, 2])
        het_bamf = het.fetch('chr' + str(chrom), xmin, xmax)
        color = plt.rcParams['axes.color_cycle'][i]
        for read in het_bamf:
            coverage_hist(read, hist, xmin)
        het_bamf = het.fetch('chr' + str(chrom), xmin, xmax)
        het_bamf = het.fetch('chr' + str(chrom), xmin, xmax)
        hist = 1e3 * hist/coverages[b_example]
        ax2[i + 1].plot(np.linspace(xmin, xmax, num=len(hist)),
                    hist, color)  
        ax2[i + 1].fill_between(np.arange(xmin, xmax),0, hist, facecolor=color)
        example_ylims.append(np.max(hist))
        # Need this to draw between every single one
        
        for tran_i in de.tid:
            junc_norm = 1e3 * out_frame.ix[b_example, tran_i]/coverages[b_example]
            ax2[i + 1] = draw_junction_arcs(plot_dict[tran_i], hist, xmin, ax2[i+1], 
                    color=color, text=junc_norm, y_buffer=np.max(hist) *0.20)
        ax2[i + 1].set_ylabel('{0} Genotype\n RPKM'.format(geno_string))
        #from lin_test import _temp_plot
        # Resave the pickeld file with correct int type
    example_ylim = max(example_ylims) * 1.2
    for i in range(3):
        ax2[i + 1].set_ylim((0, example_ylim))
    #dfmean = (out_frame - out_frame.mean())/(out_frame.max() - out_frame.min())
    #pcafit = pca.fit(dfmean)
    ax2[0].text((xmax-xmin)/2, ax2[0].get_ylim()[1]- 1, str(min(pvalues)))
    ax2[0].axvline(sann.ix[best_snp, 'pos'],color='r', linestyle='solid')
    ax2[0].set_title('{0!s}'.format(gene.symbol))
    ax2[-1].set_xlabel('Position')
    embed()
    fig, ax = plt.subplots(nrows=3, sharex=True)
    ax[0] = plot_eQTL(pvalues, 'FAM13B', sann, goi, ax=ax[0], 
            focus_snp='rs17171731') 
    plt.tight_layout()
    fig.savefig(base_path +\
            'eQTL/graphs/{0!s}_cis_eqtl_transcript.png'.format(gene.symbol))
    gr = gene_reference(chrom=5, gene=de.tid[0],
            rsID = best_snp)
    fig, ax = plt.subplots(figsize=(6, 10),nrows=2)
    gr = gene_reference(chrom=5, gene=de.tid[0],
            rsID = best_snp)
    ax[0], pv_1 = plot_dosage_by_rsID(gr, goi, X, out_frame.ix[X.index,:].T, ax=ax[0])
    ax[0].set_title(de.tid[0])
    gr = gene_reference(chrom=5, gene=de.tid[1],
            rsID = best_snp)
    ax[1], pv_2 = plot_dosage_by_rsID(gr, goi, X, out_frame.ix[X.index,:].T, ax=ax[1])
    ax[1].set_title(de.tid[1])
    plt.tight_layout()
    fig.savefig(base_path +\
            'eQTL/graphs/CAV1_bestfig.png')
    fig2.savefig(base_path + 'eQTL/graphs/{0!s}_transcript.png'.format(gene.symbol))
    '''
    fig, ax = plt.subplots()
    ax = plot_dosage_by_rsID(gr, goi, X2, prop3[X2.index].T, ax=ax)
    fig.savefig(base_path +\
            'eQTL/graphs/prop3.png')
    '''
    embed()
    return(propi)


if __name__ == '__main__':
    # Need to split by chromosome for cluster and iterate 
    # across all genes
    config = configparser.RawConfigParser() 
    config.read(sys.argv[1])  
    gtf_path = config.get('annotations', 'gene_annot')
    gtf = pysam.Tabixfile(gtf_path)
    expr = pd.read_csv('/home/hsuj/Afib/eQTL/Exons/mmseq_counts_normed.csv',
            sep="\t", index_col=0)
    goi_full = pd.read_csv('/home/hsuj/Afib/eQTL/Exons/best_snps.csv', sep=",",
            index_col=0)
    #maybe use mpi?
    temp = goi_full.query('ensembl == "{0}"'.format('ENSG00000031003'))
    #gg = goi_full.groupby('ensembl')
    gg = temp.groupby('ensembl')
    out_frames = []
    for i, j in gg:
        gene = j.ensembl[0]
        tlist, symbol, chrom, pos = get_transcript_ids(gene, ref37=True)
        expressed_transcripts = j.index
        if len(expressed_transcripts) >= 2:
            gi = Gene(gene, chrom=chrom, 
                    start=int(pos[0]), end=int(pos[1]), symbol=symbol)
            gi.get_transcripts(gtf, buffer=500000)
            matching_exons, de = compare_two_transcripts(tlist[0], tlist[6],
                    gi.transcripts)
            # MISO output
            #diff_events = re.sub("['()\ ]", '', j['diff_event'])
            #diff_events = diff_events.split(",")
            #gi.filter_transcripts(expressed_transcripts)
            #of_interest = pairwise_transcript_comparison(gi.transcripts)
            out_frames.append(main(gi, de[0], 'rs17171731', expr))
