import ConfigParser
import sys
from decimal import Decimal
#from interval import Interval, IntervalSet
from collections import defaultdict
import Queue
import multiprocessing

import pandas as pd
import matplotlib
import statsmodels.api as sm
import numpy as np


matplotlib.use('Agg')
from IPython import embed
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.pyplot import colorbar

plt.style.use('ggplot')
from statsmodels.formula.api import ols


#import emcee

from genda.AEI import remove_tr_spines
from genda import calculate_minor_allele_frequency, calculate_ld
from genda.plotting import snp_arrow, add_size_legends
from genda.transcripts import (compare_two_transcripts, 
        pairwise_transcript_comparison, 
        get_transcript_ids, 
        Gene, DiffEvent, generate_to_plot_tuple)


from genda.eQTL import (plot_eQTL, plot_dosage_by_rsID,
        gene_reference)
from genda.plotting import (create_path,
        coverage_hist, plot_transcript, get_path_max_and_min,
        draw_arrows)


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



def main(gene, de, rsid, expr, cov=None):
    """rsid is simply a SNP within the region

    # Refactor transcript order shouldn't matter 

    Arguments
    ---------
    gene : genda.transcripts.gene object
    covariates : add covriates to the fit
    """

    '''
    pheno = pd.read_csv('/home/hsuj/Afib/eQTL/pheno_for_miso.txt',
            sep="\t", header=None, skiprows=1, index_col=1)
    '''
    pheno2 = pd.read_csv('gene_pheno_eQTL_april.txt', sep=",", index_col=0)
    new_col = [i.replace('.', '-') for i in pheno2.columns]
    pheno2.columns = new_col
    pheno2 = pheno2.T
    srs = 'ENST00000453840'
    # Normalize
    ## PCA covariates
    base_path = '/home/hsuj/Afib/'
    # Need at least 2 transcripts to compare.   
    path_dict = gene.transcripts
    chrom = gene.chrom
    sann = pd.read_pickle(base_path + 'ref/snp_annot/{0!s}.pkl'.format(chrom))
    sann['pos'] = sann['pos'].astype(int)
    i = get_genotype(chrom, rsid)
    geno = i.ix[:, expr.columns]
    gaf = calculate_minor_allele_frequency(geno)
    gaf = ((gaf >= 0.10) & (gaf <= 0.90)).values
    geno = geno.ix[gaf,:]
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
    cpath = path_dict[de.transcript_ids[sea]]
    eoi_intron_lengths = [i[1] for i in cigar_skipped if i[0] == 3]
    eoi2 = getattr(de,'cigar{0!s}'.format(2 - 1*sea))[0][1]
    # Move generation of to plot into diffevents?
    # rough size normalization factor
    norm_fact1 = (float(sum([i[1] - i[0] for i in to_plot1])) / sum([i[1] - i[0] for i
        in to_plot2]))
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    xmin = min(to_plot1[0][0], to_plot2[0][0]) 
    xmax = max(to_plot1[2][1], to_plot2[1][1])
    sample_mappings = sample_mappings.ix[
            np.logical_not(sample_mappings.index.duplicated()),:]
    nsamples =  sample_mappings.shape[0]
    out_frame = pd.DataFrame({de.transcript_ids[0]:np.zeros(nsamples),
        de.transcript_ids[1]: np.zeros(sample_mappings.shape[0])}, 
        index = sample_mappings.index)
    buffer_bp = 0 
    t1 = de.transcript_ids[0]
    t2 = de.transcript_ids[1]
    for bami in sample_mappings.index:
        # :TODO Make this more generalizable
        fname = sample_mappings.ix[bami,2].split("/")[-1]
        bi = '/home/hsuj/lustre/AFdata/{0!s}'.format(fname)
        bamf = pysam.Samfile(bi)
        bamiter = bamf.fetch('chr' + chrom,  xmin-buffer_bp,
                xmax + buffer_bp)
        # Conver this to cython
        c0 = 0
        c1 = 0
        cnot = 0
        for i in bamiter:
            #start = i.pos
            exons = [j[1] for j in i.cigar if j[0] == 0]
            introns = [j[1] for j in i.cigar if j[0] == 3]
            # Probably need to grab exact positions even though we are fetching
            # in small region
            matches1 = [zi for zi in eoi_intron_lengths if zi in introns]  
            try:
                matches2 = [zi for zi in eoi2 if zi in introns]
            except TypeError:
                matches2 = [zi for zi in [eoi2] if zi in introns]
            if (len(matches1) > 0): c0 += 1
            elif len(matches2) > 0: c1 += 1
            else: cnot +=1
        out_frame.ix[bami, t1] = c0 
        out_frame.ix[bami, t2] = c1
    # Filter low counts
    count_threshold = 100
    out_frame = out_frame.ix[out_frame.sum(axis=1) > count_threshold, :]
    read_length = 100
    if to_plot2 > 2:
        intron_factor = 3
    else: intron_factor = 1
    propi = (out_frame.ix[:, t1]+1)/intron_factor/(out_frame.ix[:,t1]/intron_factor +
                out_frame.ix[:, t2])
    '''
    bii = (sann['pos'] > xmin) &\
            (sann['pos'] < xmax) 
    bii = sann.index[bii.values]
    goi = geno.ix[bii, :]
    '''
    X = test.ix[srs, out_frame.index]
    X = sm.add_constant(X)
    X2 = X.copy()
    X2 = X2.join(pheno2.ix[:,0:5], how='inner')
    X2['prop'] = propi[X2.index]
    X2['fullsum'] = out_frame.ix[X2.index].sum(axis=1)
    try:
        prop_model =\
                ols('prop~sexFemale+ENST00000453840+MDS4+MDS1+MDS3+MDS2+fullsum', 
                data=X2,
                missing='drop').fit()
    except ValueError:
        embed()
    fig, ax = plt.subplots(figsize=(6,6))
    fig = sm.graphics.plot_partregress("prop", "ENST00000453840", ['sexFemale',
        'MDS4', 'fullsum'], data=X2, ax=ax, obs_labels=False)
    ax.text(0.5, ax.get_ylim()[1] - 0.02, 
            'p-value: %.2E' % Decimal(prop_model.pvalues['ENST00000453840']),
            size=12)
    ax.set_xlabel('SRSF10 expression')
    ax.set_ylabel('{0!s} included exon / skipped exon proportion'.format(gene.symbol))
    ax.set_title('')
    plt.tight_layout()
    fig.savefig(base_path + 'eQTL/graphs/{0!s}_srsf10_fit.png'.format(gene.symbol))
    print(prop_model.pvalues)
    from lin_test import test_transcript
    # Let's get all the SNPs that fall within a certain region



    if gene.symbol == 'CAST' or gene.symbol == 'GDAP1L1':
        for i, j in enumerate(de.transcript_ids):
            ax2[0] = plot_transcript(j, plot_dict, ax2[0], y=i*2.5, 
                    height=2.)
            ax2[0].hlines((i*2.5 + 2) - 1, xmin, xmax, colors='darkgrey', lw=2)
            ax2[0].xaxis.set_major_formatter(x_formatter)
        ax2[0] = remove_tr_spines(ax2[0])
        ax2[0].set_xlim((xmin, xmax))
        ax2[0].set_ylim((-0.5, 2*2.5 + 0.5))
        goi = geno
        goi = goi.ix[:, out_frame.index]
        gfits = goi.apply(test_transcript, axis=1, args=(X, propi))
        pvalues = [i.pvalues['geno'] for i in gfits]
        best_snp = geno.index[np.nanargmin(pvalues)]
        pvalues = pd.Series(pvalues, index=geno.index)
        print(gfits[np.nanargmin(pvalues)].pvalues)
        color = plt.rcParams['axes.color_cycle'][0]
        embed()
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
            het = pysam.Samfile(sample_mappings.ix[b_example,2])
            het_bamf = het.fetch('chr' + str(chrom), xmin, xmax)
            color = plt.rcParams['axes.color_cycle'][i]
            for read in het_bamf:
                coverage_hist(read, hist, xmin)

            ax2[i + 1].plot(np.linspace(xmin, xmax, num=len(hist)),
                        hist, color)  
            ax2[i + 1].fill_between(np.arange(xmin, xmax),0, hist, facecolor=color)
            ax2[i + 1].set_ylim((np.min(hist), np.max(hist) + 0.2 * np.max(hist)))
            try:
                ax2[i + 1].text((xmax + xmin)/2, np.max(hist),
                        str(out_frame.ix[b_example, 1]))
            except KeyError:
                pass
            ax2[i + 1].set_ylabel('{0} Genotype'.format(geno_string))
            #from lin_test import _temp_plot
            # Resave the pickeld file with correct int type
        ax2[0].text((xmax-xmin)/2, ax2[0].get_ylim()[1]- 1, str(min(pvalues)))
        ax2[0].axvline(sann.ix[best_snp, 'pos'],color='r', linestyle='solid')
        ax2[0].set_title('{0!s}'.format(gene.symbol))
        ax2[-1].set_xlabel('Position')
        fig2.savefig(base_path + 'eQTL/graphs/{0!s}_transcript.png'.format(gene.symbol))
        out_frame.columns = ['{0} IE'.format(gene.symbol), 
                            '{0} SE'.format(gene.symbol)]
    return(propi)



if __name__ == '__main__':
    # Need to split by chromosome for cluster and iterate 
    # across all genes
    config = ConfigParser.RawConfigParser() 
    config.read(sys.argv[1])  
    gtf_path = config.get('annotations', 'gene_annot')
    gtf = pysam.Tabixfile(gtf_path)
    '''
    test, _, uh, prop = get_mmseq(redo=False)
    y = apply_normalization(test, uh)
    npca = 15
    pca = PCA(n_components = npca)
    pcafit = pca.fit(test)
    test = y; del y
    '''
    sample_mappings = pd.read_csv('/home/hsuj/lustre/sample_mapping.txt',
            sep="\t", index_col=0,)
    debug = 0

    bd = {}
    gtf_iter = gtf.fetch('21')
    cur_gene = None
    paths = {}
    for i in gtf_iter:
        i = i.split("\t")
        if i[6] == 'miRNA': pass
        elif i[7] == 'exon':
            geni = i[9].split(";")
            gene_id = geni[0].lstrip(' gene_id "').rstrip('"')
            itrans = geni[1].lstrip(' transcript_id "').rstrip('"')
            exon_number = (geni[2].lstrip(' exon_number "').
                    rstrip('"'))
            if cur_gene == gene_id:
                try:
                    bd[gene_id][itrans].append([int(i[1]), int(i[2]),
                        int(exon_number)])
                except KeyError:
                    pass
            elif cur_gene == None:
                cur_gene = gene_id
                bd[cur_gene] = defaultdict(list)
            else:
                # Get rid of duplicates
                debug +=1
                if debug > 100: 
                    break
                else: pass
                cur_gene = gene_id
                bd[cur_gene] = defaultdict(list)
                try:
                    bd[gene_id][itrans].append([int(i[1]), int(i[2]),
                        int(exon_number)])
                except KeyError:
                    pass
        else: pass
    buffer = 5000
    for i, j in bd.iteritems():
        if len(j) <= 1:
            pass
        else:
            test = pairwise_transcript_comparison(j)
            if len(test) > 0:
                for event_c in test:
                    for de in event_c:
                        fetch_region = (de.start , de.end)
                        tp1, tp2, eoi, eoi2, = generate_to_plot_tuple(de, j)
                        embed()
    '''
    for i, j in goi_full.iterrows():
        tlist, symbol, chrom, pos = get_transcript_ids(i)
        print(symbol, chrom, pos)
        expressed_transcripts = test.ix[tlist,:].mean(axis=1) > -8
        expressed_transcripts =\
                expressed_transcripts.index[expressed_transcripts.values]
        print(pos[0], pos[1])
        if len(expressed_transcripts) >= 2:
            gi = Gene(i, chrom=chrom, 
                    start=int(pos[0]), end=int(pos[1]), symbol=symbol)
            gi.get_transcripts(gtf, buffer=5000000)
            #gi.filter_transcripts(expressed_transcripts)
            #of_interest = pairwise_transcript_comparison(gi.transcripts)
            embed()
            _, de = compare_two_transcripts(j.diff_event[0], j.diff_event[1],
                    gi.transcripts)
            de = de.filter('skipped_exon')
            if symbol == 'LRRFIP1':
                de = de[3]
            else:
                de = de[0]
            out_frames.append(main(gi, de, goi_full.ix[i, 'snps'], test))
            new_header.append(symbol)
        else: 
            print('{0} has no highly expressed transcripts'.format(symbol))
    full = pd.concat(out_frames, axis=1)
    full.columns = new_header
    full = full.fillna(full.mean())
    import seaborn as sns; sns.set()
    from sklearn import preprocessing
    nfull =  preprocessing.scale(full)
    nfull = pd.DataFrame(nfull)
    nfull.columns = full.columns
    g = sns.clustermap(nfull)
    g.savefig( '/home/hsuj/Afib/eQTL/graphs/clustermap.png')
    '''
