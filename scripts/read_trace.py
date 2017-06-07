#!/usr/bin/env
import pandas as pd
import numpy as np
import sys
from IPython import embed
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
import random
from collections import defaultdict
# Set seed

try:
    import configparser 
except ImportError:
    import ConfigParser as configparser

# CHUNKSIZE
CHUNKSIZE = 10000
BURN_IN = 400

# COV_mat



def get_genotype(chrom, rsid):
    """:TODO switch to reading HDF5
    """
    geno_path = ('/home/hsuj/lustre/geno/'
            'CCF_1000G_Aug2013_Chr{0}.dose.double.ATB.RNASeq_MEQTL.txt')

    geno_gen = pd.read_csv(geno_path.format(str(chrom)), 
            sep=" ", chunksize = 10000)
    for i in geno_gen:
        if rsid in i.index:
            break
        else: pass
    return(i.ix[rsid, :])


def main(chunk):
    """
    nsamps - number to sample from MCMC posterior
    """
    try:
        chunk = int(chunk)
    except TypeError:
        print('Script requires chunk to be an int')
    files = glob.glob('/home/hsuj/lustre/mmseq_analysis/remapped/*/*.mmseq.trace_gibbs.gz')
    # Split files by genes into different slurm jobs
    # number of individuals to test
    nindv = len(files)
    # 
    nsamps = 10
    best_snps_file = '/home/hsuj/lustre/matrixeQTL_mmseq/best_snps.csv'
    best_snps_columns = pd.read_csv(best_snps_file, header=0, sep=",",
            index_col=0, nrows=2)
    best_snps = pd.read_csv(best_snps_file,
            sep=",", index_col=0, header=None)
    best_snps.columns = best_snps_columns.columns
    phen_file = pd.read_csv('/home/hsuj/Afib/eQTL/Exons/pca_final.csv', 
            sep=",")
    #### Read in genotypes
    by_mcmc_sampling = []
    for i in range(nindv):
        print('Reading in file: {0!s}'.format(str(i)))
        tnames = pd.read_csv(open(files[i], 'rb'), compression='gzip', 
            sep=" ", nrows=1)
        # Sigh this should not be needed
        collided_indexes = []
        new_best_snp_index = []
        for transcript in best_snps.index:
            try:
                collided_indexes.append(tnames.columns.get_loc(transcript))
                new_best_snp_index.append(transcript)
            except KeyError:
                pass
        test = pd.read_csv(open(files[i], 'rb'), compression='gzip', sep=" ",
                usecols=collided_indexes,
                engine='python', encoding='utf-8',
                nrows=1024)
        rows = random.sample(list(test.index[BURN_IN:]), nsamps)
        temp = test.iloc[rows, :].copy()
        temp = temp.iloc[chunk, :]
        temp.name = files[i].split("/")[-1].split(".")[0]
        by_mcmc#_sampling.append(temp)
        """
        temp.columns = [files[i].split("/")[-1].split(".")[0] for trans_i in\
                range(len(collided_indexes))]
        temp.index = pd.Index(np.arange(0, nsamps))
        for k, j in enumerate(new_best_snp_index):
            tdict[j].append(temp.iloc[:, k])
        """
    fig, ax = plt.subplots()
    """
    for i, j in tdict.items():
        j = pd.concat(j, axis=1)
        tdict[i] = j
        ax.plot(j.mean(axis=1), j.std(axis=1) , 'o')
    """
    fig.savefig('/home/hsuj/lustre/output/sd_mean_mcmc.png', dpi=300)
    all_genes_single_sampling = pd.concat(by_mcmc_sampling, axis=1, join='inner')
    all_genes_single_sampling.to_csv('/home/hsuj/lustre/matrixeQTL_mmseq/temp/sampling_{0}.txt'.format(chunk))


    

if __name__ == '__main__':
    # argument 1 the chunk to parse
    main(sys.argv[1])
