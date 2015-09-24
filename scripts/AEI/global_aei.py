import sys
import requests

import ConfigParser
#from multiprocessing.dummy import Pool as ThreadPool

from IPython import embed
import pandas as pd
import numpy as np

from genda import calculate_minor_allele_frequency, calculate_ld
from genda.AEI import AEI, dosage_round
pd.options.mode.chained_assignment = None


def get_aei(aei_path):  
    # Currently 
    global sample_mix_ups
    aei = pd.read_pickle(aei_path)
    new_columns =\
            [i[0].split("/")[-1].rstrip('.MD.bam').rstrip('\
            LA').upper().replace("DN", "D").replace("D ", "D")\
            for i in aei.columns][::4]
    new_columns = zip(np.repeat(new_columns, repeats=4),
            [0,1,2,3]*len(new_columns))
    aei.columns = pd.MultiIndex.from_tuples(new_columns)
    return(aei)


def get_snp_annot(tabix_iter):
    """
    """
    snp_id = []
    pos = []
    a0 = []
    a1 = []
    for i in tabix_iter:
        i = i.split("\t")
        snp_id.append(i[3])
        pos.append(i[1])
        a0.append(i[6])
        a1.append(i[7])
    annot_out = pd.DataFrame({'pos': pos, 'a0': a0, 'a1':a1},
                             index = snp_id)
    return(annot_out)


def get_snp_annot2(tabix_file):
    outf = pd.read_csv(tabix_file, sep="\t", usecols=[3,0,1,6,7], index_col=2,
            compression='gzip')
    return(outf)


def get_dosage_ann():
    """ Convenience function for getting the dosage and snp annotation
    """
    dos = {}
    s_ann = {}
    dos_path =\
            ("/export/home/barnarj/CCF_1000G_Aug2013_DatABEL/CCF_1000G_Aug2013_Chr"
            "{0}.dose.double.ATB.RNASeq_MEQTL.txt")
    SNP_ANNOT =\
            ("/proj/genetics/Projects/shared/Studies/Impute_CCF_Arrythmia/"
             "Projects/CCF/Projects/ATB/Projects/ATB_RNASeq/OutputData/"
             "ATB.RNASeq_Variant_Ann.bed.gz")
    return(dos, s_ann)


def get_annot_table(reference, annot):
    annoti = annot.fetch(reference=reference)
    snpid = []
    pos = []
    a0 = []
    a1 = []
    for i in annoti:
        i = i.split("\t")
        snpid.append(i[3])
        pos.append(i[1])
        a0.append(i[6])
        a1.append(i[7])
    annot_out = pd.DataFrame({'pos': pos,
        'a0': a0, 'a1': a1},
        index=pd.Index(snpid))
    return(annot_out)


def get_gene_boundries(gff_df, ensid):
    """
    col index 0 - 'aggregate gene'
    col index 1 - start
    col index 2 - end
    """
    gff_df = gff_df.ix[(gff_df.ensid == ensid) &\
             (gff_df.iloc[:,0] == 'aggregate_gene')]
    return((int(gff_df.iloc[0,1]), int(gff_df.iloc[0,2])))


def get_ratios(gaei, cissnp):
    pass


def get_symbol_ensembl(ensembl):
    """Get symbol from Ensembl REST API
    """
    server = "http://rest.ensembl.org"
    ext = "/lookup/id/{0}?expand=1"
    r = requests.get(server + ext.format(str(ensembl)),
            headers={ "Content-Type" : "application/json"})
    if not r.ok:
        return('NA')
    else:
        decoded = r.json()
        try:
            return(decoded['display_name'])
        except KeyError:
            return('NA')



def _all_eQTL_comp(chrom, base_path, debug=None, count_threshold=200):
    chrom_dict = {'chrom': str(chrom)}
    config = ConfigParser.ConfigParser()
    config.read('test.cfg')
    aei_base = base_path + config.get('data', 'aei_prefix')
    dos_path = base_path + config.get('data', 'dosage_prefix', 0, chrom_dict)
    snp_annot = base_path + config.get('annotation', 'snp_annot_prefix', 0, chrom_dict)
    gene_snps_path = base_path + config.get('annotation', 'gene_snps', 0, chrom_dict)
    eqtl_path = base_path + config.get('data', 'eqtl_prefix', 0, chrom_dict)

    ###### Loading the data ##########
    aei = get_aei((aei_base + 'chr{chrom!s}.pkl'.format(chrom=chrom)))
    dos = pd.read_csv(dos_path, sep=" ", index_col=0, header=0)
    s_ann = pd.read_pickle(snp_annot)
    gene_snps = pd.read_pickle(gene_snps_path)
    eqtl_matrix = pd.read_pickle(eqtl_path)


    ##### Restrict to only Europeans
    vsd_counts = pd.read_csv(base_path + 'eQTL/gene_vsd_eQTL_april.txt', 
            sep=",", index_col=0)
    af_euro = calculate_minor_allele_frequency(dos.ix[:,vsd_counts.columns])
    maf_euro = af_euro.index[(af_euro >= 0.01) & (af_euro <= 0.99)]
    eqtl_matrix =\
            eqtl_matrix.ix[eqtl_matrix.index.get_level_values('SNP').isin(maf_euro)]  
    print(eqtl_matrix.shape)
    count_threshold = 200
    outfile = open(base_path +\
            'eQTL/tables/global_aei/{0}_aei_rep.txt'.format(chrom), 'w+')
    eqtl_matrix = eqtl_matrix.swaplevel(0, 1)
    eqtl = eqtl_matrix.groupby(level=0)

    pvalues_fdr_calc = []
    idx = eqtl.apply(lambda x: x['p-value'].argmin())
    print(eqtl.groups.keys()[0:10])
    header = ['Symbol', 'ensid', 'Chrom', 'Indicator',
            'Best AEI SNP', '']

    for i, j in gene_snps.iteritems():
        if i == 'ENSG00000054654':
            pass
        else:
            print('going')
            continue
        symbol = get_symbol_ensembl(i)
        print(i, symbol)
        empty_out = [str(symbol), str(i), str(chrom), 'NA', 'NA', 'NA', 
                'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA']
        aei_t = aei.ix[j, :]
        aei_t = aei_t.ix[aei_t.sum(axis=1) >= count_threshold, :]
        if aei_t.shape[0] == 0:
            outfile.write("\t".join(empty_out) + "\n")
        #si = 
        # Or grab from eQTL.matrix?
        snps_cis = [eqtl_i[1] for eqtl_i in eqtl.groups[i]]
        gaei2 = AEI(aei_t, dos.ix[snps_cis, :], s_ann.ix[snps_cis,:], i)
        gaei2.calc_aei(num_threshold=20)

        bt = gaei2.pvalues
        bt = bt.ix[:, np.logical_not(bt.min().isnull())]


        cur_best = None
        cur_best_pvalue = 1

        for ind in bt.columns:
            cissnp = bt.ix[:, ind].idxmin()
            outliers_g = gaei2.outliers.ix[gaei2.hets_dict[ind], ind].values
            tgeno = dosage_round(gaei2.geno.ix[cissnp,
                gaei2.hets_dict[ind]][np.logical_not(outliers_g)])
            sum_hets_homos = np.sum((tgeno == 0) | (tgeno == 2)) 
            if (sum_hets_homos > 15) and (bt.ix[cissnp, ind] < cur_best_pvalue):
                cur_best = ind
                cur_best_pvalue = bt.ix[cissnp, ind]

        if cur_best:
            good = cur_best
            pvalue_good = cur_best_pvalue
            cissnp = gaei2.pvalues.ix[:,good].idxmin()
        else:
            outfile.write("\t".join(empty_out) + "\n")
            continue


        # :TODO get beta estimate
        pvalues_fdr_calc.extend(gaei2.pvalues.ix[:,good].values)
        #indsnp = gaei2.pvalues.columns[good]
        indsnp = good

        # Beta estimates 
        try:
            outliers_g = gaei2.outliers.ix[gaei2.hets_dict[indsnp] ,indsnp].values
        except KeyError:
            continue
        tgeno = dosage_round(gaei2.geno.ix[cissnp,
            gaei2.hets_dict[indsnp]][np.logical_not(outliers_g)])
        ar = gaei2.ratios.ix[gaei2.hets_dict[indsnp],
                indsnp][np.logical_not(outliers_g)]
        ar[ar > 1]  = 1/ar[ar > 1]
        het_combined = ar[np.array(tgeno == 1)]
        homo_combined = ar[np.array((tgeno == 0) | (tgeno == 2))]
        beta_best = np.mean(het_combined)/np.mean(homo_combined)
        try:
            aei_eqtl_best = gaei2.pvalues.ix[idx[i][1], good]
            tgeno = dosage_round(gaei2.geno.ix[idx[i][1],
                gaei2.hets_dict[indsnp]][np.logical_not(outliers_g)])
            het_combined_e = ar[np.array(tgeno == 1)]
            homo_combined_e = ar[np.array((tgeno == 0) | (tgeno == 2))]
            eqtl_best_aei_beta =\
                    np.mean(het_combined_e)/np.mean(homo_combined_e)
        except KeyError:
            aei_eqtl_best = 'NA'
            eqtl_best_aei_beta = 'NA'
        if not cissnp == idx[i][1]:
            ldbest = calculate_ld(dos.ix[[cissnp, idx[i][1]],:], 
                    cissnp)[idx[i][1]] ** 2
        else:
            ldbest = 1

        out_l = [symbol, i, 
                chrom,
                indsnp,
                cissnp,
                beta_best,
                pvalue_good,
                eqtl_matrix.ix[(i, cissnp), 'beta'],
                eqtl_matrix.ix[(i, cissnp), 'p-value'],
                idx[i][1],
                eqtl_best_aei_beta,
                aei_eqtl_best,
                eqtl_matrix.ix[(i, idx[i][1]), 'beta'],
                eqtl_matrix.ix[(i, idx[i][1]), 'p-value'],
                float(np.sum(outliers_g))/len(gaei2.hets_dict[indsnp]),
                ldbest
                ]
        if (symbol == str(debug)) or (i==str(debug)):
            embed()
            break
        else:
            pass
        out_l = [str(out_column) for out_column in out_l]
        outfile.write("\t".join(out_l) + "\n")
        '''
        if debug >= 30:
            embed()
            break
        else:
            debug += 1
        '''
    '''
    outfile_pvalues.write("\n".join([str(i) for i in pvalues_fdr_calc]))
    '''
    outfile.close()
    

if __name__ == '__main__':


    base_path = '/proj/GenomicsHD/Atrial_RNASeq/'

    _all_eQTL_comp(int(sys.argv[1]), base_path, 
            debug='SYNE2')
    #p = Process(target=_all_eQTL_comp, )
