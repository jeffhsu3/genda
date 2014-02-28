""" Functions for working with tabix dosages in pandas dataframes
"""
import gzip
import numpy as np
import pandas as pd
import pysam
import pickle as pkl
import statsmodels.api as sm

class Dosage(object):

    def __init__(self, dosages, annotations, gene_name):
        # Match up the annotation dataframe with the dosage dataframe
        mindex = np.intersect1d(np.asarray(dosages.index, dtype=str), 
                np.asarray(annotations.index, dtype=str))
        self.annot = annotations.loc[mindex, :]
        ordering = self.annot.ix[:, 'pos'].argsort()
        self.annot = self.annot.iloc[ordering, :]
        self.dosages = dosages.ix[mindex, :]
        self.dosages = self.dosages.iloc[ordering, :]
        self.gene_name = gene_name

    def run_eQTL(self, count_matrix, covariates, extra_snps=None):
        #self.pvalues = self.dosages.apply()
        pvalues = self.dosages.apply(eQTL_func, axis=1, args=(covariates,
            count_matrix.ix[self.gene_name, :]))
        self.pvalues = pvalues



def get_dosages_by_range(chrm, start, end, gene_name, annotation_file, 
        dosage_df, mapping=None):
    """
    Fuzzy mapping between annotation and genotypes
    Returns Dosage instance.
    """
    ann_file = pysam.Tabixfile(annotation_file)
    ann_v = ann_file.fetch(chrm, start, end)
    rsIDs = []
    pos = []
    ref = []
    alt = []
    for i in ann_v:
        i = i.split("\t")
        rsIDs.append(i[3])
        pos.append(int(i[1]))
        ref.append(i[6])
        alt.append(i[7])
    annot = pd.DataFrame({'pos': pos, 'ref': ref, 'alt': alt}, index=pd.Index(rsIDs))
    comb_iter = []
    for dos in dosage_df:
        mindex = np.intersect1d(np.asarray(dos.index, dtype=str),
                np.asarray(annot.index, dtype=str))
        if len(mindex) > 0:
            comb_iter.append(dos.ix[mindex, :])
        else:
            pass
    out_dos = pd.concat(comb_iter)
    '''
    dosages = pd.read_csv(dosage_path + path, sep=" ", header=None,
            index_col = 0, skiprows=roughly_first, 
            nrows=roughly_end-roughly_first, names=col_names.columns)
    '''
    print(annot.shape, out_dos.shape, gene_name)
    return Dosage(out_dos, annot, gene_name)


def grab_gene_location(hgnc, cis_buffer=0, ensid = False):
    """
    """
    ann_file = "/proj/genetics/Projects/shared/Subject_Sources/" +\
            "External/Ensembl/OutputData/EnsemblAnnotationAllHumanGenes.bed.gz"
    with gzip.open(ann_file) as annot:
        for line in annot:
            if hgnc in line:
                line = line.split("\t")
                return(line[0], int(line[1]) - cis_buffer, 
                        int(line[2]) + cis_buffer, str(line[4]))


def generate_dosage_mapping(dosage_file, mapping_file = None, interval=50):
    """
    Returns dictionary of rsIDs: fileposition from a dosage file
    """
    if not mapping_file:
        with open(dosage_file) as fh:
            fh.next()
            t = 0
            debug = 0
            f_i = {}
            for i, j in enumerate(fh):
                if i % 50 == 0:
                    f_i[j.split(" ")[0]] = i - 1
                else: pass
    return(f_i)


def eQTL_func(snps, cov, expression):
    """
    """
    cov = cov.T
    cov['snps'] = snps
    cov = sm.add_constant(cov)
    model = sm.OLS(expression, cov)
    return(model.fit().pvalues['snps'])


class eQTL(object):
    """ Python class for completing eQTLs.  Does lazy loading of all large
    files. 
    """
    def __init__(self, dosages_path, expression, vannotation):
        self.dosage = dosages_path
        self.expression = expression
        self.vannotations = vannotations


    def generate_mapping():
        pass

"""
    if mapping:
        for i in ann_v:
            rsID = i.split("\t")[3]
            try:
                roughly_first = mapping[rsID]
                rsIDs.append(rsID)
                pos.append(int(i.split("\t")[1]))
                break
            except KeyError:
                pass
        for i in ann_v:
            i = i.split("\t")
            try: 
                roughly_end = mapping[i[3]]
            except KeyError:
                pass
            pos.append(int(i[1]))
            rsIDs.append(i[3])
"""
