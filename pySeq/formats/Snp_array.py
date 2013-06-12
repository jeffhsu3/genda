"""
Jeffrey Hsu 2012
"""
import gzip
import Genotype
import pandas as pd
import numpy as np

class SNP_array(Genotype.Genotype):
    def __init__(self, zipfile, fileformat = 'two column', delim = ',', encoding = None, header_lines = 0,
            startatline = 0, readnrows = None):
        
        try:
            self.df = pd.read_csv(gzip.GzipFile(zipfile))
        except:
            self.df = pd.read_csv(zipfile, delimiter = delim, header = header_lines, skiprows = startatline-1, nrows = readnrows)


        if fileformat == 'two column':
            self.df.index = self.df.ix[:,'Snp.ID']
            if encoding == None:
                self.encoder = pd.Series(['A/G']*self.df.shape[0], index=self.df.index)
            else:
                self.encoder = encoding
            self.geno = self.df.ix[:,4:].apply(_two_columns_per_individual_conversion ,encoder = self.encoder, axis=1)

        if fileformat == 'one column':
            self.df.index = self.df.ix[:,0]
            if encoding == None:
                self.encoder = pd.Series(['A/G']*self.df.shape[0], index=self.df.index)
            else:
                self.encoder = encoding
            self.geno = self.df.ix[:,3:].apply(_single_column_allele, encoder = self.encoder, axis = 1)


def combine_afib_with_hapmap(afib, hapmap):
    """
    """
    geno = hapmap.join(afib, how='inner')
    print(geno.shape)
    afib_index = [0]
    afib_index.extend(range(187, afib.shape[1]))
    afib_new = geno.ix[:, afib_index].copy()
    afib_geno = convert_illumina(geno.ix[:, afib_index])
    print(len(afib_index))
    print(afib_geno.shape)
    hapmap_index = [0]
    hapmap_index.extend(range(11,185))
    hapmap_geno = convert_hapmap(geno.ix[:,hapmap_index], recode = False)
    out_geno = geno.ix[:,0:3].join([hapmap_geno, afib_geno], how='inner')
    return(out_geno)

def flat_convert(afib):
    encoding = pd.Series(['A\\G'] * afib.shape[0], index=afib.index)
    afib['enc'] = encoding
    return(afib)


def simp_mode(x):
    major_allele = mode(x)[0][0]
    return(major_allele)


def create_encoding_dict(allele1, allele2, no_ambig=True):
    comp = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
    enc_dict = {
            allele1 + allele1 : 0,
            allele1 + allele2 : 1,
            allele2 + allele1 : 1,
            allele2 + allele2 : 2,
            'NN' : np.nan,
            '--' : np.nan,
            'NA' : np.nan,
            np.nan : np.nan,
            }

    if no_ambig:
        enc_dict.update({
                comp[allele1] + comp[allele1] : 0,
                comp[allele1] + comp[allele2] : 1,
                comp[allele2] + comp[allele1] : 1,
                comp[allele2] + comp[allele2] : 2,
                })


    return enc_dict

def _single_column_allele(genotype_array, encoder):
    """ Converts dataframes containing values the format described below into a integer dataframe. The first
    element in genotype_array should be the encoding.

    eg:  Ind01    Inde02
          AA       AG
    """
    e = encoder[genotype_array.name]
    encoding_dict = {e[0]+e[0] : 0, e[-1]+e[0] : 1, e[0]+e[-1] : 1, e[-1]+e[-1] : 2}
    new_array = genotype_array[:].map(encoding_dict)

    return new_array

def _two_columns_per_individual_conversion(genotype_array, encoder):
    """ Convert a dataframe containing a column with two alleles

    eg:  Ind1.A1    Ind1.A2    Ind2.A1 ...
            A          G          A
    """
    #genotype_array = pd.Series(['A/G'].extend(list(genotype_array)))
    #encoding_dict = {genotype_array[0][0] : 0, genotype_array[0][2] : 1}
    temp = encoder[genotype_array.name]
    encoding_dict = {temp[0] : 0, temp[2] : 1}
    try:
        recoded_1 = genotype_array[0::2].map(encoding_dict)
        recoded_2 = genotype_array[1::2].map(encoding_dict)
        new_index = recoded_1.index.map(lambda x: x[:-3])
        recoded_1.index = new_index
        recoded_2.index = new_index
        recoded_array = recoded_1 + recoded_2
    except IndexError:
        print('Welp')

    return recoded_array


def convert_hapmap(input_dataframe, recode=False, index_col=0):
    """ Specifically deals with hapmap and 23anMe Output
    """
    complement = {'G/T': 'C/A', 'C/T': 'G/A', "G/A" : "G/A", "C/A": "C/A", "A/G" : "A/G",
            "A/C": "A/C"}
    dataframe = input_dataframe.copy()
    if recode:
        recode = dataframe.ix[:, index_col].apply(lambda x: complement[x])
        dataframe.ix[:,0] = recode
    new_dataframe = dataframe.apply(_single_column_allele, axis=1)

    return new_dataframe


def convert_illumina(dataframe, index_col = 0):
    """ Since the Illumina platform contains no ambigous snps ('A/T' or 'G/C') only
    'G/T' and 'G/C' need to be converted.

    Returns a datarame as 0,1 and 2 encoding for the genotype.

    This is about thirty times faster than the non-vectorized solution.
    """
    print(dataframe.shape)
    complement = {'G/T': 'C/A', 'C/T': 'G/A', "G/A" : "G/A", "C/A": "C/A", "A/G" : "A/G",
            "A/C": "A/C"}
    recode = dataframe.ix[:, index_col].apply(lambda x: complement[x])
    dataframe.ix[:,0] = recode
    new_dataframe = dataframe.apply(_two_columns_per_individual_conversion, axis = 1)
    print(new_dataframe.shape)

    return new_dataframe


def convert_plink(x):
    """
    """
    pass


def convert_ped(x):
    """
    """
    pass


