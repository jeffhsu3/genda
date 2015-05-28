"""
Jeffrey Hsu 2012
"""
import gzip
import Genotype
import pandas as pd
import numpy as np


class SNP_array(Genotype.Genotype):
    __doc__ = """
    Class that specifically deals with genotyping from microarrays.

    %(param)s
    zipfile : location of output

    %(extra_param)s
    fileformat : {'one column', 'two column', 'illumina'} 
    delim : Data delimiter
    sample_col : The columns in which 
    encoding : Encoding data for the SNPs.
    header_lines : Number of lines to be read in as the header.
    startatline : The line to start reading data at.

    Attributes
    ----------
    df : original dataframe loaded
    geno : converted genotype matrix
    file_format : file format
    apply_encoder :  


    Related
    -------
    Genotype.Genotype
    """
    def __init__(self, zipfile, file_format = 'one column', 
            delim = ',', samp_col = None, encoding = None, header_lines = 0,
            startatline = 0, readnrows = None):
        
        try:
            self.df = pd.read_csv(gzip.GzipFile(zipfile), 
                    delimiter = delim, header = header_lines,
                    skiprows = startatline-1, nrows = readnrows)
        # :TODO change exception
        except:
            self.df = pd.read_csv(zipfile, delimiter = delim, header = header_lines, skiprows = startatline-1, nrows = readnrows)


        if file_format == 'two_column':
            if samp_col == None:
                samp_col = 4
            self.samp_col = samp_col
            try:
                self.df.index = self.df.ix[:,'Snp.ID']
            except KeyError:
                self.df.index = self.df.ix[:,1]
            if isinstance(encoding, (int, long)):
                encoding = self.df.ix[:,encoding]
            self.encoder = encoding
            if type(self.encoder) == type(None):
                self.geno = None
            else:
                self.geno = self.df.ix[:,self.samp_col:].apply(
                        _two_columns_per_individual_conversion, 
                        encoder = self.encoder, axis=1)

        elif file_format == 'one column':
            if samp_col == None:
                samp_col = 3
            self.samp_col = samp_col
            self.df.index = self.df.ix[:,0]
            if isinstance(encoding, (int, long)):
                encoding = self.df.ix[:,encoding]
            self.encoder = encoding
            if type(self.encoder) == type(None):
                self.geno = None
            else:
                self.geno = self.df.ix[:,self.samp_col:].apply(
                        _single_column_allele, 
                        encoder = self.encoder, axis = 1)

        elif file_format == 'illumina':
            pass

        else:
            raise NotImplementedError("Only 'one_column', 'two_column' and\
                    'illumina' formats are currently supported")

        self.filetype = file_format


    def apply_encoder(self, encoder):
        """
        To apply the encoder to the SNP_array object after the fact if one was sot 
        supplied with creation
        """
        if self.filetype == 'one column':
            return self.df.ix[:,self.samp_col:].apply(_single_column_allele, 
                    encoder = encoder, axis = 1)
        elif self.filetype == 'two column':
            return self.df.ix[:,self.samp_col:].apply(_two_columns_per_individual_coversion, 
                    encoder = encoder, axis = 1)
        else:
            return None


def combine_afib_with_hapmap(geno , hapmap):
    """ Convenience function for combining
    """
    geno = hapmap.join(geno, how='inner')
    afib_index = [0]
    afib_index.extend(range(187, geno.shape[1]))
    geno = geno.ix[:, afib_index].copy()
    afib_geno = convert_illumina(geno.ix[:, afib_index])
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
    """ Converts dataframes containing values the format described 
    below into a integer dataframe. The first element in 
    genotype_array should be the encoding.

    eg:  Ind01    Inde02
          AA       AG
    """
    e = encoder[genotype_array.name]
    encoding_dict = {e[0]+e[0] : 0, e[-1]+e[0] : 1, e[0]+e[-1] : 1, e[-1]+e[-1] : 2}
    new_array = genotype_array[:].map(encoding_dict)

    return new_array


def _two_columns_per_individual_conversion(genotype_array, encoder):
    """ Convert a dataframe containing a column with two alleles

    Parameters
    ----------
    genotype_array : 
    encoder : 

    eg:  Ind1.A1    Ind1.A2    Ind2.A1 ...
            A          G          A
    """

    try:
        temp = encoder[genotype_array.name]
    except AttributeError:
        temp = encoder
    except TypeError:
        temp = encoder
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


def convert_illumina(dataframe, encoding):
    """ Since the Illumina platform contains no ambigous snps ('A/T' or 'G/C') only
    'G/T' and 'G/C' need to be converted.

    Returns a datarame as 0,1 and 2 encoding for the genotype.

    This is about thirty times faster than the non-vectorized solution.
    """
    complement = {'G/T': 'C/A', 'C/T': 'G/A', "G/A" : "G/A", "C/A": "C/A", "A/G" : "A/G",
            "A/C": "A/C"}
    recode = encoding.apply(lambda x: complement[x])
    dataframe.ix[:,0] = recode
    new_dataframe = dataframe.apply(_two_columns_per_individual_conversion,
            args=(encoding,), axis = 1)

    return new_dataframe


def convert_plink(x):
    """
    """
    pass


def convert_ped(x):
    """
    """
    pass


