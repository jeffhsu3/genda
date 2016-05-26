import sys
import numpy as np
from scipy.stats import pearsonr


def calculate_minor_allele_frequency(genotypes):
    """
    Returns a dataframe of the allele frequencies

    Parameters
    ----------
    genotypes : a genotype dataframe with polymorphisms going row-wise
    """
    try:
        maf = genotypes.sum(axis=1) / (2 * genotypes.shape[1])
    except IndexError:
        maf = genotypes.sum() / (2 * len(genotypes))
    return(maf)


def calculate_ld(genotypes, snp):
    """ Calculates SNP linkage disquilibrium

    Parameters
    ----------
    genotypes : a genotype dataframe with polymorphisms going rowise and samples
    columnwise.
    snp: snp-identifier

    Returns
    -------
    A vector of linkage disequilibrium

    """
    ld = genotypes.corrwith(genotypes.ix[snp, :], axis=1)
    ld[np.isnan(ld.values)] = 0
    ld = np.sqrt(ld ** 2)
    return(ld)


def regionParse(string):
    """
    Parses region:start_position-end_position, returns a tuple
    >>> regionParse('chr4:200000-200100')
    ('chr4', 200000, 200100)
    """
    start = string[string.find(':') + 1:string.find('-')]
    end = string[string.find('-') + 1:]
    start = start.replace(',', '')
    end = end.replace(',', '')

    try:
        start = int(start)
        end = int(end)
    except ValueError:
        sys.exit()
    return(string[0:string.find(':')], start, end)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
