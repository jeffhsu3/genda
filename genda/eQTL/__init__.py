import numpy as np
import pandas as pd

from .plotting import *
from .multiple_eQTL import *


def var_boxplot(ax, x, colors=None):
    """
    Arguments:
    ----------
    ax - axis object
    """
    for i in x:
        pass


def get_meQTL(filename, gene):
    """ Get a subset of the matrix-eQTL object
    """
    meqtl_iter = pd.read_csv(filename, sep="\t", chunksize=30000,
            index_col=(1,0))

    meqtl = pd.concat([
        chunk.ix[chunk.index.get_level_values('gene') == gene,:]\
                for chunk in meqtl_iter])
    return(meqtl)


def subset_meQTL(meQTL, gene_name):
    """
    """
    try:
        index = meQTL.gene == gene_name
        subset = meQTL.ix[index, :]
    except AttributeError:
        # Case where eQTL matrix is already 1 gene
        subset = meQTL
    return(subset)
