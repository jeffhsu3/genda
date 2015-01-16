""" Functions for calculating multiple eQTLS
"""

import pandas as pd
import numpy as np
import statsmodels.api as sm


def ceQTL(counts, dos, cov_mat, rsid):
    """ Convinence function to calculate the single SNP comparison
    """
    acov_mat = cov_mat.copy(deep=True)
    acov_mat['soi'] = dos.ix[rsid, acov_mat.index]
    res = sm.OLS(counts, acov_mat).fit()
    return(res)

def calculate_top_cond(counts, dos, cov_mat, rsid1, 
        return_rsid1=False):
    """ Need to mask rsid1 due to perfect multicolinneraity
    """
    acov_mat = cov_mat.copy(deep=True)
    acov_mat[rsid1] = dos.ix[rsid1, acov_mat.index]

    def _calc_model(geno, counts, cov_mat):
        t = cov_mat.copy(deep=True)
        t['ts'] = geno
        return(sm.OLS(counts, t).fit().pvalues['ts'])

    def _calc_model_two(geno, counts, cov_mat, rsid1):
        t = cov_mat.copy(deep=True)
        t['ts'] = geno
        fit = sm.OLS(counts, t).fit()
        return(pd.Series({'aic': fit.aic,
            'rsid2': fit.pvalues['ts']}))

    '''
    pvalues = dos.apply(_calc_model, axis=1, args=(counts, 
        acov_mat))
    '''
    pvalues = dos.apply(_calc_model_two, axis=1, args=(counts, 
        acov_mat, rsid1))
    # Mask collinearity
    pvalues.ix[rsid1, 'rsid2'] = 1
    return(pvalues)


def calculate_two_snp_model(counts, dos, cov_mat, rsid1, rsid2):
    """ 
    """
    pass
