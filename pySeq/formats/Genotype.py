import pandas as pd
import numpy as np
#from numba import autojit
import functools

class Genotype:
    def __init__(self,data):
        self.df = data
        self.ix = self.df.ix

    def dendrogram(self):
        from scipy.cluster.hierarchy import linkage, dendrogram
        g = self.geno.copy()
        X = g.as_matrix()
        Z = linkage(X,'single')
        dendrogram(Z)

    #@autojit
    def hardyweinberg(self, snp, excludeNan = True):
        from scipy.stats import chisquare
        if excludeNan:
            n=sum([0 if np.isnan(x) else 1 for x in self.geno.ix[snp,:]])
        else:
            n=self.geno.ix[snp,:].size[1]
        q = float(sum([0 if np.isnan(x) else x for x in self.geno.ix[snp,:]]))/(2*n)
        p = 1-q
        probs=[p**2,2*p*q,q**2]
        exp=np.array([probs[0]*n,probs[1]*n,probs[2]*n])
        if excludeNan:
            obs=np.array([sum([1 if x == 0 else 0 for x in self.geno.ix[snp,:]]),sum([1 if x == 1 else 0 for x in self.geno.ix[snp,:]]),\
                sum([1 if x == 2 else 0 for x in self.geno.ix[snp,:]])])
        else:
            obs=np.array([sum([1 if x == 0 or np.isnan(x) else 0 for x in self.geno.ix[snp,:]]),\
                    sum([1 if x == 1 else 0 for x in self.geno.ix[snp,:]]), sum([1 if x == 2 else 0 for x in self.geno.ix[snp,:]])])
        if chisquare(obs,exp)[1] > 0.05:
            return True
        else:
            return False

def comparable(data1, data2):
    common_indexes = data1.index.join(data2.index, how = 'inner')
    data1 = data1.ix[common_indexes, :]
    data2 = data2.ix[common_indexes, :]
    return (data1, data2)

def chi2_association(control, case, excludeNan = True):
    from scipy.stats import chi2_contingency
    probabilities = {}
    probs_in_order = []
    for snp in case.index:
        if excludeNan:
            exp = [sum([1 if x == 0 else 0 for x in control.ix[snp,:]]),\
                    sum([x if x == 1 or x == 2 else 0 for x in control.ix[snp,:]])]
        else:
            exp = [sum([1 if x == 0 or np.isnan(x) else 0 for x in control.ix[snp,:]]),\
                    sum([x if x == 1 or x == 2 else 0 for x in control.ix[snp,:]])]
        if excludeNan:
            obs = [sum([1 if x == 0 else 0 for x in case.ix[snp,:]]),\
                    sum([x if x == 1 or x == 2 else 0 for x in case.ix[snp,:]])]
        else:
            obs = [sum([1 if x == 0 or np.isnan(x) else 0 for x in case.ix[snp,:]]),\
                    sum([x if x == 1 or x == 2 else 0 for x in case.ix[snp,:]])]
        try:
            p = chi2_contingency([exp, obs])[1]
            probabilities[snp] = p
            probs_in_order.append(p)
        except:
            pass
    return (probabilities, probs_in_order)
