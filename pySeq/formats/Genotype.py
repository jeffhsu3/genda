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
