import pandas as pd
import numpy as np
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


def Fst(subpopulations, loci, method = 'WC', excludeNan = True):
    Fsts = []
    N = []
    D = []
    for snp in loci:
        pbar = 0
        total_size=0
        sub_ps=[]
        h = []
        ns=[]
        for pop in subpopulations:
            if excludeNan:
                n=sum([0 if np.isnan(x) else 1 for x in pop.ix[snp,:]])
            else:
                n=pop.ix[snp,:].size[1]
            ns.append(n)
            total_size+=n
            allele_counts = sum([0 if np.isnan(x) else x for x in pop.ix[snp,:]])
            pbar+=allele_counts
            p = float(allele_counts)/(2*n)
            sub_ps.append(p)
            if method in ['W', 'R']:
                h.append(sum([1 if x==1 else 0 for x in pop.ix[snp,:]])/float(n))
        pbar /= float(total_size)
        r = len(subpopulations)
        s2 = 0
        for p in sub_ps:
            s2 += ((p - pbar)**2) / (r - 1)
        if method == 'WC':
            Fsts.append(s2/(pbar * (1 - pbar)))
        if method == 'W':
            N.append(s2-((1/(2*float(total_size)/r-1))*(pbar*(1-pbar)-s2*(r-1)/r-sum(h)/r)))
            D.append(pbar*(1-pbar) + s2/r)
        if method == 'R':
            N.append((sub_ps[0] - sub_ps[1])**2 - (1-(sub_ps[0]**2 + (1-sub_ps[0])**2))/n[0] \
                    - (1-(sub_ps[1]**2 + (1-sub_ps[1])**2))/n[1])
            D.append(N[-1] + (1-(sub_ps[0]**2 + (1-sub_ps[0])**2)) + (1-(sub_ps[1]**2 + (1-sub_ps[1])**2)))
    if method == 'WC':
        Fst = sum(Fsts)/len(Fsts)
    if method in ['W', 'R']:
        Fst = sum(N)/sum(D)
    return Fst
