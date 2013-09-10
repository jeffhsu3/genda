import pandas as pd
import numpy as np
import functools
from pickle import dump, loads

class Genotype(object):
    __doc__ = """ 
    Base class for genotypes. Not meant to be used directly.

    %{param)s
    data : genotype dataframe

    %(extra_param)s
    """
    def __init__(self, data):
        self.df = data
        self.ix = self.df.ix

    def dendrogram(self):
        """ Plots dendrogram
        """
        from scipy.cluster.hierarchy import linkage, dendrogram
        g = self.geno.copy()
        X = g.as_matrix()
        Z = linkage(X,'single')
        dendrogram(Z)

    def hardyweinberg(self, snp, excludeNan = True):
        """ Calculate hardyweinberg

        Parameters
        ----------
        snp : identifier to run the test at 
        excludeNan : exlude nas

        Returns
        -------
        boolean
        """
        from scipy.stats import chisquare
        if excludeNan:
            n=sum([0 if np.isnan(x) else 1 for x in self.geno.ix[snp,:]])
        else:
            n=self.geno.ix[snp,:].size[1]
        q = float(sum([0 if np.isnan(x) else x for x in self.geno.ix[snp,:]]))/(2*n)
        p = 1-q
        probs=[p**2, 2*p*q, q**2]
        exp=np.array([probs[0]*n, probs[1]*n, probs[2]*n])
        if excludeNan:
            obs=np.array([sum([1 if x == 0 else 0 for x in self.geno.ix[snp,:]]), 
                sum([1 if x == 1 else 0 for x in self.geno.ix[snp,:]]),\
                sum([1 if x == 2 else 0 for x in self.geno.ix[snp,:]])])
        else:
            obs=np.array([sum([1 if x == 0 or np.isnan(x) else\
                    0 for x in self.geno.ix[snp,:]]),\
                    sum([1 if x == 1 else 0 for x in self.geno.ix[snp,:]]), 
                    sum([1 if x == 2 else 0 for x in self.geno.ix[snp,:]])])
        if chisquare(obs,exp)[1] > 0.05:
            return True
        else:
            return False


    def geno_save(self, f):
        """ Save genotype object

        Paremeters
        ----------
        f : file location
        """
        # Hmm how to remove the raw df temporarily before saving?
        del self.df
        with open(f, 'w') as temp:
            dump(self, temp)
            temp.close()



def load_geno(self, f):
    """ Load genotype object from file

    Parameters
    ----------
    f : 

    Returns
    -------
    Genotype object

    """
    return loads(open(f, 'r').read())




def comparable(data1, data2):
    """
    """
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


def Fst(subpopulations, loci, method = 'W', excludeNan = True, disable_warning = False):
    """
    Fixation index (Fst calcualtions)


    Parameters
    ----------
    subpopulations :  
    loci : 
    method : {'W', 'R', 'WC'}

    Returns
    -------


    """
    if method != 'W' and disable_warning == False:
        print('In our limited testing, the method you selected gave \
                questionable answers. We are not exactly sure right now \
                if that is simply randomness and differences in the methods \
                or an incorrect implamentation. Unless you need to use the \
                method you selected, it is suggested you use the W method,  \
                which is selected by default.')
    method = method.upper()
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
            if n > 0:
                total_size+=n    
                ns.append(n)
                allele_counts = sum([0 if np.isnan(x) else x for x in pop.ix[snp,:]])
                pbar+=allele_counts
                p = float(allele_counts)/(2*n)
                sub_ps.append(p)
                if method in ['W', 'R']:
                    h.append(sum([1 if x==1 else 0 for x in pop.ix[snp,:]])/float(n))
        if total_size > 0:
            pbar /= float(2*total_size)
            r = len(subpopulations)
            s2 = 0
            for p in sub_ps:
                s2 += ((p - pbar)**2) / (r - 1)
            if method == 'WC':
                if pbar == 1 or pbar == 0:
                    Fsts.append(0)
                else:
                    Fsts.append(s2/(pbar * (1 - pbar)))
            if method == 'W':
                N.append(s2-((1/(2*float(total_size)/r-1))*(pbar*(1-pbar)-s2*(r-1)/r-sum(h)/r)))
                D.append(pbar*(1-pbar) + s2/r)
            if method == 'R':
                N.append((sub_ps[0] - sub_ps[1])**2 - (1-(sub_ps[0]**2  +\
                        (1-sub_ps[0])**2))/ns[0] - (1-(sub_ps[1]**2 + \
                        (1-sub_ps[1])**2))/ns[1])
                D.append(N[-1] + (1-(sub_ps[0]**2 + (1-sub_ps[0])**2)) + \
                        (1-(sub_ps[1]**2 + (1-sub_ps[1])**2)))
    if method == 'WC':
        Fst = sum(Fsts)/len(Fsts)
    if method in ['W', 'R']:
        Fst = sum(N)/sum(D)
    return Fst
