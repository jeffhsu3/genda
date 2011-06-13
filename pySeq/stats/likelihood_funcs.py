"""
An implementation of the model of allelic imbalance found in:
    Statistical Inference of Allelic Imbalance from Transcriptome Data.
    Nothangel, Michael et. Al. Human Mutation Volume 32, Issue 1 December 2010
    
    The empirical error numbers were taken from the same publication.  

    Need to asses bias of aligning.  Make assumption that most loci show
    no-allelic imbalance.  Fit slope to curve and adjust.  

Written by Jeff Hsu 2011

"""

from scipy import optimize
from scipy.misc import comb
from scipy import stats
import numpy as np
from itertools import izip
from math import log, exp
# See cited publication to how these values were obtained from 5 HapMap samples
transMatrix = np.array([[0,0.40210, 0.48006, 0.11784],[0.22213, 0, 0.25456,
    0.52330], [0.53012, 0.26919, 0, 0.20070],[0.13354,0.44217,0.42429,0]])

argmax = lambda array: max(izip(array, xrange(len(array))))[1]

def logLik(f,counts, pi=2.17e-3, piTrans=transMatrix):
    """
    f = real allelic imbalance
    counts = [A counts, C counts, G counts, T counts] 
    pi = error rate
    piTrans = conditional transistion error matrix, rows are the actual nuclear
    base and columns are the miscall (ACGT ordering)
    """
    ind = counts.argsort()[::-1]  
    # Formula 1 pi: is the mis the error rate
    
    loglik = counts[ind[0]]*log(f*(1-pi)+(1-f)*pi*piTrans[ind[1],ind[0]])+\
	    counts[ind[1]]*log((1-f)*(1-pi)+f*pi*piTrans[ind[0],ind[1]])+\
	    counts[ind[2]]*log(pi*(f*piTrans[ind[0],ind[2]]+(1-f)*piTrans[ind[1],ind[2]]))+\
	    counts[ind[3]]*log(pi*(f*piTrans[ind[0],ind[3]]+(1-f)*piTrans[ind[1],ind[3]]))
    return(loglik)


def lik(f, counts, pi=2.17e-3, piTrans=transMatrix):
    ind = counts.argsort()[::-1]
    lik = ((f*(1-pi)+(1-f)*pi*piTrans[ind[1],ind[0]])**counts[ind[0]])*\
	    (((1-f)*(1-pi)+f*pi*piTrans[ind[0],ind[1]])**counts[ind[1]])*\
	    ((pi*(f*piTrans[ind[0],ind[2]]+(1-f)*piTrans[ind[1],ind[2]]))**counts[ind[2]])*\
	    ((pi*(f*piTrans[ind[0],ind[3]]+(1-f)*piTrans[ind[1],ind[3]]))**counts[ind[3]])
    return(lik)

def logLikHomo(counts, pi=2.17e-3, piTrans=transMatrix):
    """
    The log Likelihood that the data is homozygous
    """
    ind = counts.argsort()[::-1]  
    loglik = counts[ind[0]]*log(1-pi)+counts[ind[1]]*\
	    log(pi*piTrans[ind[0],ind[1]])+counts[ind[2]]*log(pi*piTrans[ind[0],ind[2]])+\
	    counts[ind[3]]*log(pi*piTrans[ind[0], ind[3]])
    return(loglik)


def likHomo(counts, pi=2.17e-3, piTrans=transMatrix):
    ind = counts.argsort()[::-1]
    lik = ((1-pi)**counts[ind[0]])*\
	    ((pi*piTrans[ind[0],ind[1]])**counts[ind[1]])*\
	    ((pi*piTrans[ind[0],ind[2]])**counts[ind[2]])*\
	    ((pi*piTrans[ind[0],ind[3]])**counts[ind[3]])
    return(lik)


def ratioLik(counts):
    f_max = optimize.fminbound(lambda x: -logLik(x, counts), 0.5,1, xtol=1e-8,
	    disp=0)
    L_null=logLik(0.5, counts)
    L_max = logLik(f_max, counts)    
    test_stat = -2*(L_null - L_max)
    return(stats.chi2.sf(test_stat, 1)) 



def isHet(counts, prior=1):
    """ 
	Infers if the genotype is a Het or not at that particular base position
	not accurate at really imbalanced alleles and high coverage
    """
    # If statement is necessary for really high read coverage
    if sum(counts) >= 150:
	if float(max(counts))/sum(counts)>=0.95: 
	    return False
	else:
	    return True
    else:
	logLik_het = logLik(0.5, counts)
	p_D =likHomo(counts)+lik(0.5,counts)
	p_D = log(p_D)
	post_prob = exp(logLik_het+log(prior)-p_D)
	if post_prob >= 0.5: return True
	else: return False


def err_handler(type, flag):
    pass
    

    
