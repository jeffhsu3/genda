""" Population Genetic based calculations
"""

from genda import calculate_minor_allele_frequency

def calculate_D_prime(geno1, geno2):
    """ Calculates pairwise D_prime of genotype with the

    geno1 - an array or series  
    geno2 - an array or series
    """
    assert len(geno1) == len(geno2)
    p1 = calculate_minor_allele_frequency(geno1)
    p2 = calculate_minor_allele_frequency(geno2)


    p12 = (geno1 == 2) & (geno2 == 2)
    p21 = (geno1 == 0) & (geno2 == 0)
    return(p1 * p2)



