*****************
Working With VCFs
*****************

Genotypes
=========
    VCF.geno
        Pandas dataframe of genotype data


Hardy-Weinberg
==============
    VCF.hardyweinberg(snp, excludeNan=True)
        parameters:
            snp - name of the snp or locus to test

            excludeNans - When true any calls of Nan, will be excluded, when False, they will be treated as zeroes
        output:
            Boolean - True if snp is within Hardy-Weinberg, False, if it is not.
