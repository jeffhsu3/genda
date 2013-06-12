Genotype
========

chi2_association
----------------
    Uses the chi-squared statistic to determine the likelihood that each SNP determines a certain trait.
    
    chi2_association(control, case, excludeNan = True)
        parameters:
            control - A pandas dataframe (SNP_array.geno for example) of samples that do not exhibit the desired trait.
        
            case - A pandas dataframe (SNP_array.geno for example) of samples that do exhibit the desired trait.

            excludeNan - Defaults to True. If True, means that Nans are not counted, otherwise they are counted as zeros.
        output:
            (SNP_dict, ordered_p_values)
            SNP_dict - A dictionary of each SNP ID to its probability.

            ordered_p_values - A list of all of the p-values in the order they were inputted. Useful for graphing.

Hardy-Weinberg
--------------
    VCF.hardyweinberg(snp, excludeNan=True)
        parameters:
            snp - Name of the snp or locus to test.

            excludeNans - When true any calls of Nan, will be excluded, when False, they will be treated as zeroes.
        output:
            Boolean - True if snp is within Hardy-Weinberg, False, if it is not.
