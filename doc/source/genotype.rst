Genotype
========
    Most of the analytic tools are built into the Genotype class, for which VCFs, SNP arrrays, and PED files are all subclasses of.

chi2_association
----------------
    Uses the chi-squared statistic to determine the likelihood that each SNP determines a certain trait. `Example <http://nbviewer.ipython.org/be0590cd0cb37cc58a96>`_.

    chi2_association(control, case, excludeNan = True)
        parameters:
            control - A pandas dataframe (SNP_array.geno for example) of samples that do not exhibit the desired trait.
        
            case - A pandas dataframe (SNP_array.geno for example) of samples that do exhibit the desired trait.

            excludeNan - Defaults to True. If True, means that Nans are not counted, otherwise they are counted as zeros.
        output:
            (SNP_dict, ordered_p_values)
            SNP_dict - A dictionary of each SNP ID to its probability.

            ordered_p_values - A list of all of the p-values in the order they were inputted. Useful for graphing.

Hierarcical Clustering
----------------------
    Creates a dendrogram showing the hierarchical clustering of the data. `Example <http://nbviewer.ipython.org/90a548316eeae6bfb476>`_.

    object.dendrogram()

Hardy-Weinberg
--------------
    Determines if a given SNP is in Hardy-Weinberg equilibrium. `Example <http://nbviewer.ipython.org/90a548316eeae6bfb476>`_.

    object.hardyweinberg(snp, excludeNan=True)
        parameters:
            snp - Name of the snp or locus to test.

            excludeNans - When True, any calls of Nan will be excluded, when False, they will be treated as zeroes.
        output:
            Boolean - True if snp is within Hardy-Weinberg, False, if it is not.

Fst
---
    Fst(subpopulations, loci, method = 'W', excludeNan = True, disable_warning = False)
        parameters:
            subpopulations - A list of genotype matricies (ex. VCF.geno, PED.geno, etc.) to be analyzed.

            loci - A list of snps to be included in the analysis.
            
            method - Which method of analysis to use. Options are 'WC' for the Weir & Cockerham method, 'W' for the Wright method, or 'R' for the Reich method (pairwise comparison). The default method is the Wright method because while testing it gave us by far the most accurate answers. For now, if you attempt to use the other two methods, a warning message will be printed to remind you of the possibility of an incorrect implimentation. There is more information about each method `here <http://www.plosone.org/article/info:doi/10.1371/journal.pone.0042649>`_.

            excludeNan - When True, any calls of Nan will be excluded, when False, they will be treated as zeroes.


            disable_warning - Defaults to False. Set to True to make the warning that is printed when using the WC or R methods disappear.
        Output:
            Float - The Fst statistic.
