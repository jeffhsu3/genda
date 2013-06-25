**********************
Working with PED files
**********************

The PED class
===================
    PED(PED, MAP, encoder = None):
    
        parameters:
            PED - The ped file containing the genotype data.
            
            MAP - The map file which contains labels for the snps in the ped file.
            
            encoder - A pandas series which contains the reference/alternate alleles (eg. 'A/G')  for each position with ids that correspond to the file. Defaults to None. Must provide one in order to get a genotype matrix.

Data
----
    PED.PED
        Data that was passed in parsed as a pandas dataframe.

    PED.MAP
        The map file parsed as a pandas dataframe.

    PED.geno
        A pandas dataframe of the genotype data. 0 represents homozygous for the reference allele, 1 is heterozygous, and 2 is homozygous for the alternate allele.

Hardy-Weinberg
--------------
    PED.hardyweinberg(snp, excludeNan=True)
        parameters:
            snp - name of the snp or locus to test

            excludeNans - When true any calls of Nan, will be excluded, when False, they will be treated as zeroes
        output:
            Boolean - True if snp is within Hardy-Weinberg, False, if it is not
