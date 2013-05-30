***********************
Working with SNP Arrays
***********************

The SNP Array class
===================
    SNP_array(zipfile, fileformat = 'two column', delim = ',', encoding = None):
    
        parameters:
            zipfile - File with data in it. Will try gunzip, but if this fails, it will read it as plain text.
            
            fileformat - Either 'one column' or 'two column' depending on wheter each individual's alleles are shown in one or two columsn. Defaults to 'two column'.
            
            delim - What seperates the values in your data? Defaults to ',' but other common options are tabs ('\t') or spaces (' ').
            
            encoding - A pandas series which contains the reference/alternate alleles (eg. 'A/G')  for each position with ids that correspond to the file.

Data
----
    SNP_array.df
        Data that was passed in parsed as a pandas dataframe.

    SNP_array.encoder
        The encoder passed in on creation of the object.

    SNP_array.geno
        A pandas dataframe of the genotype data. 0 represents homozygous for the reference allele, 1 is heterozygous, and 2 is homozygous for the alternate allele.

Hardy-Weinberg
--------------
    SNP_array.hardyweinberg(snp, excludeNan=True)
        parameters:
            snp - name of the snp or locus to test

            excludeNans - When true any calls of Nan, will be excluded, when False, they will be treated as zeroes
        output:
            Boolean - True if snp is within Hardy-Weinberg, False, if it is not
