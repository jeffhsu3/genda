***********************
Working with SNP Arrays
***********************

The SNP Array class
===================
    SNP_array(zipfile, fileformat = 'one column', delim = ',', samp_col = None, encoding = None, heade_lines = 0, startatline = 0, readnrows = None):
    
        parameters:
            zipfile - File with data in it. Will try gunzip, but if this fails, it will read it as plain text.
            
            file_format - Either 'one column' or 'two column' depending on wheter each individual's alleles are shown in one or two columsn. Defaults to 'one column'.
            
            delim - What seperates the values in your data? Defaults to ',' but other common options are tabs ('\t') or spaces (' ').
            
            samp_col - What is the first column where the data for each sample starts. If no value is supplied, it defaults to 3 for one column data and 4 for two column data.

            encoding - This arguement will take either an integer or a pandas series. If an integer is supplied, the encoder will be taken from the corresponding column (0 based of course). Alternatively a pandas series which contains the reference/alternate alleles (eg. 'A/G')  for each position with ids that correspond to the file can be supplied. If no encoder is supplied, no genotype matrix will be created.

            header_lines - The number of lines to be read as the reader.

            startatline - What line to start reading the data at.

            readnrows - How many lines of data to read. When None, the file is read until the end.

Data
----
    SNP_array.df
        Data that was passed in parsed as a pandas dataframe.

    SNP_array.encoder
        The encoder passed in on creation of the object.

    SNP_array.geno
        A pandas dataframe of the genotype data. 0 represents homozygous for the reference allele, 1 is heterozygous, and 2 is homozygous for the alternate allele.

    SNP_array.apply_encoder(encoder)
        If an encoder was not supplied upon creating the object, this is how you could still get a genotype matrix after the fact.    

        parameters:
            encoder - An integer that corresponds to the appropriate column or a pandas series which contains the reference/alternate alleles (eg. 'A/G')  for each position with ids that correspond to the file.

        output:
            A pandas dataframe of genotype data. 0 represents homozygous for the reference allele, 1 is heterozygous, and 2 is homozygous f    or the alternate allele.

Hardy-Weinberg
--------------
    `SNP_array.hardyweinberg(snp, excludeNan=True) <https://pyseq.readthedocs.org/en/latest/genotype.html#hardy-weinberg>`_
