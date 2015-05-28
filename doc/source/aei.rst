.. _aei

*********************
Running Aellic Expression Imbalance
*********************

Running aellic expression imbalance analysis on RNA-Sequencing data.  Counts the 
alleles in the specified bam files from 'file_to_sample_mapping.txt' at
locations specified in the genotype file.  The script returns 
a pandas dataframe pickled to a file specified by OUTPUT.  The columns of
the aei dataframe is a multi-index with the first index being the sample
id and the second index composing of base pairs in the order as specificed:
[A, G, C, T].

:TODO Need to handle indels tag SNPs.

.. code-block
    python genda_path/scripts/allele_counter.py [genotype_data] file_to_sample_mapping.txt -o OUTPUT


file_to_sample_mapping.txt should have the following tab-delimited format.

'''
sample_name    path_to_file
Sample_1    /proj/data/sample_1.bam
'''

genotype_data is a file specificying which SNPs to count at in the 1000G VCF 
type format.

At the moment the script has a lot of limitations.  Tested on STAR and tophat
aligned bams.
