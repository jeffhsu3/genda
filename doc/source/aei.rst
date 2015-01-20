.. _aei

*********************
Running Aellic Expression Imbalance
*********************

Running aellic expression imbalance analysis on RNA-Sequencing data.  Returns 
a pandas dataframe pickled to a file specified by OUTPUT.

.. code-block
    python genda_path/scripts/aei_count.py genotype_data file_to_sample_mapping.txt -o OUTPUT

file_to_sample_mapping.txt should have the following tab-delimited format.

.. code-block
sample_name    path_to_file
Sample_1    /proj/data/sample_1.bam

genotype_data is a file specificying which SNPs to count at in the 1000G VCF 
type format.

At the moment the script has a lot of limitations.  Tested on STAR and tophat
aligned bams.
