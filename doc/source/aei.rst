.. _aei

*********************
Running Aellic Expression Imbalance
*********************

Running aellic expression imbalance analysis on RNA-Sequencing data. 

python AEI_2.py genotype_data file_to_sample_mapping.txt -o OUTPUT

file_to_sample_mapping.txt should have the following tab-delimited format.

sample_name    path_to_file
Sample_1    /proj/data/sample_1.bam

At the moment the script has a lot of limitations.  Tested on STAR and tophat
aligned data.
