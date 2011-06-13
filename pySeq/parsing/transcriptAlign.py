"""
Series of functions that deal with collapsing bowtie BAM alignment to a
transcriptome where a read can align to many transcripts.  Bowtie should be
run as follows:
    bowtie -a --best --strata -S -X 400 -p 8 transcriptome.ref -1 1.fastq -2
    2.fastq > out.sam

There is a script in scripts that runs the whole pipeline.
"""

class readSet(object):

    def __init__(self):
