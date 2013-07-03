"""
Class to count genes that have multiple exons using pysam
"""

class geneCounter():
    """
	Callback class
    """

    def __init__(self, gene):
        self.gene = gene
        # span_slice contains a list of the reads that span splice junctions
        self.splicereads = []
        self.mCounts = 0
        self.transcripts = []
        #self.transcript = transcript

    def __call__(self, alignment):
        if alignment.is_proper_pair and alignment.is_read1 and alignment.is_proper_pair:
            cigarTypes = [ x[0] for x in alignment.cigar]
            if 3 in cigarTypes:
            # 3 corresponds to skipping reference ie it spans a splice
            # junction
                if alignment.qname in self.splicereads:
                    # Don't double count spliced reads
                    pass
                else:
                    self.splicereads.append(alignment.qname)
                    self.mCounts +=  1
                    # self.allReads.append(alignment.qname)
            else:
                self.mCounts += 1
