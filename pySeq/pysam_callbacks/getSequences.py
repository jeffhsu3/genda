"""
A callback class that gets the sequences that overlap a genomic interval
"""

class seqRetriever(object):
    """
	Callback class
    """
    
    def __init__(self, gene):
	self.gene = gene
	# span_slice contains a list of the reads that span splice junctions
	self.splicereads = []
	self.mCounts = 0

	#self.transcript = transcript
    def __call__(self, alignment):
	if alignment.is_proper_pair:
	    # Right now Tophat requires perfect alignment of the splice
	    # junction, therefore, the 2nd CIGAR entry is always the splice
	    # junction.  

	    # Compare the splice reads for differences.  Need to make sure the
	    # starts are the same.  In addition this means I have to check the
	    # other reads as well
	    cigarTypes = [ x[0] for x in alignment.cigar]
	    if 3 in cigarTypes:
		# 3 corresponds to skipping reference ie it spans a splice
		# junction
		if alignment.qname in self.splicereads:
		    pass
		else: 
		    self.splicereads.append(alignment.qname)
		    print(">" + self.gene)
		    print(alignment.seq) 
		    self.mCounts +=  1 
		    # self.allReads.append(alignment.qname)
	    else:
		print(">" + self.gene)
		print(alignment.seq)
		self.mCounts += 1
