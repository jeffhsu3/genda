class FastA(object):
    """
    Flexible wrapper for a FastA file and its index.  
    """

    def __init__(self, fastAfile):
	self.fh = open(fastAfile, 'rU')
	self.index = self.__indexFA(fastAfile)

    def __indexFA(self, fastAfile):
	"""
	Try and incorporate multiple types of indexes including various binary
	and tree based ones.  
	"""
	index = {}
	with open(fastAfile+".fai", 'rU') as indexfile:
	    for line in indexfile:
		line = line.strip().split('\t')
		index[line[0]]=int(line[2])
	return(index)

    def grabSequence(self, region, start, end):
	"""
	Grabs the sequence in region:start-end from self.fh.  Need to adjust
	for multiple line lengths.  Also need to alter so it doesn't seek to
	the beginning of the file each time.    
	"""
	try:
	    current = self.fh.tell()
	    length = end - start + 1
	    newlines = start/50  
	    # Calculates how many newlines are in the query
	    mod = start%50
	    if mod == 0:
		q_nl = length/50 + 1  
		start = start - 2
		length = length + q_nl
	    else:
		q_nl = length/50  
		start = start - 1
		length = length + q_nl
	    self.fh.seek(self.index[region]+start+newlines)
	    seq = (self.fh.read(length))
	    seq = seq.replace("\n", "")
	    self.fh.seek(0)
	except IndexError:
	    print("Region '%s' was not in the index file" % (region)) 
	return(seq)
	 

def main():
    import sys
    fa = FastA(sys.argv[1])
    what = fa.grabSequence('chr10', 1500003, 1500017) 
    assert(what=="cctgacaatgtcaat")
    what = fa.grabSequence('chr7', 116201200, 116201230) 
    assert(what.lower()=='caagcctgcacaataaaaatgtttaacggtt')

if __name__ == '__main__':
    main()
