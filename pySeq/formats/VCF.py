"""Written by Jeff Hsu 2010
"""

import collections, os
class VCFfile(object):
    """ A flexible container for VCF files.  
    
    http://www.1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcf4.0
    for VCF 4.0 specification.

    :TODO better parse the metadata information

    """
 
    def __init__(self, filename):
	""" Removes the header
	"""
	if isinstance(filename, str):
	    try:
		self.handle = open(filename, "rU")
	    except IOError as (errno, strerror):
		print "I/O error({0}): {1}".format(errno,strerror)
	elif isinstance(filename,file):
	    self.handle = filename
	else: 
	    print("VCFfile needs a filename or a file")
	    sys.exit()
	# Header processing.  Makes it easy since the metadata starts with
	# '##' and the header starts with '#'
	self.metadata = []
	line = self.handle.readline().strip('\n')
	self.metadata.append(line[2:])
	while line[0:2]=="##":
	    self.metadata.append(line[2:])
	    line = self.handle.readline().strip('\n')
	# The header
	self.header = line[1:].strip('\n').split('\t')
	
    
    def __iter__(self):
	return self
    
    def header(self):
	"""
	Returns the VCF header
	"""
	print(self.header)
    # Iterator 
    def next(self):
	if self.handle.closed: 
	    raise StopIteration
	else: pass
	line = self.handle.readline()
	if not line:
	    self.handle.close()
	    raise StopIteration
	line = line.strip('\n').split('\t')
	propString = 'chr, pos, id, ref, alt, qual, filter, info, format'
	# Need to convert the propString into appropiat types
	conv = ''
	variantCall = collections.namedtuple('vcfline', propString)
	genotypes = line[9:]
	line = line[0:8]
	line.append(genotypes)
	modline = variantCall(*line)
	return modline

    def chunk(self, n, l):
	pass
    
    def __parse(self,line):
	    """Parses a VCF line and returns a collection
	    """
	    line = variant.strip('\n').split('\t') 
	    return line
