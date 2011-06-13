"""
Script to ascertain splicing varients

Usage:
    python geneCounts.py [options] GTFfile BAMFile


Author: Jeffrey Hsu
2010
"""

import optparse, pysam, sys, csv
from pySeq.pysam_callbacks.geneCounter import geneCounter

def convertEnsembl(chrom):
    if chrom == "MT":
	return("chrM")
    else:
	return("chr"+chrom)

def main():
    ###########################
    #
    # Option Parsing
    #
    ###########################
    
    # :TODO the contigs should be grabbed form the .bam header.  
    contigs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
	    "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
	    "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
	    "chrX", "chrY", "chrM"]
    p = optparse.OptionParser(__doc__)
    p.add_option("-D", "--debug", action="store_true", dest="D", help="debug")
    p.add_option("-r", "--refSeq", action="store_true", dest="r", help=\
	    "Input GTF is from refSeq Flat File")

    options, args = p.parse_args()
    debug = 0 
    fh = open(args[0], "rU")
    reader = csv.reader(fh, delimiter="\t")
    bamfile = pysam.Samfile(args[1], "rb")

    # Ensure that exons are not double counted
    seen_exon = []
    
    c = geneCounter("Gene")
    for line in reader:
	if line[2] == 'exon':
	    # Grab the identifiers
	    if options.r: 
		chrom = line[0]
		gene_name = line[8].split(";")[0][8:].strip('"')
		transcript_id = line[8].split(";")[1][15:].strip('"')
	    else: 
		chrom = convertEnsembl(line[0])
		gene_name = line[8].split(";")[3][11:].strip('"')
		transcript_id = line[8].split(";")[4][17:].strip('"')
	    if chrom in contigs:
		try:
		    start = int(line[3])
		    end = int(line[4])
		except ValueError:
		    print("Input GTF has non interger value where start or end\
			    position should be at line %s" % fh.tell())
    		    sys.exit()
		# :TODO Handle the case in which genes overlap, ie one exon.
		# Shouldn't be a problem if the GTF file clusters by gene

		# :TODO Doesn't handle imperfect start end matches  :(. So
		# near exon matches to the same gene (different isoforms)
		# overcount the number of reads. Although a read should be
		# placed in seen reads, and that should mitigate it?

		if gene_name == c.gene:
		    if (start,end) in seen_exon:
			pass
		    else:
			bamfile.fetch(chrom, start, end, callback= c)
			seen_exon.append((start,end))  
		else:
		    print("\t".join([c.gene, str(c.mCounts)]))
		    c = geneCounter(gene=gene_name)
		    bamfile.fetch(chrom, start, end, callback= c)
		    # Reset seen exons
		    # Change this a tree for faster look ups
		    seen_exon = []
		    seen_exon.append((start,end))

	    else:
		pass    
	    if options.D:
		debug += 1
		if debug > 1000:
		    break
		else:pass
	    else: pass
    print("\t".join([c.gene, str(c.mCounts)]))



if __name__=='__main__':
    main()
    #import doctest
    #doctest.testmod()
