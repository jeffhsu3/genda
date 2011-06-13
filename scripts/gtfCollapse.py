"""
Script to ascertain splicing varients

Usage:
    python geneCounts.py [options] GTFfile BAMFile


Author: Jeffrey Hsu
2010
"""

import optparse, pysam, sys, csv
import numpy as np
from bx.intervals.intersection import IntervalTree, Interval
from pySeq.pysam_callbacks.geneCounter import geneCounter

def convertEnsembl(chrom):
    if chrom == "MT":
	return("chrM")
    else:
	return("chr"+chrom)

def compareInt(int1, int2):
    """
    Compare two intervals and returns the whole span if two of them are equal. 
    >>> compareInt((1,10),(3,25))
    (1, 25)
    """
    assert int1[0] <= int1[1]
    assert int2[0] <= int2[1]
    if int1[0] <= int2[0] <= int1[1]:
	return(int1[0], max(int2[1], int1[1]))
    elif int1[0] <= int2[1] <= int1[1]:
	return(min(int1[0], int2[0]), int1[1])
    else: 
	return(None) 

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

    options, args = p.parse_args()
    debug = 0 
    fh = csv.reader(open(args[0], "rU"), delimiter="\t")
    bamfile = pysam.Samfile(args[1], "rb")
    trees = {}
    for i in contigs:
	trees[i] = IntervalTree()
    
    
    for line in fh:
	#Generates Tree
	if line[2] == 'exon':
	    gene_name = line[8].split(";")[3][11:].strip('"')
	    transcript_name = line[8].split(";")[4][17:].strip('"')
	    start, end = (int(line[3]), int(line[4]))
	    chrom=convertEnsembl(line[0])
	    intersects = trees[chrom].find(start,end)
	    if len(intersects) > 0:
		for i in intersects:
		    i.value.append(transcript_name)
	    else:
		trees[chrom].insert_interval(Interval(start, end, [transcript_name]))

    # Once tree is generated for a GTF file should save it.
    # Requries a dict lookup for each one. Probably still faster than going to
    # disk
    gene_count = {} 

    for read in bamfile.fetch():
	if read.is_proper_pair and read.is_read1 and read.mapq > 10:
	    # Probably a faster way to do this
	    cigar_types = [x[0] for x in read.cigar]
	    intersect_1 = trees[bamfile.getrname(read.rname)].find(read.pos, read.pos)
	    if len(intersect_1) >= 0:
		for i in intersect_1:
		    for transcript in i.value:
			gene_count[transcript] = gene_count.get(transcript, 0)\
				+ 1
	    else: pass
	if options.D:
	    debug += 1
	    if debug > 2000:
		break
	    else:pass
	else: pass

    for i in gene_count.keys():
	print("\t".join([i, str(gene_count[i])]))


if __name__=='__main__':
    main()
    #import doctest
    #doctest.testmod()
