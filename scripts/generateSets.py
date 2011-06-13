"""
Takes a BAM file and generates counts
Usage:
"""
import re, sys, optparse, csv
import pysam

def main():
    p = optparse.OptionParser(__doc__)
    p.add_option("-D", "--debug", action="store_true", dest="D", help="debug")
    opts, args = p.parse_args()
    #fh = csv.reader(open(args[0], 'rU'), delimiter="\t") 
    bamfile = pysam.Samfile(args[0], 'r')
    bam_iter = bamfile.fetch()

    isoform_A = re.compile(r'*_A$')
    isoform_B = re.compile(r'*_B$')

    rawCounts = {}

    # Check for each queryname?
    for read in bam_iter:
	# EC: likely a read maps twice into the same transcript?
	rawCounts[read.rname] = rawCounts.get(read.rname, 0) + 1
	


    
    # Build Matrix M

if __name__ == '__main__':
    main()
