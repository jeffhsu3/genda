"""
    python walker.py [OPTIONS] <bed/vcf/gff/tab> <bamfile1> [bamfile2, ...]

    Right now takes a diff tab-delimited produced by cufflinks output from cuffdiff and outputs a
    BEDFILE sequences.


"""
import optparse, sys, csv, os, re 
from pySeq.common_funcs import regionParse

def main():
    p = optparse.OptionParser(__doc__)
    p.add_option("-D", "--debug", action="store_true", dest="D", help="debug")
    p.add_option("-N", "--novel", action="store_true", dest="novel",\
	    help="Only check at novel loci")
    options, args = p.parse_args()

    fh = csv.reader(open(args[0], "rU"), delimiter="\t")
    
    debug = 0 
    unanno = re.compile(r"^RP[0-9]")
    for line in fh:
	if line[11] == "yes" and (line[1] =="-" or unanno.match(line[1])):
	    region, start, end = regionParse(line[2])
	    print("\t".join([region, str(start-1), str(end), line[0], line[1]]))


	if options.D:
	    debug += 1
	    if debug > 10:
		break


if __name__=='__main__':
    main()

