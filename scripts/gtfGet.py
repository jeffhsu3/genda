"""
From a list of identifiers extracts the GTF file with only those entries
Usage:
    python gtfGet [options] <ID.list> <GTF_FILE>
"""
import optparse, sys, csv 
from pySeq.parsing.gtfTools import parseID

def main():
    p = optparse.OptionParser(__doc__)
    p.add_option("-D", "--debug", action="store_true", dest="D", help="debug")
    p.add_option("-d", "--database", dest="d",\
	    help="database")
    p.add_option("-T", "--transcript", dest="t",\
	    help="Transcript_id")
    p.add_option("-C", "--class", dest="C",\
	    help="Look only at a particular Cufflinks class of transcript")
    options, args = p.parse_args()
    debug = 0

    IDs = open(args[0], "rU").read().rstrip("\n").split("\n")
    fh = csv.reader(open(args[1], "rU"), delimiter="\t")
    

    ident = "transcript_id"

    for line in fh:
	transcript = parseID(line, ident)
	if transcript in IDs and line[2]=="exon":
	    print("\t".join(line))
	else: pass
	if options.D:
	    debug += 1
	    if debug > 40:
		break


if __name__=='__main__':
    main()
