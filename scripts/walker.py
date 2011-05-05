"""
    python walker.py [OPTIONS] <bed/vcf/gff/tab> <bamfile1> [bamfile2, ...]

    Right now takes a diff tab-delimited output from cuffdiff and grabs the
    sequences from bamfiles.  The goal of which is to take the sequences and
    Blast them to verify the mapping and ensure it is not simply an alignment
    artifact.

    :TODO
	Right now the gene.diff file only contains the bounding region, not the
	particular exons and such, where to get this information for the novel
	genes?  Is it in the cuffcompare output? gene.diff has XLOC

"""
import optparse, sys, csv, pysam
from pySeq.pysam_callbacks.seq_call import Seq_call
from pySeq.common_funcs import regionParse

def fetchSeq(diff_line, bamFile):
    region = regionParse(diff_line[2])
    seq_call = Seq_call()
    bamFile.fetch(region[0], region[1], region[2], callback=seq_call)
    return(seq_call.sequences)


def main():
    p = optparse.OptionParser(__doc__)
    p.add_option("-D", "--debug", action="store_true", dest="D", help="debug")
    p.add_option("-N", "--novel", action="store_true", dest="novel",\
	    help="Only check at novel loci")
    options, args = p.parse_args()

    fh = csv.reader(open(args[0], "rU"), delimiter="\t")
    bamFile = pysam.Samfile(args[1], "rb")
    header = fh.next()
    
    debug = 0 
    for line in fh:
	if options.novel:
	    if line[1] == "-":
		print(line[1])
		seqs=fetchSeq(line, bamFile) 
		print("\n".join(seqs))
	    else: pass
	else:
	    print(line[1])
	    seqs = fetchSeq(line, bamFile)
	    print("\n".join(seqs))


	if options.D:
	    debug += 1
	    if debug > 10:
		break


if __name__=='__main__':
    main()
