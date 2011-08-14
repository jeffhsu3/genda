"""
Takes a bigwig file and annotates a bed file
Usage:
    python BigWigAnnotateBed.py [OPTIONS] bed bigwig bed.out
"""

from bx.bbi.bigwig_file import BigWigFile
import optparse

def main():
    p = optparse.OptionParser(__doc__)
    p.add_option('-A', '--absolute', action='store_true',dest='A',default=False,\
	help='absolute threshold')
    p.add_option('-s','--standard_background', action='store_true', dest='stdbg')
    p.add_option('-D', '--debug', action='store_true', dest='debug')
    options, args = p.parse_args()
    debug_c = 0

    BEDFILE = open(args[0], 'rU')
    BW = BigWigFile(file=open(args[1]))
    BEDout = open(args[2], 'w')

    for line in BEDFILE:
	line = line.strip().split('\t')
	x = BW.query(line[0], int(line[1]), int(line[2]),1)
	line.append(str(round(x[0]['mean'], 5)))
	BEDout.write("\t".join(line)+"\n")
	"""
	for i in x:
	    print i['mean']
	"""

	if options.debug:
	    debug_c +=1
	    if debug_c >= 10:
		break


if __name__ == '__main__':
    main()
