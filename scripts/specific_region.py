import sys, pysam
import numpy as np
import rpy2
from pySeq.common_funcs import regionParse 

def main():
    import optparse
    p = optparse.OptionParser(__doc__)
    p.add_option("-D", "--debug", action="store_true", dest="D", help="debug")
    p.add_option("-S", "--stam", action="store_true", dest="S", help="DNAseI\
	    is generated from STAM's group")
    p.add_option("-s", "--shift", action="store", dest="shift", help="Amount to\
	    shift the negative strand", default = 36)

    options, args = p.parse_args()

    options.shift = int(options.shift)
    bamfile = pysam.Samfile(args[0], 'rb')
    region_str = regionParse(args[1])

    chrom = region_str[0]
    start = region_str[1]
    end = region_str[2]
    diff = end-start
    a = np.zeros(2*(diff), dtype=np.int)
    try:
	for alignment in bamfile.fetch(chrom, start, end):
	    if alignment.pos-start < 0: pass
	    else:
		if alignment.is_reverse:
		    try:
			a[alignment.pos-start+diff+options.shift] += 1
		    except IndexError:
			pass
		else:
		    a[alignment.pos-start] += 1
    except ValueError:
	pass
    print("\t".join(map(str,a)))

    if options.D:
	debug += 1
	if debug >= 400:
	    break
	    

if __name__ == '__main__':
    main()

