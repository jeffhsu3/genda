"""
    Generates cut-site matrix for use in CENTIPEDE
Usage:
    python cutsite_matrix.py something.bam annotation.bed
"""
import pysam
import numpy as np	

def get_cutsite():
    pass

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
    PWM_bed = open(args[1], 'rU')
    debug = 0
    for line in PWM_bed:
	line = line.split('\t')
	chrom = line[0]
	start = int(line[1]) - 100 - 1
	end = int(line[2]) + 100
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
