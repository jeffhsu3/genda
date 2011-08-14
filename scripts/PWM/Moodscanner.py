"""
Creates a BEDFile with each row being a matrix match

Usage:
    python Moodscanenr.py [OPTIONS] pwmfile sequence.fa

Kai Smith 
Jeff Hsu
"""

import MOODS, pySeq.parsing.PWMparser, Bio.Seq, Bio.SeqIO, sys, optparse

def strand_adjust(position, match_size):
    """
    Gets the position and adjusts the start, end appropiately if the strand is
    negative and also returns the strand
    :TODO Doctest
    """
    if position > 0:
	strand = '+'
	start = position
	end = position + match_size
    else: 
	strand = '-'
	start = position*-1
	end = position*-1 + match_size
    return(start, end, strand)

def main():
    p = optparse.OptionParser(__doc__)
    p.add_option('-t', '--thresh', action='store', dest='threshold', default=0.0,\
	    help='determines threshold')
    p.add_option('-a', '--append', action='store', dest='name',default='resultsfor',\
	    help='appends pwm name to this when creating files')
    p.add_option('-A', '--absolute', action='store_true',dest='A',default=False,\
	    help='absolute threshold')
    p.add_option('-s','--standard_background', action='store_true', dest='stdbg')
    p.add_option('-M', '--specific_Matrix', action='store', dest='specific')
    options, args = p.parse_args()

    pwm = open(args[0], 'rU')
    fa = open(args[1], 'rU')
    pfa = list(Bio.SeqIO.parse(fa, 'fasta'))
    index, matricies, sizes = pySeq.parsing.PWMparser.parse(pwm)

    underorequal20 = []
    over20 = []
    under20names = []
    over20names = []
    pwmdata={}
    fileout = {}
    bgt = False
    if options.stdbg:
	bgt = [0.25,0.25,0.25,0.25]
 

    # Construct Matrices to search and files to write to.  
    for k in index.keys():
	if options.specific:
	    if k == options.specific:
		filename = options.name + k + '.bed'
		fileout[k] = open(filename, 'w')
		if sizes[k] <= 20:
		    underorequal20.append(matricies[k])
		    under20names.append(k)
		else:
		    over20.append(matricies[k])
		    over20names.append(k)
	else:
	    filename = options.name + k + '.bed'
	    fileout[k] = open(filename, 'w')
	    if sizes[k] <= 20:
		underorequal20.append(matricies[k])
		under20names.append(k)
	    else:
		over20.append(matricies[k])
		over20names.append(k)
    

    for chrom in pfa:
	print(chrom.name)
	#Run under 20s
	# Should we sort the results as all downstream applications require a
	# sort first

	res = MOODS.search(chrom.seq, underorequal20, float(options.threshold),
		absolute_threshold=options.A , both_strands =
		True,bg=bgt,
		algorithm='lf')


	for n,r in enumerate(res):
	    for position,score in r:
		start, end, strand = strand_adjust(position,
			sizes[under20names[n]])
		fileout[under20names[n]].write('\t'.join([chrom.name,
		    str(start), str(end),under20names[n],str(score), strand])+'\n')

	#Run over 20s
	res = MOODS.search(chrom.seq, over20, float(options.threshold),
		absolute_threshold=options.A , both_strands =
		True,bg=bgt,
		algorithm='supera')

	for n,r in enumerate(res):
	    for position,score in r:
		start, end, strand = strand_adjust(position,
			sizes[over20names[n]])
		fileout[over20names[n]].write('\t'.join([chrom.name, str(start),
		    str(end),over20names[n],str(score)],strand)+'\n')


if __name__ == '__main__':
    main()
