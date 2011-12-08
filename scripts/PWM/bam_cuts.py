"""
    Splits bam by whether they are positive or negative strand 
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
    p.add_option("-x", "--file1", dest="positive", help="Write the\
            positive reads to this file")
    p.add_option("-y", "--file2", dest="negative", help="Write the\
            positive reads to this file")
    options, args = p.parse_args()


    bamfile = pysam.Samfile(args[0], 'rb')
    bamout1 = pysam.Samfile(options.positive, "wb", template=bamfile)
    bamout2 = pysam.Samfile(options.negative, "wb", template=bamfile)
    PWM_bed = open(args[1], 'rU')
    debug = 0
    for line in PWM_bed:
        line = line.split('\t')
        chrom = line[0]
        start = int(line[1]) - 100 - 1
        end = int(line[2]) + 100
        diff = end-start
        try:
            for alignment in bamfile.fetch(chrom, start, end):
                if alignment.is_reverse:
                    alignment.pos = (alignment.pos-start) + 200000 
                    alignment.rname= 1
                    if alignment.pos < 0:
                        pass
                    else:
                        bamout2.write(alignment)
                else:
                    alignment.pos = (alignment.pos-start) + 200000 
                    alignment.rname = 1
                    if alignment.pos < 0:
                        pass
                    else:
                        bamout1.write(alignment)
        except ValueError:
            pass

        if options.D:
            debug += 1
            if debug >= 400:
                break


if __name__ == '__main__':
    main()
