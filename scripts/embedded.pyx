#cython: embedsignature=True
#cython: profile=True
"""
Script to make all sequence reads match the reference sequence

Usage:
    python bamtoreference [options] BAMFile reference_file

Author: Jeff Hsu
2012
"""


import optparse

cimport csamtools as csam

from pysam import Fastafile
# For some reason csamtools.fastaFile isn't working for me

#import bcolors

def main():
    p = optparse.OptionParser(__doc__)
    p.add_option("-D", "--debug", action="store_true", dest="debug", help="debug")
    p.add_option("-r", "--refSeq", action="store_true", dest="r", help=\
            "Input GTF is from refSeq Flat File")
    p.add_option("-o", "--output", action="store", dest="o", help=\
            "Output file, if not set will default to INPUT.an.bam")
    options, args = p.parse_args()



    # File Inputs ####################################################

    samfile = csam.Samfile(args[0], 'rb')
    fasta = Fastafile(args[1])

    # :TODO add this program into the header
    if not options.o:
        outfile = csam.Samfile(args[0].rstrip(".bam") + "anon.bam", "wb", template=samfile)
    else:
        outfile = csam.Samfile(options.o, "wb", template=samfile)




    ##################################################################
    cdef csam.AlignedRead read
    cdef int debug_num
    cdef int ref_id
    cdef int first_exon_len, second_exon_len, junction_len
    cdef bytes new_string, first_exon, second_exon, chrom
    debug_num = 0

    for read in samfile:
        try:
            # For a sorted file this should remain the same for most
            chrom = samfile.getrname(read.tid)
        except ValueError:
            # I think tophat puts unmapped reads as tid = -1 and that thses occur
            # at the end
            break
        try:
            # :TODO what are the other possibilities here, such as a read spanning
            # multiple exons
            first_exon_len = read.cigar[0][1]
            junction_len = read.cigar[1][1]
            second_exon_len = read.cigar[2][1]
            first_exon = fasta.fetch(chrom, read.pos,
                    read.pos + first_exon_len).upper()
            second_exon = fasta.fetch(chrom,
                    read.pos + first_exon_len + junction_len,
                    read.pos + first_exon_len + junction_len + second_exon_len).upper()
            new_seq = first_exon + second_exon
        except IndexError:
            new_seq = fasta.fetch(chrom, read.pos, read.pos + read.qlen).upper()

        read.seq = new_seq
        outfile.write(read)



        if debug_num > 10000 and options.debug:
            break
        else:
            debug_num += 1

    samfile.close()
    outfile.close()
    fasta.close()



if __name__ == '__main__':
    main()
