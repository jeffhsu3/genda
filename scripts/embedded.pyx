#cython: embedsignature=True
#cython: profile=True
"""
Script to make all sequence reads match the reference sequence

Usage:
    python bamtoreference [options] BAMFile reference_file

Author: Jeff Hsu
2012
"""


import optparse, sys
#from csamtools cimport *
cimport csamtools as csam
from pySeq.parsing.fastA import FastA

import bcolors

cdef handle_exons():
    """
    """
    pass

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
    fasta = FastA(args[1])

    # :TODO add this program into the header
    if not options.o:
        outfile = csam.Samfile(args[0].rstrip(".bam"), "wb", template=samfile)
    else:
        outfile = csam.Samfile(options.o, "wb", template=samfile)




    ##################################################################
    cdef csam.AlignedRead read
    cdef int debug_num
    cdef int ref_id
    cdef int first_exon_len, second_exon_len, junction_len
    cdef bytes new_string
    debug_num = 0

    for read in samfile:
        try:
            # :TODO what are the other possibilities here
            first_exon_len = read.cigar[0][1]
            junction_len = read.cigar[1][1]
            second_exon_len = read.cigar[2][1]
            print(first_exon_len, junction_len, second_exon_len)
            print('Exon junction')
            # :TODO C string handling
            first_exon = fasta.grabSequence(samfile.getrname(read.tid), read.pos+1,
                    read.pos + first_exon_len).upper()
            second_exon = fasta.grabSequence(samfile.getrname(read.tid),
                    read.pos + first_exon_len + junction_len,
                    read.pos + first_exon_len + junction_len + second_exon_len).upper()
            new_seq = first_exon + second_exon
            print(read.seq)
            print(bcolors.bcolors.WARNING + first_exon + bcolors.bcolors.ENDC\
                    + bcolors.bcolors.OKBLUE + second_exon + bcolors.bcolors.ENDC)
        except IndexError:
            new_seq = fasta.grabSequence(samfile.getrname(read.tid), read.pos+1, read.pos + read.qlen).upper()
        #print(read.seq)
        #print(bcolors.bcolors.WARNING + str(new_seq)+ bcolors.bcolors.ENDC)

        read.seq = new_seq
        outfile.write(read)



        if debug_num > 20000 and options.debug:
            break
        else:
            debug_num += 1

    samfile.close()
    outfile.close()



if __name__ == '__main__':
    main()
