"""Script to annotate SNPs with missense mutations
"""

import optparse
import itertools
import GTF
from pySeq.formats.VCF import VCFfile
from pySeq.parsing.fastA import FastA
from Bio.Seq import Seq

def build_sequence():
    pass

def main():
    p = optparse.OptionParser(__doc__)
    p.add_option("-D", "--debug", action="store_true", dest="debug", help="debug")
    p.add_option("-r", "--refSeq", action="store_true", dest="r", help=\
            "Input GTF is from refSeq Flat File")
    p.add_option("-o", "--output", action="store", dest="o", help=\
            "Output file, if not set will default to INPUT.an.bam")
    options, args = p.parse_args()

    gtf = GTF.GTF(args[0])
    gtf.load_gtf()

    vcf = VCFfile(args[1])
    fa = FastA(args[2])

    pfa = Seq(open("HCN4.pfa.txt", "rU").read())
    print("Length of the protein sequence:%i" % len(pfa))

    bleh = open("test.txt", "w+")

    sequence = ""
    first = True
    for x in gtf.data:
        print(x.start, x.end)
        if first:
            # Get last codon
            gstart = x.start - 3
            gend = x.end - 1
        else:
            gstart = x.start
            gend = x.end
        temp = Seq(fa.grabSequence("chr15", gstart+1, gend+1).rstrip("\n"))
        bleh.write(str(temp))
        temp = temp.reverse_complement()
        bleh.write(str(temp))
        bleh.write("\n")
        bleh.write("\n")
        sequence = temp + sequence
        first = False
    print("Our sequence:%i " % len(sequence))
    print("pp: %i " % len(sequence.translate()))
    """
    bleh.write(str(sequence).upper())
    bleh.write("\n")
    bleh.write("\n")
    bleh.write(str(sequence.translate()).upper())
    bleh.write("\n")
    bleh.write(str(pfa))
    """
    print(str(sequence.translate().upper()) == str(pfa).upper())
    print "\t".join(["rsID", "chrom", "pos", "bp position","ref", "alt", "AA position", "AA Ref"])
    for line in vcf:
        try:
            chrom = int(line.chr[3:])
            pos = int(line.pos)
        except TypeError:
            chrom = ord(line.chr[3:])
        if chrom == 15:
            hit = itertools.takewhile(lambda x: x.start < pos and x.end > pos, gtf.data)
            for i in hit:
                position_within_gene = i.end - pos + i.dist_tss + 2

                aa_num = position_within_gene/3+1
                out_list = [line.id, str(chrom), str(pos), str(position_within_gene), line.ref, line.alt, str(aa_num), pfa[aa_num-1]]
                print("\t".join(out_list))
        else: pass


if __name__ == '__main__':
    main()
