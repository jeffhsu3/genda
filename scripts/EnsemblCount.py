""" Count reads from an ensembl file
usage:
    python EnsemblCount.py GTFfile bamfile
"""

import HTSeq, optparse, sys

def main():
    p = optparse.OptionParser(__doc__)
    p.add_option("-D", "--debug", action="store_true", dest="D", help="debug")
    p.add_option("-r", "--refSeq", action="store_true", dest="r", help=\
            "Input GTF is from refSeq Flat File")

    options, args = p.parse_args()
    gtf_file = HTSeq.GFF_Reader(args[0])

    for feature in gtf_file:
        if feature.type == "exon":
            exons[ feature.iv ] = feature

if __name__ == '__main__':
    main()
