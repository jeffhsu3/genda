""" Converts chr1 -> 1 or 1 -> chr1 for use in other downstream applications,
ie BEDTOOLS
"""

import sys

def main():
    fh = open(sys.argv[1], "rU")

    for line in fh:
        line = line.rstrip("\n").split("\t")
        line[0] = line[0].lstrip("chr")
        print("\t".join(line))


if __name__ == '__main__':
    main()
