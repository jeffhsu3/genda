#!/usr/bin/python
""" Creates a new fasta file by masking out.  Requries a decent amount of
memory to load the whole annotation file.  Doesn't assume the files are sorted
in the same manner (Actually it does for the moment). 

Written by Jeffrey Hsu

usage:
    python maskSNPs hg19.fa dbSNP-131.vcf
:TODO Need to deal with areas where there are more than 2 SNPS within a line,
recursion?
"""
_BASES = ["A", "G", "C", "T"]

import random

def create_keys(annot):
    # Parse the Header
    line = annot.readline().strip('\n')
    while line[0:2] == "##":
        line = annot.readline().strip('\n')
    l_dict = {} 
    def _parse(entry):
        """
        """
        _BASES = ["A", "G", "C", "T"]
        t = entry.split("\t")
        try:
            _BASES.remove(t[3])
            alts = t[4].split(",")
            _BASES.remove(alts[0])
            new_base = random.choice(_BASES)
            return(t[0], int(t[1]), new_base)
        except ValueError:
            return None
        except IndexError:
            return None

    data = annot.read()
    # Load the dictionaries in 
    for i in data.split('\n'):
        k = _parse(i)
        if k:
            try:
                l_dict[k[0]].append((k[1],k[2]))
            except KeyError:
                l_dict[k[0]] = []
                l_dict[k[0]].append((k[1],k[2]))
        else: pass

    for i in l_dict:
        l_dict[i].sort(key = lambda entry: entry[0]) 

    return(l_dict)


def main():
    import sys 

    annot = open(sys.argv[1], 'rU')
    fasta = open(sys.argv[2], 'rU')
    fileout = open(sys.argv[3], 'w')
    l_dict = create_keys(annot)


    debug = 0

    for line in fasta:
        if line[0] == ">":
            position = 0
            try:
                fileout.write(line)
                chrom = line.lstrip(">chr").rstrip('\n')
                print(chrom)
                k = l_dict[chrom].pop(0)
                print(k)
            except ValueError:
                pass
        else: 
            temp = list(line)
            while position < k[0] < position + 50:
                qbase = k[0] - position - 1   
                print("Replacement", position,k[0], qbase, k[1])
                temp[qbase] = k[1]
                k = l_dict[chrom].pop(0)
                print("".join(temp))
            out = "".join(temp) + "\n"
            fileout.write(out)
            position += 50

        debug += 1
        if debug >= 11000:
            break
        else: pass

if __name__ == '__main__':
    main()
        
