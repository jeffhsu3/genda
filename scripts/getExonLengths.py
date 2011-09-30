"""
usage:
    python getExonLengths.py GTF file 
"""

from exonUtils import uniqufy, overlapCases, collideIntervals, swapIntervals,\
        collapseIntervals, calcOverlap 


def calcOverlap(intervals):
    """
    Calculates the number of base pairs that overlap in a list of intervals

    All intervals must be unique
    :TODO Full enclosed case fails.  Probably should just use interval trees,
    but constructing trees might take longer, if matches are sparse? or would

    """
    bp = 0 
    for i in intervals:
        bp += sum([overlapCases(i, j) for j in intervals])
    return(bp)



def main():
    import optparse, sys
    p = optparse.OptionParser(__doc__)
    opts, args = p.parse_args()
    fh = open(args[0], 'rU')
    debug=0
    r = fh.read()
    t = r.replace('\n', '\t').split('\t')
    r = r.split('\n')
    r = [i.split('\t') for i in r]
    gene_info = t[8::9]
    genes = [i[9:i.find('";')] for i in gene_info]
    u_gene = uniqufy(genes)
    gene_dict = {}
    for gene in u_gene:
        ints = [(int(r[gr][3]),int(r[gr][4]))\
                for gr in [i for i,x in enumerate(genes) if x == gene] if\
                    r[gr][2]=='exon']
        ints = uniqufy(ints)
        gene_dict[gene] = str(sum([i[1]-i[0]+1 for i in\
                collapseIntervals(ints)]))
    for i in gene_dict:
        print("\t".join([i,gene_dict[i]])) 



    """
    dict_length = {}
    for line in fh:
        if line[2] == "exon":
            line = line.split('\t')
            gene = line[8][9:line[8].find('";')]
            dict_length[gene] = dict_length.get(gene, 0) +\
                    int(line[4])-int(line[3])
        else:
            pass

        debug += 1
        if debug > 30:
            break
    print(dict_length)
    """


if __name__ == '__main__':
    main()
