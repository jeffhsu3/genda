""" Python Implementation of Jim Kent's genomic binning
http://genome.cshlp.org/content/12/6/996.full
http://genomewiki.ucsc.edu/index.php/Bin_indexing_system
"""

from collections import defaultdict

binOffsets = [512+64+8+1, 64+8+1, 8+1, 1, 0]
_BINFIRSTSHIFT = 17
_BINNEXTSHIFT = 3

class GTab(object):
    """ Object with start and end and optional info for holding row information
    from tab delimited genomic files
    """
    __slots__ = ('start', 'end', 'info')
    def __init__(self, start, end, info=None):
        self.start = start
        self.end = end
        self.info = info

def binFromRangeStandard(start, end):
    global binOffsets
    global _BINFIRSTSHIFT
    global _BINNEXTSHIFT
    startBin = start
    endBin = end - 1
    startBin >>= _BINFIRSTSHIFT
    endBin >>= _BINFIRSTSHIFT
    for i in binOffsets:
        if startBin == endBin:
            return i + startBin
        else: pass
        startBin >>= _BINNEXTSHIFT
        endBin >>= _BINNEXTSHIFT
    print("start %d, end %d out of range in findBin (max is 512M)", start,
          end)
    return 0

def read_and_Kent_index(filename):
    """ Create an index for a BED or VCF file
    """
    chr_dict = defaultdict(lambda : defaultdict(list))
    debug = 0
    with open(filename, 'rU') as fh:
        # Skip comment lines
        # :TODO Fix this and make more general
        fh.next()
        fh.next()
        for line in fh:
            p_line = line[:-1].split("\t")
            try:
                start = int(p_line[1])
                end = int(p_line[2])
                kent_bin = binFromRangeStandard(start, end)
            except ValueError:
                # Case for VCF files
                start = int(p_line[1]) - 1
                end = int(p_line[1])
                kent_bin = binFromRangeStandard(start, end)
            chr_dict[p_line[0]][kent_bin].append(GTab(start, end))
            """
            if debug > 100000:
                break
            else:
                debug += 1
            """
    return(chr_dict)

def get_overlaps_counts(region, start, end, tab_map):
    """
    """
    global binOffsets, _BINFIRSTSHIFT, _BINNEXTSHIFT
    startbin = start >> _BINFIRSTSHIFT
    endbin = end >> _BINFIRSTSHIFT
    counts = 0
    for offset in binOffsets:
        for i in range(startbin+offset, endbin+offset+1):
            potential = tab_map[region][i]
            #print("%s : %s " % (i, len(potential)))
            for k in potential:
                max_start = max(k.start, start)
                min_end = min(k.end, end)
                overlap = min_end - max_start
                if overlap > 0:
                    counts += 1
                else:
                    pass
        startbin >>= _BINNEXTSHIFT
        endbin >>= _BINNEXTSHIFT
    return(counts)



if __name__ == '__main__':
    import sys
    t = read_and_Kent_index(sys.argv[1])
    interestBed = open(sys.argv[2], 'rU')
    debug = 0
    for line in interestBed:
        line = line[:-1].split('\t')
        counts = get_overlaps_counts(line[0], int(line[3]) - 10000,
                                     int(line[4]) + 10000, t)
        if counts > 0:
            info = line[8].split(";")
            transcript_id = info[1][15:].replace('"', '')
            gene_symbol = info[3][10:].replace('"', '')
            print("\t".join([transcript_id, gene_symbol, str(counts)]))
