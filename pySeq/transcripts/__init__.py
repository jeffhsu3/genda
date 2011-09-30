import numpy as np


class SNPInfo(object):
    """ Docstring

    """

    def __init__(self, position):
        self.position = position

    def __call__(self, position):
        pass

class Exon(object):
    """ A region
    """

    def __init__(self, region, start, end):
        self.region = region
        try:
            self.start = int(start)
            self.end = int(end)
        except ValueError:
            print("Start and End positions need to be integers")

class Transcript(object):
    """ A collection of exons
    """

    def __init__(self, exons=[]):
       self.exons = exons


class Gene(object):
    """ A collection of transcripts
    """

    def __init__(self, transcripts=[]):
        self.transcripts = []

    def _unique_junctions(self):
        pass

    def hmm(self):
        pass

def break_exons(exon1, exon2):
    """ Returns a list of new exons

    Breaks a interval into 2 intervals if it overlaps, upstream exon always
    comes first.

    """
    if exon1.region == exon2.region:
        if exon2.start > exon1.start and exon2.end < exon1.end:
            # Case where exon1 is wholly enclosed in exon2
            return([Exon(exon1.region, exon1.start, exon2.start),
                    Exon(exon1.region, exon2.start, exon2.end),
                    Exon(exon1.region, exon2.end, exon1.end)])
        elif exon2.end < exon1.end:
            return([Exon(exon1.region, exon1.start, exon2.start),
                    Exon(exon1.region, exon2.start, exon1.end),
                    Exon(exon1.region, exon1.end, exon2.end)])
        else:
            # Do nothing if there is no overlap
            return(exon1, exon2)
    else:
        pass

def unique_sets(set_list):
    """ Returns a list of sets
    Calculates the unique regions defined such that only the same elements
    are bound to it.
    """

    # Since most exons are shared, remove the most common elements first
    common_set = reduce( lambda x, y: x & y, list_of_sets)
    rs = [i - common_set for i in list_of_sets]
    out_sets = []
    all_unique = True
    n = len(set_list)
    for i in set_list:
        intersection = [t & i for t in set_list]
        # uniqufy the intersections
        u_s = set(intersection)
        for t in u_s:
                

        set_list = [i
    while len(rs) > 0:
        intersection = [t & i for i in set_list]
        out_sets.append(

    return(out_sets)

