

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
            pass


        set_list = [i]
    while len(rs) > 0:
        intersection = [t & i for i in set_list]
        out_sets.append(i)

    return(out_sets)


def compare_two_transcripts(trans1, trans2, transcript_dict):
    """
    Returns the splice differences between two transcripts.
    Note this ignores TSS and just looks for splice junctions

    Parameters
    ----------
    :TODO refactor this
    :TODO still not finished
    """

    t1 = transcript_dict[trans1]
    t2 = transcript_dict[trans2]

    starts1 = [i[0] for i in t1]
    starts2 = [i[0] for i in t2]
    # Assume sorted
    reverse = False
    if min(starts1) >= min(starts2):
        s1 = t2 
        s2 = t1
        reverse = True
    else:
        s1 = t1
        s2 = t2
    # Assume sorted order, ordering in reverse by construction
    if reverse:
        torder = (trans2, trans1)
    else:
        torder = (trans1, trans2)
    oc = 0
    #start_matches = False
    exclusive_juncs = []
    matching_exons = []
    skipped_exons = []
    for i in s1:
        if (i[0] < s2[0][0]) and (i[1] <= s2[0][0]): 
            pass
        elif (i[0] < s2[oc][0]) and (i[1] > s2[oc][0]):
            # 5' difference
            print('5 primt diff')
            sdiff = i[0] - s2[oc][0]
            ediff = i[1] - s2[oc][1]
            print(s2[oc])
            exclusive_juncs.append((i[0], s2[oc][0], 
                sdiff, ediff, (i[2], s2[oc][2])))
            oc += 1
        elif (i[0] > s2[oc][0]) and (i[0] < s2[oc][1]):
            print('3 prime diff')
            print(i[0] - s2[oc][0])
            print(i)
            print(s2[oc])
            #:TODO FIX this
            oc += 1
        elif (i[0] == s2[oc][0]) and (i[1] == s2[oc][1]):
            matching_exons.append((i[2], s2[oc][2]))
            oc += 1
            pass
    return(exclusive_juncs, torder)
