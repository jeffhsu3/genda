from itertools import tee, izip
import requests
# Requires bx python
from bx.intervals.intersection import Intersecter, Interval
from HTSeq import GenomicInterval as GI



def pairwise(iterable):
    a, b = tee(iterable)
    next(b, None)
    return(izip(a,b))


class EventCollection(object):
    """ A collection of diffevents

    Arguments
    ---------
    """

    def __init__(self, events=[]):
        self.events = events




class DiffEvent(object):
    """A differential splicing event between
    two transcripts

    Arguments
    ---------
    event_type : ['skipped_exon', 'mxe', 'A5SE', 'ATS']
    start : start position of differential event
    end : end position of differential event
    transcript_id : 
    chrom : optional
    exon_num : exon number 
    cigar : a tuple of cigar differences 
    """

    def __init__(self, event_type, start, end,
            transcript_ids, cigar=None, chrom=None, 
            exon_num=None):
        self.event_type = event_type
        self.start = start
        self.end = end
        self.transcript_ids = transcript_ids
        self.exon_num = exon_num
        self.chrom = chrom
        self.cigar = cigar

    def __repr__(self):
        return(str(self.start) + '-' +  str(self.end) +\
                ':' + str(self.cigar))

    def __eq__(self, other):
        if self.start == other.start and\
                (self.end == other.end):
            return(True)
        else:
            return(False)

    def _extend(self, new_diff_event):
        """Extends a diffevent for example
        a mutual exclusive exon or two skipped
        exons and other complex events
        """
        self.cigar = new_cigar
        raise NotImplementedError



class Exon(object):
    """ A region
    """
    def __init__(self, region, start, end):
        self.region = region
        try:
            self.start = int(start)
            self.end = int(end)
        except ValueError:
            print("start and end positions need to be integers")


def get_transcript_ids(gene, gff=None):
    """ Get transcript ids from Ensembl REST api
    Returns a list of transcripts,

    :TODO make out_transcripts a list of transcript class
    :TODO try GTF first and except to Ensembl otherwise
    """
    server = "http://rest.ensembl.org"
    ext = "/lookup/id/{0}?species=homo_sapiens;expand=1"
    r = requests.get(server+ext.format(gene), headers={ "Content-Type" :
        "application/json"})
    out_transcripts = []
    decoded = r.json()
    for i in decoded['Transcript']:
        out_transcripts.append(i['id'])
    return(out_transcripts, 
            decoded['display_name'],
            decoded['seq_region_name'], 
            (int(decoded['start']), int(decoded['end'])))


class Transcript(object):
    """ A collection of exons
    """
    def __init__(self, exons=[]):
       self.exons = exons






def break_exons(exon1, exon2):
    """ 
    Returns a list of new exons by breaking an interval into 2 intervals if it overlaps, 
    upstream exon always
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
    """ 
    Returns a list of sets
    Calculates the unique regions defined such that only the same elements
    are bound to it.
    """

    # Since most exons are shared, remove the most common elements first
    common_set = reduce( lambda x, y: x & y, set_list)
    rs = [i - common_set for i in set_list]
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
    trans1 : string of transcript of interest
    trans2 : string of second transcript of interest
    transcript_dict : a dictionary of transcript names with 
    values being a list of exons

    Returns
    -------
    Exclusive Junctions : 
    5' upstream exons : 
    3' downstram exons : 
    Skipped Exons : Diffevent 

    """
    # TODO refactor this
    t1 = transcript_dict[trans1]
    t2 = transcript_dict[trans2]
    tree = Intersecter()
    starts1 = [i[0] for i in t1]
    starts2 = [i[0] for i in t2]
    reverse = False
    if min(starts1) <= min(starts2):
        s1 = t1
        s2 = t2
        s2_beg = min(starts2)
    else:
        s1 = t2 
        s2 = t1
        reverse = True
        s2_beg = min(starts1)
    if reverse: torder = (trans2, trans1)
    else: torder = (trans1, trans2)
    for i in s1:
        tree.add_interval(Interval(int(i[0]), int(i[1]), 
            value={'anno':i[2]}))
    matching_exons = []
    exclusive_juncs = []
    skipped_exons = []
    # Perform the query
    start_of_exons = None
    s1.sort(key=lambda x: x[0])
    s2.sort(key=lambda x: x[0])
    max_exon_1 =  s1[-1][2]
    max_exon_2 = s2[-1][2]
    end_position_s2 = max([i[1] for i in s2])
    s1_end = max([i[1] for i in s1])
    pcurr = 0
    for start, end, exon_n in s2:
        start = int(start)
        end = int(end)
        overlap = tree.find(int(start), int(end))
        if len(overlap) == 0:
            if start_of_exons and (start < s1_end):
                pint = s2[pcurr-1][0:2]
                nint = s2[pcurr + 1][0:2]
                c1 = max([start - pint[1], start-nint[1]])
                c2 = max([nint[0]-end, pint[0] - end])
                cigar = [(3, c1), 
                        (0, end-start), 
                        (3, c2)]
                # :TODO handle multiple skips
                skipped_exons.append(DiffEvent('skipped_exon', start, end,
                        torder, cigar=cigar, exon_num = (None, exon_n)))
            elif start > s1_end:
                break
            else: pass
        elif len(overlap) == 1:
            if start_of_exons: pass
            else: start_of_exons = overlap[0].value['anno']
            if start == overlap[0].start and end == overlap[0].end:
                matching_exons.append((start, end, (overlap[0].value['anno'], 
                    exon_n), (0, 0)))
            else:
                sstart = min(start, overlap[0].start)
                ssend = max(end, overlap[0].end)
                # Ignore 5' or 3' differences
                if (exon_n == max_exon_2 and
                        overlap[0].value['anno'] == max_exon_1):
                    pass
                else:
                    exclusive_juncs.append(
                            (sstart, ssend,
                            (overlap[0].value['anno'], exon_n), 
                            (overlap[0].start - start, overlap[0].end - end) ))
        else:
            if start_of_exons:
                pass
            else: start_of_exons = overlap[0].value['anno']
        pcurr += 1
    # Exons in s1 that are hit
    hit_exon = [i[2][0] for i in exclusive_juncs] 
    hit_exon.extend([i[2][0] for i in matching_exons])
    # Case where there is no overlap
    if s2_beg > s1_end: pass
    else:
        pcurr = 0
        for start, end, exon_n in s1:
            if exon_n <= start_of_exons: pass
            elif exon_n in hit_exon: pass
            elif start > end_position_s2: break
            else:
                pint = s1[pcurr - 1][0:2]
                nint = s1[pcurr + 1][0:2]
                c1 = max([start - pint[1], start- nint[1]])
                c2 = max([nint[0]-end, pint[0] - end])
                cigar = [(3, c1), 
                        (0, end-start), 
                        (3, c2)]
                diffevent = DiffEvent('skipped_exon', start, end,
                        torder, cigar = cigar, exon_num = (exon_n, None))
                skipped_exons.append(diffevent) 
            pcurr += 1
    return(exclusive_juncs, torder, matching_exons, skipped_exons)


def pairwise_transcript_comparison(transcript_dict):
    """ Returns pairwise transcripts where there are skipped exons 
    and the exon that is skipped.
    :TODO add this to the gene class?
    """
    skipped_exons_out = []
    for key1, key2 in pairwise(transcript_dict.keys()):
        _, _, _, skipped_exons = compare_two_transcripts(
                key1, key2, transcript_dict)
        if len(skipped_exons) >=1:
            skipped_exons_out.append(skipped_exons)
    return(skipped_exons_out)
