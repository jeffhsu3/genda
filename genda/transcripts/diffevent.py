class EventCollection(object):
    """ A collection of diffevents

    Arguments
    ---------
    events - a list of genda.transcripts.DiffEvents
    transcript_ids - a list of transcripts tested
    """
    equivalent_events = {
            'SE' : 'skipped_exon',
            'AS' : 'alternative_start',
            'MXE' : 'mutually_exclusive_exon',
            }

    def __init__(self, events, transcript_ids = None):
        self.events = events
        self.transcript_ids = transcript_ids

    def __getitem__(self, key):
        return(self.events[key])

    def __len__(self):
        return(len(self.events))

    def filter(self, event_type):
        out_events = []
        for i in self.events:
            if i.event_type == event_type or\
            i.event_type == self.equivalent_events[event_type]:
                out_events.append(i)
        return(out_events)

    def collapse(self):
        """
        """
        raise NotImplemented


class DiffEvent(object):
    """A differential splicing event between
    two or more transcripts

    Arguments
    ---------
    event_type : ['skipped_exon', 'mxe', 'A5SE', 'ATS', 'AFE']
    start : start position of differential event
    end : end position of differential event
    transcript_id : 
    chrom : optional
    exon_num : exon number 
    cigar1 :  the cigar tuple generated across an event
    cigar2 : the cigar tuple generated across the other transcript
    event
    """
    def __init__(self, event_type, start, end,
            transcript_ids, cigar1=None, cigar2=None, 
            chrom=None, 
            exon_num=None, exon2 = None, exons = None):
        self.event_type = event_type
        self.start = start
        self.end = end
        # Change to list
        # Make this less chunky
        self.tid = transcript_ids
        self.exon_num = exon_num
        self.chrom = chrom
        self.cigar1 = cigar1
        self.cigar2 = cigar2
        self.exon2 = exon2
        self.exons = exons 

    def __repr__(self):
        return(str(self.start) + '-' +  str(self.end) +\
                ':' + str(self.cigar1) +  ':' + str(self.cigar2))

    def __eq__(self, other):
        if self.start == other.start and\
                (self.end == other.end) and\
                ((self.cigar1 == other.cigar1) |\
                (self.cigar1 == other.cigar2)) and\
                self.cigar2 in [other.cigar1, other.cigar2]:
            return(True)
        else:
            return(False)

    def _extend(self, new_cigar, cig=1):
        """Extends a diffevent for example
        a mutual exclusive exon or two skipped
        exons and other complex events
        """
        old_cig = getattr(self, 'cigar{0!s}'.format(cig)) 
        old_cig.extend(new_cigar)
        setattr(self, 'cigar{0!s}'.format(cig), old_cig)

    def add_transcript(self, transcripts=[]):
        """ Add other transcript pairs that share 
        the same differential event
        """
        raise NotImplementedError


    def cig_to_string(self):
        raise NotImplementedError


    def calc_size_norm(self):
        """ Calculate rough normalization factor
        """
        if self.event_type == 'skipped_exon':
            return()
