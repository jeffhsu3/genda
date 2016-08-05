#from pymc import Categorical
#from pymc import MCMC
from bx.intervals.intersection import Intersecter, Interval


def count_reads(transcripts, bam_iter):
    ''' Count the reads in a given transcript
    Arguments
    ---------
    transcripts : list
        list of exons
    bam_iter : pysam.BamFileIterator
        gotton after fetching from a pysam.BamFile
    '''
    # Convert this to Cython
    tree = Intersecter()
    intron_lengths = []
    read_vector = []
    # :TODO
    for j, i in enumerate(transcripts):
        tree.add_interval(Interval(int(i[0]), int(i[1]), 
            value={'anno':i[2]}))
        if j != 0:
            intron_lengths.append(transcripts[j-1][1] - transcripts[j][0])
    counter = 0
    for read in bam_iter:
        blocks = read.get_blocks()
        junction_lengths = []
        for i,j in enumerate(blocks):
            if i != 0:
                junction_lengths.append(blocks[i - 1][1] - j[0])
            else: pass
        tc = 0
        for k in blocks:
            overlap = tree.find(k[0], k[1]) 
            if len(overlap) == 0:
                break
            else:
                if (k[0] >= overlap[0].start) and\
                        (k[1] <= overlap[0].end):
                    tc += 1
                    # :TODO check length
        junctions = [True if i in intron_lengths else False  for i in junction_lengths]
        if tc == len(blocks) and all(junctions):
            counter += 1
            read_vector.append(True)
        else:
            read_vector.append(False)
    return(counter)


def bias_estimation():
    """
    """
    pass


def precalculate_lr(transcripts):
    """ Precalculate P(p|r).  Ignores Phred score at the moment
    """
    pass


def prob_m_read(read):
    """
    """
    transcript_length = 2
    prob_position_given_m = 1/(read.query_length + transcript_length)
    return prob_position_given_m * prob_read_sequence



def prob_m_pair_read():
    """
    """
    pass
