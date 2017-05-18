#from pymc import Categorical
#from pymc import MCMC
from bx.intervals.intersection import Intersecter, Interval
from numpy import zeros
from numpy import repeat as nrepeat
from numpy import logical_and
from numpy import array
from numpy import sum as sum_

from collections import defaultdict


def msample_count_reads():
    """
    """


def count_reads(transcripts, bam_iter, number_of_counts=1):
    """ Count the reads in a given transcript

    :TODO rename 
    :TODO change to cython
    Arguments
    ---------
    transcripts : list
        list of exons
    bam_iter : pysam.BamFileIterator
        gotton after pysam.BamFile.fetch() call
    """
    # Convert this to Cython
    out_counts = zeros(len(transcripts))
    intron_lengths = []
    read_vector = []
    tree = Intersecter()
    # Assume exons are position sorted
    for ti, transcript in enumerate(transcripts):
        ex_list = []
        for j, i in enumerate(transcript):
            tree.add_interval(Interval(int(i[0]), int(i[1]), 
                value={'anno':ti}))
            if j != 0:
                ex_list.append(transcript[j-1][1]\
                        - transcript[j][0])
        intron_lengths.append(ex_list)
    for read in bam_iter:
        block_counter = zeros((len(transcripts),))
        intron_match =  zeros((len(transcripts),))
        blocks = read.get_blocks()
        junction_lengths = []
        for i,j in enumerate(blocks):
            if i != 0:
                junction_lengths.append(blocks[i - 1][1] - j[0])
            else: pass
        junction_lengths = set(junction_lengths)
        for i, k in enumerate(blocks):
            overlap = tree.find(k[0], k[1]) 
            if len(overlap) == 0:
                break
            else:
                for s in overlap:
                    if (k[0] >= s.start) and\
                            (k[1] <= s.end):
                        block_counter[s.value['anno']] += 1
        for ij, il in enumerate(intron_lengths):
            if set(junction_lengths).issubset(set(il)):
                intron_match[ij] = 1
            else: pass
        smatch = nrepeat(len(blocks), len(transcripts))
        gg = logical_and(block_counter == smatch, intron_match)
        read_vector.append(gg)
        out_counts += gg
    read_matrix = array(read_vector)
    uniq_r = sum_(read_matrix, axis=1) == 1
    #normalization_constant = [for i in transcripts]
    return(out_counts)


def bias_estimation():
    """Estimate 5' and 3' bias from unambigously assigned transcript reads
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
