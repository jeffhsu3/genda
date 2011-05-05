"""
Callback to retrieve sequences from a bam file over a region
"""


class Seq_call(object):

    def __init__(self):
	self.sequences = []

    def __call__(self, alignment):
	self.sequences.append(alignment.seq)
