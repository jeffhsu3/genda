from collections import defaultdict


class Gene(object):
    """ A collection of transcripts

    Arguments
    ---------
    gene - a unique identifier for the gene
    chrom - chromosome identifier
    start - start position of the gene
    end - end position of the gene 
    transcripts - a transcript class
    """

    def __init__(self, gene, chrom = None, start=0, end=3e8, 
            transcripts=[], strand=None, 
            symbol=None):
        self.gene = gene
        self.chrom = str(chrom)
        self.start = start
        self.end = end
        self.transcripts = []
        self.symbol = symbol

    def get_transcripts(self, gtf_tabix, 
            buffer = 500000):
        """
        Generates transcripts from gtf_tabix.  Replaces
        self.transcripts.

        Arguments
        ---------
        gtf_tabix - A pysam.Tabixfile gtf
        """
        try:
            pls = gtf_tabix.fetch(self.chrom, self.start - buffer,
                    self.end + buffer)
        except ValueError:
            pls = gtf_tabix.fetch(self.chrom, 0,
                    self.end + buffer)
        path_dict = defaultdict(list)
        # Convert to numba or cython
        for i in pls:
            i = i.split("\t")
            geni = i[9].split(";")
            gene_id = geni[0].lstrip(' gene_id "').rstrip('"')
            if gene_id == self.gene:
                if i[7] == 'exon':
                    try:
                        itrans = geni[1].lstrip(' transcript_id "').rstrip('"')
                        exon_number = (geni[2].
                                lstrip(' exon_number "').rstrip('"'))
                        path_dict[itrans].append([int(i[1]), int(i[2]),
                            int(exon_number)])
                    except KeyError:
                        pass
                else: pass
        self.transcripts = path_dict

    def filter_transcripts(self, transcripts):
        ''' Filter transcripts that are not in transcripts 

        Arguments
        ---------
        transcripts : an iterator of transcripts that are to be kept
        '''
        if len(self.transcripts) == 0:
            self.transcripts = []
        else:
            for key in self.transcripts.keys():
                if key in transcripts:
                    pass
                else:
                    del self.transcripts[key]

    def get_common_parts(self):
        raise NotImplementedError


