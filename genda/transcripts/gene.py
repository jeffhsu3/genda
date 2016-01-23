from collections import defaultdict


class Gene(object):
    """ A collection of transcripts

    """

    def __init__(self, gene, chrom = None, start=0, end=3e8, 
            transcripts=[], strand=None, 
            symbol=None):
        self.gene = gene
        self.chrom = chrom.encode('utf8')
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
        pls = gtf_tabix.fetch(self.chrom, self.start - buffer,
                self.end + buffer)
        path_dict = defaultdict(list)
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
