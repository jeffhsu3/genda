import Genotype
import pandas as pd
import numpy as np

class PED(Genotype.Genotype):

    def __init__(self, PED, MAP, encoder = None):
        self.PED = pd.read_table(PED, delim_whitespace = True, header = None)
        self.MAP = pd.read_table(MAP, delim_whitespace = True, header = None)
        if encoder != None:
            self.geno = self.genotype(self.PED, self.MAP, encoder)
        else:
            self.geno = None

    def genotype(self, PED, MAP, encoder):
        geno = PED.ix[:, 6:]
        snps = []
        for s in MAP.ix[:,1]:
            snps.extend([s]*2)
        geno.columns = snps
        geno = geno.apply(self.apply_encoder, encoder = encoder)
        geno = pd.DataFrame([[geno.ix[r,c*2] + geno.ix[r,(c*2)+1]\
                for c in range(geno.shape[1]/2)] for r in range(geno.shape[0])])
        geno = geno.transpose()
        geno.index = MAP.ix[:,1]
        geno.columns = PED.ix[:,1]
        return geno
        
    def apply_encoder(self, PED, encoder):
        snp = encoder[PED.name]
        encoding_dict = {snp[0] : 0, snp[2] : 1, '0': np.nan}
        return PED.map(encoding_dict)
