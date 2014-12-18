"""
Tools for analyzing allelic expression imbalance (AEI) from RNA-sequencing
==========================================================================

=========
AEI funcs
=========
"""

from .AEI_plot import *
from .AEI import *



def get_aei(aei_path, dosage_matrix):
    """
    """
    aei = pd.read_pickle(aei_path)
    new_columns = [i[0].split("\t")]
    return(aei)
