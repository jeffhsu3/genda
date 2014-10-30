from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import re
#from scipy.stats import linregress, pearsonr
#from matplotlib import figure
import pandas as pd
import numpy as np
import pysam

from matplotlib.path import Path
from genda.formats.dosages import grab_gene_location
from .plotting_utils import *
#from genda.plotting.plotting_utils import create_path


def make_rectangle(start, end, y1, y2):
    verts = [(start, y1),
            (start, y2),
            (end, y2),
            (end, y1),
            (start, y1),
            ]
    codes = [
            Path.MOVETO,
            Path.LINETO,
            Path.LINETO,
            Path.LINETO,
            Path.CLOSEPOLY,
            ]
    return (Path(verts, codes))


def create_path(gtf_iterator, gene, height=0.5, x_scale=1):
    """ Create a path for a given transcript
    """
    transcripts = defaultdict(list)
    t_rxe = re.compile('transcript_id "(ENST[0-9]+)"')
    t_add = 0
    for i in gtf_iterator:
        mverts = []
        i = i.split("\t")
        if (i[3] == gene or gene in i[9]) and i[7] == 'exon':
            verts = [(int(i[1])/x_scale, t_add),
                     (int(i[1])/x_scale, height + t_add),
                     (int(i[2])/x_scale, height + t_add),
                     (int(i[2])/x_scale, t_add),
                     (int(i[1])/x_scale, t_add),
                    ] 
           
	    codes = [
		    Path.MOVETO,
		    Path.LINETO,
		    Path.LINETO,
		    Path.LINETO,
		    Path.CLOSEPOLY,
		    ]
	    transcript = t_rxe.search(i[9]).group(1)
	    if transcript not in transcripts.keys():
	        t_add += 1
	    else: pass
	    transcripts[transcript].append(Path(verts, codes))
        else: pass
    return(transcripts)
