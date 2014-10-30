""" Functions for adding to an axis.  Most functions here require passing in an
axis object from matplotlib
"""

#from collections import defaultdict
import matplotlib.patches as patches
from matplotlib.path import Path
import numpy as np
import pandas as pd

from genda.formats.dosages import grab_gene_location


def should_not_plot(x):
    """
    """
    if x is None:
        return True
    elif isinstance(x, np.ndarray):
        return x.size==0
    elif isinstance(x, pd.DataFrame):
        return True
    else:
        return(bool(x))


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


def plot_transcript(transcript_id, paths, ax, y=0, height=2., 
        full_gene=False):
    """
    """
    xmin = None
    xmax = None
    ps = []
    for i in paths[transcript_id]:
        i.vertices[:, 1] = np.asarray([y, y + height, y + height, y, y])
        if not xmin or not xmax:
            if not xmin or not xmax:
                xmin = np.min(i.vertices[:, 0])
                xmax = np.max(i.vertices[:, 0])
            else:
                if xmin > np.min(i.vertices[:, 0]):
                    xmin = np.min(i.vertices[:, 0])
                elif xmax < np.max(i.vertices[:, 0]):
                    xmax = np.max(i.vertices[:, 0])
                else: pass
        ps.append(patches.PathPatch(i, facecolor='darkgray', lw=0))
    for patch in ps:
        ax.add_patch(patch)
    ax.hlines(y+height/2., xmin, xmax, colors='darkgray', lw=2)


def add_gene_bounderies(ax, gene_annot, gene_name, x_scale):
    """ Add gene bounderies to an axis
    """
    gene_bounds = grab_gene_location(gene_name, cis_buffer=0)
    path_a = make_rectangle(float(gene_bounds[1])/x_scale, 
            float(gene_bounds[2])/x_scale,
            0, 200)
    patch = patches.PathPatch(path_a, facecolor='grey', alpha=0.25)
    return(patch)
