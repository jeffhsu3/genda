""" Functions for adding to an axis.  Most functions here require passing in an
axis object from matplotlib
"""
from collections import defaultdict
import re
import matplotlib.patches as patches
from matplotlib.path import Path
import numpy as np
import pandas as pd

from genda.formats import grab_gene_location


def create_path(gtf_iterator, gene, height=0.5, x_scale=1):
    """ Create a path for a given transcript

    :TODO change name
    """
    transcripts = defaultdict(list)
    t_rxe = re.compile('transcript_id "(ENST[0-9]+)"')
    t_add = 0
    for i in gtf_iterator:
        mverts = []
        i = i.split("\t")
        if (i[3] == gene or gene in i[9]) and i[7] == 'exon':
            verts = [(int(i[1]) / x_scale, t_add),
                     (int(i[1]) / x_scale, height + t_add),
                     (int(i[2]) / x_scale, height + t_add),
                     (int(i[2]) / x_scale, t_add),
                     (int(i[1]) / x_scale, t_add),
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
            else:
                pass
            transcripts[transcript].append(Path(verts, codes))
        else:
            pass
    return(transcripts)


def should_not_plot(x):
    """Returns true if not None
    x : pandas.DataFrame or numpy array
    """
    if x is None:
        return True
    elif isinstance(x, np.ndarray):
        return x.size == 0
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


def plot_box_(transcript, ax):

    path = make_rectangle()
    return(path)


def plot_transcript(transcript_id, paths, ax, y=0, height=2.,
                    full_gene=False):
    """ Plot a transcript
    """
    xmin = None
    xmax = None
    ps = []
    ymin, ymax = ax.get_ylim()
    try:
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
                    else:
                        pass
            ps.append(patches.PathPatch(i, facecolor='darkgray', lw=0))
        for patch in ps:
            ax.add_patch(patch)
    except AttributeError:
        for i in paths[transcript_id]:
            try:
                path = make_rectangle(i[0], i[1], y, y + height)
            except IndexError:
                path = make_rectangle(i.start, i.end, y, y + height)
            patch = patches.PathPatch(path, facecolor='grey', alpha=0.6)
            ax.add_patch(patch)
    # hlines doesn't in this function but outside it?
    #draw_arrows(xmin, xmax, y+height/2.0, ax)
    #ax.hlines(y+height/2.0, xmin, xmax, colors='darkgray', lw=2)
    return(ax)


def add_gene_bounderies(ax, gene_annot, gene_name, x_scale):
    """ Add gene bounderies to an axis
    """
    gene_bounds = grab_gene_location(gene_name,
                                     ann_file=gene_annot, cis_buffer=0)
    path_a = make_rectangle(float(gene_bounds[1]) / x_scale,
                            float(gene_bounds[2]) / x_scale,
                            0, 200)
    patch = patches.PathPatch(path_a, facecolor='grey', alpha=0.25)
    return(patch)


def add_snp_arrow(adj_pv, x, snp, ax):
    """Adds an arrow with snp text to plot 
    :TODO The best snp should be some sort of object
    """
    print(type(ax.get_ylim()))
    print(type(adj_pv))
    arrow_start = ax.get_ylim()[1] - adj_pv / 6.0 / 2
    arrow_length = adj_pv / 6.0 / 2 - adj_pv / 6.0 / 2 / 8
    ax.arrow(x, arrow_start, 0, -arrow_length, head_width=0.01,
             head_length=0.1, fc='k', ec='k')
    ax.text(x - 0.05, arrow_start + adj_pv / 6.0 / 5.0,
            snp, style='italic')
    return(ax)


def plot_read2(read, y_start, height=1, ymax=30, scale=1e6):
    """
    """
    gout = {'patches': [],
            'hlines': [],
            }
    cont = 0
    for i, j in read.cigar:
        if i == 0:
            path = make_rectangle(read.pos / scale + cont / scale,
                                  read.pos / scale + cont / scale + j / scale,
                                  ymax - y_start,
                                  ymax - y_start - float(height))
            patch = patches.PathPatch(path, facecolor='grey',
                                      alpha=0.4)
            gout['patches'].append(patch)
            cont += j
        elif i == 3:
            gout['hlines'].append((ymax - y_start - float(height) / 2.0,
                                   read.pos / scale + cont / scale,
                                   read.pos / scale + cont / scale + j / scale))
    return(gout)


def plot_read(read, y_start, ax, height=1, ymax=30, scale=1e6):
    """ Plot a sequencing read

    """
    cont = 0
    for i, j in read.cigar:
        if i == 0:
            path = make_rectangle(read.pos / scale + cont / scale,
                                  read.pos / scale + cont / scale + j / scale,
                                  ymax - y_start,
                                  ymax - y_start - float(height))
            patch = patches.PathPatch(path, facecolor='grey',
                                      alpha=0.4)
            ax.add_patch(patch)
            cont += j
        elif i == 3:
            ax.hlines(ymax - y_start - float(height) / 2.0,
                      read.pos / scale + cont / scale,
                      read.pos / scale + cont / scale + j / scale,
                      colors='grey')
            cont += j


def coverage_hist(read, hist_array, start):
    """ updates hist_array with coverage counts at a given position
    # Move to cython
    """
    ii = read.pos - start
    cont = 0
    for i, j in read.cigar:
        if i == 0:
            if ii + cont > 0:
                hist_array[ii + cont:ii + cont + j] += 1
            elif ii + cont + j > 0:
                hist_array[0: ii + cont + j] += 1
            else:
                pass
            cont += j
        elif i == 3:
            # cigar is always ordered
            cont += j


def create_gene_path(gtf_iter, gene, x_scale):
    """
    Returns
    -------
    None
    """
    raise NotImplementedError


def draw_arrows(xmin, xmax, y, ax, spacing=None):
    """
    xmin : int or float
        start of where to draw arrows
    xmax : int or float
        where to stop drawing arrows
    y : y position
        y position where to center arrows
    """
    ylim = ax.get_ylim()
    xlim = ax.get_xlim()
    # :TODO make sensible defaults here
    height = (ylim[1] - ylim[0])/5.0
    length = (xlim[1] - xlim[0])/25.0
    
    if spacing:
        pass
    else:
        spacing = (xlim[1] - xlim[0])/10.0
    x_range = np.arange(xmin, xmax, spacing)[1:-1]
    for i in x_range:
        draw_arrow(i, y, ax, x_adjust=length, y_adjust=height)


def draw_arrow(x, y, ax, x_adjust=0.05, y_adjust=0.05):
    # :TODO x_adjust and y_adjust scale based on axis limits
    ax.plot([x, x + x_adjust], [y, y + y_adjust], 'k-', color='darkgrey')
    ax.plot([x, x + x_adjust], [y, y - y_adjust], 'k-', color='darkgrey')


def get_path_max_and_min(gene_dict):
    """ Return min, max of a path
    """
    points = []
    try:
        for i in gene_dict.values():
            for j in i:
                for k in j.vertices:
                    points.append(k[0])
    except AttributeError:
        # Case where only 1 transcript is chosen
        for j in gene_dict:
            for k in j.vertices:
                points.append(k[0])
    return(min(points), max(points))


def plot_transcript_(transcript, ax, y=0, height=2.,
                     intron_scale=0.10, scale=1):
    '''
    Arguments
    ---------  
    transcript - a list of exons (start, end, exon_number) or
    a genda.transcripts.Transcript object
    ax - matplotlib.axes
    y - where to place the base of the transcript on the y-axis
    height - height of the boxes of the transcript
    ax - matplotlib axes object

    :TODO transcript could be an object.  
    :TODO write tests for this

    Returns
    ------- 
    matplotlib.axis
    new_scaled_ylim and xlim
    '''
    xmin, xmax = 1e12, 0
    ymin, ymax = ax.get_ylim()
    beg_exon = transcript[0]
    path = make_rectangle(beg_exon[0], beg_exon[1], y, y + height)
    patch = patches.PathPatch(path, facecolor='grey', alpha=0.6)
    ax.add_patch(patch)
    rx = transcript[0][1]
    for i, exon in enumerate(transcript[1:]):
        intron_l = (exon[0] - transcript[i - 1][1]) * intron_scale
        exon_length = exon[1] - exon[0]
        ns = rx + intron_l
        ne = rx + exon_length
        path = make_rectangle(ns, ne, y, y + height)
        patch = patches.PathPatch(path, facecolor='grey', alpha=0.6)
        ax.add_patch(patch)
        rx = rx + intron_l + exon_length
        if ns < xmin:
            xmin = ns
        if ne > xmax:
            xmax = ne
    return(ax, (beg_exon[0], xmax))


def draw_arc_label(start, end, ax, color = 'black',
        mirror=False):
    """ dra 
    Primarily used
    """
    codes = [Path.MOVETO,
             Path.CURVE4,
             Path.CURVE4,
             Path.CURVE4, ]
    if mirror:
        P1 = (end[0], start[1])
    else:
        P1 = (start[0], end[1])
    P2 = end
    verts = [start,
             P1,
             P2,
             end]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, lw=2, ec=color, facecolor='none', 
            linestyle='--')
    ax.add_patch(patch)
    return(ax)


def draw_junction_arcs(to_plot, hist,  xmin, ax, color='black', 
        text='hmm', y_buffer=0):
    """
    :TODO refactor so that to_plot is a diff event
    to_plot - tuple
    """
    x1 = to_plot[0][1]
    x2 = to_plot[1][0]
    start = (x1, hist[x1 - xmin - 1])
    xmid = (to_plot[1][0] - to_plot[0][1])/2 + to_plot[0][1]
    ymid =  max(hist[x1 - xmin - 1], hist[x2 - xmin]) + y_buffer
    ax = draw_arc_label(start, 
            (xmid, ymid),
            ax, color=color)
    ax = draw_arc_label((xmid, ymid),
            (x2, hist[x2-xmin]),
            ax, color = color, mirror = True)
    ax.text(xmid, ymid, '{0:.2g}'.format(text), size='x-large')
    return(ax)


def add_size_legends(ax, size_base=200):
    sizes = [((size_base * i) + 20) for i in np.linspace(0, 0.5, num=4)]
    l1 = ax.scatter([], [], s=sizes[0])
    l2 = ax.scatter([], [], s=sizes[1])
    l3 = ax.scatter([], [], s=sizes[2])
    l4 = ax.scatter([], [], s=sizes[3])
    leg = ax.legend(
        [l1, l2, l3, l4],
        labels=np.around(np.linspace(0, 0.5, num=4), decimals=2),
        fontsize=14, loc=0,
        frameon=True, title='maf', scatterpoints=1)
    return(ax)
