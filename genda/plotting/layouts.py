""" Common layouts
"""
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import numpy as np


def two_eQTL_layout():
    """ Two eQTL layouts and a axes from genes/annotations
    """
    fig, ax = plt.subplots(nrows=3, figsize=(20, 14),
            sharex=True,)
    gs = gridspec.GridSpec(5, 3)
    ax[0] = plt.subplot(gs[0:2, :])
    ax[1] = plt.subplot(gs[2:4, :])
    ax[2] = plt.subplot(gs[4:, :])
    ss = np.linspace(0.1, 0.5, num=3)
    sizes = [((200 * i ) + 20) for i in ss]
    l1 = ax[0].scatter([],[], s=sizes[0], c='b')
    l2 = ax[0].scatter([],[], s=sizes[1], c='b')
    l3 = ax[0].scatter([],[], s=sizes[2], c='b')
    leg = ax[0].legend([l1, l2, l3], 
        labels = np.around(ss, decimals=2),
        fontsize=14, loc=0, 
        frameon=True, title='maf', scatterpoints=1
        )
    return(fig, ax)
