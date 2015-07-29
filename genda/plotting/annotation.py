

def snp_arrow(pos, height, snpid, ax, 
        hadjust=0.0, xadjust=0.0,
        **kwargs):
    """ Draws an arrow pointing at a SNP in 
    a Manhattan plot


    Arguments
    =========
    pos - an object with pos and height
    height - height
    snpid - text label for the arrow
    ax - matplotlib axes
    """

    rel_width = 15
    #p_width = ax.get_xlim()
    p_height = ax.get_ylim()
    spacer = (p_height[1] - p_height[0])/rel_width
    ax.annotate(snpid, xy=(pos, height),
            xytext = (pos + xadjust, p_height[1] - spacer + hadjust),
            arrowprops=dict(facecolor='black', shrink=0.05), 
            **kwargs)
    return ax
