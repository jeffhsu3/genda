

def snp_arrow(pos, height, snpid, ax):
    """ Draws an arrow pointing at a SNP in 
    a Manhattan plot

    pos - an object with pos and height
    height - height
    snpid - text label for the arrow
    ax - matplotlib axes

    """

    rel_width = 15
    p_width = ax.get_xlim()
    p_height = ax.get_ylim()
    spacer = (p_height[1] - p_height[0])/rel_width
    ax.annotate(snpid, xy=(pos, height),
            xytext = (pos, p_height[1] - spacer),
            arrowprops=dict(facecolor='black', shrink=0.05))


    return ax


