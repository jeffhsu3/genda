

def snp_arrow(pos, height, snpid, ax):
    """ Draws an arrow pointing at a SNP in 
    a Manhattan plot

    pos - an object with pos and height
    height - height
    ax - matplotlib axes

    """

    rel_width = 25

    p_width = ax.get_xlim()
    p_height = ax.get_ylim()
    '''
    spacer = (p_height[1] - p_height[0])/rel_width
    a_head_width = (p_width[1] - p_width[0])/rel_width
    a_head_length = p_height[1] - p_height[0]/rel_width
    arrow_length = p_height - height - spacer
    ax.arrow(pos, p_height[1], 0, 
            arrow_length, 
            head_width = a_head_width,
            head_length = a_head_length, 
            fc = 'k')
    '''

    ax.annotate(snpid, xy=(pos, height),
            xytext = (pos, p_height[1]),
            arrowprops=dict(facecolor='black', shrink=0.05))


    return ax


