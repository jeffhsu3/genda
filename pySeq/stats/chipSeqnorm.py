"""
Norms DNAse data from a BAM file. 

Citation:

Boyle, A.P.; Song, L.; Lee High-resolution genome-wide in vivo footprinting of
diverse transcription factors in human cells. Genome Research. 2011
"""

import numpy as np

def window_norm(y, window_size):
    """
    Normalizes read counts by averaging reads/bp across the window size and
    dividing. 
    """

    # Handles the edge cases
    pass

def savitzky_golay(y, window_size, order, deriv=0):
    """
    From scipy recipies (www.scipy.org/Cookbook/SavitzkyGolay 
    array : array of counts at a position

    Citations:

    """
    try: 
	window_size = np.abs(np.int(window_size))
	order = np.abs(np.int(order))
    except ValueError, msg:
	raise ValueError("window_size and order have to be type int")
    if window_size % 2 != 1 or window_size < 1:
	raise TypeError("window_size must be a positive odd number")
    if window_size < order + 2:
	raise TypeError("window_size is too small for polynomials order")
    order_range = range(order +1)
    half_window = (window_size -1)//2
    b = np.mat([[k**i for i in order_range] for k in range(-half_window,
	half_window+1)])
    m = np.linalg.pinv(b).A[deriv]
    # padding 
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1]-y[0])
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1]-y[-1]) 
    y = np.concatenate((firstvals, y, lastvals))
    return(np.convolve(m, y, mode='valid'))


