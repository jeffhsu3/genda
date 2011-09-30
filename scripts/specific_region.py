import sys, pysam
import numpy as np
import rpy2.robjects as ro
import rpy2.robjects.lib.ggplot2 as ggplot2
from rpy2.robjects.packages import importr
from pySeq.common_funcs import regionParse 
from rpy2.robjects.numpy2ri import numpy2ri
import matplotlib.pyplot as plt
from pySeq.stats.chipSeqnorm import savitzky_golay

def main():
    import optparse
    p = optparse.OptionParser(__doc__)
    p.add_option("-D", "--debug", action="store_true", dest="D", help="debug")
    p.add_option("-S", "--stam", action="store_true", dest="S", help="DNAseI\
	    is generated from STAM's group")
    p.add_option("-s", "--shift", action="store", dest="shift", help="Amount to\
	    shift the negative strand", default = 36)

    
    ro.conversion.py2ri = numpy2ri
    options, args = p.parse_args()
    options.shift = int(options.shift)
    bamfile = pysam.Samfile(args[0], 'rb')
    region_str = regionParse(args[1])
    chrom = region_str[0]
    start = region_str[1]
    end = region_str[2]
    diff = end-start
    a = np.zeros(2*(diff), dtype=np.int)
    try:
	for alignment in bamfile.fetch(chrom, start, end):
	    if alignment.pos-start < 0: pass
	    else:
		if alignment.is_reverse:
		    try:
			a[alignment.pos-start+diff+options.shift] += 1
		    except IndexError:
			pass
		else:
		    a[alignment.pos-start] += 1
    except ValueError:
	pass

    #########################################
    # R plotting set-up
    #########################################
    r = ro.r
    graphics = importr('graphics')
    grdevices = importr('grDevices')
    grdevices.png(file="/home/hsuj/projects/dnaseq/Stam/test.png", width=1000,
            height=1000)
    """
    base = importr('base')
    datasets = importr('datasets')
    mtcars = datasets.mtcars
    pp = ggplot2.ggplot(mtcars) + \
            ggplot2.aes_string(x='wt', y='mpg', col='factor(cyl)') + \
            ggplot2.geom_point() + \
            ggplot2.geom_smooth(ggplot2.aes_string(group='cyl'),
                    method='lm')
    pp.plot()
    """	    
    
    graphics.plot(a[1:diff], main="test", pch=19, cex=0.5)
    graphics.points(a[diff+1:2*diff], col="red", pch=19, cex=0.5)
    grdevices.dev_off() 
    x = np.arange(1, diff, 1)
    data2 =savitzky_golay(a, 9, 2)
    plt.plot(x, data2[1:diff], 'b-',x, data2[diff+1:2*diff], 'r-')
    plt.show()


if __name__ == '__main__':
    main()

