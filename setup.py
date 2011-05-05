import os, sys, glob

name = 'pySeq'
version = '0.1'

from distutils.core import setup
from Cython.Distutils import build_ext
metadata = {'name':name,
	'version':version,
	'description':'pySeq',
	'author':'Jeff Hsu',
	'packages':['pySeq', 'pySeq.pysam_callbacks', 'pySeq.parsing',
	    'pySeq.stats', 'pySeq.formats'],
}



if __name__ == '__main__':
    dist = setup(**metadata)
