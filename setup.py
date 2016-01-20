import os, sys, glob
import numpy
import pysam

name = 'genda'
version = '0.1'

from setuptools import setup
from distutils.extension import Extension

try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
    print("Cython not found")
else:
    use_cython = True

cmdclass = {}
ext_modules = []
includes = pysam.get_include()
includes.append(numpy.get_include())


if use_cython:
    print('**********using cython*********')
    ext_modules += [
        Extension("genda.transcripts.exon_utils",
                  ["genda/transcripts/exon_utils.pyx" ], 
                  include_dirs=includes),

        Extension("genda.pysam_callbacks.allele_counter",
            ["genda/pysam_callbacks/allele_counter.pyx"],
            include_dirs=includes,
            define_macros=pysam.get_defines()),

        Extension("genda.stats.aei_count_samples",
            ["genda/stats/aei_count_samples.pyx"],
            include_dirs=includes,
            define_macros=pysam.get_defines()),
    ]
    print(ext_modules)
    cmdclass.update({'build_ext': build_ext})

else:
    ext_modules += [
        Extension("genda.transcripts.exon_utils",
                  ["genda/transcripts/exon_utils.c"]),
        Extension("genda.pysam_callbacks.allele_counter",
            ["genda/pysam_callbacks/allele_counter.c"]),
        Extension("genda.pysam_callbacks.gene_counter",
            ["genda/stats/aei_count_samples.c"]),
    ]


metadata = {'name':name,
            'version': version,
            'cmdclass': cmdclass,
            'ext_modules': ext_modules,
            'scripts': glob.glob('scripts/*.py'),
            'description':'genda',
            'author':'Jeffrey Hsu',
            'packages':['genda', 
                        'genda.stats',
                        'genda.parsing',
                        'genda.formats',
                        'genda.pysam_callbacks',
                        'genda.transcripts', 
                        'genda.AEI', 
                        'genda.plotting',
                        'genda.eQTL'],
}


if __name__ == '__main__':
    dist = setup(**metadata)
    """
        Extension("genda.pysam_callbacks.gene_counter",
            ["genda/pysam_callbacks/gene_counter.pyx"],
            include_dirs=pysam.get_include()),
    """
