
import os, sys, glob
import pysam

name = 'genda'
version = '0.1'

from distutils.core import setup
from distutils.extension import Extension

def make_ext(modname, pyxfilename):
    import pysam, os
    dirname = os.path.dirname( pysam.__file__)[:-len("pysam")]
    return Extension(name = modname, sources = [pyxfilename],
                    extra_link_args = [os.pth.join(dirname, "csamtools.so")],
                    include_dirs = ['../samtools',
                                    '../pysam',
                                    ])

try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
    print("Cython not found")
else:
    use_cython = True

cmdclass = {}
ext_modules = []

if use_cython:
    print('using cython')
    ext_modules += [
        Extension("genda.transcripts.exon_utils",
                  ["genda/transcripts/exon_utils.pyx" ],),
        Extension("genda.pysam_callbacks.allele_counter",
            ["genda/pysam_callbacks/allele_counter.pyx"],
            include_dirs=pysam.get_include(),
            define_macros=pysam.get_defines()),
        Extension("genda.pysam_callbacks.gene_counter",
            ["genda/pysam_callbacks/gene_counter.pyx"],
            include_dirs=pysam.get_include()),
    ]
    print(ext_modules)
    cmdclass.update({'build_ext': build_ext})
else:
    ext_modules += [
        Extension("genda.transcripts.exon_utils",
                  ["genda/transcripts/exon_utils.c"]),
        Extension("genda.pysam_callbacks.allele_counter",
            ["genda/pysam_callbacks.allele_counter.c"]),
        Extension("genda.pysam_callbacks.gene_counter",
            ["genda/pysam_callbacks.allele_counter.c"]),
    ]


metadata = {'name':name,
            'version': version,
            'cmdclass': cmdclass,
            'ext_modules': ext_modules,
            'scripts': glob.glob('scripts/*.py'),
            'description':'genda',
            'author':'Jeffrey Hsu',
            'packages':['genda', 'genda.stats',
                        'genda.parsing','genda.formats',
                        'genda.pysam_callbacks',
                       'genda.transcripts', 'genda.AEI'],
}


if __name__ == '__main__':
    dist = setup(**metadata)
    """
        Extension("genda.pysam_callbacks.gene_counter",
            ["genda/pysam_callbacks/gene_counter.pyx"],
            include_dirs=pysam.get_include()),
    """
