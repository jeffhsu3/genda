import os, sys, glob

name = 'pySeq'
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
print(use_cython)

if use_cython:
    ext_modules += [
        Extension("pySeq.transcripts.exon_utils",
                  ["pySeq/transcripts/exon_utils.pyx" ]),
        Extension("pySeq.pysam_callbacks.allele_counter",
            ["pySeq/pysam_callbacks/allele_counter.pyx"]),
    ]
    cmdclass.update({'build_ext': build_ext})
else:
    ext_modules += [
        Extension("pySeq.transcripts.exon_utils",
                  ["pySeq/transcripts/exon_utils.c"]),
        Extension("pySeq.pysam_callbacks.allele_counter",
            ["pySeq/pysam_callbacks.allele_counter.c"]),
    ]

metadata = {'name':name,
            'version': version,
            'cmdclass': cmdclass,
            'ext_modules': ext_modules,
            'description':'pySeq',
            'author':'Jeffrey Hsu',
            'packages':['pySeq', 'pySeq.stats',
                        'pySeq.parsing','pySeq.formats',
                        'pySeq.pysam_callbacks', 
                       'pySeq.transcripts', 'pySeq.AEI'],
}


if __name__ == '__main__':
    dist = setup(**metadata)
