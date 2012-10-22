from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import subprocess, os, pysam

dirname = "/gen_local/hsuj/src/pysam/pysam"
dirname2 = "/gen_local/hsuj/src/pysam/samtools"


def make_ext(modname, pyxfilename):
    import pysam, os
    dirname = os.path.dirname( pysam.__file__)[:-len("pysam")]
    return Extension(name = modname, sources = [pyxfilename],
                    extra_link_args = [os.path.join(dirname, "csamtools.so")],
                    include_dirs = [dirname,
                                    dirname2,
                                    ])


ext_modules = [make_ext("bleh", "test_embed.pyx")]

setup(
        name = 'Bleh',
        cmdclass = {'build_ext': build_ext},
        ext_modules = ext_modules,
        )
"""
dirname = os.path.dirname(pysam.__file__)[:-len("pysam")]
dirname = "/gen_local/hsuj/src/pysam/pysam"
dirname2 = "/gen_local/hsuj/src/pysam/samtools"
print(dirname)

subprocess.call(["ls", "%s" % dirname])

subprocess.call(["cython", "-I %s" % dirname, "--embed",
                "-I %s" % dirname2, "test_embed.pyx"])
"""
