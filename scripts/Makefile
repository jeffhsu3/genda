SHELL = /bin/sh
.SUFFIXES: .c .o
# Makefile
program_PYX_SRC:= $(wildcard *.pyx)
PYTHON=python
PYVERSION=$(shell $(PYTHON) -c "import sys; print(sys.version[:3])")
PYPREFIX=/usr
INCLUDES=-I/$(PYPREFIX)/include/python$(PYVERSION)
#INCLUDES=-I/$(PYPREFIX)/include/python$(PYVERSION)

INCDIR=$(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_python_inc())")
LIBDIR1=$(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
LIBDIR2=$(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_config_var('LIBPL'))")

LINKCC=$(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_config_var('LINKCC'))")
LINKFORSHARED=$(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_config_var('LINKFORSHARED'))")
LIBS=$(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_config_var('LIBS'))")
SYSLIB=$(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_config_var('SYSLIBS'))")
PYSAM=$(shell $(PYTHON) -c "import pysam; print(pysam.get_include()[0])")
SAMTOOLS=$(shell $(PYTHON) -c "import pysam; print(pysam.get_include()[1])")


# How to correctly find the source files for pysam?
# SAMTOOLS=/gen_local/hsuj/src/pysam/samtools/
# PYSAM=/gen_local/hsuj/src/pysam/pysam/

gene_counter: gene_counter.o
	$(LINKCC) -o $@ $^ -L$(LIBDIR1) -L$(LIBDIR2) -lpython$(PYVERSION) $(LIBS) $(SYSLIB) $(LINKFORSHARED)

gene_counter.o: gene_counter.c
	gcc -c $^ -I$(SAMTOOLS) -I$(PYSAM) -I$(INCDIR)

gene_counter.c: gene_counter.pyx
	cython --embed -I$(PYSAM) gene_counter.pyx 

embedded: embedded.o
	$(LINKCC) -o $@ $^ -L$(LIBDIR1) -L$(LIBDIR2) -lpython$(PYVERSION) $(LIBS) $(SYSLIB) $(LINKFORSHARED)

embedded.o: embedded.c
	gcc -c $^ -I$(SAMTOOLS) -I$(PYSAM) -I$(INCDIR)

embedded.c: embedded.pyx
	cython --embed -I$(PYSAM) embedded.pyx 

htseqcount: htseqcount.o
	$(LINKCC) -o $@ $^ -L$(LIBDIR1) -L$(LIBDIR2) -lpython$(PYVERSION) $(LIBS) $(SYSLIB) $(LINKFORSHARED)

htseqcount.o: htseqcount.c
	gcc -c $^ -I$(SAMTOOLS) -I$(PYSAM) -I$(INCDIR)

htseqcount.c: htseqcount.pyx
	cython --embed -I$(PYSAM) htseqcount.pyx

all: embedded gene_counter htseqcount

clean:
	@echo Cleaning Demos/embeded
	@rm -f *~ *.o *.so core core.* *.c
