# Makefile
PYTHON=python
PYVERSION=$(shell $(PYTHON) -c "import sys; print(sys.version[:3])")
PYPREFIX=/usr
INCLUDES=-I/$(PYPREFIX)/include/python$(PYVERSION)

INCDIR=$(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_python_inc())")
LIBDIR1=$(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
LIBDIR2=$(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_config_var('LIBPL'))")

LINKCC=$(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_config_var('LINKCC'))")
LINKFORSHARED=$(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_config_var('LINKFORSHARED'))")
LIBS=$(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_config_var('LIBS'))")
SYSLIB=$(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_config_var('SYSLIBS'))")


SAMTOOLS=/gen_local/hsuj/src/pysam/samtools/
PYSAM=/gen_local/hsuj/src/pysam/pysam/

embedded: embedded.o
	$(LINKCC) -o $@ $^ -L$(LIBDIR1) -L$(LIBDIR2) -lpython$(PYVERSION) $(LIBS) $(SYSLIB) $(LINKFORSHARED)

embedded.o: embedded.c
	gcc -c $^ -I$(SAMTOOLS) -I$(PYSAM) -I$(INCDIR)

embedded.c: embedded.pyx
	cython --embed -I$(PYSAM) embedded.pyx 

all: embedded

clean:
	@echo Cleaning Demos/embeded
	@rm -f *~ *.o *.so core core.* *.c
