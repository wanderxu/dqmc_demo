.SUFFIXES: .f90

#include ../build/make.sys
include ./make.sys

mod = m_constants.o m_spring.o m_variance.o mkl_dfti.o
sub = s_fft.o s_matrix.o s_util.o
objects = $(mod) $(sub)

default: all

all: build-mod build-sub build-lib

build-mod: $(mod)
build-sub: $(sub)
build-lib: $(objects)
	$(ARCHIVER) libMM.a $(objects)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f libMM.a

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
