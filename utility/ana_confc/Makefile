# first delete all default suffixes
.SUFFIXES:
# then add your suffixes
.SUFFIXES: .f90 .F90 .o

include ./make.sys

default: all

modu1 = 
modu2 = 
dqmcf1 = generate_neighbor.o
main = main.o

objects = $(modu1) $(modu2) $(dqmcf1) $(main)

#LIBS+= ../common/libMM.a

all: cpmod x_confc

x_confc: $(objects)
	$(LINKER) $(objects) -o x_confc $(LFLAGS) $(LIBS)

cpmod:
	cp ../../common/*.mod .

.f90.o:
	$(F90) $(FFLAGS) $*.f90

.F90.o:
	$(F90) $(FFLAGS) $*.F90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f x_confc

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
