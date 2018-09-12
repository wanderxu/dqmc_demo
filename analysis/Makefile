FC= ifort -O3
#FLAGS= -c  -w   -r8   -O3   
FLAGS= -c  -w   -r8    -check bounds
#LF= -pg
LIBS= ../common/libMM.a \
       -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm

all:
	cp ../common/*.mod .;\
	(make -f Compile_tau  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) ;\
	(make -f Compile_tau_flux  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) ;\
	(make -f Compile_tau_flux_l4  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) ;\
	(make -f Compile_tau_flux_l8  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) ;\
	
clean: 	
	(make -f Compile_tau  clean );\
	(make -f Compile_tau_flux  clean );\
	(make -f Compile_tau_flux_l4  clean );\
	rm *.mod *~ \#* *.out
