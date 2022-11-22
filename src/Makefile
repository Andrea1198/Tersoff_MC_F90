F90=gfortran
FFLAGS=-O1 -fopenmp
.SUFFIXES= .o .f90
src=calcener.o\
	eqmc.o\
	library.o\
	potential.o\
	system.o

eqmc: calcener.o
calcener.o: library.o potential.o system.o
potential.o: library.o
system.o: library.o

%.o: %.f90
	$(F90) $(FFLAGS) -c $<

all: run_test.x

run_test.x: $(src) run_test.o 
	$(F90) $(FFLAGS) -o $@ run_test.o $(src)
	ln -sf $@ ./bin/$@

clean:
	rm -rf *.o *.mod

-include make.deps