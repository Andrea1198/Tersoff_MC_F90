F90=gfortran
FFLAGS=-O1 -fopenmp
.SUFFIXES= .o .f90
src=library.o\
	potential.o\
	system.o\
	calcener.o\
	eqmc.o\

exe=sim.x

%.o: %.f90
	$(F90) $(FFLAGS) -c $<

all: sim.x

$(exe): $(src) run_test.o 
	$(F90) $(FFLAGS) -o $@ run_test.o $(src)
	ln -sf $@ ../bin/$@

run:
	./$(exe)
	
clean:
	rm -rf *.o *.mod 

-include make.deps
