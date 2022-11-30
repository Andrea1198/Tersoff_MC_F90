F90=gfortran
FFLAGS=-O3 -fopenmp -pg 
.SUFFIXES= .o .f90
src=library.o\
	potential.o\
	system.o\
	compute_gdr.o\
	calcener.o\
	calcener_MC.o\
	eqmc.o\

exe=sim.x

%.o: %.f90
	$(F90) $(FFLAGS) -c $<

all: sim.x

$(exe): $(src) run_test.o 
	$(F90) $(FFLAGS) -o $@ run_test.o $(src)
	cp $@ ../bin/$@

run:
	./$(exe)
	
clean:
	rm -rf *.o
	rm ./files/output/*

-include make.deps
