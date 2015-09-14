# STATUS : Tested on 9/11/15. 
# Result : Passed

# Notes  : Runs with mpiexec -n 2 but reads the file twice.

FC      = gfortran
FFLAGS  = -O3
OBJECTS = fourier.o fileio.o main.o 
MODULES = fourier.mod fileio.mod

.PHONY : ALES clean

ALES : ALES.exe
	./ALES.exe

ALES.exe: $(MODULES) $(OBJECTS)
	$(FC) $(OBJECTS) -o ALES.exe

%.o : %.f90
	$(FC) $(FFLAGS) -c $< 
%.mod: %.f90
	$(FC) $(FFLAGS) -c $<
clean:
	rm -f *.o *.exe  *.mod *# *~ 


