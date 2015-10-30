# STATUS : Integrating LAPACK links
# Notes : Make the links more concise 

# external library directories

HDF5_DIR = /opt/hdf5-1.8.14


# HDF5
HDF5_INC = -I$(HDF5_DIR)/include
HDF5_LIB = -L$(HDF5_DIR)/lib -lhdf5 -lhdf5_fortran

FC      = gfortran
LFLAGS  = -lblas -llapack
FFLAGS  = -O3

.PHONY : testCut clean

testCut : testCut.exe
	./testCut.exe

testCut.exe: fileio.o linsolve.o testCut.o
	$(FC) -o testCut.exe testCut.o fileio.o $(HDF5_LIB) linsolve.o $(LFLAGS)

%.o : %.f90
	$(FC) $(FFLAGS) -c $< $(HDF5_INC)

clean:
	rm -f *.o *.exe  *.mod *# *~ 

