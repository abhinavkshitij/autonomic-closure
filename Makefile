# STATUS : Operational 
# Notes : Make the links more concise 

# external library directories
FFTW_DIR = /opt/fftw-3.3.4
HDF5_DIR = /opt/hdf5-1.8.14

# FFTW 
FFTW_INC = -I$(FFTW_DIR)/include
FFTW_LIB = -L$(FFTW_DIR)/lib -lfftw3 -lm

# HDF5
HDF5_INC = -I$(HDF5_DIR)/include
HDF5_LIB = -L$(HDF5_DIR)/lib -lhdf5 -lhdf5_fortran

FC      = gfortran
FFLAGS  = -O3

.PHONY : main clean

main : main.exe
	./main.exe

main.exe: fileio.o fourier.o main.o
	$(FC) -o main.exe main.o fileio.o $(HDF5_LIB) fourier.o $(FFTW_LIB)

%.o : %.f90
	$(FC) $(FFLAGS) -c $< $(HDF5_INC) $(FFTW_INC)

clean:
	rm -f *.o *.exe  *.mod *# *~ 


