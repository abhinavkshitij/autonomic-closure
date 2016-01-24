# USE THIS MAKEFILE TO RUN main.exe THAT WILL TAKE THE FFT AND SAVE ON AN EXTERNAL FILE.

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
LFLAGS  = -lblas -llapack
FFLAGS  = -O3

.PHONY : main clean rmpic

main : main.exe
	./main.exe

main.exe: fileio.o fourier.o linsolve.o main.o
	$(FC) -o main.exe main.o fileio.o $(HDF5_LIB) fourier.o $(FFTW_LIB) linsolve.o $(LFLAGS)

%.o : %.f90
	$(FC) $(FFLAGS) -c $< $(HDF5_INC) $(FFTW_INC)

clean:
	rm -f *.o *.exe  *.mod *# *~ 

rmpic:
	rm -f *.jpg *.png

