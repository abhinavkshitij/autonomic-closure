# STATUS : Tested on 9/11/15. 
# Result : Passed

# external library directories
FFTW_DIR = /opt/fftw-3.3.4
HDF5_DIR = /opt/hdf5-1.8.14

# external library link and include commands
FFTW_INC = -I$(FFTW_DIR)/include
FFTW_LIB = -L$(FFTW_DIR)/lib -lfftw3 -lm

HDF5_INC = -I$(HDF5_DIR)/include
HDF5_LIB = -L$(HDF5_DIR)/lib -lhdf5 -lhdf5_fortran

# Link with external libraries
FFTW_DIR = /opt/fftw-3.3.4

FFTW_INC = -I$(FFTW_DIR)/include
FFTW_LIB = -L$(FFTW_DIR)/lib -lfftw3 -lm 

FC      = gfortran

FLIBS = $(FFTW_LIB) $(HDF5_LIB)
FFLAGS  = -O3 $(HDF5_INC) $(FFTW_INC)
OBJECTS_FFTW = fftw.o fourier.o tryfftw1.o
OBJECTS_HDF5 = h5_crtdat.o


.PHONY : ALES h5test clean



# MAIN:
main : main.exe
	./main.exe


ALES.exe: $(OBJECTS_FFTW)
	$(FC) $(OBJECTS_FFTW) -o ALES.exe $(FLIBS)

h5test : h5test.exe
	./h5test.exe

h5test.exe: $(OBJECTS_HDF5)
	$(FC) $(OBJECTS_HDF5) -o h5test.exe $(FLIBS)

%.o : %.f90
	$(FC) $(FFLAGS) -c $<  
#%.mod: %.f90
#	$(FC) $(FFLAGS) -c $<


clean:
	rm -f *.o *.exe  *.mod *# *~ 


