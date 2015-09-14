# STATUS : Tested on 9/11/15. 
# Result : Passed

# Notes  : Runs with mpiexec -n 2 but reads the file twice.

# Link with external libraries
FFTW_DIR = /opt/fftw-3.3.4

FFTW_INC = -I$(FFTW_DIR)/include
FFTW_LIB = -L$(FFTW_DIR)/lib -lfftw3 -lm 

FC      = gfortran
FFLAGS  = -O3

#OBJECTS = fourier.o fileio.o main.o 
OBJECTS = fourier.o tryfftw.o
#MODULES = fourier.mod fileio.mod
MODULES = fourier.mod

.PHONY : main fft clean

# MAIN:
main : main.exe
	./main.exe

main.exe: $(MODULES) $(OBJECTS)
	$(FC) $(OBJECTS) -o main.exe

%.o : %.f90
	$(FC) $(FFLAGS) -c $< 
%.mod: %.f90
	$(FC) $(FFLAGS) -c $<

# FFT:

fft : fft.exe
	./fft.exe

fft.exe:  $(MODULES) $(OBJECTS)
	$(FC) $(OBJECTS) -o fft.exe $(FFTW_LIB)

%.o : %.f90
	$(FC) $(FFLAGS) -c $<
%.mod:%.f90
	$(FC) $(FFLAGS) -c $<


# CLEAN:
clean:
	rm -f *.o *.exe  *.mod *# *~ 


