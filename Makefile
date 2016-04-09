include ../config/config.mk

.PHONY : main optimize apples clean rmpic

# MAIN : takes in .bin files and gives cutout
main : main.exe clean
	./$<
main.exe: global.o fileio.o fourier.o linsolve.o main.o
	$(FC) -o $@ main.o global.o fileio.o $(HDF5_LIB) fourier.o $(FFTW_LIB) linsolve.o $(LFLAGS)


# OPTIMIZE : runs the optimization problem on the cutout data.
optimize : optimize.exe
	time ./$<
	rm -f *.o *.mod *# *~
optimize.exe: fileio.o linsolve.o optimize.o
	$(FC) -o $@ optimize.o fileio.o $(HDF5_LIB) linsolve.o $(LFLAGS)



# APPLES : Validation case.
apples : apples.exe
	./$<
apples.exe: fileio.o fourier.o  apples.o
	$(FC) -o $@ apples.o fileio.o $(HDF5_LIB) fourier.o $(FFTW_LIB) $(LFLAGS)


apples.o : apples.f90
	$(FC) $(FFLAGS) -c $<
global.o : global.f90
	$(FC) $(FFLAGS) -c $< 
fileio.o : fileio.f90
	$(FC) $(FFLAGS) -c $< $(HDF5_INC)
fourier.o : fourier.f90
	$(FC) $(FFLAGS) -c $< $(FFTW_INC)
linsolve.o : linsolve.f90
	$(FC) $(FFLAGS) -c $< 
main.o : main.f90
	$(FC) $(FFLAGS) -c $<
optimize.o : optimize.f90
	$(FC) $(FFLAGS) -c $<

# CLEAN:
clean:
	rm -f *.o  *.mod *[#~] 

