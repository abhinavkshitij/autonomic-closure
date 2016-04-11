include ../config/config.mk

.PHONY : main optimize apples clean rmpic


# MAIN : takes in .bin files and gives cutout
OBJECTS_MAIN = global.o fileio.o fourier.o linsolve.o main.o
main  : $(OBJECTS_MAIN)
	$(FC) -o $@ $^ $(LFLAGS)
	./$@



# OPTIMIZE : runs the optimization problem on the cutout data.
OBJECTS_OPTIMIZE = global.o fileio.o linsolve.o optimize.o
optimize : LFLAGS = -lblas -llapack $(HDF5_LIB)
optimize : $(OBJECTS_OPTIMIZE)
	$(FC) -o $@ $^ $(LFLAGS)
	time ./$@



# APPLES : Validation case.
OBJECTS_APPLES = global.o fileio.o fourier.o apples.o
apples :$(OBJECTS_APPLES)
	$(FC) -o $@ $(LFLAGS) $^
	./$@


%.o : %.f90
	$(FC) -c $(FFLAGS) $< 

fileio.o  : FFLAGS += $(HDF5_INC)
fourier.o : FFLAGS += $(FFTW_INC)

# apples.o : apples.f90
# 	$(FC) $(FFLAGS) -c $<
# global.o : global.f90
# 	$(FC) $(FFLAGS) -c $< 
# fileio.o : fileio.f90
# 	$(FC) $(FFLAGS) $(HDF5_INC) -c $< 
# fourier.o : fourier.f90
# 	$(FC) $(FFLAGS) $(FFTW_INC) -c $< 
# linsolve.o : linsolve.f90
# 	$(FC) $(FFLAGS) -c $< 
# main.o : main.f90
# 	$(FC) $(FFLAGS) -c $<
# optimize.o : optimize.f90
# 	$(FC) $(FFLAGS) -c $<

# CLEAN:
clean:
	rm -f *.o  *.mod *[#~] 

