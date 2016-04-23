#!----------------------------------------------------------------
#!                          MAKEFILE
#!----------------------------------------------------------------

# config.mk contains all system based configuration settings.
# macros.mk contains user-defined macros - create-exe.

# DEVELOPMENT NOTES:
# The next thing will be to build a library with global, fileio, fourier and solver.
# Archive library and link with -l flag while making the executable.

#
# INCLUDE FILES:
#----------------------------------------------------------------
include $(WORK)/../config/config.mk
include $(WORK)/../config/macros.mk
#----------------------------------------------------------------


#
# AUTONOMIC : Validation case.
OBJECTS_AUTONOMIC = global.o fileio.o fourier.o actools.o solver.o autonomic.o
autonomic : $(OBJECTS_AUTONOMIC)
	 $(build-exe)


#
# MAIN : takes in .bin files and gives cutout
OBJECTS_MAIN = global.o fileio.o fourier.o solver.o main.o
main :  $(OBJECTS_MAIN) 
	$(build-exe)


#
# OPTIMIZE : runs the optimization problem on the cutout data.
OBJECTS_OPTIMIZE = global.o fileio.o solver.o optimize.o
optimize : LFLAGS = -lblas -llapack $(HDF5_LIB)
optimize : $(OBJECTS_OPTIMIZE)
	   $(build-exe)


#
# TESTCUT : test run for colocated formuation.
OBJECT_TESTCUT = global.o fileio.o solver.o testCut.o
testCut : $(OBJECTS_TESTCUT)
	$(build-exe)

#----------------------------------------------------------------

#
# BUILD OBJECTS:
fileio.o  : FFLAGS += $(HDF5_INC)
fourier.o : FFLAGS += $(FFTW_INC)
%.o : %.f90
	@$(FC) -c $(FFLAGS) $< 


#
# PHONY TARGETS:
.PHONY :  clean rmpic
clean  :
	rm -f *.o  *.mod *[#~] 
rmpic  :
	rm *.png



