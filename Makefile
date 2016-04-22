#!----------------------------------------------------------------
#!                          MAKEFILE
#!----------------------------------------------------------------

# config.mk contains all system based configuration settings.
# macros.mk contains user-defined macros - create-exe.

# DEVELOPMENT NOTES:
# The next thing will be to build a library with global, fileio, fourier and linsolve.
# Archive library and link with -l flag while making the executable.

#
# INCLUDE FILES:
#----------------------------------------------------------------
include $(WORK)/../config/config.mk
include $(WORK)/../config/macros.mk
#----------------------------------------------------------------


#
# MAIN : takes in .bin files and gives cutout
OBJECTS_MAIN = global.o fileio.o fourier.o linsolve.o main.o
main :  $(OBJECTS_MAIN) 
	$(build-exe)


#
# OPTIMIZE : runs the optimization problem on the cutout data.
OBJECTS_OPTIMIZE = global.o fileio.o linsolve.o optimize.o
optimize : LFLAGS = -lblas -llapack $(HDF5_LIB)
optimize : $(OBJECTS_OPTIMIZE)
	   $(build-exe)


#
# APPLES : Validation case.
OBJECTS_APPLES = global.o fileio.o fourier.o actools.o linsolve.o apples.o
apples : $(OBJECTS_APPLES)
	 $(build-exe)


#
# ORDER : Test matrix order for performance
OBJECTS_ORDER = global.o order.o
order : $(OBJECTS_ORDER)
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



