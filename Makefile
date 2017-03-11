#!----------------------------------------------------------------
#!                          MAKEFILE
#!----------------------------------------------------------------

# config.mk contains all system based configuration settings.
# macros.mk contains user-defined macros - create-exe.

# DEVELOPMENT NOTES:
# 	* Build a library with global, fileio, fourier and solver.
# 	* Archive library and link with -l flag while making the executable.
# 	* Create .lst file. Send filenames to .lst files. Refer ../makefile/makefile_BLAS

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
#	 $(plot-exe)


#
# MAIN : takes in .bin files and gives cutout
OBJECTS_MAIN = global.o fileio.o fourier.o actools.o solver.o main.o
main :  $(OBJECTS_MAIN) 
	$(build-exe)


#
# OPTIMIZE : runs the optimization problem on the cutout data.
OBJECTS_OPTIMIZE = global.o fileio.o solver.o optimize.o
optimize : LFLAGS = -lblas -llapack $(HDF5_LIB)
optimize : $(OBJECTS_OPTIMIZE)
	   $(build-exe)

#
# JHU1024 : Dealias jhu1024 DATASET
OBJECTS_JHU1024 = global.o fileio.o fourier.o  jhu_dealias.o
dealias_jhu1024 : $(OBJECTS_JHU1024)
	 $(build-exe)
#	 $(plot-exe)


#
# TESTCUT : test run for colocated formulation.
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


#----------------------------------------------------------------


OUTPUT_DIR=$()
DAT_FILES=$(wildcard $(OUTPUT_DIR)/dat/*.dat)
PNG_FILES=$(patsubst %.dat, %.png, $(DAT_FILES))
vpath %.dat run/dat
vpath %.png run/plots

#
# BUILD PLOTS:
#plots : %.png : %.dat
#	$(PLOT_EXE) $< $*.dat

plots :
	$(plot-exe)


#
# PHONY TARGETS:
.PHONY :  clean rmpic
clean  :
	rm -f *.o  *.mod *[#~] 

veryclean  :
	rm -f *.o  *.mod *[#~] *.err *.out

rmpic  :
	rm *.png



