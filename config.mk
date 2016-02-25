# Settings and configurations:

# external library directories
FFTW_DIR := /opt/fftw-3.3.4
HDF5_DIR := /opt/hdf5-1.8.14

# FFTW 
#FFTW_INC := -I$(FFTW_DIR)/include
FFTW_LIB := -L$(FFTW_DIR)/lib -lfftw3 -lm

# HDF5
HDF5_INC := -I$(HDF5_DIR)/include
HDF5_LIB := -L$(HDF5_DIR)/lib -lhdf5 -lhdf5_fortran

#vpath %.f90 include

FC      := gfortran
LFLAGS  := -lblas -llapack
FFLAGS  := -O3 -fdefault-real-8


