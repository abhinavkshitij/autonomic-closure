VERSION CONTROL NOTES
======================

BRANCH : test/validation  
-------------------------

* Validation passed with ij=11 to 33
* Results compared with MATLAB - agrees
* Added binary file write in main.f90; dat uses too much space. 
* Compute deviatoric stress component
* Compute tau_ij^F globally 
* Use non-colocated formulation
* Projection on all points
* Compare with MATLAB code - at z=128 plane for ij = 11


PURPOSE
-------

1. Validate code
2. Run code on different lambda
3. Test parallel using OpenMP	

CHANGES FROM PARENT COMMIT:
----------------------------

1. Single point computation moved to test/1 branch.
2. No single point based debug case in this version.
    *  Removed DEGCON().No condition number computation in DGESV(). 
    *  Removed bubblesort(). 
    *  Only one lambda value. Removed LAMBDA and ENSEMBLE loops.

ISSUES:
-------

1. u_t calculated from u_f was earlier done only upto single precision. The cmplx() function
   was working in single precision for this opertion. This has been fixed by using 
   the -fdefault-real-8 switch in the gfortran compiler. 
----------------------------------------------------------------

REQUIRED DEPENDENCIES:
1. HDF5 
2. FFTW

Since it takes a long time to run codes, the whole process has been broken into two parts. 
Each phase has its own Makefile associated with it. 

1. In the first phase, the dataset is filtered at LES and test scales by main.f90.
   Next, based on the parameters, main.f90 writes the appropriate size of the cutout fields on external files.
   Run Makefile to complete Phase I


2. In the next phase the saved cutout data is read by optimize.f90
   This code then computes the tau_ij at LES scale by solving the damped least squares problem.
   optimize.f90 uses the module linsolve to carry out linear algebra computations.
   The key parameters are computed at the start of the linsolve module.
   Run Makefile2.mk to complete Phase II

DESCRIPTION OF MODULES:

1. LINSOLVE:

   Linsolve contains routines for:
   a) Taking cutout of the test field				-- cutout()
   b) Random number generation 			 		-- init_random_seed() 
   c) Randomly marking stencil-center points 	 		-- randTrainingSet()
   d) Computing tau_ij after solving least squares problem. 	-- synStress()
   e) Solving the inverse problem by direct iversion		-- choleskyLU() 
   f) Computing the singular values  	    			-- SVD()
   g) Sorting a list	     					-- bubblesort()
  	    			
2) FOURIER:

   Fourier contains routines for:
   a) Creating a spectrally sharp filter			-- createFilter()
   b) Performing the fftshift operation (its a MATLAB op)	-- fftshift()
   c) Applying the  spectrally sharp filter     		-- sharpFilter()
  
3) FILEIO:

   Fileio conducts read/write operations on external files either during preprocessing 
   or postprocessing. It has routines for:
   a) Reading the NRL (256^3) dataset (in binary format)	-- binRead()
   b) Reading JHU (1024^3) dataset (in HDF5 format)		-- hdf5Read()
   c) Writing planes slices to plot in MATLAB			-- matrixview()
   d) Print an array in the form of a matrix on the console	-- printmatrix2()

4) GLOBAL:
   
   Contains global definitions and variables that are common in all other modules and procedures
   	    * Includes GRID size [=256]
	    * Includes parameters
	    * Has the value of pi





