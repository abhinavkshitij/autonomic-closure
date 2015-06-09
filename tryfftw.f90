program testfft

!!$ By default the INCLUDE statement is used. Here we stick with the USE function.
!!$  use, intrinsic :: iso_c_binding
!!$  include 'fftw3.f03'
  
  use fftw
 
! Example to call 1-D real FFT routine of FFTW


  integer,parameter::    N=16        !! sample size
  integer*8:: PLAN_FOR, PLAN_BAC     !! forward and backward plans
  real*8,dimension(N)::IN,OUT,IN2    !! IN for input vector , OUT for output vector after FFT
  real*8::twopi,xj
  integer::j,k,mode
  
  twopi = 2.*acos(-1.)
  
  !! Discrete data of function f(x)=cos(x)+0.2*sin(2x)
  !! Here an explicit function has been taken for the test case. The data will be read from directly fed to IN

  do j=0,N-1
     xj = twopi*real(j) / real(N)
     IN(j) = cos(xj) + 0.2*sin(2.*xj)
  end do

  write(*,*) "Original data"

  do j=1,N
     write(*,100) j,IN(j)
  end do

100 format(i4,f12.5)

  !! Forward transform
  !! Plans are not saved on any variable!!
  call dfftw_plan_r2r_1d(PLAN_FOR,N,IN,OUT,FFTW_R2HC,FFTW_ESTIMATE)
  call dfftw_execute_r2r(PLAN_FOR,IN,OUT)

  OUT = OUT / real(N,KIND=8) ! Normalize

  write(*,*) "Fourier coefficient after forward FFT"

  do k=1,N
     mode=k-1
     if(k > N/2+1) mode=N-k+1
     write(*,100) mode,OUT(k)
  end do

  ! Backward transform
  call dfftw_plan_r2r_1d(PLAN_BAC,N,OUT,IN2,FFTW_HC2R,FFTW_ESTIMATE)
  call dfftw_execute_r2r(PLAN_BAC,OUT,IN2)

  write(*,*) "Data after backward FFT"

  do j=1,N
     write(*,100) j,IN2(j)
  end do

  ! Destroy the plans
  call dfftw_destroy_plan(PLAN_FOR)
  call dfftw_destroy_plan(PLAN_BAC)
  
end program testfft

!!$I assume that you have compiled FFTW to include F77 wrappers and that you have read Chapter 6 "Calling FFTW from Fortran" in the FFTW manual.
!!$
!!$If you look inside the file fftw3.f, you will see that it contains only integer parameter declarations, so it is a file to be "included"; you cannot compile fftw3.f into a module file and "use" it.
!!$
!!$The simplest way to compile a Fortran program and statically link the FFTW library is to issue a command like this:
!!$
!!$ifort -o test-fftw3 test-fftw3.f90 -L/usr/local/lib -lfftw3
!!$
!!$where -L gives ifort the path that allows it to find libfftw3.a (change this if your FFTW is not installed in /usr/local/lib), and -lfftw3 tells ifort to link libfftw3.a into your program. Note that libfftw3.a is the double-precision FFTW library; if you need the single-precision FFTW library, specify -lfftw3f for libfftw3f.a . Both single- and double-precision FFTW libraries can be linked into the same program if both are needed. You can add other compiler flags; the above is the minimum.
!!$
!!$It seems to me in your post that you might also have problems calling FFTW in your Fortran program. You will need to read and understand Chapter 6 in the FFTW manual, and if you know how to program in C, read Chapter 2 before you read Chapter 6. Post back if you have a specific question.
!!$
!!$I have used FFTW 3.1 with ifort 9.1.040 and it does work. I haven't used FFTW 3.1.2 with ifort 9.1.040 (the latest versions for both) but I would expect it to work too.
