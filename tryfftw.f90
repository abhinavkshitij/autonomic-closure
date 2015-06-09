program testfft

!!$ By default the INCLUDE statement is used. Here we stick with the USE function.
!!$  use, intrinsic :: iso_c_binding
!!$  include 'fftw3.f03'
  
  use fftw
 
! Example to call 1-D real FFT routine of FFTW


  integer(C_INT),parameter::    N=16                 !! sample size
  type(C_PTR):: PLAN_FOR, PLAN_BAC            !! forward and backward plans
  real(C_DOUBLE),dimension(N)::IN,OUT,IN2     !! IN for input vector , OUT for output vector after FFT
  real(C_DOUBLE)::twopi,xj
  integer(C_INT)::j,k,mode
  
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


!!$The simplest way to compile a Fortran program and statically link the FFTW library is to issue a command like this:
!!$
!!$gfortran -o test-fftw3 test-fftw3.f90 -L/usr/local/lib -lfftw3
!!$
!!$where -L gives ifort the path that allows it to find libfftw3.a (change this if your FFTW is not installed in /usr/local/lib), and -lfftw3 tells ifort to link libfftw3.a into your program. Note that libfftw3.a is the double-precision FFTW library; if you need the single-precision FFTW library, specify -lfftw3f for libfftw3f.a . Both single- and double-precision FFTW libraries can be linked into the same program if both are needed. You can add other compiler flags; the above is the minimum.
