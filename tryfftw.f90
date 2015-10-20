program tryfft

!!$ By default the INCLUDE statement is used. Here we stick with the USE function.
!!$  use, intrinsic :: iso_c_binding
!!$  include 'fftw3.f03'
  
  use fftw
 
! Example to call 1-D real FFT routine of FFTW


  integer(C_INT),parameter::    N=16                      !! sample size
  type(C_PTR)::                 PLAN_FWD, PLAN_BCK        !! forward and backward plans
  real(C_DOUBLE_COMPLEX),dimension(N):: IN,OUT,IN2                !! IN for input vector , OUT for output vector after FFT
  real(C_DOUBLE)::              twopi,xj
  integer(C_INT)::              j,k,mode
  
  twopi = 2.*acos(-1.)

  call system ('clear')
  
  !! Discrete data of function f(x)=sin(2x)
  !! Here an explicit function has been taken for the test case. The data will be read from directly fed to IN

  do i=1,N
  in(i) = cmplx(i,0.0)
  write(*,*) in(i)
enddo

  do j=0,N-1
     xj = twopi*real(j) / real(N)
     IN(j) = sin(xj)
     write(1,100) j,IN(j) 
  end do

 !! write(*,*) "Original data"

!!$  do j=1,N
!!$     write(*,100) j,IN(j)
!!$  end do



  !! Forward transform
  !! Plans are not saved on any variable!!
  
  call dfftw_plan_r2r_1d(PLAN_FWD,N,IN,OUT,FFTW_R2HC,FFTW_ESTIMATE)
  call dfftw_execute_r2r(PLAN_FWD,IN,OUT)

  OUT = OUT / real(N,KIND=8) ! Normalize

 !! write(*,*) "Fourier coefficient after forward FFT"

  do k=1,N
     mode=k-1
     if(k > N/2+1) mode=N-k+1
     write(2,100) mode,OUT(k)
  end do

 !! write (*,*) "applying filter at mode > 7"
  do k=1,N
     mode=k-1
     if(k > N/2+1) mode=N-k+1
     if (mode > 7) OUT(k) =0.
     write(3,100) mode,OUT(k)
  end do
  
  ! Backward transform
  call dfftw_plan_r2r_1d(PLAN_BCK,N,OUT,IN2,FFTW_HC2R,FFTW_ESTIMATE)
  call dfftw_execute_r2r(PLAN_BCK,OUT,IN2)

 !! write(*,*) "Data after backward FFT"

  do j=1,N
     write(4,100) j,IN2(j)
  end do

  ! Destroy the plans
  call dfftw_destroy_plan(PLAN_FWD)
  call dfftw_destroy_plan(PLAN_BCK)

  close(1)
  close(2)
  close(3)
  close(4)
  
end program tryfft


!!$The simplest way to compile a Fortran program and statically link the FFTW library is to issue a command like this:
!!$
!!$gfortran -o test-fftw3 test-fftw3.f90 -L/usr/local/lib -lfftw3
!!$
!!$where -L gives ifort the path that allows it to find libfftw3.a (change this if your FFTW is not installed in /usr/local/lib), and -lfftw3 tells ifort to link libfftw3.a into your program. Note that libfftw3.a is the double-precision FFTW library; if you need the single-precision FFTW library, specify -lfftw3f for libfftw3f.a . Both single- and double-precision FFTW libraries can be linked into the same program if both are needed. You can add other compiler flags; the above is the minimum.
