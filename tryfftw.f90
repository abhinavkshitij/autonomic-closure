program tryfft

  !! CODE AND CONVENTIONS FOR MODERN FORTRAN
  
  use fftw
  

  integer, parameter                      :: GRID = 4          
  integer(C_INT), parameter               :: M = GRID, N = GRID
  
  real(C_DOUBLE),dimension(0:M-1,0:N-1)           :: IN, IN2
  complex(C_DOUBLE_COMPLEX),dimension(M/2+1,N):: OUT
  type(C_PTR)                             :: plan_f, plan_b, data
        
             
  real(KIND=8)                            :: twopi,x,y
  integer                                 :: i,j,k,mode
  
  twopi = 2.*acos(-1.)

  call system ('clear')
   
  100 format(i4,',',f12.5)

!!$ ALLOCATE MEMORY:

!!$  data = fftw_alloc_complex(int((M/2)+1 * N, C_SIZE_T))
!!$  call c_f_pointer(data, IN, [2*(M/2+1),N])
!!$  call c_f_pointer(data, OUT, [M/2+1, N])
!!$
  
!!$ DEFINE FUNCTION:
  
do j=0,N-1
   do i=0,M-1
        x = twopi*real(j) / real(M)
        y = twopi*real(i) / real(N)
        
        IN(i,j) = sin(x) + sin(y)
        
     write(*,*) i,j,IN(i,j) 
  end do
end do

!!$ FORWARD FFT:

  plan_f = fftw_plan_dft_r2c_2d(N,M,IN,OUT,FFTW_ESTIMATE)
  call fftw_execute_dft_r2c(plan_f,IN,OUT)

  OUT = OUT / real(M*N,KIND=8) ! Normalize by MxN

 !! write(*,*) "Fourier coefficients after forward FFT"

  do j=1,N
     do i=1,M
     mode=j-1
     write(*,*) mode,OUT(1,k)
  end do
  end do

!!$ !! write (*,*) "applying filter at mode 1,2 to test "
!!$  do k=1,N
!!$     mode=k-1
!!$     if(k > N/2+1) mode=N-k+1
!!$     if (mode == 2) OUT(k) =0.
!!$     write(*,*) mode,OUT(1,k)
!!$  end do
!!$  
  ! Backward transform


  plan_b = fftw_plan_dft_c2r_2d(N,M,OUT,IN2,FFTW_ESTIMATE)
  call fftw_execute_dft_c2r(plan_b,OUT,IN2)

 !! write(*,*) "Data after backward FFT"

  do j=0,N-1
     do i=0,M-1
        write(*,*) i,j,IN2(i,j)
     end do    
  end do

  ! Destroy plans:
  call fftw_destroy_plan(plan_f)
  call fftw_destroy_plan(plan_b)

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
