program tryfft1

!!$ By default the INCLUDE statement is used. Here we stick with the USE function.
!!$  use, intrinsic :: iso_c_binding
!!$  include 'fftw3.f03'

  !use, intrinsic :: iso_c_binding
 ! use fftw
  use fourier
 
implicit none
integer,parameter           :: GRID=256
integer,parameter           :: LES_scale=4, test_scale=16
integer(C_INT) :: n,i,j,k
type(C_PTR)::plan

double complex, allocatable :: in(:,:,:), out(:,:,:)
double complex, allocatable :: LES_cmplx(:,:,:)
real(C_DOUBLE)::              twopi,xj
real,dimension(1:GRID,1:GRID,1:GRID) :: LES


n = 256
twopi = 2.*acos(-1.)
allocate(in(n,n,n))
allocate(out(n,n,n))
allocate(LES_cmplx(GRID,GRID,GRID))
in =0.; out=0.
write(*,*) 'Input data:'
do k=1,n
do j=1,n
   do i=1,n
      in(i,j,k) = cmplx(2.*sin(real(2.*twopi*i/n)) + &
                           sin(real(twopi*j/n)) + &
                       4.*sin(real(8.*twopi*k/n)),0.0)
enddo
enddo
enddo

!!$do j=1,n
!!$    print*,''
!!$do i=1,n
!!$   write(*,*) real(in(i,j,:))
!!$enddo
!!$enddo
print*, real(in(10,10,10))
print*, twopi*0.5d0
! Forward Fourier transform
call dfftw_plan_dft_3d(plan,n,n,n,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
write(*,*) 'Fourier transform of the input data:'

call window(LES,GRID,LES_scale)
call matrixview(LES)
call fftshift (LES,GRID)
call matrixview(LES)

LES_cmplx = cmplx(LES)
print*, LES_cmplx(1,1,1)
!call matrixview(real(LES_cmplx))

out = out*LES_cmplx

!!$do k=1,n
!!$do j=1,n
!!$do i=1,n
!!$  write(*,*) i,j,k,out(i,j,k)/n**3
!!$enddo
!!$enddo
!!$enddo
 

! Inverse Fourier transform
call dfftw_plan_dft_3d(plan,n,n,n,out,in,FFTW_BACKWARD,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
write(*,*) 'Recovered data from inverse Fourier transform:'

!!$do j=1,n
!!$   print*,''
!!$do i=1,n
!!$  write(*,*) real(in(i,j,:)/n**3)
!!$enddo
!!$enddo
print*, real(in(10,10,10)/n**3)
end program tryfft1
