program main
! STATUS : Integrated with readfile.f90 as a subroutine. 9/11/2015 (Passed)
!          Integrated with two .mod files. 9/13/2015.(Passed)
! Result : Passed 
! Notes  : 

use fourier
use fileio
implicit none

! Define global parameters:

integer,parameter           :: LES_scale=16, test_scale=4
integer                     :: i,j,k,DIM

! Define velocities:
real(kind=4),dimension(1:3,GRID,GRID,GRID) :: u_s
real(kind=8),dimension(1:3,GRID,GRID,GRID) :: u
real(kind=8),dimension(GRID,GRID,GRID) :: u1_t, u2_t, u3_t
real(kind=8),dimension(GRID,GRID,GRID) :: u1_th, u2_th, u3_th
real(kind=8),dimension(GRID,GRID,GRID) :: LES,test


call binRead(u_s,DIM=1)
u = u_s/1.d2

print *, u_s(1,1,1,1) , u(1,1,1,1)
print *, u_s(1,1,1,10) , u(1,1,1,10)

call hdf5Read()
stop
!!$call createFilter(LES,LES_scale)
!!$call createFilter(test,test_scale)
!!$call fftshift(LES)
!!$call fftshift(test)
!!$
!!$u1_t = sharpFilter(u,LES)
!!$u1_th = sharpFilter(u1_t,test)
!!$
!!$
!!$print *, '' ! Blank Line
!!$print *, u1(1,1,1,10)
!!$print *, u1_t(1,1,1,10)
!!$print *, u1_th(1,1,1,10)
 
end program main
