program main

! STATUS : Integrated with readfile.f90 as a subroutine. 9/11/2015 (Passed)
!          Integrated with two .mod files. 9/13/2015.
! Result : Testing
 
! Notes  : 

use fourier
use fileio

implicit none

! Define global parameters:
integer,parameter           :: GRID=32
integer,parameter           :: LES_scale=8, test_scale=16

! Loop indices:
integer                     :: i,j,k

! Define velocities:
real,dimension(1:GRID,1:GRID,1:GRID) :: u1_DNS, u2_DNS, u3_DNS
real,dimension(1:GRID,1:GRID,1:GRID) :: LES

call readfile(GRID,u1_DNS,u2_DNS,u3_DNS)

call window(GRID,LES_scale,LES)

call matrixview(LES) ! Change the matrix size in the format!
print *, '' ! Blank Line
call fftshift (LES,GRID)

call matrixview(LES)

!print *, u1_DNS(1,1,1)
!print *, u2_DNS(128,128,128)
!print *, u3_DNS(1,128,32)
 
end program main
