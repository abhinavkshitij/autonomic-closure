program optimize
  
! STATUS : Testing for damped least squares formulation. This program will read the filtered fields and compute the inverse.
!          
!          
! Result : 
! Notes  : 

use fileio
use linsolve
implicit none

! Define velocities:

real(kind=8), dimension(3,testcutSize,testcutSize,testcutSize):: u_f    , u_t
real(kind=8), dimension(6,testcutSize,testcutSize,testcutSize):: tau_ij , T_ij

integer :: n_u=3, n_uu=6
integer :: i,j,k,d=1

integer,dimension(4) :: debug=(/0,1,1,1/)
real :: tic, toc

call system('clear')
!call printParams()


!! Select file to read:
open(1,file='./testOpt/bin4020/u_f.dat') ! 3 components
open(2,file='./testOpt/bin4020/u_t.dat') ! 3 components
open(3,file='./testOpt/bin4020/tau_ij.dat') ! 6 components
open(4,file='./testOpt/bin4020/T_ij.dat')  ! 6 components

 read(1,*) u_f
 read(2,*) u_t
 read(3,*) tau_ij
 read(4,*) T_ij

 !! For testing correct read/write

! print*,"T_ij (1,11,11,11):"
! print*,T_ij(1,11,11,11)
! print*,"tau_ij (1,11,11,11):"
! print*,tau_ij(1,11,11,11)

close(1)
close(2)
close(3)
close(4)



call synStress(u_f, u_t, tau_ij, T_ij, n_u, n_uu)




end program optimize
