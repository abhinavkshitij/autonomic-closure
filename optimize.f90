program optimize

! /Users/Kshitij/Desktop/ALES/f90/src/optimize.f90  
! STATUS : Testing for damped least squares formulation. This program will read the filtered fields and compute the inverse.
!          
!          
! Result : 
! Notes  : 

use fileio
use linsolve
implicit none


real(8), dimension(3,testcutSize,testcutSize,testcutSize):: u_f    , u_t
real(8), dimension(6,testcutSize,testcutSize,testcutSize):: tau_ij , T_ij !Testing for _11 ij component

character(50):: CUT_DATA = '../derived_data/cutout/bin4020/' !Change bin4020 by 'sed' in shell script
integer :: n_u=3, n_uu=6

!!$ DEBUG SWITCHES:
integer,dimension(4) :: debug=(/0,0,0,0/)


call system('clear')
call printParams()

open(1,file= trim(CUT_DATA)//'u_f.dat') ! 3 components
open(2,file= trim(CUT_DATA)//'u_t.dat') ! 3 components
open(3,file= trim(CUT_DATA)//'tau_ij.dat') ! 6 components
open(4,file= trim(CUT_DATA)//'T_ij.dat')  ! 6 components

print *, "*** Reading files ***"
 
read(1,*) u_f
read(2,*) u_t
read(3,*) tau_ij
read(4,*) T_ij

close(1)
close(2)
close(3)
close(4)

!!$ For testing correct read/write:
if(debug(1).eq.1) then
 print*,"T_ij (1,11,11,11):", T_ij(1,11,11,11)
 print*,"tau_ij (1,11,11,11):",tau_ij(1,11,11,11)
end if

call synStress(u_f, u_t, tau_ij, T_ij, n_u, n_uu)



end program optimize
