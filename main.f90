program main
  
! STATUS : Integrated with readfile.f90 as a subroutine. 9/11/2015 (Passed)
!          Integrated with two .mod files. 9/13/2015.(Passed)
!          Incorporated HDF5 read capability. 10/22/2015(Passed)
! Result : Passed 
! Notes  : 

use fourier
use fileio
use linsolve

implicit none
integer,parameter           :: LES_scale=40, test_scale=20
integer                     :: i,j,k,d=0

! Define velocities:
real(8),allocatable,dimension(:,:,:,:) :: u, u_f, u_t
real(8),allocatable,dimension(:,:,:,:) :: tau_ij,T_ij
real(8),allocatable,dimension(:,:,:) :: LES,test

integer :: n_u, n_uu

character(50):: CUT_DATA = '../derived_data/cutout64/jhu/' !Change bin4020 by 'sed' in shell script
character(10):: f_CUT 
character(3) :: d_set = 'jhu'   ! for binRead()

! DEBUG FLAGS:
integer,dimension(2) :: debug=(/0,1/)
! 1- To select only one velocity component and 3 velocity products.
! 2- Choose between NRL/JHU256 database[1] or JHU(HDF5) database[0]

real :: tic, toc

call system('clear')
call printParams()

!! Set debug flags for velocity components:
n_u=3; n_uu = 6
if (debug(1).eq.1) then
   n_u = 1
   n_uu = 3
   print*, 'Debug mode for velocity components...'
end if


!! Select file to read:
allocate(u(n_u,GRID,GRID,GRID))

fileSelect:if (debug(2).eq.1) then
   call binRead(u,  d_set,  DIM=n_u)
   write(f_CUT,'(a3,2(i2),a1)') 'bin',LES_scale,test_scale,'/' !Write dirname
   !print*, trim(CUT_DATA)//trim(f_CUT) !Test pathname
else
   call hdf5Read() 
end if fileSelect

print*,'u(1,15,24,10):',u(1,15,24,10)



!! Create LES and test scale sharp filters:
allocate(LES(GRID,GRID,GRID))
allocate(test(GRID,GRID,GRID))
call createFilter(LES,LES_scale)
call createFilter(test,test_scale)
call fftshift(LES)
call fftshift(test)


print *, 'Applying filters:'
!! Take LES and test filtered fields:
allocate(u_f(n_u,GRID,GRID,GRID))
allocate(u_t(n_u,GRID,GRID,GRID))
filter:do i=1,n_u
   u_f(i,:,:,:) = sharpFilter(u(i,:,:,:),LES) ! Speed up this part -- Bottleneck
   u_t(i,:,:,:) = sharpFilter(u(i,:,:,:),test)
end do filter
print*,'Check sharpFilter():'
print*,'u_f(1,15,24,10):',u_f(1,15,24,10)
print*,'u_t(1,15,24,10):',u_t(1,15,24,10)


print*, '' ! Blank Line
print*,  'u(1,1,1)'   , u  (1,lBound+testCutsize-1,lBound+testCutsize-1,lBound+testCutsize-1)
print*,  'u_f(1,1,1)' , u_f(1,lBound+testCutsize-1,lBound+testCutsize-1,lBound+testCutsize-1)
print*,  'u_t(1,1,1)' , u_t(1,lBound+testCutsize-1,lBound+testCutsize-1,lBound+testCutsize-1)
print*,''



allocate(tau_ij(n_uu,GRID,GRID,GRID))
allocate(T_ij(n_uu,GRID,GRID,GRID))
call computeStress(u,u_f,u_t,tau_ij,T_ij,n_u,n_uu,LES,test,stress='abs')
deallocate(LES,test)

stop

print*,'Shape before cutout:',shape(u_t)
print*, u_t(1,testLower+lBound-1,testLower+lBound-1,testLower+lBound-1)
call matrixview(u_t(1,:,:,:),frameLim=5,z=testLower+lBound-1)

call cutout(u_t,n_u)
call cutout(u_f,n_u)
call cutout(T_ij,n_uu)
call cutout(tau_ij,n_uu)

print*, "Done cutout..."
print*, u_t(1,testLower,testLower,testLower)
print*, u_t(1,testLower+6,testLower+6,testLower)
call matrixview(u_t(1,:,:,:),frameLim=17,z=testLower) !The bottom-right element should match with the value above
stop

! Write velocities & stresses to files
open(1,file=trim(CUT_DATA)//trim(f_CUT)//'u_f.dat',status='replace')
open(2,file=trim(CUT_DATA)//trim(f_CUT)//'u_t.dat',status='replace')
open(3,file=trim(CUT_DATA)//trim(f_CUT)//'tau_ij.dat',status='replace')
open(4,file=trim(CUT_DATA)//trim(f_CUT)//'T_ij.dat',status='replace')

 write(1,*) u_f
 write(2,*) u_t
 write(3,*) tau_ij
 write(4,*) T_ij

 !! For testing correct read/write

 print*,"Testing read/write at T_ij (3,2,1,4):"
 print*,T_ij(3,2,1,4)

close(1)
close(2)
close(3)
close(4)
stop


call synStress(u_f,u_t,tau_ij,T_ij,n_u,n_uu)
!print*,'tau_ij',tau_ij(1,testLower,testLower,testLower)


end program main
