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

real(kind=8),allocatable,dimension(:,:,:,:) :: u, u_f, u_t
real(kind=8),allocatable,dimension(:,:,:,:) :: tau_ij,T_ij,temp
real(kind=8),allocatable,dimension(:,:,:) :: LES,test
real(kind=8):: dev_t
integer :: n_u, n_uu
integer,dimension(4) :: debug=(/0,1,1,1/)




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
   call binRead(u,DIM=n_u)
else
   call hdf5Read() !! Under testing to read multiple files
end if fileSelect



!! Create LES and test scale sharp filters:
allocate(LES(GRID,GRID,GRID))
allocate(test(GRID,GRID,GRID))
call createFilter(LES,LES_scale)
call createFilter(test,test_scale)
call fftshift(LES)
call fftshift(test)




!! Take LES and test filtered fields:
allocate(u_f(n_u,GRID,GRID,GRID))
allocate(u_t(n_u,GRID,GRID,GRID))
print *, 'FFT- u_ij:'
filter:do i=1,n_u
   u_f(i,:,:,:) = sharpFilter(u(i,:,:,:),LES) ! Speed up this part -- Bottleneck
   u_t(i,:,:,:) = sharpFilter(u_f(i,:,:,:),test)
end do filter
print *, '' ! Blank Line
print *,   'u(1,1,1)'  , u (1,lBound+testCutsize-1,lBound+testCutsize-1,lBound+testCutsize-1)
print *,  'u_f(1,1,1)' , u_f(1,lBound+testCutsize-1,lBound+testCutsize-1,lBound+testCutsize-1)
print *,  'u_t(1,1,1)' , u_t(1,lBound+testCutsize-1,lBound+testCutsize-1,lBound+testCutsize-1)
print*,''


if (d.eq.1) then

!! Compute tau_ij and T_ij:
allocate(tau_ij(n_uu,GRID,GRID,GRID))
allocate(T_ij(n_uu,GRID,GRID,GRID))
k = 1
do j=1,n_u
   do i=1,n_u
      if (i.ge.j) then
         Print *, 'tau(', i, ',', j, ')'
         tau_ij(k,:,:,:) = sharpFilter(u(i,:,:,:)*u(j,:,:,:),LES)     &
              - u_f(i,:,:,:)*u_f(j,:,:,:)
         Print *, 'T(', i, ',', j, ')'
         T_ij(k,:,:,:) = sharpFilter(u_f(i,:,:,:)*u_f(j,:,:,:),test)    &
              - u_t(i,:,:,:)*u_t(j,:,:,:)
         k = k+1
      end if
   end do
end do
deallocate(LES,test)


!! Print tau_ij and T_ij to check:
if (debug(3).eq.0)then
   dev_t = (tau_ij(1,200,156,129)+tau_ij(4,200,156,129)+tau_ij(6,200,156,129))/3.d0
   print*, 'w/ deviatoric:',tau_ij(4,200,156,129)-dev_t
else
  ! print*, tau_ij(4,200,156,129)
end if
!print*, T_ij(4,200,156,129)

!print*,'T_11(11,11,11)',T_ij(1,112,112,112)

! Take a cutout of the field(32x32x32)
! This is will result in a 16x16x16 test scale field
! This can be then safely used for a 8x8x8 field to find the h's
end if
print*,'Shape before cutout:',shape(u_t)
print*, u_t(1,testLower+lBound-1,testLower+lBound-1,testLower+lBound-1)
call cutout(u_t,n_u)
print*, u_t(1,testLower,testLower,testLower)
!call cutout(T_ij,n_u)

call synStress(u_t,n_u,n_uu)








end program main
