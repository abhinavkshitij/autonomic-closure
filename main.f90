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

integer,parameter           :: LES_scale=64, test_scale=32
integer                     :: i,j,k,DIM

! Define velocities:

real(kind=8),allocatable,dimension(:,:,:,:) :: u, u_f, u_t
real(kind=8),allocatable,dimension(:,:,:,:) :: tau_ij,T_ij,temp
real(kind=8),allocatable,dimension(:,:,:) :: LES,test
real(kind=8):: dev_t
integer :: n_u, n_uu
integer,dimension(4) :: debug=(/0,1,1,1/)


!call randAlloc()
!stop


!! Set debug flags for velocity components:
n_u=3; n_uu = 6
if (debug(1).eq.1) n_u = 1
if (debug(1).eq.1) n_uu = 3



!! Select file to read:
allocate(u(n_u,GRID,GRID,GRID))
fileSelect:if (debug(2).eq.1) then
   call binRead(u,DIM=n_u)
else
   call hdf5Read() !! Under testing to read multiple files
end if fileSelect

print*,shape(u)
call cutout(u)

print *, u(3,1,1,1)
print *,shape(u)
stop

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
print *,'u', u(1,15,24,10)
print *,'u_f', u_f(1,15,24,10)
print *, 'u_t',u_t(1,15,24,10)
print*,''





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
if (debug(1).eq.0)then
   dev_t = (tau_ij(1,15,24,10)+tau_ij(4,15,24,10)+tau_ij(6,15,24,10))/3.d0
   print*, tau_ij(1,15,24,10)-dev_t
else
   print*, tau_ij(1,15,24,10)
end if
print*, T_ij(1,15,24,10)

! Take a cutout of the field(32x32x32)
! This is will result in a 16x16x16 test scale field
! This can be then safely used for a 8x8x8 field to find the h's










end program main
