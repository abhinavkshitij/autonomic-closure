program testCut
  
use fileio
use linsolve
implicit none

integer                     :: i,j,k,DIM

! Define velocities:

real(kind=8),allocatable,dimension(:,:,:,:) :: u, u_f, u_t
real(kind=8),allocatable,dimension(:,:,:,:) :: tau_ij,T_ij,temp
real(kind=8),allocatable,dimension(:,:,:) :: LES,test
real(kind=8):: dev_t
integer :: n_u, n_uu
integer,dimension(4) :: debug=(/0,1,1,1/)


!! Set debug flags for velocity components:
n_u=3; n_uu = 6
if (debug(1).eq.1) n_u = 1
if (debug(1).eq.1) n_uu = 3


!! Select file to read:
allocate(u(n_u,GRID,GRID,GRID))


! Take a cutout of the field(32x32x32)
! This is will result in a 16x16x16 test scale field
! This can be then safely used for a 8x8x8 field to find the h's

call binRead(u,DIM=n_u)
print*,shape(u)
call cutout(u)

print *, u(3,1,1,1)
print *,shape(u)

!call randAlloc()
!stop


end program testCut

