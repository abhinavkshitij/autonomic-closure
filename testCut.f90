program testCut

use fileio
use linsolve
implicit none

integer                     :: i,j,k,DIM

! Define velocities:

real(kind=8),allocatable,dimension(:,:,:,:) :: u
integer :: n_u, n_uu
integer ::randMask(boxSize-coloc2)
integer,dimension(4) :: debug=(/1,0,1,1/)
integer,dimension(testcutSize,testcutSize,testcutSize):: test_cut
integer :: i_test, j_test, k_test
integer :: i_box, j_box, k_box
integer :: i_stencil,j_stencil,k_stencil
integer :: rand_count

!! Set debug flags for velocity components:
n_u=3; n_uu = 6
if (debug(1).eq.1) n_u = 1
if (debug(1).eq.1) n_uu = 3


!! Select file to read:
allocate(u(n_u,GRID,GRID,GRID))


! Take a cutout of the field(32x32x32)
! This is will result in a 16x16x16 test scale field
! This can be then safely used for a 8x8x8 field to find the h's

if (debug(2).eq.1) then
   call binRead(u,DIM=n_u)
   print*,shape(u)
   call cutout(u)

   print *, u(3,1,1,1)
   print *,shape(u)
   
end if

print*,boxCenter,smallHalf,bigHalf







! 2*(5+1)-1 for 8 = 11 ; 5 for box,1 for stencil
! As a modification, it should be independent of stride=2 , as well


!!$do k_test = 2*bigHalf+1, 2*(bigHalf-1+testSize)+1, 2 
!!$   do j_test = 2*bigHalf+1, 2*(bigHalf-1+testSize)+1, 2
!!$      do i_test = 2*bigHalf+1, 2*(bigHalf-1+testSize)+1, 2 ! i_test = 11,43,2

 do k_test = 64,64 
   do j_test = 64,64
      do i_test = 64,64
         
         test_cut = 0
         rand_count = 0
         
         do k_box = k_test-(2*(bigHalf-1)), k_test+(2*(box-bigHalf)),2
            do j_box = j_test-(2*(bigHalf-1)), j_test+(2*(box-bigHalf)),2
               do i_box = i_test-(2*(bigHalf-1)),i_test+(2*(box-bigHalf)),2 ! i_box = 3,49,2
                  call randAlloc(randMask)

                  rand_count = rand_count+1
                  if (any(randMask.eq.rand_count)) cycle
                    

                  test_cut(i_box,j_box,k_box) = 1
                  
                  do k_stencil = k_box-2,k_box+2,2
                     do j_stencil = j_box-2,j_box+2,2
                        do i_stencil = i_box-2,i_box+2,2
                           
                           test_cut(i_stencil,j_stencil,k_stencil) = 2
                          
                        end do
                     end do
                  end do
                  
                  
                  
               end do
            end do
         end do
         

        
         
      end do
   end do
end do

print*, randMask
print*, rand_count
print*,test_cut(51,51,51), test_cut(53,53,53)

         
print*,''

print*,'Last limit of test array:',i_test-2 
print*,'Last limit of bounding box:',i_box-2
print*,'Last limit of stencil:',i_stencil-2






end program testCut

