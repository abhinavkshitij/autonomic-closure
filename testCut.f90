program testCut

use fileio
use linsolve
implicit none

integer                     :: i,j,k,p,DIM

! Define velocities:

real(kind=8),allocatable,dimension(:,:,:,:) :: u,uu
real(kind=8),dimension(coloc2,coloc2) :: A

integer :: n_u
integer :: n_uu
integer :: randMask(boxSize-coloc2) !512-378=134; 512-270=240

integer :: debug(4)=(/0,1,1,0/) 
integer :: test_cut(testcutSize,testcutSize,testcutSize)
integer :: i_test,   j_test,   k_test
integer :: i_box,    j_box,    k_box
integer :: i_stencil,j_stencil,k_stencil
integer :: rand_count,lim
integer :: row_index, col_index
integer :: u_comp, uu_comp

!! Set debug flags for velocity components:
n_u=3; n_uu = 6
if (debug(1).eq.1) n_u = 1
if (debug(1).eq.1) n_uu = 1


!! Select file to read:
allocate(u(n_u,GRID,GRID,GRID))



! Take a cutout of the field(32x32x32)
! This is will result in a 16x16x16 test scale field
! This can be then safely used for a 8x8x8 field to find the h's

if (debug(2).eq.1) then
   
   call binRead(u,DIM=n_u)
   print*,shape(u)
   call cutout(u,n_u)

   print *, u(n_u,1,1,1)
   print *,'shape u cutout: ',shape(u)

   allocate(uu(n_uu,testcutSize,testcutSize,testcutSize))
   print*,'shape uu cutout:',shape(uu)

   ! Compute velocity products:
   k=0
   do j=1,n_u
      do i=1,n_u
         if (i.ge.j) then
            k=k+1
            uu(k,:,:,:) = u(i,:,:,:)*u(j,:,:,:)
         end if       
      end do
   end do
      
end if

print*,'boxCenter index:',boxCenter
print*,'smallHalf,bighalf:',smallHalf,bigHalf



! 2*(5+1)-1 for 8 = 11 ; 5 for box,1 for stencil
! As a modification, it should be independent of stride=2 , as well

if(debug(3).eq.1)then
   lim=testUpper
else
   lim=testLower
end if

!!$do k_test = testLower, testUpper, stride 
!!$   do j_test = testLower, testUpper, stride
!!$      do i_test = testLower, testUpper, stride ! i_test = 11,43,2


 do k_test = lim,lim,stride
   do j_test = lim,lim,stride
      do i_test = lim,lim,stride
         
         A   = 0   
         rand_count = 0 
         row_index  = 0 ! Reset pointer to the first position after moving
                        ! to the next cell.

         call randAlloc(randMask)
         do k_box = k_test-boxLower, k_test+boxUpper,stride
            do j_box = j_test-boxLower, j_test+boxUpper,stride
               do i_box = i_test-boxLower,i_test+boxUpper,stride ! i_box = 3,49,2
                  
                  rand_count = rand_count+1
                  if (any(randMask.eq.rand_count)) cycle   ! Skip if the point is listed in randMask
                    

                  ! test_cut(i_box,j_box,k_box) = 1
                  
                  col_index = 0 ! Reset pointer to the first column after each stencil operation
                  row_index = row_index + 1
                  
                  do k_stencil = k_box-stride,k_box+stride,stride
                     do j_stencil = j_box-stride,j_box+stride,stride
                        do i_stencil = i_box-stride,i_box+stride,stride
                           
                           test_cut(i_stencil,j_stencil,k_stencil) = 2 ! Flag in stencil-->check:first cell
                           
                           
                           do u_comp = 1,n_u ! 1 to 3
                              col_index = col_index+1
                              A(row_index,col_index)=u(u_comp,i_stencil,j_stencil,k_stencil)                           
                           end do
                           do uu_comp = 1,n_uu ! 1 to 6
                              col_index = col_index+1
                              A(row_index,col_index)=uu(uu_comp,i_stencil,j_stencil,k_stencil)
                           end do
                           
                        
                        end do
                     end do
                  end do
                  
                  
                  
                  
               end do
            end do
         end do
         

        
         
      end do
   end do
end do

print*,''
print*, 'Exclusion list:'
print*, randMask
print*,''
print *, 'Number of random cells sampled:',rand_count
print*, 'Check for first cell:',test_cut(1,1,1)
print*, 'Check for last cell:',test_cut(testcutSize,testcutSize,testcutSize)
        
print*,''

print*,'Last limit of test array:',i_test-stride 
print*,'Last limit of bounding box:',i_box-stride
print*,'Last limit of stencil:',i_stencil-stride
print*,''

print*,uu(6,51,51,51),A(coloc2,coloc2) !--> check for the last cell 
print*,uu(6,51,49,51),A(coloc2,coloc2) !--> check for the last cell

if(debug(4).eq.1) then
   
print*,'Testing for A:'

k=0 ; p=0
do k_stencil = 3-stride,3+stride,stride
   do j_stencil = 3-stride,3+stride,stride
      do i_stencil = 3-stride,3+stride,stride
         p=p+1
         print*,p

do i=1,n_u
      k=k+1
      print*, 'u',i,i_stencil,j_stencil,k_stencil,u(i,i_stencil,j_stencil,k_stencil),'A:',k,A(1,k)
   end do

   do i=1,n_uu
      k=k+1
      print*, 'uu',i,i_stencil,j_stencil,k_stencil,uu(i,i_stencil,j_stencil,k_stencil),'A:',k,A(1,k)
   end do
   print*,''
         
         
      end do
   end do
end do

end if

   call condition(A)                     
                           
                           



end program testCut

