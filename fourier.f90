module fourier

integer,parameter :: GRID=32
integer, parameter :: LES_scale=16

contains

subroutine window(GRID,scale,filter)


! STATUS : > Test for FFTW layout.
!          > Passed test for center in spectral space 
!          > Check for the axes -- should be in alignment with FFTW output.
! Result :
! Notes  : 1) Create and define variables in a module. 
!          2) Needs to be parallelized
!      
! Layout: The center should lie at an unit offset from the geometric center.
!         This is because the DC component occupies the first index. 
! 
!  1             129          256 
!1  |-------------|------------|
!   |             |            |
!   |             |            |
!   |      A      |     B      |
!128|             |            |
!129--------------o------------|
!   |             |            |
!   |      D      |     C      |
!   |             |            |
!256--------------|------------|

    
implicit none

! Define arguments:
 integer, intent(in)  :: GRID,scale
 real,dimension(1:GRID,1:GRID,1:GRID),intent(inout):: filter

!! Define local variables: 
integer           :: center 
real              :: distance

real,allocatable,dimension(:,:,:) :: temp,A,B,C,D,E,F,G,H 

! Loop Indices:
  integer              :: i,j,k

! Initialize window :
  filter = 1.d0
  center = 0.5*GRID+1.d0

print *, center

! Create window:
  do k = 1,GRID
    do j = 1,GRID
      do i = 1,GRID
      
       distance = sqrt( real((i-center)**2) &
                   +real((j-center)**2)&
                   +real((k-center)**2) )
   
       if (distance.gt.scale) filter(i,j,k) = 0.d0
         
      end do
    end do
  end do

! FORMAT TO PRINT MATRIX:
!dims = shape(filter)
!    do i=1,dims(1)
!       write(*,20) int(filter(i,17,:)) ! Print the ith row
!    end do
! 20  format (33(i2))  

!Sanity checks:
print*, 'Sanity Checks:'
write(*,*) filter(center+scale+1,center,center) ! should be 0
write(*,*) filter(center,center-scale,center)   ! should be 1
write(*,*) filter(center+scale,center,center)   ! should be 1
write(*,*) filter(center,center,center+scale)   ! should be 1
write(*,*) filter(center,center-scale-1,center) ! should be 0

return
end subroutine window



subroutine fftshift(filter,GRID)

! fftshift:
! 
! Layout:                     k 
!             E       F       |
!           A       B         |
!                             /----> j
!             H       G      /
!           D       C       i
!  
!
! In clockwise direction from the topleft, 
! the front half(i = 129,256) has A,B,C,D
! The rear half (i = 1,128) is labelled as E,F,G,H
! For fftshift, the swaps needed are:

! Define arguments:
integer :: GRID
real,dimension(1:GRID,1:GRID,1:GRID),intent(inout):: filter

! Define local variables: 
real,allocatable,dimension(:,:,:) :: temp,A,B,C,D,E,F,G,H 
integer                           :: center


center = 0.5*GRID+1.d0

allocate(temp( 1:(center-1),1:(center-1),1:(center-1) ))
allocate(A( 1:(center-1),1:(center-1),1:(center-1) ) )
allocate(B( 1:(center-1),1:(center-1),1:(center-1) ) )
allocate(C( 1:(center-1),1:(center-1),1:(center-1) ) )
allocate(D( 1:(center-1),1:(center-1),1:(center-1) ) )
allocate(E( 1:(center-1),1:(center-1),1:(center-1) ) )
allocate(F( 1:(center-1),1:(center-1),1:(center-1) ) )
allocate(G( 1:(center-1),1:(center-1),1:(center-1) ) )
allocate(H( 1:(center-1),1:(center-1),1:(center-1) ) )

!Define subarrays:

A = filter( center:GRID, 1:(center-1), center:GRID )
B = filter( center:GRID, center:GRID, center:GRID )
C = filter( center:GRID, center:GRID, 1:(center-1))
D = filter( center:GRID, 1:(center-1), 1:(center-1))
E = filter( 1:(center-1), 1:(center-1), center:GRID)
F = filter( 1:(center-1), center:GRID, center:GRID)
G = filter( 1:(center-1), center:GRID , 1:(center-1) )
H = filter( 1:(center-1), 1:(center-1), 1:(center-1) )


! A - G:
temp = A
A = G
G = temp

! C - E:
temp = C
C = E
E = temp

! D - F:
temp = D
D = F
F = temp

! B - H:
temp = B
B = H
H = temp

! Recreate the filter matrix:

filter( center:GRID, 1:(center-1), center:GRID )  = A
filter( center:GRID, center:GRID, center:GRID )   = B
filter( center:GRID, center:GRID, 1:(center-1))   = C
filter( center:GRID, 1:(center-1), 1:(center-1))  = D
filter( 1:(center-1), 1:(center-1), center:GRID)  = E
filter( 1:(center-1), center:GRID, center:GRID)   = F
filter( 1:(center-1), center:GRID , 1:(center-1) )= G
filter( 1:(center-1), 1:(center-1), 1:(center-1) )= H


deallocate (temp,A,B,C,D,E,F,G,H)
return
end subroutine fftshift




! To print results in matrix view:
! Change format each time!!

subroutine matrixview(array)
implicit none

real, dimension(:,:,:), intent(inout)::array
integer :: i, dims(3)

dims = shape(array)
    do i=1,dims(1)
       write(*,20) int(array(i,1,:)) ! Print the ith row                                                              
    end do
20  format (33(i2))
return
end subroutine matrixview

end module fourier
