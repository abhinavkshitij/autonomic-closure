program DNS

  use FFTW

  implicit none
  
  INTEGER,parameter:: GRID=256
  INTEGER,parameter :: m = GRID, n = GRID, p = GRID
  INTEGER ::i,j,k,count
  REAL,DIMENSION (1:m, 1:n ,1:p):: x

!!     Open direct access unformatted file with records sized for a single row
!!     You just need to make recl the smallest thing (4) you will read and adjust
!!     your read statements accordingly.

  
!! CAUTION: The value of x is NOT divided by 100

      
open (unit=1,file='../Velocity1_0460.bin',form='unformatted', access='direct',&
           convert='big_endian',recl=4)

count =0

!! READ DATA FROM BINARY FILES

do k=1,p
   do j=1,n
      do i=1,m
         count = count+1
         read (1,rec=count) x(i,j,k)
!! DEBUG MODE: write(*,*) i ,'', j , '' , k, '' ,x(i,j,k)/100
      end do
   end do
end do
   
close(1)
           
end program DNS
 
