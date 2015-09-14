program dataread      
implicit none

integer                         	:: npos, i, j, k
integer,dimension(3)            	:: nx
real,dimension(:,:,:),allocatable 	:: var

character*150 varn,time

!Set variable, time, and grid dimensions
varn='Velocity1'
time='0300'
nx(1)=512
nx(2)=512
nx(3)=512

!Allocate variable array
allocate(var(nx(1),nx(2),nx(3)))

!Open data file and read using stream access
open(1,file=trim(varn)//'_'//trim(time)//'.bin', &
       status='old',access='stream',form='unformatted',convert='big_endian')
npos=1
do k=1,nx(3)
   do j=1,nx(2)
      do i=1,nx(1)
      	 read(1,pos=npos) var(i,j,k)
         npos=npos+4
      enddo !i
   end do !j
end do !k
close(1)
write(*,*)'Opened file:',trim(varn)//'/'//trim(varn)//'_'//trim(time)//'.bin'

!Write some data to check
write(*,*) var(3,15,30)

end program dataread