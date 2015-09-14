module fileio

contains

subroutine readfile(GRID,u1_DNS,u2_DNS,u3_DNS)
  
! STATUS : Tested on 9/11/2015
! Result : Passed with full precision of digits.

  implicit none
  
! Define arguments:
  integer,intent(in)  :: GRID
  real,dimension(1:GRID,1:GRID,1:GRID),intent(inout):: u1_DNS ,u2_DNS, u3_DNS

! Define local variables:

  integer                             :: i,j,k,  fID, position
  real,allocatable,dimension(:,:,:)   :: temp 

  ! fID -- is the counter for the file number.                                                                           ! fIndex -- char cast from int(fID).
        
  ! File variables:
  character*150 variableName, time, PATH, fIndex

  variableName='Velocity'
  time = '0460'
  PATH = '/Users/Kshitij/Desktop/ALES/DNS_Data/'


!! Allocate temp array
  allocate (temp(1:GRID,1:GRID,1:GRID))
  
do fID = 1,3

   print 10, fID
10 format("Reading... u",i1,"_DNS")
 
   ! Convert from int to char
   write(fIndex,'(i1)') fID

open(unit=fID, file = trim(PATH)//trim(variableName)//trim(fIndex)//'_'//trim(time)//'.bin', &
       status='old',form='unformatted',access='direct',convert='big_endian',recl=4)

    position = 0
    do k=1,GRID
       do j=1,GRID
          do i=1,GRID
             position = position+1
             read (fID,rec=position) temp(i,j,k)
          end do
       end do
    end do

    if (fID.eq.1) u1_DNS = temp/1.d2
    if (fID.eq.2) u2_DNS = temp/1.d2
    if (fID.eq.3) u3_DNS = temp/1.d2
  
close(fID)
end do

deallocate (temp)


end subroutine readfile
 
end module fileio
