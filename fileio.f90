module fileio  
  use HDF5  
  integer, parameter :: GRID=256
contains
  
subroutine binRead(u_d,set,DIM) 
 ! STATUS : Added capablity in matrixView() , hdf5read()
 ! Result : 
implicit none  

integer,intent(in)                         :: DIM
character(3),intent(in)                    :: set
real(8),dimension(DIM,GRID,GRID,GRID),intent(out)      :: u_d ! Data read in will be double precision kind.

real(4),dimension(:,:,:),allocatable       :: u_s ! Data is stored in single precision kind.
real(8),dimension(:,:,:),allocatable       :: u

integer                                    :: fID, position
character(LEN=32)                             PATH, variableName, time, fIndex

integer                                    :: i,j,k
logical                                    :: test = .TRUE.

! fID -- is the counter for the file number.
! fIndex -- char cast from int(fID).      

! CHANGE PATH HERE TO READ IN DATA:
! %% NRL:
if (set.eq.'nrl') then
   print('(a9,a5,/)'),'Dataset:',set

   PATH = '../data/nrl/'
   variableName ='Velocity'; time = '0460'

   allocate(u_s(GRID,GRID,GRID))
   u_s = 0.
   
   do fID = 1,DIM
      print 10,fID
10    format("Reading... u",i1,"_DNS")
      write(fIndex,'(i1)') fID

      open(unit=fID, file = trim(PATH)//trim(variableName)//trim(fIndex)//'_'//trim(time)//'.bin', &
                     status = 'old',                                                               &
                     form   = 'unformatted',                                                       &
                     access = 'direct',                                                            &
                     convert= 'big_endian',                                                        &
                     recl   =  4                                                                   )
      position = 0
      do k=1,GRID
      do j=1,GRID
      do i=1,GRID
         position = position + 1
         read (fID,rec = position) u_s(i,j,k)
      end do
      end do
      end do
      u_d(fID,:,:,:) = u_s*1.d-2 ! scale  original data by 1/100
      close(fID)
   end do

   deallocate(u_s)

! %% JHU:
elseif (set.eq.'jhu')then
   print('(a9,a5,/)'),'Dataset:',set

   variableName='Velocity'; time = '256'
   PATH = '../data/jhu256/bin/'

   allocate(u(GRID,GRID,GRID))
   u=0.d0

   do fID = 1,DIM
      print 20, fID
20    format("Reading... u",i1,"_DNS")
      write(fIndex,'(i1)') fID

      open(unit=fID, file = trim(PATH)//trim(variableName)//trim(fIndex)//'_'//trim(time)//'.bin', &
                     status = 'old',                                                               &
                     form   = 'unformatted',                                                       &
                     access = 'direct',                                                            &
                     recl   =  8                                                                   )
      position = 0
      do k=1,GRID
      do j=1,GRID
      do i=1,GRID
         position = position+1
         read (fID,rec=position) u(i,j,k)
      end do
      end do
      end do
      u_d(fID,:,:,:) = u
      close(fID)
   end do
   deallocate(u)

   !  TEST:
   if (test) then
      if (u_d(1,15,24,10).ne.-0.99597495794296265) then
            print*, 'Error reading data!'
            print*, u_d(1,15,24,10)
            print*, 'precision', precision(u_d(1,15,24,10))
           stop
         end if
   end if

! %% SIN:
elseif (set.eq.'sin')then
   print('(a9,a5,/)'),'Dataset:',set

   variableName='Velocity'; time = '256'
   PATH = '../data/sin3D/bin/'

   allocate(u(GRID,GRID,GRID))
   u=0.d0

   do fID = 1,DIM
      print 20, fID
      write(fIndex,'(i1)') fID

      open(unit=fID, file = trim(PATH)//trim(variableName)//trim(fIndex)//'_'//trim(time)//'.bin', &
                     status = 'old',                                                               &
                     form   = 'unformatted',                                                       &
                     access = 'direct',                                                            &
                     recl   =  8                                                                   )
      position = 0
      do k=1,GRID
      do j=1,GRID
      do i=1,GRID
         position = position+1
         read (fID,rec=position) u(i,j,k)
      end do
      end do
      end do
      u_d(fID,:,:,:) = u
      close(fID)
   end do
   deallocate(u)
print*, u_d(1,20,20,20)


! %% DEFAULT:
else
   print*,"Dataset must be either nrl or jhu"
   stop
end if
return
end subroutine binRead


subroutine hdf5Read()
implicit none

integer , parameter :: n_files=2 ! n is the number of files to be loaded.

character(LEN=20),parameter :: PATH = "../data/"
character(LEN=30),parameter :: fName = "isotropic1024coarse"
character(LEN=100) :: filename
character(LEN=4) :: L_Index, U_Index
character(LEN=10)  :: dset_u = "u10240", dset_p = "p10240", fL_Index, fU_Index  

integer(HID_T) :: file_id       ! File identifier
integer(HID_T) :: dset_id       ! Dataset identifier

!  real(kind=8), dimension(1,1024,1024,96) :: dset_data
real(kind=8), dimension(3,1024,1024,96*n_files) :: u ! Data buffers
integer(HSIZE_T), dimension(4) :: data_dims

integer     ::   error 
integer     ::   i, j, fCount
integer     ::   lowerIndex, upperIndex

!
! Initilize file indices:
!
lowerIndex = 0
upperIndex = 95
 
 
call h5open_f(error)   ! Begin HDF5

  do fCount=1,n_files
     if (fCount.eq.1) then
        write(L_Index,'(i1)') lowerIndex
        write(U_Index,'(i2)') upperIndex
     else if (fCount.eq.2) then
        write(L_Index,'(i2)') lowerIndex
        write(U_Index,'(i3)') upperIndex
     else if(fCount.ge.3.and.fCount.lt.10) then
        write(L_Index,'(i3)') lowerIndex
        write(U_Index,'(i3)') upperIndex
     else
        write(L_Index,'(i3)') lowerIndex
        write(U_Index,'(i4)') upperIndex
     end if
   
     filename = trim(PATH)//trim(fName)//trim(L_Index)//'_'//trim(U_Index)//'_10240.h5'
     
     call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error) ! Access file
     call h5dopen_f(file_id, dset_u, dset_id, error) ! Access dataset
     
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dset_data, data_dims, error) ! Read into dset_data
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, u(:,:,:,lowerIndex+1:upperIndex+1), data_dims, error) ! Read into dset_data 
     !u(1,:,:,lowerIndex+1:upperIndex+1) = dset_data
     
     call h5dclose_f(dset_id, error)    
     call h5fclose_f(file_id, error)

     lowerIndex = lowerIndex + 96
     upperIndex = upperIndex + 96
  end do
  call h5close_f(error)

  ! Check for bin files:
     print*,u(2,1,5,89) , u(3,5,5,89)
     print*,u(2,1,5,185), u(3,5,5,185)
  return
end subroutine hdf5Read


subroutine matrixview(array,frameLim,x,y,z,fID)
implicit none
  real(8), dimension(:,:,:), intent(in) :: array
  integer, intent(in), optional :: frameLim, fID, x, y, z
  integer :: i, lim, fileID
! Define flexible length format:                                                                                             
  character(4) :: form2

  ! Check for frameLimits:                                                                                                   
  lim = size(array,dim=1)
  if (present(frameLim))lim = frameLim
  write(form2,'(i4)') lim       ! Int to char                                                                                

  ! Check for fID:                                                                                                           
  fileID = 6 ! Default to dump results on the screen                                                                         
  if (present(fID)) fileID = fID


  ! Check for plane slice input:                                                                                       
  if(present(z)) then
     write(fileID,'('//trim(form2)//'f10.4)') , (array(i,1:lim,z),i=1,lim);print *, '' !z-plane                         
  elseif(present(x)) then
     write(fileID,'('//trim(form2)//'f10.4)') , (array(x,1:lim,i),i=1,lim);print *, '' !x-plane                        
  elseif(present(y)) then
     write(fileID,'('//trim(form2)//'f10.4)') , (array(1:lim,y,i),i=1,lim);print *, '' !y-plane                  
  else
     write(fileID,'('//trim(form2)//'f10.4)') , (array(i,1:lim,1),i=1,lim);print *, '' !z=1 default                     
  end if

  return
end subroutine matrixview

subroutine printmatrix2(matrix,M,N)
implicit none

integer,intent(in) :: M,N
real(8),dimension(M,N),intent(inout) :: matrix
character(4) :: form2
integer :: i

write(form2,'(i4)') N
write(*, '('//trim(form2)//'f30.15)'), ( matrix(i,1:N), i=1,M ); print*, ''
print*,''
end subroutine printmatrix2



end module fileio
