module fileio  
  use HDF5  
  integer, parameter :: GRID=256

contains
  
subroutine binRead(u_dpk,DIM) 
 ! STATUS : Added capablity in matrixView() , hdf5read()
 ! Result : 
  implicit none
  
  
  integer,intent(in)                         :: DIM
  real,dimension(:,:,:),allocatable          :: u_spk ! Data is stored in single precision kind.
  real(kind=8),dimension(DIM,GRID,GRID,GRID) :: u_dpk ! Data read in will be double precision kind.
  integer                                    :: fID, position
  integer :: i,j,k
  character(LEN=16) variableName, time, fIndex
  character(LEN=100)PATH
  
  ! fID -- is the counter for the file number.
  ! fIndex -- char cast from int(fID).      

   
  variableName='Velocity'; time = '0460'
  PATH = '/Users/Kshitij/Desktop/ALES/DNS_Data/'

  allocate(u_spk(GRID,GRID,GRID))
  
do fID = 1,DIM
   print 10, fID
10 format("Reading... u",i1,"_DNS")
   write(fIndex,'(i1)') fID
   u_spk = 0.
   open(unit=fID, file = trim(PATH)//trim(variableName)//trim(fIndex)//'_'//trim(time)//'.bin', &
       status='old',form='unformatted',access='direct',convert='big_endian',recl=4)

    position = 0
    do k=1,GRID
       do j=1,GRID
          do i=1,GRID
             position = position+1
             read (fID,rec=position) u_spk(i,j,k)
          end do
       end do
    end do
  u_dpk(fID,:,:,:) = u_spk*1d-2 ! Multiply original data by 0.01.
close(fID)
end do
deallocate(u_spk)

return
end subroutine binRead


subroutine hdf5Read()
  implicit none

  integer , parameter :: n=3 ! n is the number of files to be loaded.

  character(LEN=20),parameter :: PATH = "./h5/"
  character(LEN=30),parameter :: fName = "isotropic1024coarse"
  character(LEN=100) :: filename
  character(LEN=4) :: L_Index, U_Index
  character(LEN=10)  :: dset_u = "u10240", dset_p = "p10240", fL_Index, fU_Index  

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id       ! Dataset identifier

!  real(kind=8), dimension(1,1024,1024,96) :: dset_data
  real(kind=8), dimension(1:3,1024,1024,96*n) :: u ! Data buffers
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

  ! I cannot find a way to read multiple files right now.
  ! So I am fixing it by using an IF condition.

  do fCount=1,n
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
     print*,u(1,1,1,6)
     
     call h5dclose_f(dset_id, error)    
     call h5fclose_f(file_id, error)

     lowerIndex = lowerIndex + 96
     upperIndex = upperIndex + 96
  end do
  call h5close_f(error)
  return
end subroutine hdf5Read

subroutine matrixview(array,frameLim,x,y,z,fID)
implicit none
  real(8), dimension(:,:,:), intent(inout) :: array
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

end module fileio
