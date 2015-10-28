module fileio  
  use HDF5  
  integer, parameter :: GRID=256
contains
  
subroutine binRead(u_dpk,DIM) 
 ! STATUS : Tested on 10/23/2015
 ! Result : Passed with full precision of digits. Added DIM keyword for debugging.  
  implicit none
  
  
  integer,intent(in)                         :: DIM
  real,dimension(:,:,:),allocatable          :: u_spk
  real(kind=8),dimension(DIM,GRID,GRID,GRID) :: u_dpk
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
  u_dpk(fID,:,:,:) = u_spk*1d-2
close(fID)
end do
deallocate(u_spk)

return
end subroutine binRead


subroutine hdf5Read()
  implicit none

  integer , parameter :: n=2
  character(LEN=100),parameter :: filename="./h5/isotropic1024coarse0_95_10240.h5"
  character(LEN=10)  :: dsetname = "u10240" , fL_Index, fU_Index  

  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id       ! Dataset identifier

  integer     ::   error ! Error flag
  integer     ::   i, j
  integer     :: lowerIndex, upperIndex

  real(kind=8), dimension(1:3,1024,1024,96) :: dset_data
  !real(kind=8), dimension(1:3,1024,1024,96*n) ::data_out ! Data buffers
  integer(HSIZE_T), dimension(4) :: data_dims
 
  call h5open_f(error)   
     call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
     call h5dopen_f(file_id, dsetname, dset_id, error)
     
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dset_data, data_dims, error)
     print*,dset_data(1,1,1,6)
     
     call h5dclose_f(dset_id, error)    
     call h5fclose_f(file_id, error)
  call h5close_f(error)
  return
end subroutine hdf5Read


subroutine matrixview(array,frameLim,fID)
implicit none

real(kind=8), dimension(:,:,:), intent(inout)::array
integer, intent(in), optional :: frameLim, fID
integer :: i, lim, fileID

! Define flexible length format:
  character(1) :: form1
  character(4) :: form2
  character(6) :: form3

  ! Check for frameLimits:
  lim = size(array,dim=1)
  if (present(frameLim)) lim = frameLim

  ! Check for fID:
  fileID = 6 ! Default to dump results on the screen
  if (present(fID)) fileID = fID
  
  form1 = '('; write(form2,'(i4)') lim ; form3 = 'f10.4)'
  write(fileID,form1//trim(form2)//form3) , (array(i,1:lim,1),i=1,lim);print *, ''
  
  return
  end subroutine matrixview

  
  

end module fileio
