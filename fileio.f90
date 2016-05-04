!****************************************************************
!                              FILEIO
!****************************************************************

!----------------------------------------------------------------
! USE:  1) Generate data(u_i, p) from DATA dir          : SOURCE  
!       2) Save intermediate results in TEMP dir        : BUFFER
!       3) [SAVE] final results in RUN dir              : SINK 
!
! FORM: module fileio
!          contains
!       subroutine readData            [SOURCE]
!           subroutine readBinSingle   [SOURCE]
!           subroutine readBinDouble   [SOURCE]
!       subroutine readHDF5            [SOURCE]
!       subroutine writeBin            [BUFFER]
!       subroutine check_dataRead      [TEST] 
!       subroutine contour             [SINK]
!
!
! BEHAVIOR: readData can read binary format dataset in single or
!           double precision. Each read operation is checked at 
!           probe point P(15,24,10).
!          
!
! STATUS : Refactoring this unit.
! 
!----------------------------------------------------------------


module fileio  
  use HDF5  
  use global

  interface readBin
     module procedure readBinSingle, readBinDouble
  end interface

contains
  
  !****************************************************************
  !                              READ_DATA
  !****************************************************************

  !----------------------------------------------------------------
  ! USE : Read binary datasets - NRL, JHU, TEST(sinusoidal). This
  !      is a source procudure and generates data -- u_i
  !
  ! FORM: 
  !
  ! BEHAVIOR:
  !
  ! STATUS : Refactoring this unit.
  ! 
  !----------------------------------------------------------------


  subroutine readData(u_i, DIM) 
    implicit none  
    !
    !    ..SCALAR ARGUMENTS..
    integer,intent(in)                         :: DIM
    !
    !    ..ARRAY ARGUMENTS..
    real(8),dimension(:,:,:,:),intent(out)     :: u_i
    !
    !    ..WORK ARRAY..
    real(4),dimension(:,:,:),allocatable       :: u_s ! Data is stored in single precision kind.
    real(8),dimension(:,:,:),allocatable       :: u_d
    !
    !    ..INTERNAL VARIABLES..
    integer                                    :: fID
    character(32)                              :: PATH, variableName, time, fIndex
    character(64)                              :: filename
    character(16)                              :: endian

    !    ..DEFAULT VALUES..
    PATH = trim(DATA_DIR) // trim(d_set) // '/' // trim(ext) // '/'
    variableName = trim(var(1)%name)! var(1) --> 'Velocity'; Expand for Pressure dataset.
    time = '256'
    endian = 'little_endian'

10  format("Reading... u",i1,"_DNS")

     !  READ SINGLE PRECISION DATA - NRL
     if (d_set.eq.'nrl') then
        time = '0460'
        endian = 'big_endian'

        allocate(u_s(i_GRID,j_GRID,k_GRID))
        do fID = 1,DIM
           write(fileID, 10) fID
           write(fIndex,'(i1)') fID
           filename = trim(PATH)//trim(variableName)//trim(fIndex)//'_'//trim(time)//'.'//trim(ext)
           call readBin(u_i,u_s,fID,filename,endian)
         end do
        deallocate(u_s)

        u_i = u_i * 1.d-2 ! scale  original data by 1/100

        !  READ DOUBLE PRECISION DATA - SIN3D, JHU256
     elseif (d_set.eq.'jhu256'.or.d_set.eq.'sin3D') then

        allocate(u_d(i_GRID,j_GRID,k_GRID))
        do fID = 1,DIM
           write(fileID, 10) fID
           write(fIndex,'(i1)') fID
           filename = trim(PATH)//trim(variableName)//trim(fIndex)//'_'//trim(time)//'.'//trim(ext)
           call readBin(u_i,u_d,fID,filename,endian)
        end do
        deallocate(u_d)

      elseif (d_set.eq.'jhu'.and.ext.eq.'h5') then
         call readHDF5()
      
        ! DEFAULT:
     else
        print*,"No dataset found by the name", d_set
        stop
     end if
        
     ! CHECK DATA READ:
     call check_dataRead(u_i(1,15,24,10))

     return
   end subroutine readData

   !****************************************************************
   !                           READ_BIN_SINGLE
   !****************************************************************

   !----------------------------------------------------------------
   ! USE : Reads single precision data for NRL 256 dataset
   !      
   ! FORM: readBinSingle(u_i, u_s, fID, filename, endian)
   !
   ! BEHAVIOR: Convert to big endian while reading in data.
   !
   ! STATUS : Interfaced with readBin()
   ! 
   !----------------------------------------------------------------

   subroutine readBinSingle(u_i,u_s,fID,filename,endian)
     implicit none
     !
     !    ..ARRAY ARGUMENTS..
     real(8), dimension(:,:,:,:),intent(inout) :: u_i
     real(4), dimension(:,:,:),intent(inout)   :: u_s
     !
     !    ..SCALAR ARGUMENTS..
     character(*), intent(in) :: filename
     character(*), intent(in),optional :: endian
     !
     !    ..LOCAL VARIABLES..
     integer :: fID
     integer :: position


     open(unit=fID, file = filename, &
          status = 'old',            &
          form   = 'unformatted',    &
          access = 'direct',         & 
          convert = endian,          & 
          recl   =  4                 )

     position = 0
     do k=1,k_GRID
        do j=1,j_GRID
           do i=1,i_GRID
              position = position + 1
              read (fID,rec=position) u_s(i,j,k)
           end do
        end do
     end do
     u_i(fID,:,:,:) = u_s
     close(fID)

   end subroutine readBinSingle

   !****************************************************************
   !                           READ_BIN_DOUBLE
   !****************************************************************

   !----------------------------------------------------------------
   ! USE : Reads double precision binary dataset.
   !      
   ! FORM: readBinDouble(u_i, u_d, fID, filename, endian)
   !
   ! BEHAVIOR: 
   !
   ! STATUS : Interfaced with readBin()
   ! 
   !----------------------------------------------------------------

   subroutine readBinDouble(u_i,u_d,fID,filename,endian)
     implicit none
     !
     !    ..ARRAY ARGUMENTS..
     real(8), dimension(:,:,:,:),intent(inout) :: u_i
     real(8), dimension(:,:,:),intent(inout)   :: u_d
     !
     !    ..SCALAR ARGUMENTS..
     character(*), intent(in) :: filename
     character(*), intent(in),optional :: endian
     !
     !    ..LOCAL VARIABLES..
     integer :: fID
     integer :: position


     open(unit = fID, file = filename, &
          status = 'old',            &
          form   = 'unformatted',    &
          access = 'direct',         & 
          convert = endian,          & 
          recl   =  8                 )

     position = 0
     do k=1,k_GRID
        do j=1,j_GRID
           do i=1,i_GRID
              position = position + 1
              read (fID,rec=position) u_d(i,j,k)
           end do
        end do
     end do
     u_i(fID,:,:,:) = u_d
     close(fID)

   end subroutine readBinDouble


   !****************************************************************
   !                            WRITEBIN
   !****************************************************************

   !----------------------------------------------------------------
   ! USE : Writes a 3D array into binary files in [TEMP] dir
   !      
   !
   ! FORM: 
   !
   ! BEHAVIOR: Needs allocated, defined arrays.
   !
   ! STATUS : 
   ! 
   !----------------------------------------------------------------

   subroutine writeBin(var, var_name, PATH)
     implicit none
     !
     !    ..ARRAY ARGUMENTS..
     real(8), dimension (:,:,:), intent(in) :: var
     !
     !    ..SCALAR ARGUMENTS..
     character(*), intent(in) :: var_name
     character(*), intent(in) :: PATH
     !
     !    ..LOCAL VARIABLES..
     character(16) :: scale
     character(32) :: filename
     integer :: fID

!     write(scale,'a3,2(i2),a1')

!     PATH = trim(TEMP_DIR) // trim(d_set) // '/' // trim(ext) // '/'
     filename = trim(PATH)


   end subroutine writeBin




   !****************************************************************
   !                           READHDF5
   !****************************************************************
   !----------------------------------------------------------------
   ! USE : Read HDF5 dataset JHU
   !      * source procudure, generates data (u_i)
   !
   ! FORM: 
   !
   ! BEHAVIOR:
   !
   ! STATUS : HDF5 files do not exist. They are saved in Ocotillo
   !          Need to create a test section.
   !          To handle large dataset - implement parallel I/O (MPI)
   !----------------------------------------------------------------

   subroutine readHDF5()
     implicit none

     integer, parameter :: n_files = 2 ! n is the number of files to be loaded.

     character(32) :: PATH 
     character(30) :: fName = "isotropic1024coarse"
     character(100) :: filename
     character(4) :: L_Index, U_Index
     character(10)  :: dset_u = "u10240", dset_p = "p10240", fL_Index, fU_Index  

     integer(HID_T) :: file_id       ! File identifier
     integer(HID_T) :: dset_id       ! Dataset identifier
     integer(HSIZE_T), dimension(4) :: data_dims
     !  real(kind=8), dimension(1,1024,1024,96) :: dset_data
     real(8), dimension(3,1024,1024,96*n_files) :: u ! Data buffers


     integer     ::   error 
     integer     ::   fCount
     integer     ::   lowerIndex, upperIndex

     !
     ! Initilize file indices:
     !
     lowerIndex = 0
     upperIndex = 95

     PATH = trim(DATA_DIR)//trim(fName)

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
        print*,filename
        stop
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
   end subroutine readHDF5

   !****************************************************************
   !                           CHECK_DATA_READ
   !****************************************************************

   !----------------------------------------------------------------
   ! USE : Writes a 3D array into binary files in [TEMP] dir
   !      
   !
   ! FORM: 
   !
   ! BEHAVIOR: Needs allocated, defined arrays.
   !
   ! STATUS : 
   ! 
   !----------------------------------------------------------------

  subroutine check_dataRead(value)
    implicit none
    real(8), intent(in) :: value
    
    ! NRL
    if (d_set.eq.'nrl') then
       if (value.ne.-11.070811767578125d0) then
          print*, 'Error reading data!'
          print*, value
          stop
       end if

    ! JHU256
    elseif (d_set.eq.'jhu256') then
       if (value.ne.-0.99597495794296265d0) then
          print*, 'Error reading data!'
          print*, value
          stop
       end if

    ! SIN 3D
    elseif (d_set.eq.'sin3D') then 
       if (value.ne.2.2787249947361117d0) then
          print*, 'Error reading data!'
          print*, value
          stop
       end if


    end if
    return
  end subroutine check_dataRead
  
  
  !****************************************************************
  !                           CONTOUR
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE : Saves a plane in .dat format (delimiter: space)
  !      
  !
  ! FORM: 
  !
  ! BEHAVIOR: Calls printplane() with fID to save in [RESULT] dir.
  !
  !
  ! STATUS : 
  ! 
  !----------------------------------------------------------------
  
  
  subroutine contour(var,var_name)
    implicit none

    !    ..ARRAY ARGUMENTS..
    real(8),dimension(:,:), intent(in) :: var
    character,intent(in) :: var_name

  end subroutine contour



  !****************************************************************
  !                           XYPLOT
  !****************************************************************

  !----------------------------------------------------------------
  ! USE : Saves XY data with option of multiple columns. Can be used
  !       for line plot, scatter plots and bar plots.
  !      
  !
  ! FORM: 
  !
  ! BEHAVIOR: Needs allocated, defined arrays.
  !
  ! STATUS : 
  ! 
  !----------------------------------------------------------------

  subroutine xyplot(var,var_name)
    implicit none

    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:), intent(in) :: var
    character, intent(in) :: var_name

  end subroutine xyplot


end module fileio
