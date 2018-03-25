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
!       subroutine check_dataRead      [TEST]
!       subroutine plotVelocities      [SINK]
!       subroutine loadFFT_data        [BUFFER]
!       subroutine saveFFT_data        [BUFFER]
!       subroutine plotOriginalStress  [SINK]
!       subroutine plotProductionTerm  [SINK]
!       subroutine plotPDF             [SINK]
!       subroutine contour             [SINK]
!       subroutine xyplot              [SINK]
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
  !      is a source procedure and generates data -- u_i
  !
  ! FORM: 
  !
  ! BEHAVIOR: u_s: for binary read for single precision files.
  !
  ! STATUS : Refactoring this unit.
  !       ** TEMPORARY CHANGE BASED ON PROF HAMLINGTON'S MAIL
  !----------------------------------------------------------------


  subroutine readData(DIM) 
    implicit none  
    !
    !    ..SCALAR ARGUMENTS..
    integer,intent(in)                         :: DIM
    !
    !    ..WORK ARRAY..
    real(4),dimension(:,:,:), allocatable      :: u_s !  single precision kind
    !
    !    ..INTERNAL VARIABLES..
    integer                                    :: fID
    character(32)                              :: PATH, variableName, fIndex
    character(64)                              :: filename
    character(16)                              :: endian

    !    ..DEFAULT VALUES..
    DATA_PATH = trim(DATA_DIR) // trim(dataset) // '/' // trim(ext) // '/'
    variableName = trim(var(1)%name)! var(1) --> 'Velocity'; Expand for Pressure dataset.
    endian = 'little_endian'

10  format("Reading... u",i1,"_DNS")
11  format("Reading... p_DNS")

     !  READ SINGLE PRECISION DATA - NRL
     if (dataset.eq.'nrl') then
        endian = 'big_endian'

        allocate(u_s(i_GRID,j_GRID,k_GRID))
        
        do fID = 1,DIM
           write(fileID, 10) fID
           write(fIndex,'(i1)') fID
           filename = trim(DATA_PATH)//trim(variableName)//trim(fIndex)//'_'//trim(time)//'.'//trim(ext)
           call readBin(u_s,fID,filename,endian)
         end do
        deallocate(u_s)

        ! scale  original data by 1/100
        u = u * 1.d-2 


     !  READ SINGLE PRECISION DATA - HST
     elseif (dataset.eq.'hst') then
        DATA_PATH = trim(DATA_DIR) // trim(dataset) // '/' // trim(ext) // '/' // trim(hst_set) // '/'

        allocate(u_s(i_GRID, bigHalf(j_GRID),k_GRID))
        do fID = 1,DIM
           write(fileID, 10) fID
           write(fIndex,'(i1)') fID
           filename = trim(DATA_PATH)//trim(variableName)//trim(fIndex)//'_'//trim(time)//'.'//trim(ext)
!           print*, filename
           call readBin(u_s,fID,filename,endian)
         end do
        deallocate(u_s)

        
        !  READ DOUBLE PRECISION DATA - SIN3D, JHU256
     elseif (dataset.eq.'jhu256'.or.dataset.eq.'sin3D') then

        if (allocated(u).eqv..false.) allocate(u(n_u, i_GRID, j_GRID, k_GRID))
        do fID = 1,DIM
           write(fIndex,'(i1)') fID
           if (fID.eq.4) then
              write(fileID, 11) 
              variableName = trim(var(2)%name) ! Pressure term
           else
              write(fileID, 10) fID
           end if
           filename = trim(DATA_PATH)//trim(variableName)//trim(fIndex)//'_'//trim(time)//'.'//trim(ext)
           call readBin(fID,filename,endian)
        end do

        ! READ DOUBLE PRECISION DATA - JHU1024 - HDF5
      elseif (dataset.eq.'jhu1024'.and.ext.eq.'h5') then
         call readHDF5(n_files=1)
      
        ! DEFAULT:
     else
        print*,"No dataset found by the name", dataset
        stop
     end if
        
     ! CHECK DATA READ:
     call check_dataRead(u(1,15,24,10))

     return
   end subroutine readData

   !****************************************************************
   !                           READ_BIN_SINGLE
   !****************************************************************

   !----------------------------------------------------------------
   ! USE : Reads single precision data for NRL 256 dataset
   !      
   ! FORM: readBinSingle(u, u_s, fID, filename, endian)
   !
   ! BEHAVIOR: Convert to big endian while reading in data.
   !
   ! STATUS : Interfaced with readBin()
   !          ## Commented out j_GRID; replaced with 129 for HST.
   !         ### Move back to (:,:,:)
   !----------------------------------------------------------------

   subroutine readBinSingle(u_s,fID,filename,endian)
     implicit none
     !
     !    ..ARRAY ARGUMENTS..
     real(4), dimension(:,:,:),intent(inout)   :: u_s
     !
     !    ..SCALAR ARGUMENTS..
     character(*), intent(in) :: filename
     character(*), intent(in), optional :: endian
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
        do j=1,129! ##  j_GRID
           do i=1,i_GRID
              position = position + 1
              read (fID,rec=position) u_s(i,j,k)
           end do
        end do
     end do
     u(fID,1:256,1:129,1:256) = u_s(1:256,1:129,1:256) ! ### move back to (:,:,:)
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

   subroutine readBinDouble(fID,filename,endian)
     implicit none
     !
     !    ..SCALAR ARGUMENTS..
     integer,intent(in) :: fID
     character(*), intent(in) :: filename
     character(*), intent(in),optional :: endian
     !
     !    ..LOCAL VARIABLES..
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
              read (fID,rec=position) u(fID,i,j,k)
           end do
        end do
     end do
     close(fID)

   end subroutine readBinDouble




   !****************************************************************
   !                           READHDF5
   !****************************************************************
   !----------------------------------------------------------------
   ! USE : Read HDF5 dataset JHU
   !      * source procudure, generates data (u_i)
   !
   ! FORM:    subroutine readHDF5(u, [n_files])
   !
   ! BEHAVIOR: If no n_files is specified then it reads all .h5 files
   !           Allocates u(:,:,:,:) accessed from GLOBAL.
   !
   !
   ! STATUS : HDF5 files do not exist. They are saved in Ocotillo
   !          Need to create a test section.
   !          To handle large dataset - implement parallel I/O (MPI)
   ! STASHED LINES:
   !     # real(8), dimension(3,1024,1024,96) :: dset_data
   !     ...
   !    + call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dset_data, data_dims, error) 
   !   ++ u(1,:,:,lowerIndex+1:upperIndex+1) = dset_data(1,:,:,:)
   !----------------------------------------------------------------

   subroutine readHDF5(n_files)
     implicit none
     !
     !    ..SCALAR ARGUMENTS..
     integer, optional, intent(in)  :: n_files  
     !
     !    ..LOCAL VARS..
     character(30) :: fName = "isotropic1024coarse"
     !character(32) :: PATH 
     character(100) :: filename
     !#
     !    ..HDF5 VARS..
     character(10)  :: dset_u = "u10240" 
     character(10)  :: dset_p = "p10240"
     integer(HID_T) :: file_id       ! POINTER TO FILE
     integer(HID_T) :: dset_id       ! POINTER TO DATASET
     integer(HSIZE_T), dimension(4) :: data_dims
     !
     !    ..INDICES..
     integer     ::   error 
     integer     ::   fCount
     integer     ::   lowerIndex, upperIndex
     character(4)::   L_Index,    U_Index
     !
     !SET PATH:
     DATA_PATH = trim(DATA_DIR) // trim(dataset) // '/' // trim(ext) // '/'
     !
     ! Initialize file indices:
     !
     lowerIndex = 0
     upperIndex = 95
     !
     ! ALLOCATE ARRAY SIZE:
     if (present(n_files)) then

        if (n_files.le.10)then
           allocate(u(n_u, i_GRID,j_GRID, 96*n_files))
        end if
     else
        allocate(u(n_u, i_GRID, j_GRID, k_GRID))
     end if

        ! HDF5_INIT
        call h5open_f(error)  

        ! SET INDICES
        do fCount=1,n_files
           if (fCount.eq.1) then
              write(L_Index,'(i1)') lowerIndex
              write(U_Index,'(i2)') upperIndex
           else if (fCount.eq.2) then
              write(L_Index,'(i2)') lowerIndex
              write(U_Index,'(i3)') upperIndex
           else if(fCount.ge.3.and.fCount.le.10) then
              write(L_Index,'(i3)') lowerIndex
              write(U_Index,'(i3)') upperIndex
           else
              write(L_Index,'(i3)') lowerIndex
              write(U_Index,'(i4)') upperIndex
           end if

           filename = trim(DATA_PATH)//trim(fName)//trim(L_Index)//'_'//trim(U_Index)//'_10240'//'.'//trim(ext)
           print*, filename

           ! ACCESS [F]ILE
           call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error) 

           ! ACCESS [D]ATASET
           call h5dopen_f(file_id, dset_u, dset_id, error) 

           ! LOAD DATA
           ! +     
           call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, u(:,:,:,lowerIndex+1:upperIndex+1), data_dims, error) 
           ! ++ 
           ! ALTERNATE METHOD:READ INTO DSET_DATA. STASHED ABOVE.
           call h5dclose_f(dset_id, error)    
           call h5fclose_f(file_id, error)

           lowerIndex = lowerIndex + 96
           upperIndex = upperIndex + 96
           if (fCount == 10) upperIndex = 1023
        end do

        ! HDF5_FINALIZE
        call h5close_f(error)

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
    if (dataset.eq.'nrl') then
       if (value.ne.-11.070811767578125d0) then
          print*, 'Error reading data!'
          print*, value
          stop
       end if

    ! JHU256
    elseif (dataset.eq.'jhu256') then
       if (value.ne.-0.99597495794296265d0) then
          print*, 'Error reading data!'
          print*, value
          stop
       end if

    ! SIN 3D
    elseif (dataset.eq.'sin3D') then 
!       if (value.ne.2.2787249947361117d0) then
       if (value.ne.5.9589822233083432d0) then
          print*, 'Error reading data!'
          print*, value
!          stop
       end if

       !
    elseif (dataset.eq.'jhu1024') then 
       if (value.ne.0.41232076287269592d0) then
          print*, 'Error reading data!'
          print*, value
!          stop
       end if

    end if
    return
  end subroutine check_dataRead
  
   !****************************************************************
   !                            PLOT VELOCITIES
   !****************************************************************

   !----------------------------------------------------------------
   ! USE : Saves [z-midplane] in RESULTS directory
   !      
   !
   ! FORM:    subroutine plotVelocities()
   !
   ! BEHAVIOR: Needs allocated, defined arrays.
   !
   ! STATUS :
   !      
   !   
   !----------------------------------------------------------------
   
   subroutine plotVelocities(plotOption)
     implicit none
     !
     !    ..SCALAR ARGUMENTS..
     character(*), optional, intent(in) :: plotOption
     !
     !    ..LOCAL VARIABLES..
     character(64) :: filename
     integer :: i, n_ij
     character(1) :: ij


     n_ij = 1
     if (plotOption.eq.'All') n_ij = 3
     !
     ! SAVE VELOCITIES:
     print*
     print*,'Saving filtered velocities in', RES_PATH

     do i=1, n_ij
        write(ij, '(i0)') i
        open(20,file=trim(RES_PATH)//'u_i'//trim(ij)//'.dat')
        open(21,file=trim(RES_PATH)//'u_f'//trim(ij)//'.dat')
        open(22,file=trim(RES_PATH)//'u_t'//trim(ij)//'.dat')
        write(20,*) u  (i,:,:,z_plane)
        write(21,*) u_f(i,:,:,z_plane)
        write(22,*) u_t(i,:,:,z_plane)
        close(20)
        close(21)
        close(22)
     end do

end subroutine plotVelocities


   !****************************************************************
   !                            PLOT PRESSURE
   !****************************************************************

   !----------------------------------------------------------------
   ! USE : Saves [z-midplane] in RESULTS directory
   !      
   !
   ! FORM:    subroutine plotPressure()
   !
   ! BEHAVIOR: Needs allocated, defined arrays.
   !
   ! STATUS :
   !      
   !   
   !----------------------------------------------------------------
   
   subroutine plotPressure()
     implicit none
     !
     ! SAVE PRESSURE:
     print*
     print*,'Saving filtered pressure in', RES_PATH
          
     open(30,file=trim(RES_PATH)//'p.dat')
     open(31,file=trim(RES_PATH)//'p_f.dat')
     open(32,file=trim(RES_PATH)//'p_t.dat')
     write(30,*) u  (4,:,:,z_plane)
     write(31,*) u_f(4,:,:,z_plane)
     write(32,*) u_t(4,:,:,z_plane)
     close(30)
     close(31)
     close(32)

   end subroutine plotPressure


   !****************************************************************
   !                            PLOT VORTICITY
   !****************************************************************

   !----------------------------------------------------------------
   ! USE : Saves [z-midplane] in RESULTS directory
   !      
   !
   ! FORM:    subroutine plotVorticity()
   !
   ! BEHAVIOR: Needs allocated, defined arrays.
   !
   ! STATUS :
   !      
   !   
   !----------------------------------------------------------------
   
   subroutine plotVorticity()
     implicit none
     !
     !   ..COUNTER..
     integer :: i
     !
     ! SAVE VORTICITY:
     print*
     print*,'Saving vorticity(omega) in', RES_PATH
    

     open(30,file=trim(RES_PATH)//'omega1.dat')
     open(31,file=trim(RES_PATH)//'omega2.dat')
     open(32,file=trim(RES_PATH)//'omega3.dat')

     write(30,*) omega (1,:,:,z_plane)
     write(31,*) omega (2,:,:,z_plane)
     write(32,*) omega (3,:,:,z_plane)

     close(30)
     close(31)
     close(32)
   end subroutine plotVorticity
   



   !****************************************************************
   !                          LOAD FFT DATA
   !****************************************************************

   !----------------------------------------------------------------
   ! USE : Loads a 3D array into binary files in [TEMP] dir
   !      
   !
   ! FORM:   subroutine loadFFT_data()
   !
   ! BEHAVIOR: Needs allocated, defined arrays.
   !
   ! STATUS : 
   ! 
   !----------------------------------------------------------------

  subroutine loadFFT_data()
    implicit none

     !
     !    ..LOCAL VARIABLES..
     character(64) :: filename
     integer :: i
  
     print*
     print*,'Load filtered variables ... '
     do i = 1,size(var_FFT)
        if (i.eq.1.or.i.eq.3) then
        filename = trim(TEMP_PATH)//trim(LESfilterType)//'/'//trim(var_FFT(i)%name)//'.bin'
      else
        filename = trim(TEMP_PATH)//trim(LESfilterType)//'/'//trim(TestfilterType)//'/'//trim(var_FFT(i)%name)//'.bin'
      endif
        
      print*, filename
      open(i, file = filename,form='unformatted')
      if (i.eq.1) read(i) u_f
      if (i.eq.2) read(i) u_t
      if (i.eq.3) read(i) tau_ij
      if (i.eq.4) read(i) T_ij
      close(i)
     end do

   end subroutine loadFFT_data

   !****************************************************************
   !                     LOAD FFT DATA - 3 FILTER
   !****************************************************************

   !----------------------------------------------------------------
   ! USE : Loads u_f, tau_ij from LES/ dir and u_t, T_ij from test/ 
   !      dir
   !
   ! FORM:   subroutine loadFFT_data3()
   !
   ! BEHAVIOR: Needs allocated, defined arrays. Is a brute force method 
   !           of reading a file. 
   !
   ! STATUS : 
   ! 
   !----------------------------------------------------------------

  subroutine loadFFT_data3(LESFilterType)
    implicit none
    !
    !    ..ARGUMENTS..
    character(*),intent(in), optional :: LESFilterType
    !
    !    ..LOCAL VARIABLES..
    character(64) :: filename
    character(8) :: Testfilter1, Testfilter2, Testfilter3
    integer :: i

    select case  (TestfilterType)
        case ('GaussBox')
          Testfilter1 = 'Gauss'
          Testfilter2 = 'Box'
        case ('GaussTri')
          Testfilter1 = 'Gauss'
          Testfilter2 = 'Tri'
        case ('BoxTri')
          Testfilter1 = 'Box'
          Testfilter2 = 'Tri'
        case ('All')
          Testfilter1 = 'Gauss'
          Testfilter2 = 'Box'
          Testfilter3 = 'Tri'
    end select

    print*
    print*,'Load filtered variables ... '

     ! READ LES DATA - u_f, tau_ij:
    filename = trim(TEMP_PATH)//trim(LESFilterType)//'/'//trim('u_f.bin')
    print*, filename
    open(1, file = filename,form='unformatted')
    read(1) u_f
    close(1)

    filename = trim(TEMP_PATH)//trim(LESFilterType)//'/'//trim('tau_ij.bin')
    print*, filename
    open(2, file = filename,form='unformatted')
    read(2) tau_ij
    close(2)

    ! READ TEST DATA - u_t, T_ij:
    filename = trim(TEMP_PATH)//trim(LESFilterType)//'/'//trim(Testfilter1)//'/'//trim('u_t.bin')
    print*, filename
    open(1, file = filename,form='unformatted')
    read(1) u_t
    close(1)

    filename = trim(TEMP_PATH)//trim(LESFilterType)//'/'//trim(Testfilter1)//'/'//trim('T_ij.bin') 
    print*, filename
    open(2, file = filename,form='unformatted')
    read(2) T_ij
    close(2)
    
    filename = trim(TEMP_PATH)//trim(LESFilterType)//'/'//trim(Testfilter2)//'/'//trim('u_t.bin')
    print*, filename
    open(1, file = filename,form='unformatted')
    read(1) u_tB
    close(1)

    filename = trim(TEMP_PATH)//trim(LESFilterType)//'/'//trim(Testfilter2)//'/'//trim('T_ij.bin')
    print*, filename
    open(2, file = filename,form='unformatted')
    read(2) T_ijB
    close(2)
    
    ! filename = trim(TEMP_PATH)//trim(LESFilterType)//'/'//trim(Testfilter3)//'/'//trim('u_t.bin')
    ! print*, filename
    ! open(1, file = filename,form='unformatted')
    ! read(1) u_tG
    ! close(1)

    ! filename = trim(TEMP_PATH)//trim(LESFilterType)//'/'//trim(Testfilter3)//'/'//trim('T_ij.bin')
    ! print*, filename
    ! open(2, file = filename,form='unformatted')
    ! read(2) T_ijG
    ! close(2)
              
      

   end subroutine loadFFT_data3


   !****************************************************************
   !                          SAVE FFT_DATA
   !****************************************************************

   !----------------------------------------------------------------
   ! USE : Writes a 3D array into binary files in [TEMP] dir
   !      
   !
   ! FORM:    subroutine saveFFT_data()
   !
   ! BEHAVIOR: Saves only one set of LES and TEST filtered data at 
   !            a time. 
   !
   ! STATUS : 
   !        On GIT BRANCH:
   !          1) test/validation: no arguments required.
   !          2) filter3: both arguments required. 
   ! 
   !----------------------------------------------------------------
   
   subroutine saveFFT_data(LESFilterType,TestFilterType)
     implicit none
     !
     !    ..ARGUMENTS..
     character(*),intent(in) :: LESFilterType
     character(*),intent(in) :: TestFilterType
     !
     !    ..LOCAL VARIABLES..
     character(64) :: filename 
     integer :: i


     ! save FFT_DATA : Filtered velocities and stress files: 
     print*
     print*,'Write filtered variables ... '
 
     do i = 1, size(var_FFT)
       if (i.eq.1.or.i.eq.3) then
        filename = trim(TEMP_PATH)//trim(LESfilterType)//'/'//trim(var_FFT(i)%name)//'.bin'
       ! filename = trim(TEMP_PATH)//trim(var_FFT(i)%name)//'.bin'
       else
        filename = trim(TEMP_PATH)//trim(LESfilterType)//'/'//trim(TestfilterType)//'/'//trim(var_FFT(i)%name)//'.bin'
      endif
        print*, filename
        open(i, file = filename,form='unformatted')
        if (i.eq.1) write(i) u_f
        if (i.eq.2) write(i) u_t
        if (i.eq.3) write(i) tau_ij
        if (i.eq.4) write(i) T_ij
        close(i)
     end do

   end subroutine saveFFT_data
  

   !****************************************************************
   !                          CHECK FFT_DATA
   !****************************************************************

   !----------------------------------------------------------------
   ! USE : Check FFT dataset
   !      
   !
   ! FORM:    subroutine saveFFT_data()
   !
   ! BEHAVIOR: Needs allocated, defined arrays.
   !
   ! STATUS : 
   ! 
   !----------------------------------------------------------------
   
   subroutine checkFFT_data()
     implicit none

     print*, 'check FFT_data'

        if (withPressure) then
           if (u_f(4,15,24,10).ne.0.12164158031779296d0) then
              print*, 'ERROR READING DATA'
              print*, u_f(4,15,24,10)
              stop
           end if
        end if

        if (dataset.eq.'jhu256') then
           if (u_t(1,15,24,10).ne.-0.48241021987284982d0) then
              print*, 'ERROR READING DATA'
              print*, u_t(1,15,24,10)
              stop
           elseif( ( (stress.eq.'dev') .and. (T_ij(1,15,24,10).ne.-5.2544371578038401d-3))   &
                .or. ( (stress.eq.'abs') .and. (T_ij(1,15,24,10).ne.-4.0152351790891661d-3) ) )then
              print*, 'ERROR COMPUTING T_ij'
              print*,'T_ij(1,15,24,10)',T_ij(1,15,24,10)
              stop
           else
              print*, 'Read data saved from main.f90: Passed \n'
           end if
        end if


   end subroutine checkFFT_data

  
   !****************************************************************
   !                        PLOT ORIGINAL STRESS
   !****************************************************************

   !----------------------------------------------------------------
   ! USE : Saves [z-midplane] in RESULTS directory
   !      
   !
   ! FORM:    subroutine plotOriginalStress()
   !
   ! BEHAVIOR: Needs allocated, defined arrays.
   !
   ! STATUS :
   !       do i = 1, size(var_FFT)
   !         filename = trim(RES_PATH)//trim(var_FFT(i)%name)//'.dat'
   !         print*, filename
   !         open(i, file = filename)
   !         if (i.eq.1) call printPlane (u_f(1,:,:,1),   fID=i)
   !         if (i.eq.2) call printPlane (u_t(1,:,:,1),   fID=i)
   !         if (i.eq.3) call printPlane (tau_ij(1,:,:,1),fID=i)
   !         if (i.eq.4) call printPlane (T_ij(1,:,:,1),  fID=i)
   !         close(i)
   !      end do
   !
   !   
   !----------------------------------------------------------------
   
   subroutine plotOriginalStress(plotOption)
     implicit none
     !
     !    ..SCALAR ARGUMENTS..
     character(*), optional, intent(in) :: plotOption
     !
     !    ..LOCAL VARIABLES..
     character(64) :: filename
     integer :: i, n_ij
     character(1) :: ij


     n_ij = 1
     if (plotOption.eq.'All') n_ij = 6
     
     ! SAVE ORIGINAL STRESS
     print*,'Saving original stress in', RES_PATH
         
     do i=1,n_ij
        write(ij, '(i0)') i
        open(10,file=trim(RES_PATH)//'T_ij_'//trim(stress)//trim(ij)//'.dat')
        open(11,file=trim(RES_PATH)//'tau_ij_'//trim(stress)//trim(ij)//'.dat')

        write(10,*) T_ij  (i,:,:,z_plane)
        write(11,*) tau_ij(i,:,:,z_plane)

        close(10)
        close(11)
     end do

   end subroutine plotOriginalStress



   
   !****************************************************************
   !                        PLOT DYNAMIC SMAG
   !****************************************************************
   
   !----------------------------------------------------------------
   ! USE : Saves [z-midplane] in RESULTS directory
   !      
   !
   ! FORM:    subroutine plotDynSmag()
   !
   ! BEHAVIOR: Needs allocated, defined arrays.
   !
   ! STATUS :
   !   
   !   
   !----------------------------------------------------------------
   subroutine plotDynSmag()
     implicit none
     !
     !    ..LOCAL VARIABLES..
     integer :: i, n_ij
     character(1) :: ij

     print*, 'Saving tau_DS in ', RES_PATH
     ! SAVE tau_DS:
     do i = 1,6
        write(ij,'(i0)') i
        open(87, file = trim(RES_PATH)//'43/tau_DS_'//trim(stress)//trim(ij)//'.dat')
        write(87,*), tau_DS(i,:,:,43)
        close(87)
        open(87, file = trim(RES_PATH)//'129/tau_DS_'//trim(stress)//trim(ij)//'.dat')
        write(87,*), tau_DS(i,:,:,z_plane)
        close(87)
        open(87, file = trim(RES_PATH)//'212/tau_DS_'//trim(stress)//trim(ij)//'.dat')
        write(87,*), tau_DS(i,:,:,212)
        close(87)
     end do

     ! SAVE Pij_DS:
      print*,'Saving DS production field in', RES_PATH
      open(53, file=trim(RES_PATH)//'43/Pij_DS_'//trim(stress)//'.dat')
      write(53,*) Pij_DS(:,:,43)
      close(53)
      open(53, file=trim(RES_PATH)//'129/Pij_DS_'//trim(stress)//'.dat')
      write(53,*) Pij_DS(:,:,z_plane)
      close(53)
      open(53, file=trim(RES_PATH)//'212/Pij_DS_'//trim(stress)//'.dat')
      write(53,*) Pij_DS(:,:,212)
      close(53)

   end subroutine plotDynSmag

   
   !****************************************************************
   !                        PLOT BARDINA
   !****************************************************************
   
   !----------------------------------------------------------------
   ! USE : Saves [z-midplane] in RESULTS directory
   !      
   !
   ! FORM:    subroutine plotBardina()
   !
   ! BEHAVIOR: Needs allocated, defined arrays.
   !
   ! STATUS :
   !   
   !   
   !----------------------------------------------------------------
   subroutine plotBardina(BD_CASE)
     implicit none
     !
     !    ..ARGUMENTS..
     character(*), intent(in) :: BD_CASE
     !
     !    ..LOCAL VARIABLES..
     integer :: i, n_ij
     character(1) :: ij

     print*, 'Saving tau_BD in ', RES_PATH
     ! SAVE tau_BD:
     do i = 1,6
        write(ij,'(i0)') i
        open(87, file = trim(RES_PATH)//'tau_'//trim(BD_CASE)//'_BD_'//trim(stress)//trim(ij)//'.dat')
!        open(88, file = trim(RES_PATH)//'T_a105_BD'//trim(ij)//'.dat')
        write(87,*), tau_BD(i,:,:,z_plane)
!        write(88,*), T_ij(i,:,:,z_plane)
        close(87)
!        close(88)
     end do

     ! SAVE Pij_BD:
      print*,'Saving BD production field in', RES_PATH
      open(53, file = trim(RES_PATH)//'Pij_'//trim(BD_CASE)//'_BD_'//trim(stress)//'.dat')
      write(53,*) Pij_BD(:,:,z_plane)
      close(53)

   end subroutine plotBardina
   

  
   !****************************************************************
   !                        PLOT COMPUTED STRESS
   !****************************************************************

   !----------------------------------------------------------------
   ! USE : Saves [z-midplane] in RESULTS directory
   !      
   !
   ! FORM:    subroutine plotComputedStress()
   !
   ! BEHAVIOR: Needs allocated, defined arrays.
   !
   ! STATUS :
   !   
   !   
   !----------------------------------------------------------------
   
   subroutine plotComputedStress(lambda, plotOption)
     implicit none
     !
     !    ..SCALAR ARGUMENTS..
     real(8), intent(in) :: lambda
     character(*), optional, intent(in) :: plotOption
     !
     !    ..LOCAL VARIABLES..
     integer :: i, n_ij
     character(64) :: lambda_char
     character(1) :: ij

     n_ij = 1
     if (plotOption.eq.'All') n_ij = 6

     write(lambda_char,'(ES6.0E2)') lambda

     ! SAVE COMPUTED STRESS
     print*,'Saving computed stress in', RES_PATH
 
     do i = 1, n_ij
        write(ij, '(i0)') i
        open(10,file=trim(RES_PATH)//'T_ijOpt_'//trim(stress)//trim(ij)//trim(lambda_char(4:6)) // '.dat',status='replace')
        open(11,file=trim(RES_PATH)//'tau_ijOpt_'//trim(stress)//trim(ij)//trim(lambda_char(4:6)) // '.dat',status='replace')

        write(10,*) T_ijOpt  (i,:,:,z_plane)
        write(11,*) tau_ijOpt(i,:,:,z_plane)

        close(10)
        close(11)
     end do

   end subroutine plotComputedStress
  
   !****************************************************************
   !                        LOAD COMPUTED STRESS
   !****************************************************************

   !----------------------------------------------------------------
   ! USE : Loads precomputed stresses for each lambda value
   !      
   !
   ! FORM:    subroutine loadComputedStress()
   !
   ! BEHAVIOR: Needs allocated, defined arrays.
   !
   ! STATUS :
   !   
   !   
   !----------------------------------------------------------------
   
   subroutine loadComputedStress(lambda, loadOption)
     implicit none
     !
     !    ..SCALAR ARGUMENTS..
     real(8), intent(in) :: lambda
     character(*), optional, intent(in) :: loadOption
     !
     !    ..LOCAL VARIABLES..
     integer :: i, n_ij
     character(64) :: lambda_char
     character(1) :: ij

     n_ij = 1
     if (loadOption.eq.'All') n_ij = 6

     write(lambda_char,'(ES6.0E2)') lambda

     ! SAVE COMPUTED STRESS
     print*,'Load computed stress from', RES_PATH
 
     do i = 1, n_ij
        write(ij, '(i0)') i
        open(10,file=trim(RES_PATH)//'T_ijOpt_'  //trim(stress)//trim(ij)//trim(lambda_char(4:6)) // '.dat',status='old')
        open(11,file=trim(RES_PATH)//'tau_ijOpt_'//trim(stress)//trim(ij)//trim(lambda_char(4:6)) // '.dat',status='old')

        read(10,*) T_ijOpt  (i,:,:,z_plane)
        read(11,*) tau_ijOpt(i,:,:,z_plane)

        close(10)
        close(11)
     end do

   end subroutine loadComputedStress



  !****************************************************************
  !                       PLOT PRODUCTION TERM
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE : Writes .dat file at z-midplane
  !      
  !
  ! FORM:    subroutine plotProductionTerm()
  !
  ! BEHAVIOR: Called from autonomic.f90 to plot original Pij
  !           Called from solver.f90 to plot Pij_Opt for a given lambda
  !
  !
  ! STATUS : 
  ! 
  !----------------------------------------------------------------

   subroutine plotProductionTerm(lambda)
     implicit none
     !
     ! 
     real(8), intent(in), optional :: lambda
     integer :: i
     character(64) :: lambda_char

     write(lambda_char,'(ES6.0E2)') lambda
     
     if (present(lambda)) then
        ! SAVE COMPUTED PRODUCTION TERM
        print*,'Saving computed production field in', RES_PATH
        call system ('mkdir -p '//trim(RES_PATH))

        open(1,file=trim(RES_PATH)//'Pij_fOpt_'//trim(stress)//trim(lambda_char(4:6)) // '.dat')
        open(2,file=trim(RES_PATH)//'Pij_tOpt_'//trim(stress)//trim(lambda_char(4:6)) // '.dat')

        write(1,*) Pij_fOpt(:,:,z_plane)
        write(2,*) Pij_tOpt(:,:,z_plane)
        close(1)
        close(2)

     else 
        ! SAVE ORIGINAL PRODUCTION TERM
        print*
        print*,'Saving original production field in', RES_PATH
        call system ('mkdir -p '//trim(RES_PATH))

        open(1,file = trim(RES_PATH)//'Pij_f_'//trim(stress)//'.dat')
        open(2,file = trim(RES_PATH)//'Pij_t_'//trim(stress)//'.dat')
        write(1,*) Pij_f(:,:,z_plane)
        write(2,*) Pij_t(:,:,z_plane)
        close(1)
        close(2)
     end if

   end subroutine plotProductionTerm


   !****************************************************************
   !                            PLOT PDF
   !****************************************************************

   !----------------------------------------------------------------
   ! USE : Write a csv file to generate pdf plot
   !      
   !
   ! FORM: subroutine plotPDF(x, pdf, [fieldname], [avg], [SD])
   !
   ! BEHAVIOR: 
   !
   ! STATUS :
   !      
   !   
   !----------------------------------------------------------------
   
   subroutine plotPDF(x, pdf, fieldname, avg, SD)
     implicit none
     !
     !    ..ARRAY ARGUMENTS..
     real(8), dimension(:), intent(in) :: x
     real(8), dimension(:), intent(in) :: pdf
     !
     !    ..SCALAR ARGUMENTS..
     character(*), intent(in) :: fieldname
     real(8), optional, intent(in) :: avg, SD

     print*
     print*,'Saving PDF'//trim(fieldname)//' in ', RES_PATH

     !
     ! PLOT PDF
     open(10, file = trim(RES_PATH)//trim('PDF')//trim(fieldname)//'.csv')
     write(10,"( F10.4,',',F10.4 )"), (x(i), pdf(i), i=1,n_bins)
     close(10)

     ! SAVE STAT DATA
     if (present(avg)) then
        open(11, file = trim(RES_PATH)//trim('statistics')//trim(fieldname)//'.dat')
        write(11, *) avg
        if (present(SD))  write(11, *) SD
        close(11)
     end if

   end subroutine plotPDF

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
  ! FORM:   subroutine xyplot(var,var_name)
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
