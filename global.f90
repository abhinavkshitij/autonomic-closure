!****************************************************************
!                          MODULE: GLOBAL
!****************************************************************
!  
!----------------------------------------------------------------
! USE:  
!     Contains global definitions.
!       
! FORM:
!      module global
!          contains
!          interface printplane
!              subroutine plane3
!              subroutine plane2
!        subroutine printParams       [SINK]
!
! BEHAVIOR:
!      * Check array rank while passing into printplane.
!      * 
!        
!      * 
!      * 
!  
!----------------------------------------------------------------

module global

  integer, parameter :: dp = selected_real_kind(15,307)

  real,    parameter :: eps = 1e-5
  real(8), parameter :: pi = 4.d0 * atan(1.d0)

  ! DEFINE STRING ARRAY:
  type str16
     character (16) :: name
  end type str16

  type(str16), parameter :: l_machine(2) = [str16 ('local'),     &
                                            str16 ('remote')]

  type(str16), parameter :: var(2) = [str16 ('Velocity'),      &
                                      str16 ('Pressure')]

  type(str16), parameter :: l_dataset(5) = [str16 ('nrl'),     & 
                                            str16 ('jhu256'),  &
                                            str16 ('hst'),     &
                                            str16 ('sin3D'),   &
                                            str16 ('jhu1024')]

  type(str16), parameter :: var_FFT(4) = [str16 ('u_f'),     &
                                          str16 ('u_t'),     &
                                          str16 ('tau_ij'),  &
                                          str16 ('T_ij')]  
                                      

  type(str16), parameter :: l_solutionMethod(2) = [str16 ('LU'),         &
                                                   str16 ('SVD')]

  type(str16), parameter :: l_trainingPoints(2) = [str16 ('ordered'),    &
                                                   str16 ('random')] 

  type(str16), parameter :: l_formulation(2) = [str16 ('colocated'),     &
                                                str16 ('noncolocated')]

              
  character(8) :: machine        = trim (l_machine(2) % name)
  character(8) :: dataset        = trim (l_dataset(3) % name)
  logical      :: withPressure   = 0
  character(8) :: solutionMethod = trim (l_solutionMethod(1) % name) ! [LU, SVD]
  character(2) :: hst_set = 'S1' ! [S1, S3, S6]
  character(3) :: stress = 'dev' ! [dev, abs]
  character(16):: formulation    = trim (l_formulation(2) % name)
  character(8) :: trainingPoints = trim (l_trainingPoints(1) % name)


  !----------------------------------------------------------------
  !
  !****************************************************************
  !                       ..NAMESPACE..                          !
  !****************************************************************

  !
  !    ..VELOCITIES..
  real(8), dimension(:,:,:,:), allocatable :: u
  real(8), dimension(:,:,:,:), allocatable :: u_f 
  real(8), dimension(:,:,:,:), allocatable :: u_t
  !
  !    ..STRESSES..
  real(8), dimension(:,:,:,:), allocatable :: tau_ij
  real(8), dimension(:,:,:,:), allocatable :: T_ij
  real(8), dimension(:,:,:,:), allocatable :: T_ijOpt
  real(8), dimension(:,:,:,:), allocatable :: tau_ijOpt
  !
  !    ..FILTERS..
  real(8), dimension(:,:,:), allocatable :: LES
  real(8), dimension(:,:,:), allocatable :: test
  !
  !    ..STRAIN RATES..
  real(8), dimension(:,:,:,:,:), allocatable :: Sij_f
  real(8), dimension(:,:,:,:,:), allocatable :: Sij_t
  !
  !    ..COEFFICIENTS..
  real(8), dimension(:,:), allocatable :: h_ij
  !
  !    ..PRODUCTION TERM..
  real(8), dimension(:,:,:), allocatable :: Pij_f
  real(8), dimension(:,:,:), allocatable :: Pij_t
  real(8), dimension(:,:,:), allocatable :: Pij_fOpt
  real(8), dimension(:,:,:), allocatable :: Pij_tOpt
 

  integer :: i_GRID 
  integer :: j_GRID 
  integer :: k_GRID 
  integer :: f_GRID 
  integer :: center

  real(8) :: dx ! To calculate gradient


  integer :: M               ! Number of training points 3x3x3x9
  integer :: N 
  integer :: P                  ! Number of components (1 to 6)
  integer :: stencil_size 
  integer :: LES_scale 
  integer :: test_scale 


  ! Stencil parameters:
  integer :: stride ! Is the ratio between LES(taken as 1) and test scale
  integer :: skip 
  integer :: X ! Number of realizations
  integer :: n_lambda  ! Number of lambdas
  real(8) :: lambda, lambda_0(2)
 
  ! Bounding Box parameters:  
  integer :: box(3) 
  integer :: boxSize  
  integer :: maskSize  
  integer :: boxCenter
  integer :: boxLower 
  integer :: boxUpper 

  ! Test field parameters: 
  integer :: testSize 
  integer :: testcutSize 
  integer :: testLower
  integer :: testUpper
  

  ! Cutout parameters:
  integer :: lBound 
  integer :: uBound 

  ! Statistics parameters:
  integer :: samples 

  !
  !    ..FILEIO..
  integer :: z_print
  integer :: path_txt   = 22
  integer :: params_txt = 23
  integer :: cross_csv  = 81

  !    ..TIME..
  character(32) :: time = '256' ! 256 is the initial value
  integer :: time_init
  integer :: time_incr 
  integer :: time_final
  integer :: n_time

  !
  !    .. DIRS..
  character(*), parameter :: DATA_DIR = '../data/'
  character(*), parameter :: TEMP_DIR = '../temp/'
  character(*), parameter :: RES_DIR  = '../results/'
  
  character(64) :: DATA_PATH
  character(64) :: TEMP_PATH
  character(64) :: RES_PATH

  character(4) :: ext   = 'bin'   ! Dataset extension: [bin]ary, [h5] or [txt] format. 

  !
  !    ..MPI..
  integer :: n_proc
  integer :: rank
  integer :: ierr

  !
  !    ..FORMATS..
  character(*), parameter :: long = 'es23.17'
  character(*), parameter :: short = 'f10.4'
  character(*), parameter :: csv_Table = "( 2(f10.4, ',') , / )"
  !
  !    ..INDICES..
  integer  :: fileID = 6
  integer  :: n_u 
  integer  :: n_uu 
  integer  :: i, j, k  
  integer  :: iter
  real(8)  :: tic, toc
  

  !----------------------------------------------------------------
  !
  !****************************************************************
  !                       ..INTERFACES..                          !
  !****************************************************************

  interface printplane
     module procedure plane2, plane3
  end interface

contains
  
  subroutine setEnv()
    
    RES_PATH = RES_DIR

    ! GRID
    i_GRID = 256
    j_GRID = 256
    k_GRID = 256
  

    ! TIME
    if (dataset.eq.'nrl') time = '0460'
    if (dataset.eq.'hst') then
       if (hst_set.eq.'S1') then
          time = '016'; time_init = 16; time_incr = 1; time_final = 28
       end if
       if (hst_set.eq.'S3') then
          time = '015'; time_init = 15; time_incr = 1; time_final = 24
       end if
       if (hst_set.eq.'S6') then
          time = '070'; time_init = 70; time_incr = 5; time_final = 100
       end if
    end if
    read(time,*) time_init
    if (machine.eq.'local') time_final = time_init
    

    ! COMPONENT
    n_u = 3
    if (withPressure) n_u = 4
    ! CHECK PRESSURE DATA AVAILABLIITY:
    if(dataset.ne.'jhu256'.and.n_u.eq.4) then
       print*, 'Dataset ',dataset,' has no pressure data files'
       stop
    end if
    n_uu = 6
    P = 6

    ! SCALE
    LES_scale  = 40
    test_scale = 20


    ! FFT 
    f_GRID = 256 !For FFT ops, the grid must be cubic.
    center = (0.5d0 * f_GRID) + 1.d0
    dx = 2.d0*pi/dble(i_GRID) !Only for JHU data. Change for others.    


    ! Stencil parameters:
    stride = 1              ! No subsampling on the original DNS GRID 
    skip = 10               ! For training points
    X = 1     
    stencil_size = n_u * (3*3*3) !3 * 27 = 81  

    lambda_0 = 1.d-11 * [1,3]
    n_lambda = 10

 
    ! Bounding Box parameters:  
    if (formulation.eq.'noncolocated') then
       ! Set number of features
       N = nCk(stencil_size + 1, 2) + nCk(stencil_size, 1) + 1

       if (trainingPoints.eq.'ordered') then
          ! Set bounding box size
          box = [i_GRID-2, j_GRID-2, k_GRID-2]
          ! Set number of training points
          M =    (floor((real(box(1) - 1)) / skip) + 1)    &
               * (floor((real(box(2) - 1)) / skip) + 1)    &
               * (floor((real(box(3) - 1)) / skip) + 1)  !17576 

       else if (formulation.eq.'colocated') then
          box  = [8,8,8]
          if (trainingPoints.eq.'random') then
          end if
       end if
    end if

    

    boxSize   = product(box)
    maskSize  = boxSize - M ! 512-Training points(243) = 269
    boxCenter = smallHalf(box(1)) * box(1)*(box(1) + 1) + bigHalf(box(1))

    boxLower  = stride * (bigHalf(box(1)) - 1)
    boxUpper  = stride * (box(1) - bigHalf(box(1)))


    ! Test field parameters: 
    testSize    = 1
    testLower   = stride *  bigHalf(box(1)) + 1
    testUpper   = stride * (bigHalf(box(1)) - 1 + testSize) + 1
    testcutSize = stride * (testSize + box(1)) + 1

    ! Cutout parameters:
    lBound = 0.5*(f_GRID - testcutSize)
    uBound = 0.5*(f_GRID + testcutSize) - 1

    ! Statistics parameters:
    samples = 250

  end subroutine setEnv


  !****************************************************************
  !                         PRINT PARAMETERS                      !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Display main parameters either onto the screen or writes
  !      to a file. Writes PATHS and parameters for postprocessing
  !      in MATLAB
  !      
  !      
  ! FORM:   subroutine printParams()
  ! 
  !      
  ! BEHAVIOR: 
  !           
  !
  ! STATUS : 
  !          
  !
  !----------------------------------------------------------------

  
  subroutine printParams(displayOption)

    character(*), optional :: displayOption
    !
    ! Write parameters in params.txt for MATLAB to read in data for
    ! post-processing the results.
    open(params_txt, file = trim(RES_DIR)//'params.txt')
    write(params_txt,*) i_GRID
    write(params_txt,*) j_GRID
    write(params_txt,*) k_GRID
    write(params_txt,*) dataset
    write(params_txt,*) hst_set
    close(params_txt)

    if (displayOption.eq.'display') then
     write(fileID, * ), 'Dataset:            ', dataset, '\n'

     write(fileID, * ), 'Stencil parameters: \n'

     write(fileID, * ), 'Training points :    ',M
     write(fileID, * ), 'Features :           ',N, '\n'

     write(fileID, * ), 'Box parameters:      '
     write(fileID, * ), 'Bounding box:        ',box
     write(fileID, * ), 'boxCenter:           ',boxCenter
     write(fileID, * ), 'boxLower:            ',boxLower
     write(fileID, * ), 'boxUpper:            ',boxUpper, '\n'

     write(fileID, * ), 'Test field parameters:'
     write(fileID, * ), 'testcutSize:         ',testcutSize
     write(fileID, * ), 'testLower:           ',testLower
     write(fileID, * ), 'testUpper:           ',testUpper, '\n'

     write(fileID, * ), 'Cutout parameters:'
     write(fileID, * ), 'lower bound:         ',lBound
     write(fileID, * ), 'upper bound:         ',uBound, '\n'

     write(fileID, * ),  'Number of samples'    ,samples
  end if
  end subroutine printParams


  !****************************************************************
  !                          PLANE3
  !****************************************************************
  !  
  !----------------------------------------------------------------
  ! USE:  
  !     Display [Save] [x], [y], or [z] plane of a 3D array (BOX).
  !       
  ! FORM:
  !      plane3 (box, [frameLim], [x], [y], [z], [fID])
  !         
  ! BEHAVIOR:
  !      * Dumps the entire z-plane (z=1) to the standard output (default:screen)
  !      * Looks for [frameLim] to limit output size of the plane (useful for
  !        large boxes).
  !      * Looks for specified [x-  y-  z-] plane. (default: z=1 plane)
  !      * Looks for specified file [fID] to [Save] plane in an external file.
  !  
  !----------------------------------------------------------------

  subroutine plane3 (box, frameLim, x, y, z, fID)
    implicit none
   
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:,:),intent(in) :: box
    !
    !    ..SCALAR ARGUMENTS.. 
    integer, intent(in), optional :: frameLim 
    integer, intent(in), optional :: x, y, z
    integer, intent(in), optional :: fID
    !
    !    ..LOCAL SCALARS..
    integer :: x_lim, y_lim, z_lim
    !
    !    ..FORMATTING..
    character(4) :: form 
    character(8) :: digits = short 



    ! Check [frameLimits]:
    x_lim = size (box, dim=1)
    y_lim = size (box, dim=2)
    z_lim = size (box, dim=3)

    if (present(frameLim)) then
       x_lim = frameLim
       y_lim = frameLim
       z_lim = frameLim
    end if

      
    ! Write plane to an external file with 15 digits.
    if (present(fID)) then
       fileID = fID
       digits = long 
    end if
    
    
    ! Write [ x y z ] plane to fileID (default: screen)
    if (present(x)) then
       write (form,'(i4)') y_lim     
       write (fileID,'('//trim(form)//trim(digits)//')'), (box(x,1:y_lim,i),i=1,z_lim); print*,'' ! +x-plane
    elseif (present(y)) then
       write (form,'(i4)') x_lim     
       write (fileID,'('//trim(form)//trim(digits)//')'), (box(1:x_lim,y,i),i=1,z_lim); print*,'' ! -y-plane                  
    elseif (present(z)) then
       write (form,'(i4)') y_lim     
       write (fileID,'('//trim(form)//trim(digits)//')'), (box(i,1:y_lim,z),i=1,x_lim); print*,'' ! +z-plane
    else
       write (form,'(i4)') y_lim     
       write (fileID,'('//trim(form)//trim(digits)//')'), (box(i,1:y_lim,1),i=1,x_lim); print *,'' !z=1 default
    end if

    return
  end subroutine plane3

  !----------------------------------------------------------------
  !                          PLANE2
  !----------------------------------------------------------------
  ! USE:  
  !     Display [Save] selected plane of a 2D array (matrix).
  !       
  ! FORM:
  !      plane2 (matrix, [frameLim], [fID])
  !         
  ! BEHAVIOR:
  !      * Dumps the entire plane to the standard output (screen)
  !      * Looks for [frameLim] to limit output size of the plane (useful for
  !        large arrays).
  !      * Limits frame size to 8 to avoid large outputs to screen.
  !      * Looks for specified file [fID] to [Save] plane in an external file.
  !  
  !----------------------------------------------------------------

  subroutine plane2 (matrix, frameLim, fID)
    implicit none

    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:), intent(in) :: matrix
    !
    !    ..SCALAR ARGUMENTS.. 
    integer, intent(in), optional :: frameLim 
    integer, intent(in), optional :: fID
    !
    !    ..LOCAL SCALARS..
    integer :: x_lim, y_lim
    integer :: frameLim_Max = 8
    !
    character(4) :: form 
    character(8) :: digits = short 
    
    ! Check [frameLimits]:
    x_lim = size (matrix, dim=1)
    y_lim = size (matrix, dim=2)
   

    if (present(frameLim)) then
       if (frameLim.lt.8) then
          x_lim = frameLim
          y_lim = frameLim
       else
          print*,"frameLim max set to 8 to limit large output to screen."
          x_lim = frameLim_Max
          y_lim = frameLim_Max
       end if
    end if

      
    ! Write plane to an external file with 15 digits.
    if (present(fID)) then
       fileID = fID
       digits = long 
    end if
    
    write (form,'(i4)') y_lim     
    write (fileID,'('//trim(form)//trim(digits)//')'), (matrix(i,1:y_lim),i=1,x_lim); print*,'' ! +z-plane
    
  end subroutine plane2



  !----------------------------------------------------------------
  !                          COMBINATION - nCk
  !----------------------------------------------------------------
  ! USE:  
  !     Computes nCk
  !       
  ! FORM:
  !      function combination (int n, int k)
  !         
  ! BEHAVIOR:
  !
  !----------------------------------------------------------------
  
  function nCk(n,k) result(res)
    implicit none
    !
    !    ..SCALAR ARGUMENTS..
    integer, value:: n, k
    integer :: res
    !
    !    ..LOCAL VARS..
    integer :: Nr, Dr, count

    res = 1; Nr = 1; Dr = 1

    do while (k.ne.0)
       Nr = Nr * n
       n = n - 1
       Dr = Dr * k
       k = k - 1
    end do

    res = Nr / Dr 
    
  end function nCk

  !----------------------------------------------------------------
  !                          FACTORIAL
  !----------------------------------------------------------------
  ! USE:  
  !     Computes the factorial of a positive integer
  !       
  ! FORM:
  !      function recursive factorial(int n)
  !         
  ! BEHAVIOR:
  ! 
  !  
  !----------------------------------------------------------------
  
  recursive function factorial(num) result (res)
    implicit none
    !
    !    ..SCALAR ARGUMENTS..
    integer :: num
    integer :: res
    
    if (num.eq.1.or.num.eq.0) then
       res = 1
    else
       res = num * factorial(num-1)
    end if  
  end function factorial



  !----------------------------------------------------------------
  !                          MEAN
  !----------------------------------------------------------------
  ! USE:  
  !     Computes the mean value of an array
  !     
  ! FORM:
  !      function mean(double array)
  !         
  ! BEHAVIOR:
  ! 
  !  
  !----------------------------------------------------------------
  
  function mean(array) 
    real(8),dimension(:,:,:),intent(in):: array  
    real(8) :: mean
    mean = sum(array) / size(array)
  end function mean

  
  !----------------------------------------------------------------
  !                          STANDARD DEVIATION
  !----------------------------------------------------------------
  ! USE:  
  !     Computes the standard deviation of an array based on
  !     estimated average [not true average]
  !     
  !     
  ! FORM:
  !      function stdev(double array)
  !         
  ! BEHAVIOR:
  ! 
  !  
  !----------------------------------------------------------------
  
  function stdev(array)
    real(8), dimension(:,:,:),intent(in):: array
    real(8) :: stdev
    stdev = sqrt ( sum((array - mean(array))**2) / (size(array)-1) )
  end function stdev



  !----------------------------------------------------------------
  !                          NORM
  !----------------------------------------------------------------
  ! USE:  
  !     Computes norm-1, norm-2 or norm-Inf of a 3D array     
  !     
  ! FORM:
  !      function norm(double array,['1', or '2', or 'Inf'])
  !         
  ! BEHAVIOR: default is norm-2
  ! 
  !  
  !----------------------------------------------------------------
  
  function norm(array, normType)
    real(8), dimension(:,:,:),intent(in):: array
    character(*), optional, intent(in) :: normType 
    real(8) :: norm

    if (present(normType)) then
       ! 
       ! L1 norm
       if (normType.eq.'1') then 
          norm = sum(abs(array))
       !
       ! L2 norm   
       elseif (normType.eq.'2') then
          norm = sqrt(sum(array**2))
       !
       ! L-Inf norm    
       elseif (normType.eq.'Inf') then
          norm = maxval(array)
       else
          print*, 'Norms computed are 1-norm, 2-norm or Inf-norm'
       end if
    else
       ! Default: norm-2 
       norm = sqrt(sum(array**2)) 
    end if   
  end function norm


  !----------------------------------------------------------------
  !                          BIGHALF
  !----------------------------------------------------------------
  ! USE:  
  !     Computes bigHalf 
  !     
  ! FORM:
  ! 
  !         
  ! BEHAVIOR: 
  ! 
  !  
  !----------------------------------------------------------------
  
  function bigHalf(num)
    implicit none
    !
    !    ..SCALAR ARGUMENTS..
    integer, intent(in) :: num
    integer :: bighalf
    !
    !    ..LOCAL SCALARS..
    real :: eps = 1e-5

    bigHalf   = ceiling(0.5*real(num) + eps) 

  end function bigHalf


  !----------------------------------------------------------------
  !                          SMALLHALF
  !----------------------------------------------------------------
  ! USE:  
  !     Computes smallHalf 
  !     
  ! FORM:
  ! 
  !         
  ! BEHAVIOR: 
  ! 
  !  
  !----------------------------------------------------------------
  
  function smallHalf(num)
    implicit none
    !
    !    ..SCALAR ARGUMENTS..
    integer, intent(in) :: num
    integer :: smallhalf
    !
    !    ..LOCAL SCALARS..
    real :: eps = 1e-5

    smallHalf   = floor(0.5*real(num) + eps) 

  end function smallHalf

  

end module global
