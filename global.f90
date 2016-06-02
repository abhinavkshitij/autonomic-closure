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

  real,    parameter :: eps = 1e-3 
  real(8), parameter :: pi = 4.d0 * atan(1.d0)

  ! DEFINE STRING ARRAY:
  type str16
     character (16) :: name
  end type str16


  type(str16), parameter :: var(2) = [str16 ('Velocity'),    &
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

              

  character(8) :: dataset        = trim (l_dataset(2) % name)
  character(8) :: solutionMethod = trim (l_solutionMethod(1) % name) ! LU, SVD
  character(2) :: hst_set = 'S6' ! HST datasets - S1, S3, S6
  character(3) :: stress = 'dev' ! dev,abs
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
  real(8), dimension(:,:,:), allocatable :: P_f
  real(8), dimension(:,:,:), allocatable :: P_t
 
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
  integer :: n_DAMP  ! Number of lambdas
  real(8) :: lambda
 
  ! Bounding Box parameters:  
  integer :: box(3) 
  integer :: boxSize  
  integer :: maskSize  
  integer :: bigHalf   
  integer :: smallHalf  
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

    i_GRID = 256
    j_GRID = 256
    k_GRID = 256
    if (dataset.eq.'hst') then
       j_GRID = 129
    end if

    n_u = 3
    n_uu = 6
    P = 6

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
    n_DAMP = 1 
    stencil_size = 3*(3*3*3) !3 * 27 = 81  
 
    ! Bounding Box parameters:  
    if (formulation.eq.'noncolocated') then
       ! Set number of features
       N = 3403                 ! 5995 for velocity and pressure
!       N = 5995

       if (trainingPoints.eq.'ordered') then
          ! Set bounding box size
          box = [i_GRID, j_GRID, k_GRID]
          ! Set number of training points
          M =  (floor((real(box(1) - 1)) / skip) + 1)    &
               * (floor((real(box(2) - 1)) / skip) + 1)  &
               * (floor((real(box(3) - 1)) / skip) + 1)  !17576 

       else if (formulation.eq.'colocated') then
          box  = [8,8,8]
          if (trainingPoints.eq.'random') then
          end if
       end if
    end if
    

    boxSize = product(box)


    maskSize = boxSize - M ! 512-Training points(243) = 269
    bigHalf   = ceiling(0.5*real(box(1)) + eps) ! for 8->5
    smallHalf = floor(0.5*real(box(1)) + eps)   ! for 8->4
    boxCenter = smallHalf * box(1)*(box(1) + 1) + bigHalf
    boxLower  = stride * (bigHalf - 1)
    boxUpper  = stride * (box(1) - bigHalf)

    ! Test field parameters: 
    testSize = 1
    testcutSize = stride * (testSize + box(1)) + 1
    testLower = stride * bigHalf + 1
    testUpper = stride * (bigHalf - 1 + testSize) + 1


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

  
  subroutine printParams()

    
    open(23, file = trim(RES_DIR)//'params.txt')
    write(23,*) i_GRID
    write(23,*) j_GRID
    write(23,*) k_GRID
    write(23,*) dataset
    write(23,*) hst_set
    close(23)

    
     write(fileID, * ), ''
     write(fileID, * ), 'Dataset:            ', dataset
     write(fileID, * ), ''
!     write(fileID, * ), 'Stencil parameters:'
!     write(fileID, * ), ''
!     write(fileID, * ), 'Training points :    ',M
!     write(fileID, * ), 'Features :           ',N
!     write(fileID, * ), ''
!     write(fileID, * ), 'Box parameters:      '
!     write(fileID, * ), 'Bounding box:        ',box
!     write(fileID, * ), 'bigHalf:             ',bigHalf
!     write(fileID, * ), 'smallHalf:           ',smallHalf
!     write(fileID, * ), 'boxCenter:           ',boxCenter
!     write(fileID, * ), 'boxLower:            ',boxLower
!     write(fileID, * ), 'boxUpper:            ',boxUpper
!     write(fileID, * ), ''
!     write(fileID, * ), 'Test field parameters:'
!     write(fileID, * ), 'testcutSize:         ',testcutSize
!     write(fileID, * ), 'testLower:           ',testLower
!     write(fileID, * ), 'testUpper:           ',testUpper
!     write(fileID, * ), ''
!     write(fileID, * ), 'Cutout parameters:'
!     write(fileID, * ), 'lower bound:         ',lBound
!     write(fileID, * ), 'upper bound:         ',uBound
!     write(fileID, * ), ''
!     write(fileID, * ),  'Number of samples'    ,samples
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
  !                          COMBINATION
  !----------------------------------------------------------------
  ! USE:  
  !     Computes nCk
  !       
  ! FORM:
  !      function combination (int n, int k)
  !         
  ! BEHAVIOR:
  ! 
  !  
  !----------------------------------------------------------------
  
  recursive function combination(num,k) result (res)
    implicit none
    !
    !    ..SCALAR ARGUMENTS..
    integer ::  num
    integer ::  k
    integer :: res
    !
    !    ..LOCAL VARIABLES..
    integer :: calls = 0

    calls  = calls + 1
    
    res = num
    
    if (num.eq.0.or.k.eq.0) then
       res = 1
    else if (calls.lt.k) then
       res = num * combination(num-1,k)
    end if

  end function combination
  


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


end module global
