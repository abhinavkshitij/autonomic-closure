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

  type(str16), parameter :: dataset(5) = [str16 ('nrl'),     & 
                                          str16 ('jhu256'),  &
                                          str16 ('hst'),     &
                                          str16 ('sin3D'),   &
                                          str16 ('jhu1024')]

  type(str16), parameter :: var_FFT(4) = [str16 ('u_f'),     &
                                          str16 ('u_t'),     &
                                          str16 ('tau_ij'),  &
                                          str16 ('T_ij')]  

  type(str16), parameter :: solv(2) = [str16 ('LU'),         &
                                       str16 ('SVD')]


  character(8) :: d_set     = trim (dataset(2) % name)
  character(8) :: linSolver = trim (solv(2) % name)
  character(2) :: hst_set = 'S6' ! HST datasets - S1, S3, S6
 

  integer, parameter :: M = 17576              ! Number of training points 3x3x3x9
  integer, parameter :: N = 3403
  integer, parameter :: P = 6
  integer, parameter :: stencil_size = 3*(3*3*3) !3 * 27 = 81  
  integer, parameter :: LES_scale  = 40
  integer, parameter :: test_scale = 20

  character(3) :: stress = 'dev'


  
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
  real(8), dimension(:,:,:,:), allocatable :: TijOpt
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


  ! Stencil parameters:
  integer :: stride ! Is the ratio between LES(taken as 1) and test scale
  integer :: skip 
  integer :: X ! Number of realizations
  integer :: n_DAMP  ! Number of lambdas
  
 
  ! Bounding Box parameters:  
  integer :: box 
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

  !
  !    ..INDICES..
  integer  :: i, j, k
  integer  :: fileID = 6
  integer  :: n_u = 3 , n_uu = 6
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
    
    i_GRID = 256
    j_GRID = 256
    k_GRID = 256
    if (d_set.eq.'hst') then
       j_GRID = 129
    end if

    ! FFT 
    f_GRID = 256 !For FFT ops, the grid must be cubic.
    center = (0.5d0 * f_GRID) + 1.d0
    dx = 2.d0*pi/dble(i_GRID) !Only for JHU data. Change for others.    


  ! Stencil parameters:
  stride = 1 
  skip = 10
  X = 1     
  n_DAMP = 1 
  
 
  ! Bounding Box parameters:  
  box       = 252
  boxSize   = box**3
  maskSize = boxSize - M ! 512-Training points(243) = 269
  bigHalf   = ceiling(0.5*real(box) + eps) ! for 8->5
  smallHalf = floor(0.5*real(box) + eps)   ! for 8->4
  boxCenter = smallHalf * box*(box + 1) + bigHalf
  boxLower  = stride * (bigHalf - 1)
  boxUpper  = stride * (box - bigHalf)

  ! Test field parameters: 
  testSize = 1
  testcutSize = stride * (testSize + box) + 1
  testLower = stride * bigHalf + 1
  testUpper = stride * (bigHalf - 1 + testSize) + 1


  ! Cutout parameters:
  lBound = 0.5*(f_GRID - testcutSize)
  uBound = 0.5*(f_GRID + testcutSize) - 1


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
    write(23,*) d_set
    write(23,*) hst_set
    close(23)

    
     write(fileID, * ), ''
     write(fileID, * ), 'Dataset:            ', d_set
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

  
end module global
