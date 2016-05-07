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


  type(str16), parameter :: var(2) = [ str16('Velocity'), str16('Pressure')]
  type(str16), parameter :: dataset(4) = [ str16('nrl'),  str16('jhu256'), str16('sin3D'), str16('jhu1024')]
  type(str16), parameter :: var_FFT(4) = [ str16('u_f'),  str16('u_t'),    str16('tau_ij'), str16('T_ij')]  
  type(str16), parameter :: solv(2) = [ str16('LU'), str16('SVD') ]

  
  integer, parameter :: i_GRID = 256
  integer, parameter :: j_GRID = 256
  integer, parameter :: k_GRID = 256
  integer, parameter :: f_GRID = 256 !For FFT ops, the grid must be cubic.


  integer, parameter :: M = 17576              ! Number of training points 3x3x3x9
  integer, parameter :: N = 3403
  integer, parameter :: P = 6
  integer, parameter :: stencil_size = 3*(3*3*3) !3 * 27 = 81  
  integer, parameter :: LES_scale  = 40
  integer, parameter :: test_scale = 20

  character(3) :: stress = 'dev'
  character(8) :: linSolver = trim(solv(1)%name)

  
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
  !
  !    ..FILEIO..
  !
  !    .. DIRS..
  character(*), parameter :: DATA_DIR = '../data/'
  character(*), parameter :: TEMP_DIR = '../temp/'
  character(*), parameter :: RES_DIR  = '../results/'
  

  character(8) :: d_set = trim (dataset(2) % name)
  character(4) :: ext   = 'bin'   ! Dataset extension: [bin]ary or [h5] format. 

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
    

  

  end subroutine setEnv


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
