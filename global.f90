module global

  real,    parameter :: eps = 1e-3 
  real(8), parameter :: pi = 4.d0 * atan(1.d0)


  integer,parameter :: GRID = 256
  integer,parameter :: M = 17576              ! Number of training points 3x3x3x9
  integer,parameter :: N = 3403
  integer,parameter :: P = 6
  
  integer,parameter :: LES_scale  = 40
  integer,parameter :: test_scale = 20


  character(3) :: stress = 'dev'


  integer  :: i, j, k
  integer  :: fileID = 6
  real(8)  :: tic, toc
  
  
  ! DIRS:
  character(*), parameter :: DATA_DIR = '../data/'
  character(*), parameter :: TEMP_DIR = '../temp/'
  character(*), parameter :: RUN_DIR  = '../run/'


  character(8) :: d_set = 'jhu256'!jhu256, nrl, sin3D   
  character(4) :: ext   = 'bin'
  

  interface printplane
     module procedure plane2, plane3
  end interface


  
contains

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
  !      * Dumps the entire z-plane (z=1) to the standard output (screen)
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
    character(8) :: digits = 'f10.4'



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
       digits = 'es23.17'
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
    !
    character(4) :: form 
    character(8) :: digits = 'f10.4'

    
    ! Check [frameLimits]:
    x_lim = size (matrix, dim=1)
    y_lim = size (matrix, dim=2)
   

    if (present(frameLim)) then
       x_lim = frameLim
       y_lim = frameLim
    end if

      
    ! Write plane to an external file with 15 digits.
    if (present(fID)) then
       fileID = fID
       digits = 'es23.17'
    end if
    
    write (form,'(i4)') y_lim     
    write (fileID,'('//trim(form)//trim(digits)//')'), (matrix(i,1:y_lim),i=1,x_lim); print*,'' ! +z-plane
    
  end subroutine plane2


  
end module global
