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
!          namespace
!          interface printplane
!              subroutine plane3
!              subroutine plane2
!          interface norm
!              subroutine norm3
!              subroutine norm2
!        subroutine printParams       [SINK]
!        function nCk
!        function factorial
!        function mean
!        function stdev
!        function bigHalf
!        function smallHalf
!
! BEHAVIOR:
!      * Check array rank while passing into printplane.
!      * 
!        
!      * 
!      * 
!  
! STASH:
!    boxCenterFlatIndex = smallHalf(box(1)) * box(1)*(box(1) + 1) + bigHalf(box(1)) 
!----------------------------------------------------------------

module global
  implicit none
  integer, parameter :: dp = selected_real_kind(15,307)
  real,    parameter :: eps = 1e-5
  real(8), parameter :: PI = 4.d0 * atan(1.d0)


  ! DEFINE STRING ARRAY:
  type str16
     character (16) :: name
  end type str16

  type(str16), parameter :: l_machine(2) = [str16 ('local'),     &
                                            str16 ('remote')]

  type(str16), parameter :: var(2) = [str16 ('Velocity'),      &
                                      str16 ('Pressure')]

  type(str16), parameter :: l_dataset(5) = [str16 ('nrl'),     & 
                                            str16 ('jhu256'),  & ! jhu256
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

  type(str16), parameter :: l_scheme(2) = [str16 ('local'),     &
                                           str16 ('global')]


              
  character(8) :: machine        = trim (l_machine(1) % name)
  character(8) :: dataset        = trim (l_dataset(2) % name)        ! [...,JHU, HST,...]
  logical      :: withPressure   = 0

  character(8) :: solutionMethod = trim (l_solutionMethod(1) % name) ! [LU, SVD]
  character(2) :: hst_set        = 'S3'                              ! [S1, S3, S6]
  character(3) :: stress         = 'abs'                             ! [dev, abs]
  character(16):: formulation    = trim (l_formulation(1) % name)    ! [colocated, non-colocated]
  character(8) :: trainingPoints = trim (l_trainingPoints(2) % name) ! [ordered, random]
  character(8) :: scheme         = trim (l_scheme(1) % name)         ! [local, global]
  integer      :: order          = 2                                 ! [first, second]

  !----------------------------------------------------------------
  !
  !****************************************************************
  !                       ..CONTROL SWITCHES..                    !
  !****************************************************************
  

  logical :: turbulentStats       =  0
  logical :: useTestData          =  0
  logical :: readFile             =  1
  logical :: filterVelocities     =  1
  logical :: plot_Velocities      =  1
  logical :: computeFFT_data      =  0 ! **** ALWAYS CHECK THIS ONE BEFORE A RUN **** !
  logical :: save_FFT_data        =  1

  logical :: plot_Stress          =  1
  logical :: production_Term      =  1
  logical :: save_ProductionTerm  =  1
  logical :: compute_Stress       =  0



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
  real(8), dimension(:,:,:,:), allocatable :: Sij
  real(8), dimension(:,:,:,:), allocatable :: Sij_f
  real(8), dimension(:,:,:,:), allocatable :: Sij_t
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

  real(8) :: dx                 ! To calculate gradient
  real(8) :: nu                 ! Viscosity

  integer :: M                  ! Number of training points 3x3x3x9
  integer :: N = 0
  integer :: P                  ! Number of components (1 to 6)
  integer :: stencil_size 
  integer :: LES_scale 
  integer :: test_scale 
  integer :: Freq_Nyq

  !    ..STENCIL..
  integer :: Delta_LES          ! Grid spacing at LES scale = Nyquist frequ/LES_scale  (128/40 = 3)
  integer :: Delta_test         ! Grid spacing at test scale = Nyquist freq/ test_scale (128/20 = 6)
  integer :: trainingPointSkip   ! 1 for RANDOMLY sampled training points; Some fixed value for ORDERED [REGULARIZED]set
  integer :: X                   ! Number of realizations
  integer :: n_lambda            ! Number of lambdas
  real(8) :: lambda, lambda_0(2) !lambda_0(2)
 
  !   ..BOUNDING BOX..
  integer :: box(3) 
  integer :: boxSize  
  integer :: maskSize  
 
  integer :: boxLower 
  integer :: boxUpper 
  integer :: boxFirst
  integer :: boxLast
  integer :: boxCenterSkip

  !   .. EXTENDED DOMAIN..
  integer :: extLower           ! Lower index of extended domain using ghost cells
  integer :: extUpper           ! Upper index                "

  !    ..OPT..
  integer :: optLower
  integer :: optUpper

  !   ..STATISTICS..
  integer :: n_bins            ! Number of bins
  integer :: N_cr               ! Number of cross-validation points

  !
  !    ..FILEIO..
  integer :: z_plane           
  integer :: path_txt          = 22
  integer :: params_txt        = 23
  integer :: cross_csv_T_ij    = 81
  integer :: cross_csv_tau_ij  = 82

  !
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

  character(4) :: ext   = 'h5'   ! Dataset extension: [bin]ary, [h5] or [txt] format. 

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

  !
  !    ..CASE NAME ..
  character(*), parameter :: CASE_NAME = 'scratch-col'

  !----------------------------------------------------------------
  !
  !****************************************************************
  !                       ..INTERFACES..                          !
  !****************************************************************

  interface printplane
     module procedure plane2, plane3
  end interface

  interface norm
     module procedure norm3, norm2
  end interface

contains
  
  subroutine setEnv()
    
    RES_PATH = RES_DIR

    ! GRID
    i_GRID = 256
    j_GRID = 256
    k_GRID = 256
    Freq_Nyq = i_GRID/2
    z_plane = bigHalf(k_GRID)
    
    ! TIMESTEPS:
    if (dataset.eq.'jhu256') then
       time = '256'; time_init = 256; time_incr = 1; time_final = 256
       nu = 1.85d-4
    end if

    if (dataset.eq.'nrl') time = '0460'
    if (dataset.eq.'hst') then
       nu = 5.67d-4 !6.093d-3
       if (hst_set.eq.'S1') then
!          time = '016'; time_init = 16; time_incr = 1; time_final = 28
          time = '017'; time_init = 17; time_incr = 1; time_final = 28
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
    

    n_u = 3 
    if (withPressure) n_u = n_u + 1
   
    ! CHECK PRESSURE DATA AVAILABLIITY:
    if(dataset.ne.'jhu256'.and.n_u.eq.4) then
       print*, 'Dataset ',dataset,' has no pressure data files'
       stop
    end if

    P = 6

    ! SCALE
    LES_scale  = 40
    test_scale = 20
    Delta_LES = floor(real(Freq_Nyq) / real(LES_scale))
    Delta_test = floor(real(Freq_Nyq) / real(test_scale))
    
    ! LAMBDA: !   [1.d-03, 1.d-01, 1.d+00, 1.d+01]
    lambda_0 = 1.d-01 * [1,3]
    n_lambda = 1
    

    ! FFT 
    f_GRID = 1024 !For FFT ops, the grid must be cubic.
    center = (0.5d0 * f_GRID) + 1.d0
    dx = 2.d0*PI/dble(i_GRID) !Only for JHU data. Check for others.    

    

    if (trainingPoints.eq.'ordered') then
       trainingPointSkip = 10   ! ORDERED 
    else                    
       trainingPointSkip = Delta_test    ! RANDOM   
    end if
    

    ! Bounding Box parameters:  YES ITS VERBOSE AND LOOKS DUMB BUT IT IS A NEAT LAYOUT [FOR THE TIME BEING]
    box = [1, 1, 1]             ! Initialize bounding box size
    if (formulation.eq.'noncolocated') then 
       stencil_size = n_u * (3*3*3) 
       do i = 0, order
          N = N + nCk(stencil_size - 1 + i, i) 
       end do
       if (trainingPoints.eq.'ordered') then ! NON-COLOCATED, ORDERED
          box = 256 * box
          M =    (floor((real(box(1) - 1)) / trainingPointSkip) + 1)    &
               * (floor((real(box(2) - 1)) / trainingPointSkip) + 1)    &
               * (floor((real(box(3) - 1)) / trainingPointSkip) + 1)  !17576 
       elseif (trainingPoints.eq.'random') then ! NON-COLOCATED, RANDOM
          M = 5 * N
          box = ceiling(M**(1./3.)) * trainingPointSkip * box
          boxSize   = product(box/trainingPointSkip)
          maskSize  = boxSize - M 
          X = 1     
       end if
       
    elseif  (formulation.eq.'colocated') then
       if (order == 2) then
          do i = 1, order
             n_uu = n_uu + nCk(n_u, i)    ! NUMBER OF uu, [up, pp] PRODUCTS
          end do
       end if
       N = 1 + ((n_u + n_uu) * (3*3*3))
       if (trainingPoints.eq.'ordered') then ! COLOCATED, ORDERED
          box = 256 * box
          M =    (floor((real(box(1) - 1)) / trainingPointSkip) + 1)    &
               * (floor((real(box(2) - 1)) / trainingPointSkip) + 1)    &
               * (floor((real(box(3) - 1)) / trainingPointSkip) + 1)  
       elseif (trainingPoints.eq.'random') then ! COLOCATED, RANDOM
          M = 8 * N
          box  = ceiling(M**(1./3.)) * trainingPointSkip * box
          boxSize   = product(box/trainingPointSkip)
          maskSize  = boxSize - M ! 512-Training points(243) = 269
          X = 1     
       end if
    end if
    
    boxLower  = bigHalf(box(1)) - 1
    boxUpper  = box(1) - bigHalf(box(1))

    if (scheme.eq.'local') then
       boxFirst  = 1
       boxLast   = i_GRID       ! Assuming the domain is equal in all directions

       boxCenterSkip = 1
       extLower = 1  - boxLower - Delta_test
       extUpper = i_GRID + boxLower + Delta_test ! Should be i_GRID + boxUpper + 2 [optimal but makes code complex]
       
       optLower = 0
       optUpper = 0
    else if (scheme.eq.'global') then  ! GLOBAL
       boxFirst  = bigHalf(box(1))  
       boxLast   = bigHalf(box(1)) + i_GRID - box(1)

       boxCenterSkip = box(1)
       extLower = 1  -  Delta_test 
       extUpper = i_GRID + Delta_test  

       optLower = boxLower
       optUpper = boxUpper
    end if

    ! Statistics parameters:
    n_bins = 250
    N_cr = 3   !(11x11x11)

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
    write(params_txt,*) dataset(1:3)
    write(params_txt,*) hst_set
    write(params_txt,*) LES_scale
    close(params_txt)

    if (displayOption.eq.'display') then
     write(fileID, * ), 'Dataset:             ', dataset, '\n'
     
     if (dataset.eq.'hst') then
        write(fileID, * ), 'HST set:          ', hst_set, '\n'
     end if

     write(fileID, * ), 'Scheme:              ', scheme
     write(fileID, * ), 'Formulation:         ', formulation
     write(fileID, * ), 'Pressure Term:       ', withPressure

     write(fileID, * ), 'Order:               ', order
     write(fileID, * ), 'Training points(M):  ', M, '\t', trainingPoints
     write(fileID, * ), 'Features(N):         ', N, '\n'
     write(fileID, * ), 'Filter scales:       ', LES_scale, test_scale
     write(fileID, * ), 'Delta_LES, _test:    ', Delta_LES, Delta_test, '\n'


     write(fileID, * ), 'BOX PARAMETERS:      '
     write(fileID, * ), 'Bounding box:        ', box
     write(fileID, * ), 'boxLower:            ', boxLower
     write(fileID, * ), 'boxUpper:            ', boxUpper, '\n'

     write(fileID, * ), 'boxFirst:            ', boxFirst
     write(fileID, * ), 'boxLast:             ', boxLast, '\n'
     
     write(fileID, * ), 'boxCenterSkip:       ', boxCenterSkip
     write(fileID, * ), 'trainingPointSkip:   ', trainingPointSkip
     write(fileID, '(a22,ES8.1)' ), 'lambda_0:              ', lambda_0(1)

     
     write(fileID, * ), 'extendedDomain:      ', extLower, '\t', extUpper ,'\n'     
     write(fileID, * ), 'Number of bins:      ', n_bins
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
  
  real(8) function mean(array) 
    real(8),dimension(:,:,:),intent(in):: array  
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
  
  real(8) function stdev(array)
    real(8), dimension(:,:,:),intent(in):: array
    stdev = sqrt ( sum((array - mean(array))**2) / (size(array)-1) )
  end function stdev



  !----------------------------------------------------------------
  !                          NORM3
  !----------------------------------------------------------------
  ! USE:  
  !     Computes norm-1, norm-2 or norm-Inf of a 3D array     
  !     
  ! FORM:
  !      function norm3(double array,['1', or '2', or 'Inf'])
  !         
  ! BEHAVIOR: default is norm-2
  ! 
  !  
  !----------------------------------------------------------------
  
  real(8) function norm3(array, normType)
    real(8), dimension(:,:,:),intent(in):: array
    character(*), optional, intent(in) :: normType 

    if (present(normType)) then
       ! 
       ! L1 norm
       if (normType.eq.'1') then 
          norm3 = sum(abs(array))
       !
       ! L2 norm   
       elseif (normType.eq.'2') then
          norm3 = sqrt(sum(array**2))
       !
       ! L-Inf norm    
       elseif (normType.eq.'Inf') then
          norm3 = maxval(array)
       else
          print*, 'Norms computed are 1-norm, 2-norm or Inf-norm'
       end if
    else
       ! Default: norm-2 
       norm3 = sqrt(sum(array**2)) 
    end if   
  end function norm3

  !----------------------------------------------------------------
  !                          NORM2
  !----------------------------------------------------------------
  ! USE:  
  !     Computes norm-1, norm-2 or norm-Inf of a 3D array     
  !     
  ! FORM:
  !      function norm2(double array,['1', or '2', or 'Inf'])
  !         
  ! BEHAVIOR: default is norm-2
  ! 
  !  
  !----------------------------------------------------------------
  
  real(8) function norm2(array, normType)
    real(8), dimension(:,:),intent(in):: array
    character(*), optional, intent(in) :: normType 

    if (present(normType)) then
       ! 
       ! L1 norm
       if (normType.eq.'1') then 
          norm2 = sum(abs(array))
       !
       ! L2 norm   
       elseif (normType.eq.'2') then
          norm2 = sqrt(sum(array**2))
       !
       ! L-Inf norm    
       elseif (normType.eq.'Inf') then
          norm2 = maxval(array)
       else
          print*, 'Norms computed are 1-norm, 2-norm or Inf-norm'
       end if
    else
       ! Default: norm-2 
       norm2 = sqrt(sum(array**2)) 
    end if   
  end function norm2


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
  
  integer function bigHalf(num)
    integer, intent(in) :: num
    bigHalf  = ceiling(0.5*real(num) + eps) 
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
  
  integer function smallHalf(num)
    integer, intent(in) :: num
    smallHalf  = floor(0.5*real(num) + eps) 
  end function smallHalf


  !----------------------------------------------------------------
  !                          PROGRESS BAR
  !----------------------------------------------------------------
  ! USE:  
  !     Displays progress bar at every 10th percentage point
  !     
  ! FORM:
  ! 
  !         
  ! BEHAVIOR: 
  ! 
  !  
  !----------------------------------------------------------------

  subroutine progressBar(count, count_max)
    integer, intent(in) :: count
    integer, intent(in) :: count_max
    integer, save :: progress(2)

    progress(2) = floor(real(count) / real(count_max) * 100.)
    if ((mod(progress(2),10) == 0) .and. (progress(2) >= 10) .and. (progress(2).ne.progress(1))) then
          print*, progress(2), ' %'
    end if
    progress(1) = progress(2)
  end subroutine progressBar


  !----------------------------------------------------------------
  !                          PROGRESS BAR
  !----------------------------------------------------------------
  ! USE:  
  !     Displays progress bar at every 10th percentage point
  !     
  ! FORM:
  ! 
  !         
  ! BEHAVIOR: 
  ! 
  !  
  !----------------------------------------------------------------
  
 
end module global
