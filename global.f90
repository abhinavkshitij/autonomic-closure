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

  ! DEFINE LIST: KEY[CASE ID]-CASE-CASE NAME 
  type list
     character (8)  :: id
     character (16) :: key
     character (64) :: name
  end type list

  ! CASES
  type(list), parameter :: l_case(20) = [list('1a','CL14',  'colocated_local_O1_4N'),     & ! 1
                                         list('1b','CL14''','colocated_local_O1_4N_P'),   & ! 2
                                         list('2a','CL18',  'colocated_local_O1_8N'),     & ! 3
                                         list('2b','CL18''','colocated_local_O1_8N_P'),   & ! 4
                                         list('3a','CL24',  'colocated_local_O2_4N'),     & ! 5
                                         list('3b','CL24''','colocated_local_O2_4N_P'),   & ! 6
                                         list('4a','CL28',  'colocated_local_O2_8N'),     & ! 7
                                         list('4b','CL28''','colocated_local_O2_8N_P'),   & ! 8
                                         list('5a','NG2',   'noncolocated_global'),       & ! 9
                                         list('5b','NG2''', 'noncolocated_global_P'),     & ! 10
                                         list('6a','CG24',  'colocated_global_O2_4N'),    & ! 11
                                         list('6b','CG28',  'colocated_global_O2_8N'),    & ! 12
                                         list('7a','CL1(3)','colocated_point_O1_3'),      & ! 13
                                         list('7b','CL2(3)','colocated_point_O2_3'),      & ! 14
                                         list('8a','CL1(5)','colocated_point_O1_5'),      & ! 15
                                         list('8b','CL2(5)','colocated_point_O2_5'),      & ! 16
                                         list('9a','CL1(7)','colocated_point_O1_7'),      & ! 17
                                         list('9b','CL2(7)','colocated_point_O2_7'),      & ! 18
                                         list('3a(o)','CL24o','colocated_local_O2_4N_O'), & ! 19
                                          list('3a(t)','CL24t','colocated_local_O2_4N_temp')]  ! 20


  ! MACHINE TYPE 
  type(str16), parameter :: l_machine(2) = [str16 ('local'),             &
                                            str16 ('remote')]

  ! BASE VARIABLES - VELOCITY, PRESSURE
  type(str16), parameter :: var(2) = [str16 ('Velocity'),                &
                                      str16 ('Pressure')]

  ! DATASET
  type(str16), parameter :: l_dataset(6) = [str16 ('nrl'),               & 
                                            str16 ('jhu256'),            & ! jhu256
                                            str16 ('hst'),               &
                                            str16 ('sin3D'),             &
                                            str16 ('jhu1024'),           &
                                            str16 ('jhu42')]

  ! FILTERED VARIABLES 
  type(str16), parameter :: var_FFT(4) = [str16 ('u_f'),                 &
                                          str16 ('u_t'),                 &
                                          str16 ('tau_ij'),              &
                                          str16 ('T_ij')]  
                                      
  ! SOLUTION METHOD
  type(str16), parameter :: l_solutionMethod(2) = [str16 ('LU'),         &
                                                   str16 ('SVD')]

  ! TRAINING POINTS
  type(str16), parameter :: l_trainingPoints(2) = [str16 ('ordered'),    &
                                                   str16 ('random')] 

  ! FORMULATION - COLOCATED OR NONCOLOCATED 
  type(str16), parameter :: l_formulation(2) = [str16 ('colocated'),     &
                                                str16 ('noncolocated')]

  ! SCHEME - LOCAL[OVERLAP] OR GLOBAL[NO OVERLAP] 
  type(str16), parameter :: l_scheme(2) = [str16 ('local'),              &
                                           str16 ('global')]

  ! COMPUTING DOMAIN
  type(str16), parameter :: l_compDomain(2) = [str16 ('all'),            &
                                               str16 ('plane')]

  ! ROTATION AXIS 
  type(str16), parameter :: l_rotationAxis(3) = [str16('none'),          &
                                                 str16('X'),             &
                                                 str16('Y')]


  ! ROTATION PLANE NAME 
  type(str16), parameter :: l_rotationPlane(3) = [str16('z_plane/'),      &
                                                  str16('rotateX/'),      &
                                                  str16('rotateY/')]

  ! TEST SCALE FILTER TYPE    
  type(str16), parameter :: l_filterType(8) = [str16('Sharp'),            &
                                               str16('Gauss'),            &
                                               str16('Box'),              &
                                               str16('Tri'),              &
                                               str16('GaussBox'),         &
                                               str16('GaussTri'),         &
                                               str16('BoxTri'),           &
                                               str16('All')] 

  !*****************************************************************
              
  character(8) :: machine        = trim (l_machine(1) % name)        ! [local, remote]
  character(8) :: dataset        = trim (l_dataset(6) % name)        ! [...,JHU[2], HST[3],...]
  logical      :: withPressure   = 0                                 ! [pressure[1], no pressure[0]]

  integer      :: case_idx       = 20                                 ! [1 - CL14, ...]          
  character(8) :: solutionMethod = trim (l_solutionMethod(1) % name) ! [LU, SVD]
  character(2) :: hst_set        = 'S6'                              ! [S1, S3, S6]
  character(3) :: stress         = 'abs'                             ! [dev[DS], abs[BD]]
  character(16):: formulation    = trim (l_formulation(1) % name)    ! [colocated, non-colocated]
  character(8) :: trainingPoints = trim (l_trainingPoints(2) % name) ! [ordered, random]
  character(8) :: scheme         = trim (l_scheme(1) % name)         ! [local, global]
  integer      :: order          = 2                                 ! [first, second]
  character(8) :: compDomain     = trim (l_compDomain(2) % name)     ! [all, plane]
  character(8) :: rotationAxis   = trim(l_rotationAxis(1) % name)    ! [none:z, X:y, Y:x]
  character(8) :: rotationPlane  = trim(l_rotationPlane(1) % name)   ! [none:z, X:y, Y:x]
  integer      :: M_N_ratio      = 4

  character(8) :: LESfilterType  = trim(l_filterType(1) % name)      ! [Sharp,Gauss,Box,Tri]
  character(8) :: TestfilterType = trim(l_filterType(1) % name)      ! [Sharp,Gauss,Box,Tri]
                                                                     ! [GaussBox, GaussTri, BoxTri]
                                                                     ! [All]
 
  
  real(8), parameter :: lambda_0(1) =  1.d-03
!  real(8), parameter :: lambda_0(2) =  [1.d-03, 1.d-01]!, 1.d-01,  1.d+01]


  character(48) :: CASE_NAME
!  character(*), parameter :: CASE_NAME = 'scratch'
!  character(*), parameter :: CASE_NAME = 'z_plane/43/colocated_global'
  character(16) :: z_plane_name



  
  !----------------------------------------------------------------
  !
  !****************************************************************
  !                       ..CONTROL SWITCHES..                    !
  !****************************************************************
  

  logical :: turbulentStats       =  0
  logical :: useTestData          =  0
  logical :: readFile             =  1
  logical :: filterVelocities     =  1
  logical :: plot_Velocities      =  0
  logical :: computeFFT_data      =  1! **** ALWAYS CHECK THIS ONE BEFORE A RUN **** !
  logical :: save_FFT_data        =  0

  logical :: computeDS            =  0
  logical :: compute_vorticity    =  0
  logical :: plot_Stress          =  0
  logical :: production_Term      =  0
  logical :: save_ProductionTerm  =  0
  logical :: compute_Stress       =  0

  logical :: make_Deviatoric      =  0
  logical :: multiFilter          =  0



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
  !    ..VORTICITY..
  real(8), dimension(:,:,:,:), allocatable :: omega
  !
  !    ..STRESSES..
  real(8), dimension(:,:,:,:), allocatable :: tau_ij
  real(8), dimension(:,:,:,:), allocatable :: T_ij
  real(8), dimension(:,:,:,:), allocatable :: T_ijOpt
  real(8), dimension(:,:,:,:), allocatable :: tau_ijOpt

  !     ..3-FILTER MODIFICATION..
  real(8), dimension(:,:,:,:), allocatable :: u_tB
  real(8), dimension(:,:,:,:), allocatable :: T_ijB
  real(8), dimension(:,:,:,:), allocatable :: u_tG
  real(8), dimension(:,:,:,:), allocatable :: T_ijG

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
  !
  !     ..DYNAMIC SMAGORINSKY..
  real(8), dimension(:,:,:,:), allocatable :: S_f_Sij_f_t
  real(8), dimension(:,:,:,:), allocatable :: tau_DS
  real(8), dimension(:,:,:), allocatable :: Pij_DS
  !
  !     ..BARDINA..
  real(8), dimension(:,:,:,:), allocatable :: tau_BD
  real(8), dimension(:,:,:), allocatable :: Pij_BD

 

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
  integer :: n_filter = 1       ! Number of filters

  !    ..STENCIL..
  integer :: Delta_LES          ! Grid spacing at LES scale = Nyquist frequ/LES_scale  (128/40 = 3)
  integer :: Delta_test         ! Grid spacing at test scale = Nyquist freq/ test_scale (128/20 = 6)
  integer :: trainingPointSkip   ! 1 for RANDOMLY sampled training points; Some fixed value for ORDERED [REGULARIZED]set
  integer :: X                   ! Number of realizations
  integer :: n_lambda            ! Number of lambdas
  real(8) :: lambda
 
  !   ..BOUNDING BOX..
  integer :: box(3)             ! Box size at DNS grid (256,30,...)
  integer :: boxSize            ! Number of available trainingPoints (343,216,...)
  integer :: box_effective      ! Box size on test grid
  integer :: maskSize  
 
  integer :: boxLower 
  integer :: boxUpper 
  integer :: boxFirst
  integer :: boxLast
  integer :: boxCenterSkip
  
  integer :: boxCenter


  !   .. EXTENDED DOMAIN..
  integer :: extLower           ! Lower index of extended domain using ghost cells
  integer :: extUpper           ! Upper index
  integer :: n_extendedLayers   ! Ghost cells size on each side

  !   .. BLOCK Z-DIR ..         ! FOR EXTENDED VARS - u_f, u_t, T_ij
  integer :: z_extLower         ! Block below z_plane
  integer :: z_extUpper         ! Block above z_plane
  
  !   .. PLANE Z-DIR ..         ! FOR NON-EXTENDED VARS - tau_ij. Sij_f, Sij_t, Pij_f, Pij_t, tau_ijOpt, T_ijOpt, Pij_fOpt, Pij_tOpt
  integer :: zLower             ![ALL] (1:k_GRID) OR [PLANE](z_plane:z_plane) 
  integer :: zUpper             ![ALL] (1:k_GRID) OR [PLANE](z_plane:z_plane) 

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
  
  character(96) :: DATA_PATH
  character(96) :: TEMP_PATH
  character(96) :: RES_PATH

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
  integer  :: i, j, k, l
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

  interface norm
     module procedure norm3, norm2
  end interface

contains


  !****************************************************************
  !                          SET ENVIRONMENT                      !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Set global environment based on the initial settings
  !      
  !      
  ! FORM: subroutine setEnv()
  ! 
  !      
  ! BEHAVIOR: 
  !           
  !
  ! STATUS : 
  !  **        
  !  box must be set such that its cube(boxSize) > M
  !
  ! Also include a "DILATION FACTOR" to make boxSize 
  ! larger than M by a factor f_D
  ! e.g. M = 960, box = 60 will give boxSize = 1000
  ! There are just enough points in the box available to fit in 960 
  ! training points. f_D = 1000/960 = 1.04 (will ideally remain close to 1)  
  ! Now I can set M = 960, box = 72, such that boxSize = 1728
  ! and f_D = 1728/960 = 1.8. I have almost 2x more avaiable points in the 
  ! bounding box. Information density(D_i) can be deifned the reciprocal of f_D.
  ! M = 960, box = 60; boxSize = (60/6)^3 = 1000; f_D = 1.04; D_i = 0.96
  ! M = 960, box = 72; boxSize = (72/6)^3 = 1728; f_D = 1.8;  D_i = 0.55
  ! CL14'  - 436/512   = 0.85   CL14 - 328/343   = 0.956
  ! CL18'  - 872/1000  = 0.872  CL18 - 656/729   = 0.899
  ! CL24' - 1516/1728  = 0.877  CL24 - 976/1000  = 0.976
  ! CL28' - 3032/3375  = 0.898  CL28 - 1952/2197 = 0.888
  ! NG2'  - 17576/74088 = 0.23 = NG2  
  !
  !----------------------------------------------------------------

  
  subroutine setEnv()


    ! SET PATH:    
    RES_PATH = RES_DIR

    ! GRID:
!    i_GRID = 256;    j_GRID = 256;    k_GRID = 256
    i_GRID = 42;    j_GRID = 42;    k_GRID = 42
    
    Freq_Nyq = i_GRID/2

    ! CASE_NAME:
    z_plane = 23!bigHalf(k_GRID) [43, 129, 212]
    write(z_plane_name,'(i0)'), z_plane
    if (case_idx == 0) then
      CASE_NAME = 'scratch-col'
    else
    !CASE_NAME = trim(l_rotationPlane(1)%name)//trim(z_plane_name)//'/'//trim(l_case(case_idx)%name)
      CASE_NAME = trim(l_case(case_idx)%name)
    
    end if

    ! SPACING:[EQUIDISTANT IN X,Y,Z-DIR]
    dx = 2.d0*PI/dble(i_GRID) 
    
    ! TIMESTEPS:
    if (dataset.eq.'jhu256' .or. dataset.eq.'sin3D') then
      time = '256'; time_init = 256; time_incr = 1; time_final = 256
    elseif (dataset.eq.'jhu42') then
       time = '42'; time_init = 42; time_incr = 1; time_final = 42
       nu = 1.85d-4
    end if
    if (dataset.eq.'nrl') time = '0460'
    if (dataset.eq.'hst') then
       nu = 5.67d-4 !6.093d-3
       if (hst_set.eq.'S1') then
          time = '016'; time_init = 16; time_incr = 1; time_final = 28
       end if
       if (hst_set.eq.'S3') then
          time = '015'; time_init = 15; time_incr = 1; time_final = 24
       end if
       if (hst_set.eq.'S6') then
          time = '100'; time_init = 100; time_incr = 5; time_final = 100
       end if
    end if
    read(time,*) time_init
    if (machine.eq.'local') time_final = time_init
    
    ! VELOCITY [PRESSURE] COMPONENTS: 
    n_u = 3 
    if (withPressure) n_u = n_u + 1
   
    ! CHECK PRESSURE DATA AVAILABLIITY:
    if(dataset.ne.'jhu256'.and.n_u.eq.4) then
       print*, 'Dataset ',dataset,' has no pressure data files'
       stop
    end if
    P = 6

    ! SCALE
    if (dataset.eq.'jhu256') then
       LES_scale  = 40;    test_scale = 20
    elseif (dataset.eq.'jhu42') then
       LES_scale  = 20;    test_scale = 10   
    else if (dataset.eq.'hst'.and.hst_set.eq.'S6') then
       LES_scale  = 20;    test_scale = 10
    else if (dataset.eq.'hst'.and.hst_set.eq.'S1') then
       LES_scale  = 16;    test_scale = 8
    end if

    Delta_LES = floor(real(Freq_Nyq) / real(LES_scale))
    Delta_test = floor(real(Freq_Nyq) / real(test_scale))    

    ! FFT GRID: Set f_GRID[CUBIC]
    f_GRID = i_GRID!256 
    center = (0.5d0 * f_GRID) + 1.d0

    ! STENCIL: Set N  [noncolocated/colocated]
    if (formulation.eq.'noncolocated') then 
       stencil_size = n_u * (3*3*3) 
       do i = 0, order
          N = N + nCk(stencil_size - 1 + i, i) 
       end do
    elseif  (formulation.eq.'colocated') then
       if (order == 2) then
          do i = 1, order
             n_uu = n_uu + nCk(n_u, i)    ! NUMBER OF uu, [up, pp] PRODUCTS
          end do
       end if
!       N = 1 + ((n_u + n_uu) * (3*3*3))
       N = 163
    end if

    ! BOUNDING BOX: Set trainingPointSkip, box, M [ordered/random]
    box = [1, 1, 1]             
    ! ORDERED
    if (trainingPoints.eq.'ordered') then 
       box = 256 * box ! Default = 256 [Initial value]
              ! change here for CP(*)3,5,7 = 18,30,42 
       trainingPointSkip = 10   ! Default=10 [M=17576]
              ! Change here for CP(*)3,5,7 = 8,,10like cases
       M =    (floor((real(box(1) - 1)) / trainingPointSkip) + 1)    &
            * (floor((real(box(2) - 1)) / trainingPointSkip) + 1)    &
            * (floor((real(box(3) - 1)) / trainingPointSkip) + 1)  
    ! RANDOM
    elseif (trainingPoints.eq.'random') then 
!       M = M_N_ratio * N ! Default
       M = 1000      ! change here for CP(*)3,5,7 = 27,64,125
       trainingPointSkip = Delta_test

       if (scheme.eq.'global') then
          box = 42 * box ! ** Initial value for box; change this part.
          boxSize = product(box/trainingPointSkip+1)

       elseif (scheme.eq.'local') then
! (a) Box is either inflated by M training points  [RANDOM] - default
          box = ceiling(M**(1.d0/3.d0)) * trainingPointSkip * box
! (b) Or given a predefined size (Used in CP2(3), ...  cases) [ORDERED]
!     Specify trainingPointSkip under the ORDERED section.
!          box = 3 * trainingPointSkip * box ! ADD DILATION FACTOR

          boxSize   = product(box/trainingPointSkip)
       end if
       maskSize  = boxSize - M
       X = 1     
    end if

    box_effective  = box(1)/trainingPointSkip

    boxLower  = bigHalf(box(1)) - 1
    boxUpper  = box(1) - bigHalf(box(1))

    ! BOUNDING BOX:         Set boxCenterSkip 
    ! EXTENDED DOMAIN[X,Y]: Set extUpper, extLower 
    ! COMPUTED STRESS:      Set optLower, optUpper  [local/global]
    if (scheme.eq.'local') then
       boxFirst  = 1
       boxLast   = i_GRID ! Domain equal in X,Y,Z-dir
       boxCenterSkip = 1

       extLower = 1 - boxLower - Delta_test
       extUpper = i_GRID + boxLower + Delta_test ! Should be i_GRID + boxUpper + 2 [optimal but makes code complex]
       n_extendedLayers = boxLower + Delta_test
       
       optLower = 0
       optUpper = 0
    else if (scheme.eq.'global') then  
       boxFirst  = bigHalf(box(1))  
       boxLast   = bigHalf(box(1)) + i_GRID - box(1)
       boxCenterSkip = box(1)

       extLower = 1 - Delta_test 
       extUpper = i_GRID + Delta_test  
       n_extendedLayers = Delta_test

       optLower = boxLower
       optUpper = boxUpper
    end if
!%%%%
    boxCenter= smallHalf(box(1))*box(1)*(box(1)+1) + bigHalf(box(1))

    ! EXTENDED DOMAIN[Z]: Set z_extUpper, z_extLower 
    if (compDomain.eq.'all') then
       zLower = 1
       zUpper = k_GRID
       z_extLower = extLower
       z_extUpper = extUpper
    else if (compDomain.eq.'plane') then
       zLower = z_plane
       zUpper = z_plane
       z_extLower = z_plane - boxLower - Delta_test
       z_extUpper = z_plane + boxUpper + Delta_test
    end if

    ! LAMBDA: 
    !   Cross-validation points:
    !    lambda_0 = 1.d-03 * [1,3]
    !    lambda = lambda_0(1)
    !   n_lambda = 1! size(lambda_0)
    
    ! STATISTICS:
    n_bins = 250
    N_cr = 1   ! Number of cross-validation points in each dir (11x11x11)

    ! 3-Filter MODIFICATION:
    if (multiFilter) n_filter = 2
    M = n_filter * M

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
       write(fileID, *) , '\n\n\n\n'
       write(fileID, *) , '****************************************************************','\n'
       
       call system ('date')
       write(fileID, * ), 'Dataset:             ', dataset
     
     if (dataset.eq.'hst') then
        write(fileID, * ), 'HST set:             ', hst_set
     end if
     write(fileID, * ), 'Time step:           ', time, '\n'
     write(fileID, * ), 'Computation domain:  ', compDomain, z_plane
     write(fileID, * ), 'Rotation axis:       ', rotationAxis

     write(fileID, * ), 'Scheme:              ', scheme
     write(fileID, * ), 'Formulation:         ', formulation
     write(fileID, * ), 'Pressure Term:       ', withPressure

     write(fileID, * ), 'Order:               ', order
     write(fileID, '(a22,f5.2)' ), 'M:N                  ', real(M)/real(N)/real(n_filter)
     write(fileID, * ), 'Training points(M):  ', M/n_filter, '\t', trainingPoints
     write(fileID, * ), 'Features(N):         ', N, '\n'

     write(fileID, * ), 'Filter scales:       ', LES_scale,  test_scale
     write(fileID, * ), 'Filter types:        ', LESfilterType, TestfilterType
     write(fileID, * ), 'multiFilter:         ', multiFilter

     write(fileID, * ), 'Delta_LES, _test:    ', Delta_LES, Delta_test, '\n'
     write(fileID, * ), 'Stress computation:  ', stress, '\n'
     write(fileID, * ), 'Convert to deviatoric:', make_Deviatoric, '\n'


     write(fileID, * ), 'BOX PARAMETERS:      '
     write(fileID, * ), 'Bounding box:        ', box
     write(fileID, * ), 'Box size:            ', box_effective
     write(fileID, * ), 'boxLower:            ', boxLower
     write(fileID, * ), 'boxUpper:            ', boxUpper, '\n'

     write(fileID, * ), 'boxFirst:            ', boxFirst
     write(fileID, * ), 'boxLast:             ', boxLast, '\n'

     write(fileID, * ), 'extLower:            ', extLower
     write(fileID, * ), 'extUpper:            ', extUpper, '\n'
     write(fileID, * ), 'z_extLower:          ', z_extLower
     write(fileID, * ), 'z_extUpper:          ', z_extUpper, '\n'
     

     write(fileID, * ), 'boxCenterSkip:       ', boxCenterSkip
     write(fileID, * ), 'trainingPointSkip:   ', trainingPointSkip
     write(fileID, '(a24,ES8.1)' ), 'lambda:              ', lambda_0(1)

     
     write(fileID, * ), 'extendedDomain:      ', extLower, '\t', extUpper ,'\n'      
     write(fileID, * ), 'Number of bins:      ', n_bins, '\n'
     write(fileID, *) , '****************************************************************','\n'
  end if
  end subroutine printParams

  !****************************************************************
  !                         MEMORY REQ.
  !****************************************************************

  !----------------------------------------------------------------
  ! USE: Display memory requirements
  !      
  ! FORM: subroutine memRequirement()
  !       
  !----------------------------------------------------------------

  subroutine memRequirement()
    real :: form0 ! 256^3 domain 
    real :: form1 ! resized in x,y,z-direction for temp var 
    real :: form2 ! extended domain after resizing 
    real :: plane ! planar slice (256^2)

    real :: mem = 0.0
    real :: mem_1, mem_2
    real :: mem_temp
    real :: mem_tau_ij, mem_T_ij, mem_u_f, mem_u_t
    real :: mem_Sij_f, mem_Sij_t, mem_Pij_f, mem_Pij_t
    real :: mem_Tij_Opt, mem_tau_ij_Opt, mem_Pij_f_Opt, mem_Pij_t_Opt
    real :: mem_h_ij, mem_V, mem_T
    real :: mem_uu_f, mem_uu_t
    real :: mem_A, mem_b, mem_WORK, mem_IPIV
    real, parameter  :: mem_real8 = 8.0
    real, parameter  :: inkB = 1./real(1024)
    real, parameter  :: inMB = 1./real(1024**2)
    real, parameter  :: inGB = 1./real(1024**3)

    form0 = real(i_GRID * j_GRID * k_GRID) * mem_real8 * inMB
    form1 = real((extUpper - extLower + 1)**2 * (max(k_GRID,z_extUpper) - min(1,z_extLower) + 1)) &
          & * mem_real8 * inMB
    form2 = real((extUpper - extLower + 1)**2 * (z_extUpper - z_extLower + 1)) * mem_real8 * inMB
    plane = real(i_GRID * j_GRID) * mem_real8 * inMB
    
    print*
    print*, 'Memory requirements' 

    ! Step 1: Load FFT data
    mem_u_f       = 3.0 * form0 
    mem_u_t       = 3.0 * form0 
    mem_tau_ij    = 6.0 * form0 
    mem_T_ij      = 6.0 * form0 

    mem = mem + (mem_u_f + mem_u_t + mem_tau_ij + mem_T_ij)
    print('(a24,f8.1,a4)'), 'Load FFT data:', mem, 'MB'

    ! Step 2: Allocate planar arrays 
    mem_Sij_f     = 6.0 * plane
    mem_Sij_t     = 6.0 * plane
    mem_Pij_f     = 6.0 * plane
    mem_Pij_t     = 6.0 * plane
    mem_Tij_Opt   = 6.0 * plane
    mem_tau_ij_Opt= 6.0 * plane
    mem_Pij_f_Opt = 6.0 * plane
    mem_Pij_t_Opt = 6.0 * plane

    mem_h_ij      = real(N*P) * mem_real8 * inMB

    mem = mem + (mem_Sij_f + mem_Sij_t + mem_Pij_f + mem_Pij_t + &
                 mem_Tij_Opt + mem_tau_ij_Opt + mem_Pij_f_Opt + mem_Pij_t_Opt + &
                 mem_h_ij)

    print('(a24,f8.1,a4)'), 'Allocate planar arrays:', mem, 'MB'

    ! Step 3: Extend domain
    mem_temp      = 6.0 * form1
    mem           = mem + mem_temp
    print('(a24,f8.1,a4)'), 'Extend domain(w/ temp):', mem, 'MB'

    mem = mem - mem_T_ij - mem_u_t - mem_u_f - mem_temp
    mem_T_ij      = 6.0 * form2
    mem_u_t       = 3.0 * form2
    mem_u_f       = 3.0 * form2
    mem = mem + mem_T_ij + mem_u_t + mem_u_f
    print('(a24,f8.1,a4)'), 'Extend domain:', mem, 'MB'

    ! Step 4: With colocated formulation
    if (formulation == 'colocated') then
       mem_uu_f      = 6.0 * form2
       mem_uu_t      = 6.0 * form2
       mem = mem + mem_uu_f + mem_uu_t
       print('(a24,f8.1,a4)'), 'Colocated vel-products:', mem, 'MB'
    end if

    ! Step 5: Autonomic Closure
    ! a) Build V
    mem_V         = real(M * N) * mem_real8 * inMB
    mem_T         = real(N * P) * mem_real8 * inMB
    mem_1 = mem_V + mem_T

    print('(a36)'), '*******  AUTONOMIC CLOSURE  *******'
    print('(a13)'), 'a) Build V,T:'
    print('(a8,f8.1,a4)'), 'V:', mem_V, 'MB'
    print('(a8,f8.1,a4)'), 'T:', mem_T * real(1024.0), 'kB'
    print('(a9)'), '\t ---------'
    print('(a8,f8.3,a4)'), 'a):', mem_1, 'MB'

    mem = mem + mem_1
    print('(a24,f8.1,a4)'), '', mem, 'MB'

    ! b) Invert Linear System 
    mem_A         = real(N * N) * mem_real8 * inMB
    mem_b         = real(N * P) * mem_real8 * inMB
    mem_WORK      = real(M * N) * mem_real8 * inMB
    mem_IPIV      = real(N)     * mem_real8 * inMB

    mem_2 = mem_A + mem_b + mem_WORK + mem_IPIV

    print('(a24)'), 'b) Invert Linear System:'
    print('(a8,f8.1,a4)'), 'A:',    mem_A    * real(1024.0), 'kB'
    print('(a8,f8.1,a4)'), 'b:',    mem_b    * real(1024.0), 'kB'
    print('(a8,f8.1,a4)'), 'WORK:', mem_WORK * real(1024.0), 'kB'
    print('(a8,f8.1,a4)'), 'IPIV:', mem_IPIV * real(1024.0), 'kB'
    print('(a12)'), '\t ---------'
    print('(a8,f8.3,a4)'), 'b):', mem_2, 'MB'

    mem = mem + mem_2
    print('(a24,f8.1,a4)'), '', mem, 'MB'

!    print*, l_case(case_idx)%id, l_case(case_idx)%key, l_case(case_idx)%name
!    print*, CASE_NAME

    ! Write to file:
!    open(44,file = trim(RES_PATH)//'mem.dat',action='write',access='append',status='old')
!    write(44,*) l_case(case_idx)%id, l_case(case_idx)%key, M,N,mem_1, mem_2
!    close(44)

  end subroutine memRequirement

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

    progress(2) = floor(real(count) / real(count_max) * 100.d0)
    if ((mod(progress(2),10) == 0) .and. (progress(2) >= 10) .and. (progress(2).ne.progress(1))) then
          print*, progress(2), ' %'
    end if
    progress(1) = progress(2)
  end subroutine progressBar


 
end module global
