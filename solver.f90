!****************************************************************
!                              SOLVER
!****************************************************************

!----------------------------------------------------------------
! USE: Core module: Performs linear algebra operations. 
!      to solve the inverse problem.
!
!
!    
! FORM: module solver
!          contains
!      subroutine printParams       [SINK]
!      subroutine cutout            [FILTER]
!      subroutine init_random_seed  [SOURCE]
!      subroutine randTrainingSet   [FILTER]
!      subroutine autonomicClosure  [FILTER] 
!      subroutine optimizedTij      [FILTER]
!      subroutine LU                [SOLVER]
!      subroutine SVD               [SOLVER]
!
!
! BEHAVIOR: 
!           
!
! STATUS : Refactoring this unit.
! 
!----------------------------------------------------------------

module solver

  use global

  ! Stencil parameters:
  integer,parameter :: stride = 1 ! Is the ratio between LES(taken as 1) and test scale
  integer,parameter :: skip = 10
  integer,parameter :: X = 1     ! Number of realizations
  integer,parameter :: n_DAMP = 1  ! Number of lambda's
  
 
  ! Bounding Box parameters:  
  integer,parameter :: box       = 252
  integer,parameter :: boxSize   = box**3
  integer,parameter :: maskSize = boxSize - M ! 512-Training points(243) = 269
  integer,parameter :: bigHalf   = ceiling(0.5*real(box) + eps) ! for 8->5
  integer,parameter :: smallHalf = floor(0.5*real(box) + eps)   ! for 8->4
  integer,parameter :: boxCenter = smallHalf * box*(box + 1) + bigHalf
  integer,parameter :: boxLower  = stride * (bigHalf - 1)
  integer,parameter :: boxUpper  = stride * (box - bigHalf)

  ! Test field parameters: 
  integer,parameter :: testSize = 1
  integer,parameter :: testcutSize = stride * (testSize + box) + 1
  integer,parameter :: testLower = stride * bigHalf + 1
  integer,parameter :: testUpper = stride * (bigHalf - 1 + testSize) + 1


  ! Cutout parameters:
  integer,parameter :: lBound = 0.5*(f_GRID - testcutSize)
  integer,parameter :: uBound = 0.5*(f_GRID + testcutSize) - 1


contains

  
  subroutine printParams()
    
    write(fileID, * ), ''
    write(fileID, * ), 'Dataset:            ', d_set
    write(fileID, * ), ''
    write(fileID, * ), 'Stencil parameters:'
    write(fileID, * ), ''
    write(fileID, * ), 'Training points :    ',M
    write(fileID, * ), 'Features :           ',N
    write(fileID, * ), ''
    !   write(fileID, * ), 'Box parameters:      '
    !   write(fileID, * ), 'Bounding box:        ',box
    !   write(fileID, * ), 'bigHalf:             ',bigHalf
    !   write(fileID, * ), 'smallHalf:           ',smallHalf
    !   write(fileID, * ), 'boxCenter:           ',boxCenter
    !   write(fileID, * ), 'boxLower:            ',boxLower
    !   write(fileID, * ), 'boxUpper:            ',boxUpper
    !   write(fileID, * ), ''
    !   write(fileID, * ), 'Test field parameters:'
    !   write(fileID, * ), 'testcutSize:         ',testcutSize
    !   write(fileID, * ), 'testLower:           ',testLower
    !   write(fileID, * ), 'testUpper:           ',testUpper
    !   write(fileID, * ), ''
    !   write(fileID, * ), 'Cutout parameters:'
    !   write(fileID, * ), 'lower bound:         ',lBound
    !   write(fileID, * ), 'upper bound:         ',uBound
    write(fileID, * ), ''
    return
  end subroutine printParams



  !****************************************************************
  !                             CUTOUT                            !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Resizes an array to a smaller cutout size. 
  !      
  !      
  ! FORM: cutout(array, n_comp)
  !       n_comp stands for number of components. 
  !      
  ! BEHAVIOR: 
  !           
  !
  ! STATUS : 
  !          
  !
  !----------------------------------------------------------------


  subroutine cutout(array, n_comp)
    implicit none
    !
    !    ..ARRAY ARGUMENTS..
    integer, intent(in) :: n_comp
    real(8), allocatable, dimension(:,:,:,:),intent(inout) :: array
    !
    !    ..WORK ARRAY..
    real(8), allocatable, dimension(:,:,:,:):: temp

    allocate (temp(n_comp,testcutSize,testcutSize,testcutSize))
    temp = array(:,lBound:uBound,lBound:uBound,lBound:uBound)
    deallocate(array)

    allocate(array(n_comp,testcutSize,testcutSize,testcutSize))
    array = temp   
    deallocate(temp)

  end subroutine cutout



  !****************************************************************
  !                       INITIALIZE RANDOM SEED                  !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Seeds random number generation process.
  !      
  !      
  ! FORM: subroutine init_random_seed
  !       
  !      
  ! BEHAVIOR: Uses system_clock and random_seed()to seed random numbers.  
  !           
  !
  ! STATUS : 
  !          
  !
  !----------------------------------------------------------------

  subroutine init_random_seed()
    implicit none
    !
    !    ..LOCAL ARRAY..
    integer, allocatable :: seed(:)
    !
    !    ..LOCAL VARS..
    integer :: i, n_seed, clock

    ! Generate random seeds to allocate seed(:)
    call RANDOM_seed(size = n_seed) 
    allocate(seed(n_seed))
    
    ! Use system clock to seed
    call SYSTEM_clock(COUNT=clock)
    seed = clock + 37 * [ (i-1, i=1,n_seed) ]
    call RANDOM_seed(PUT = seed)

    deallocate(seed)

  end subroutine init_random_seed


  !****************************************************************
  !                        RANDOM TRAINING SET                    !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Randomly samples stencil-center points in a bounding box
  !      using init_random_seed()
  !      
  !      
  ! FORM: subroutine randTrainingSet (randMask)
  !       
  !      
  ! BEHAVIOR: Stencil-center points are uniquely sampled.
  !           Creates a mask to skip cells randomly.
  !         * Always retains the center point ...
  !           So M-1 points are randomly selected ...
  !           but the center point is ALWAYS retained.
  !         
  !         
  !
  ! STATUS : 
  !          
  !
  !----------------------------------------------------------------
  
  subroutine randTrainingSet(randMask)
    implicit none
    !
    !    ..ARRAY ARGUMENT..
    integer, allocatable, intent(out) :: randMask(:)
    !
    !    ..LOCAL SCALARS..
    real(8) :: random_0_to_1 
    integer :: randomInt     !stores integer values  from 1 to 512
    !
    !    ..DEBUG..
    logical                 :: debugRandom = .FALSE.
    integer, dimension(:,:,:), allocatable :: randomMatrix
    integer                 :: i, j, k
    integer                 :: count

    
    ! Create random mask:
    allocate(randMask(maskSize))
    randMask=0; i=1          ! i starts at 1 because (M-1) points are randomly selected
    call init_random_seed()
    do
       call random_number(random_0_to_1) !returns a random value from 0 to 1
       randomInt = nint( (boxSize-1) * random_0_to_1 ) + 1 ! (511 * x) + 1
       if(any(randMask.eq.randomInt).or.(randomInt.eq.boxCenter)) cycle ! Enforce element uniqueness
       randMask(i) = randomInt       
       i = i+1   
       if (i.gt.maskSize) exit   
    end do
    
    ! DEBUG:
    if (debugRandom) then
       allocate (randomMatrix(box,box,box))
       call bubblesort(randMask)
       count=0;randomMatrix=0

       do k=1,box
       do j=1,box
       do i=1,box
          count = count+1
          if(any(randMask.eq.count)) cycle
          randomMatrix(i,j,k) = count
       enddo
       enddo
       enddo

       print*,'randomMatrix'        ! Print random training points indices:
       do k=1,box
       do i=1,box
          print*,randomMatrix(i,:,k)
       end do
       print*,''
       end do
    end if

  end subroutine randTrainingSet

  
  !****************************************************************
  !                        AUTONOMIC CLOSURE                      !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Computes Volterra coefficients h_ij
  !      
  !      
  !      
  ! FORM: subroutine autonomicClosure (u_f, u_t, tau_ij, T_ij, h_ij)
  !       
  !      
  ! BEHAVIOR: Enter only one DAMPING value.
  !           Performs entire domain computation with one lambda.
  !         * For colocated, computes h_ij and tau_ij_opt,T_ij_opt     
  !         * For non-colocated, computes h_ij; use another procedure...
  !           to compute optimized stresses. 
  !         * LAPACK routines destroy input vectors and matrices ...
  !           so T_ij is copied into T. But this can be relaxed later on.
  !       *** INDEX Nomenclature:
  !          non_col_ : for non-colocated combinations.
  !          _comp    : To select velocity and velocity products.
  !          _test    : for a point in the TEST-scale field
  !          _box     : for a point in bounding BOX (stencil-center)
  !          _stencil : for a point in the 3x3x3 STENCIL
  !          _opt     : for computed(OPTIMIZED)values at LES and Test scales.
  !          _proj    : to project stresses  for postprocessing
  !
  !          coloc    : 0  means use non-coloc
  !
  ! STATUS : Suppressed all printing and outputs the data files for images.
  !          Needs coarse grained parallization.
  !
  ! 
  !          
  !
  !----------------------------------------------------------------
  
  subroutine autonomicClosure(u_f, u_t, tau_ij, T_ij, h_ij)
    implicit none
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:,:,:), intent(in) :: u_f
    real(8), dimension(:,:,:,:), intent(in) :: u_t
    real(8), dimension(:,:,:,:), intent(in) :: tau_ij
    real(8), dimension(:,:,:,:), intent(in) :: T_ij
    real(8), dimension(:,:),     intent(out):: h_ij
    !
    !    ..LOCAL ARRAYS..
    real(8), dimension(:,:,:,:), allocatable :: uu_f, uu_t
    real(8), dimension(:,:,:,:), allocatable :: TijOpt, tau_ijOpt
    real(8), dimension(:,:),     allocatable :: V  
    real(8), dimension(:,:),     allocatable :: T  
    !
    !    .. NON-COLOCATED..
    real(8), dimension(:), allocatable :: u_n
    integer :: non_col_1, non_col_2
    !
    !    ..LOCAL SCALARS..
    real(8) :: lambda = 1.d-1
    !
    !    ..RANDOM TRAINING POINTS (STENCIL-CENTERS)..
    integer :: rand_count
    integer :: randMask(boxSize - M) !512-243=269
    !
    !    ..INDICES..
    integer :: i_test,    j_test,    k_test 
    integer :: i_box,     j_box,     k_box  
    integer :: i_stencil, j_stencil, k_stencil 
    integer :: i_opt,     j_opt,     k_opt     
    integer :: i_proj,    j_proj,    k_proj    
    integer :: row_index, col_index, row, col 
    integer :: u_comp, uu_comp 
    !
    !    ..DEBUG..
    logical :: coloc = 0    
    logical :: printval
    logical :: debug(4) = [1,0,0,0]
    ! 1 - Takes just one (tau_ij_11) for testing purpose. 
    ! 2 - Prints cutout shapes to check if cutout ops have gone well.
    ! 3 - Performs all computations on just the first point.
    ! 4 - To check if values from the stencil and bounding box are working well.


    print*, 'Computing SGS stress...'
    ! To determine stresses at coarse and fine stencil,
    ! velocities and their products must be known at LES
    ! scale. Though on the test scale it will skip every other point (coarse stencil)
    ! the fine stencil will require LES scale values to compute tau_ijOpt


    allocate(uu_t(n_uu,testcutSize,testcutSize,testcutSize)) !(6,51,51,51)
    allocate(uu_f(n_uu,testcutSize,testcutSize,testcutSize))


    if (debug(1)) then
       allocate(TijOpt(1,testSize,testSize,testSize)) !(1,17,17,17)
       allocate(tau_ijOpt(1,testSize,testSize,testSize))
    else
       allocate(TijOpt(n_uu,testSize,testSize,testSize)) !(6,17,17,17)
       allocate(tau_ijOpt(n_uu,testSize,testSize,testSize))
    end if

    allocate (V(M,N),T(M,P))
    
    TijOpt = 0.
    tau_ijOpt=0.

    if(debug(2)) then
       print*,'shape tau_ij cutout:    ',shape(tau_ij)
       print*,'shape uu_t cutout:      ',shape(uu_t)
       print*,'shape tau_ijOpt cutout: ',shape(tau_ijOpt)
    end if


    ! COLOCATED FORMULATION:
    if (coloc) then
       
!        ! Compute velocity products:(Lower triangle order for ij): Send to actools.f90
!        call velocityProducts(uu_f, u_f)
!        call velocityProducts(uu_t, u_t)
       
       
!        ! WHOLE DOMAIN COMPUTATION: 
!        do k_test = testLower, testUpper, stride 
!        do j_test = testLower, testUpper, stride
!        do i_test = testLower, testUpper, stride ! i_test = 11,43,2
       
!              row_index  = 0 
       
!              ! ENTER STENCIL-CENTER POINTS:
!              do k_box = k_test-boxLower, k_test+boxUpper, skip
!              do j_box = j_test-boxLower, j_test+boxUpper, skip
!              do i_box = i_test-boxLower, i_test+boxUpper, skip  ! i_box = 3,49,2
       
!                 col_index = 0 
!                 row_index = row_index + 1
       
!                 ! ENTER 3x3x3 STENCIL:
!                 do k_stencil = k_box-stride, k_box+stride, stride ! Vectorize using (PACK/UNPACK)
!                 do j_stencil = j_box-stride, j_box+stride, stride
!                 do i_stencil = i_box-stride, i_box+stride, stride
       
!                    ! ZERO ORDER TERMS:
!                    V(row_index, col_index) = 1.d0
       
!                    ! FIRST ORDER TERMS:
!                    do u_comp = 1,n_u ! 1 to 3 -> 3x(3x3x3) = 81
!                       col_index = col_index+1
!                       V(row_index,col_index) = u_t(u_comp,i_stencil,j_stencil,k_stencil)
!                    end do
       
!                    ! SECOND ORDER TERMS: 6x(3x3x3) = 162 (GIVES A TOTAL OF 243 TERMS)
!                    do uu_comp = 1,n_uu ! 1 to 6
!                       col_index = col_index+1
!                       V(row_index,col_index) = uu_t(uu_comp,i_stencil,j_stencil,k_stencil)
!                    end do
       
!                 end do 
!                 end do 
!                 end do ! STENCIL
       
!                 T(row_index) = T_ij(n,i_box,j_box,k_box) !Change 1 to (1-6) here.
       
!              end do
!              end do
!              end do ! BOUNDING BOX
       
       
!              ! 1) Until here, the mechanism to compute the coefficient matrix
!              ! is tested and verified. There are no issues with the basic mechanism.
!              ! Only variable strides, for more than 2 is to be adjusted in the code later.
       
!              ! 2) For testing, only the '_11' component is evaluated. Later extend to all 6 ij's
!              ! V matrix generated can be used to determine Volterra coefficients for all six ij components.
       
       
!              !call LU(V,T,h_ij,lambda)     ! Least squares by LU decomposition
       
       
!              ! Test scale computations will involve creation
!              ! of a strided matrix. Since an array stores contiguous
!              ! set of memory, two things will happen if the array index
!              ! is not carefully set - if there is a regular chuck of
!              ! heap memory, then it creates a (51,51,51) array and fills
!              ! unused space with 0. If not, it results in a
!              ! segmentation fault. Thus the index of the strided elements
!              ! that are used to form TijOpt and tau_ijOpt arrays should be
!              ! mapped with _test indices to give a contiguous array (17,17,17)
       
!              i_proj = (i_test-testLower)/stride + 1
!              j_proj = (j_test-testLower)/stride + 1
!              k_proj = (k_test-testLower)/stride + 1
       
!              ! COARSE STENCIL: Project T_ij back to itself to compare with the original T_ij field
!              ! VECTORIZE THIS PART FOR SPEED
!              col = 0
!              do k_opt = k_test-stride, k_test+stride, stride
!              do j_opt = j_test-stride, j_test+stride, stride
!              do i_opt = i_test-stride, i_test+stride, stride
       
!                 do u_comp = 1,n_u ! 1 to 3
!                    col = col+1
!                    TijOpt(1,i_proj,j_proj,k_proj) = TijOpt(1,i_proj,j_proj,k_proj)&
!                         + u_t(u_comp,i_opt,j_opt,k_opt) * h_ij(col)
!                 end do
!                 do uu_comp = 1,n_uu ! 1 to 6
!                    col = col+1
!                    TijOpt(1,i_proj,j_proj,k_proj) = TijOpt(1,i_proj,j_proj,k_proj)&
!                         + uu_t(uu_comp,i_opt,j_opt,k_opt) * h_ij(col)
!                 end do
       
!              end do
!              end do
!              end do ! coarse
       
!              ! FINE STENCIL : Calculate the SGS stress with h_ij;compare with the original tau_ij     
!              col = 0 
!              do k_opt = k_test-1, k_test+1
!              do j_opt = j_test-1, j_test+1
!              do i_opt = i_test-1, i_test+1
       
!                 do u_comp = 1,n_u ! 1 to 3
!                    col = col+1
!                    tau_ijOpt(1,i_proj,j_proj,k_proj) = tau_ijOpt(1,i_proj,j_proj,k_proj)&
!                         + u_f(u_comp,i_opt,j_opt,k_opt) * h_ij(col)
!                 end do
!                 do uu_comp = 1,n_uu ! 1 to 6
!                    col = col+1
!                    tau_ijOpt(1,i_proj,j_proj,k_proj) = tau_ijOpt(1,i_proj,j_proj,k_proj)&
!                         + uu_f(uu_comp,i_opt,j_opt,k_opt) * h_ij(col)
!                 end do
       
!              end do
!              end do
!              end do  ! fine                
       
!        end do
!        end do
!        end do ! test
       
    else

       allocate (u_n(stencil_size))

       ! WHOLE DOMAIN COMPUTATION: 
       do i_test = 129, 129, 2 
       do j_test = 129, 129, 2
       do k_test = 129, 129, 2 ! i_test = 11,43,2

          row_index  = 0 

          ! ENTER STENCIL-CENTER POINTS: C-ORDER
          do i_box = i_test-126, i_test+125, 10
          do j_box = j_test-126, j_test+125, 10
          do k_box = k_test-126, k_test+125, 10 ! i_box = 3,49,2
             ! Replace this loop with subroutine build_V()
             col_index = 0 
             row_index = row_index + 1
         
             ! ENTER 3x3x3 STENCIL: C-ORDER
             
             u_n = reshape(u_t (:, i_box-2 : i_box+2 : 2, & 
                                   j_box-2 : j_box+2 : 2, &
                                   k_box-2 : k_box+2 : 2), [stencil_size])

             ! ZERO ORDER TERMS:
             col_index = col_index + 1
             V(row_index, col_index) = 1.d0


             ! FIRST ORDER TERMS:          
             do non_col_1 = 1,stencil_size 
                col_index = col_index + 1
                V(row_index,col_index) = u_n(non_col_1)
             end do

             ! SECOND ORDER TERMS: 6x(3x3x3) = 162 (GIVES A TOTAL OF 243 TERMS)
             do non_col_1 = 1, stencil_size
             do non_col_2 = non_col_1, stencil_size
                col_index = col_index + 1
                V(row_index,col_index) = u_n(non_col_1) * u_n(non_col_2)
             end do
             end do

             T(row_index,:) = T_ij(:,i_box,j_box,k_box) !Change 1 to (1-6) here. !THIS ONE IS CORRECT; KEEP IT.

          end do
          end do
          end do ! BOUNDING BOX
          
          
          ! CHECK V:
          if (V(1500,2000).ne.2.0009431419772586d-2) then
             print*, "Error! Check sorting order in  V matrix!"
             print*, 'V(1500,2000)', V(1500,2000)
             stop
          else 
             print*,'V matrix check ... Passed'
          end if
          

          !CHECK T:
          if (T(3,1).ne.8.7759832493259110d-2)then
             print*, "Error! Check sorting order in  T vector!"
             print*, T(3,1)
             stop
          else 
             print*,'T vector check ... Passed'
          end if
          
          !
          ! CALL SOLVER
          if (linSolver.eq.'LU') then
             !Damped least squares
             call LU(V, T, h_ij, lambda)           
          elseif(linSolver.eq.'SVD') then
             ! TSVD
             call SVD(V, T, h_ij, lambda, printval) 
          else
             print*, 'Choose correct solver: LU, SVD'
             stop
          end if
          
          ! CHECK h_ij:
          if (h_ij(350,1).ne.-4.5121154730201521d-2)then
             print*, "Error! Check lambda, method or sorting order for h_ij computation:"
             print*,h_ij(350,1)
             !     stop
          else 
             print*,'SVD check ... Passed'
          end if

       end do
       end do
       end do ! test

    end if

  end subroutine autonomicClosure
  
  
  !****************************************************************
  !                        OPTIMIZE TIJ                           !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Computes tau_ij and T_ij from calculated h_ij
  !      
  !      
  
  !      
  ! FORM: subroutine optimizedTij (u_f, u_t, h_ij, TijOpt, tau_ijOpt)
  !       
  !      
  ! BEHAVIOR: TijOpt and tau_ijOpt should be allocated before. 
  !           ZERO ORDER terms are always 1.
  !           Indices are in C-order for direct comparison with MATLAB.
  !
  !     * * * INDEX Nomenclature:
  !           _test: for a point in the TEST-scale field
  !           _opt: where same h_ij is being used for point(s)
  !           non_col_: for non-colocated combinations.
  !          _comp: To select velocity and velocity products.
  !
  !
  ! STATUS : Hard coded only for non-colocated formulation.
  !          Local indices can be declared in the common namespace ...
  !          but their values get altered if passed in the argument ...
  !          list of a function call local to the procedure.       
  !
  !----------------------------------------------------------------
  
  
  subroutine optimizedTij(u_f, u_t, h_ij, TijOpt, tau_ijOpt)
    implicit none
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:,:,:), intent(in) :: u_f
    real(8), dimension(:,:,:,:), intent(in) :: u_t
    real(8), dimension(:,:),     intent(in) :: h_ij
    real(8), dimension(:,:,:,:), intent(out):: TijOpt
    real(8), dimension(:,:,:,:), intent(out):: tau_ijOpt 
    !
    !    ..LOCAL VARS..
    integer                :: col_index
    integer                :: u_comp, uu_comp 
    !
    !    ..NON-COLOCATED FORMULATION..
    real(8), dimension(stencil_size) :: u_f_stencil
    real(8), dimension(stencil_size) :: u_t_stencil
    !
    !    ..LOCAL INDICES.. 
    integer :: i_test,    j_test,    k_test 
    integer :: i_opt,     j_opt,     k_opt  
    integer :: i_stencil, j_stencil, k_stencil
    integer :: non_col_1, non_col_2  

  
    ! WHOLE DOMAIN COMPUTATION: 
    do i_test = 129, 129
    do j_test = 129, 129
    do k_test = 129, 129

       ! ENTER STENCIL-CENTER POINTS: C-ORDER
!       do i_opt = i_test-127, i_test+126
!       do j_opt = j_test-127, j_test+126
!       do k_opt = k_test-127, k_test+126
       do i_opt = 15,15
       do j_opt = 24,24
       do k_opt = 10,10

          col_index = 0 
         
          ! ENTER 3x3x3 STENCIL: C-ORDER                                      
          u_t_stencil = reshape(u_t (:, i_opt-2 : i_opt+2 : 2, &
                                        j_opt-2 : j_opt+2 : 2, &
                                        k_opt-2 : k_opt+2 : 2),  [stencil_size]) 

          u_f_stencil = reshape(u_f (:, i_opt-1 : i_opt+1 : 1, &
                                        j_opt-1 : j_opt+1 : 1, &
                                        k_opt-1 : k_opt+1 : 1),  [stencil_size]) 

          ! ZERO ORDER TERMS: WILL BE 1
          col_index = col_index + 1

          TijOpt   (:,i_opt, j_opt, k_opt) = h_ij(col_index,:) 
          tau_ijOpt(:,i_opt, j_opt, k_opt) = h_ij(col_index,:) 

          ! FIRST ORDER TERMS:             
          do non_col_1 = 1, stencil_size 
             col_index = col_index + 1

             TijOpt(:,i_opt,j_opt,k_opt) = TijOpt (:,i_opt,j_opt,k_opt)             &
                                         +                                          &
                                         (u_t_stencil(non_col_1) * h_ij(col_index,:))

             tau_ijOpt(:,i_opt,j_opt,k_opt) = tau_ijOpt(:,i_opt,j_opt,k_opt)         &
                                          +                                          &
                                          (u_f_stencil(non_col_1) * h_ij(col_index,:))
          end do

          ! SECOND ORDER TERMS: 6x(3x3x3) = 162 (GIVES A TOTAL OF 243 TERMS)
          do non_col_1 = 1, stencil_size
          do non_col_2 = non_col_1, stencil_size
             col_index = col_index+1

             TijOpt(:,i_opt,j_opt,k_opt) = TijOpt(:,i_opt,j_opt,k_opt) &
                                         +                             &
                                         (u_t_stencil(non_col_1) * u_t_stencil(non_col_2) * h_ij(col_index,:))

             tau_ijOpt(:,i_opt,j_opt,k_opt) = tau_ijOpt(:,i_opt,j_opt,k_opt) &
                                            +                                &
                                            (u_f_stencil(non_col_1) * u_f_stencil(non_col_2) * h_ij(col_index,:))
          end do
          end do

       end do
       end do
       end do ! BOUNDING BOX

    end do
    end do
    end do ! TEST

  end subroutine optimizedTij


  
  !****************************************************************
  !                                LU                             !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Implements the damped least squares method to compute h_ij
  !      Solves
  !             (V'V + lambda*I) * h_ij = (V' * T_ij)
  !      using
  !            DGEMM - BLAS LEVEL 3
  !            DSYSV - LAPACK DRIVER ROUTINE
  !
  !      
  ! FORM: subroutine LU (V, T_ij, h_ij, lambda)
  !       
  !      
  ! BEHAVIOR: 1) Input matrices are destroyed by the SV driver.
  !              #  Call by reference in function arguments. 
  !           2) Need array with defined size to pass into BLAS
  !              and LAPACK routines. Written in F77.
  !         * 3) nrhs is same as P.
  !           4) A is a symmetric matrix.        
  !           5) Check if forall() statement is safe to use in MPI.
  !
  !
  ! STATUS : Wrap with the name compute_h_ij as a procedure interface.
  !         
  !        + call DGESV(N, P, A, N, IPIV, b, N, INFO) 
  !          LWMAX is too large initially. 
  !          Check for solution accuracy.
  !
  !----------------------------------------------------------------
  
  subroutine LU (V,  T_ij, h_ij, lambda)
    implicit none

    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:), intent(in)  :: V 
    real(8), dimension(:,:), intent(in)  :: T_ij
    real(8), dimension(:,:), intent(out) :: h_ij
    !
    !    ..SCALAR ARGUMENTS..
    real(8), intent(in) :: lambda
    !
    !    ..LOCAL ARRAYS.. 
    real(8), dimension(:,:), allocatable  :: A
    real(8), dimension(:,:), allocatable  :: b
    integer, dimension(:),   allocatable  :: IPIV
    !
    !    ..DGEMM ARGUMENTS..
    real(8):: alpha = 1.d0, beta = 0.d0
    !
    !    ..DSYSV ARGUMENTS..
    real(8), dimension(:), allocatable :: WORK
    integer :: LWORK 
    integer :: LWMAX = M * N
    integer :: NRHS = P
    integer :: INFO

    !
    ! Allocate work arrays
    allocate (A(N,N), b(N,P), IPIV(N))
    allocate (WORK(LWMAX))

    ! 
    ! A(N,N) = V'(N,M) * V(M,N)
    call DGEMM('T', 'N', N, N, M, alpha, V, M, V, M, beta, A, N)

    ! Apply damping: A = A + lambda*I
    forall(i=1:N) A(i,i) = A(i,i) + lambda 

    ! b(N,P) = V'(N,M) * T_ij(M,P) 
    call DGEMM('T', 'N', N, P, M, alpha, V, M, T_ij, M, beta, b, N)

    !
    ! Solve Linear System: A(N,N) h_ij(N,P) = b(N,P)
    LWORK = -1
    call DSYSV('Lower', N, P, A, N, IPIV, b, N, WORK, LWORK, INFO) 
    LWORK = min(LWMAX, int(WORK(1)))
    call DSYSV('Lower', N, P, A, N, IPIV, b, N, WORK, LWORK, INFO) 
    ! + DGESV IS THE ALTERNATE OPTION - Stashed above.

    ! Convergence check:
    if (INFO.gt.0) then
       print*, 'Failed to converge'
       stop
    end if

    !
    ! Return h_ij
    h_ij = b

  end subroutine LU



  !****************************************************************
  !                             SVD                               !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Computes h_ij by Truncated SVD method.
  !      Computes pseudoinverse, treats singular values ...            
  !      using
  !            DGEMM  - BLAS LEVEL 3
  !            DGESVD - LAPACK DRIVER ROUTINE
  !            DGESDD - LAPACK DRIVER FOR FAST SVD
  !            DGEJSV - LAPACK EXPERT DRIVER 
  !      Can optionally view intermediate matrices for debugging.
  !      
  !
  ! FORM: subroutine SVD (A, T, h_ij, lambda, [printval])
  !
  !
  !       
  !
  ! BEHAVIOR: Input A matrix will be destroyed by DGESVD
  !           Call once per ensemble. 
  !           
  !
  ! STATUS : 
  !          Replaced slower matmul() with DGEMM. 
  !        T Matmul takes 1300 secs for (M=17576),(N=3403)
  !        + Vinv = matmul(matmul(transpose(VT),D),transpose(U)) 
  !       ++ h_ij = matmul(Vinv,T)
  !        T Psuedoinverse computation now takes about 280 secs
  !        * Can use VT(i,:) = VT(i,:) * S(i) / (S(i)**2 +lambda**2), i=1,N ...
  !          to compute VT'D first and then multiply VT'D with first N rows ...
  !          of U to get Vinv. Will require a DO-LOOP to implement in parallel.
  !       ** Run J loops for different lambda 
  !          
  ! 
  !            S 0 0 0 0 
  !       D =  0 S 0 0 0
  !            0 0 S 0 0
  !       D is a large sparse matrix!!        
  !
  !  +++ ALTERNATE FASTER METHOD TO COMPUTE PSEUDOINVERSE: GIVES ERROR - Vinv=0!!  
  !     allocate(Vinv(N,M), temp(N,N))
  !     forall (i=1:N) temp(i,:) = VT(i,:) * S(i) / (S(i)**2 + lambda**2)
  !     call DGEMM('N','T',N,M,N, alpha, temp,N, U(1:N,1:M),N, beta, Vinv, N)
  !----------------------------------------------------------------
  
  subroutine SVD (A, T, h_ij, lambda, printval)
    implicit none
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:), intent(in)  :: A
    real(8), dimension(:,:), intent(in)  :: T
    real(8), dimension(:,:), intent(out) :: h_ij
    !
    !    ..SCALAR ARGUMENTS..
    real(8), intent(in) :: lambda
    logical, intent(in), optional :: printval
    !
    !    ..LOCAL ARRAYS..
    real(8), dimension(:,:), allocatable :: U
    real(8), dimension(:,:), allocatable :: VT 
    real(8), dimension(:,:), allocatable :: D
    real(8), dimension(:),   allocatable :: S
    !
    !    .. DGESVD..
    real(8), dimension(:),allocatable :: WORK
    integer, dimension(:),allocatable :: IWORK
    integer :: LWMAX = M * N
    integer :: info, LWORK
    !
    !    ..DGEMM..
    real(8), dimension(:,:), allocatable :: temp
    real(8), dimension(:,:), allocatable :: Vinv
    real(8) :: alpha = 1.d0
    real(8) :: beta = 0.d0
    !
    !    .DEBUG..
    logical :: GESVD = 0

 
    ! 
    ! SVD DECOMPOSITION: A(M,N) = U(M,M) S([M],N) V'(N,N)
    ! DGESVD returns N diagonal entries to S(N)
    allocate (U(M,M), VT(N,N), S(N), work(LWMAX))
    LWORK = -1

    if (GESVD) then     
       call DGESVD('All','All', M, N, A, M, S, U, M, VT, N, WORK, LWORK, info)
       LWORK = min( LWMAX, int(WORK(1)) ) 
       call DGESVD('All','All', M, N, A, M, S, U, M, VT, N, WORK, LWORK,  info) !BOTTLENECK = 850 SECS
    else
       allocate(IWORK(40000))
       call DGESDD('S', M, N, A, M, S, U, M, VT, N, WORK, LWORK, IWORK, info)
       LWORK = min( LWMAX, int(WORK(1)) ) 
       call DGESDD('S', M, N, A, M, S, U, M, VT, N, WORK, LWORK, IWORK, info) !BOTTLENECK = 850 SECS
    end if

    ! Convergence check:
    if (info.gt.0) then
       print*, 'Failed to converge'
       stop
    end if

    !
    ! [RUN DIAGNOSTICS]:
    if (printval) then
       print*,'A'                      ! A MATRIX
       call printplane(A,frameLim=4)
       print*,'U'                      ! U MATRIX
       call printplane(U,frameLim=4)
       print*,'S'                      ! S VECTOR
       print*, S(1:3)
       print*,'VT'                     ! VT MATRIX
       call printplane(VT,frameLim=4)
    end if


    !
    ! Create D(N,M) with S(i) as damped diagonal entries.
    allocate(D(N,M))
    D = 0.d0
    forall(i=1:N) D(i,i) = S(i) / (S(i)**2 + lambda**2) 
    ! **  DEVELOPMENTAL NOTE IN THE HEADER
    !
    ! COMPUTE PSEUDOINVERSE: Vinv(N,M) = (VT'(N,N) * D(N,M)) * U'(M,M) 
    allocate (Vinv(N,M), temp(N,M))

    call DGEMM('T','N', N, M, N, alpha, VT,   N, D, N, beta, temp, N)
    call DGEMM('N','T', N, M, M, alpha, temp, N, U, M, beta, Vinv, N)

    deallocate (temp)
    ! + Alternate option: matmul(). Differs after 13 d.p. Stashed above.
    ! +++ Alternate method to compute pseudoinverse. Stashed above.
  
 

    ! Compute h_ij(N,P) = Vinv(N,M) * T(M,P)
    call DGEMM('N','N', N, P, M, alpha, Vinv, N, T, M, beta, h_ij, N)
    ! ++ Alternate option: matmul(). Stashed above.

    !
    ![RUN DIAGNOSTICS]:
    if(printval) then
       print*,'Vinv'                     ! Vinv MATRIX
       call printplane(Vinv,frameLim=4)
    end if

  end subroutine SVD



end module solver

