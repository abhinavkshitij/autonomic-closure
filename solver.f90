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
!      subroutine cutout            [FILTER]
!      subroutine init_random_seed  [SOURCE]
!      subroutine randTrainingSet   [FILTER]
!      subroutine autonomicClosure  [FILTER] 
!      subroutine computedStress      [FILTER]
!      subroutine LU                [SOLVER]
!      subroutine SVD               [SOLVER]
!
!
! BEHAVIOR: 
!           
!
! STATUS : Merge colocated and non-colocated formulation
! 
!----------------------------------------------------------------

module solver
  use global
  use actools
  implicit none
    integer :: i_boxCenter,    j_boxCenter,    k_boxCenter 
    real(8), dimension(:,:,:,:), allocatable :: uu_f, uu_t
contains



!   !****************************************************************
!   !                             CUTOUT                            !
!   !****************************************************************
  
!   !----------------------------------------------------------------
!   ! USE: Resizes an array to a smaller cutout size. 
!   !      
!   !      
!   ! FORM: cutout(array, n_comp)
!   !       n_comp stands for number of components. 
!   !      
!   ! BEHAVIOR: 
!   !           
!   !
!   ! STATUS : 
!   !          
!   !
!   !----------------------------------------------------------------


!   subroutine cutout(array, n_comp)
!     !
!     !    ..ARRAY ARGUMENTS..
!     integer, intent(in) :: n_comp
!     real(8), allocatable, dimension(:,:,:,:),intent(inout) :: array
!     !
!     !    ..WORK ARRAY..
!     real(8), allocatable, dimension(:,:,:,:):: temp

!     allocate (temp(n_comp,testcutSize,testcutSize,testcutSize))
!     temp = array(:,lBound:uBound,lBound:uBound,lBound:uBound)
!     deallocate(array)

!     allocate(array(n_comp,testcutSize,testcutSize,testcutSize))
!     array = temp   
!     deallocate(temp)

!   end subroutine cutout



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
  ! STASH: TO ENSURE THAT THE BOX CENTER IS ALWAYS RETAINED   
  !
  ! **
  !     ! CREATE RANDOM MASK:
  !     allocate(randMask(maskSize))
  !     randMask=0; i=1          ! i starts at 1 because (M-1) points are randomly selected
  !     call init_random_seed()
  !     do
  !        call random_number(random_0_to_1) !returns a random value from 0 to 1
  !        randomInt = nint( (boxSize-1) * random_0_to_1 ) + 1 ! (511 * x) + 1
  ! !       if(any(randMask.eq.randomInt).or.(randomInt.eq.boxCenter)) cycle ! Enforce element uniqueness
  !        i = i + 1   
  !        randMask(i) = randomInt       
  !        if (i.gt.maskSize) exit   
  !     end do
  !
  !----------------------------------------------------------------
  
  subroutine randTrainingSet(randMask)
    !
    !    ..ARRAY ARGUMENT..
    integer, allocatable, intent(out) :: randMask(:)
    !
    !    ..LOCAL SCALARS..
    real(8) :: random_0_to_1 
    integer :: randomInt     !stores integer values  from 1 to 512
    !
    !    ..DEBUG..
    logical                 :: debugRandom = 0
    integer, dimension(:,:,:), allocatable :: randomMatrix
    integer                 :: i, j, k
    integer                 :: count
    integer, save           :: debugCount ! To print this section only once 
    
    ! Create random mask: **

!       allocate(randMask(maskSize))
!       randMask=0; i=1          ! i starts at 1 because (M-1) points are randomly selected
!       call init_random_seed()
!       do
!          call random_number(random_0_to_1) !returns a random value from 0 to 1
!          randomInt = nint( (boxSize-1) * random_0_to_1 ) + 1 ! (511 * x) + 1
!   !       if(any(randMask.eq.randomInt).or.(randomInt.eq.boxCenter)) cycle ! Enforce element uniqueness
!          i = i + 1   
!          randMask(i) = randomInt       
!          if (i.ge.maskSize) exit   
!       end do
  

    allocate(randMask(maskSize))
    randMask=0; i=0          ! i starts at 1 because (M-1) points are randomly selected
    call init_random_seed()
    do
       call random_number(random_0_to_1) !returns a random value from 0 to 1
       randomInt = nint( (boxSize-1) * random_0_to_1 ) + 1 ! (511 * x) + 1
       if(any(randMask.eq.randomInt)) cycle ! Enforce element uniqueness
       if (i.ge.maskSize) exit   
       i = i + 1   
       !
       randMask(i) = randomInt       
    end do
    
    
    ! DEBUG:
    if (debugRandom.and.debugCount == 0) then
       allocate (randomMatrix(box(1),box(2),box(3)))
       count=0;randomMatrix=0

       do k=1,box(3)
       do j=1,box(2)
       do i=1,box(1)
          count = count+1
          if(any(randMask.eq.count)) cycle
          randomMatrix(i,j,k) = count
       enddo
       enddo
       enddo

       print*,'randomMatrix'        ! Print random training points indices:
       do k=1,box(3)
       do i=1,box(1)
          print*,randomMatrix(i,:,k)
       end do
       print*,''
       end do
       debugCount = 1
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
  !          _boxCenter     : bounding box center
  !          _train   : training point in bounding box (stencil-center)
  !          _stencil : for a point in the 3x3x3 STENCIL
  !          _opt     : for computed(OPTIMIZED) values at LES and test scales.
  !          _proj    : to project stresses  for postprocessing
  !
  !          coloc    : 0  means use non-coloc
  !
  ! STATUS : Suppressed all printing and outputs the data files for images.
  !          Needs coarse grained parallization.
  !
  ! 
  !         
  ! STASH:
  !   ##
  !     if (debug(1)) then
  !        allocate(T_ijOpt(1,testSize,testSize,testSize)) !(1,17,17,17)
  !        allocate(tau_ijOpt(1,testSize,testSize,testSize))
  !     else
  !        allocate(T_ijOpt(n_uu,testSize,testSize,testSize)) !(6,17,17,17)
  !        allocate(tau_ijOpt(n_uu,testSize,testSize,testSize))
  !     end if
  !  
  !  ###
  !   if(debug(2)) then
  !      print*,'shape tau_ij cutout:    ',shape(tau_ij)
  !      print*,'shape uu_t cutout:      ',shape(uu_t)
  !      print*,'shape tau_ijOpt cutout: ',shape(tau_ijOpt)
  !   end if
  !
  ! &&& : ALLOCATING STRIDED MATRICES
  ! Test scale computations will involve creation
  ! of a strided matrix. Since an array stores contiguous
  ! set of memory, two things will happen if the array index
  ! is not carefully set - if there is a regular chuck of
  ! heap memory, then it creates a (51,51,51) array and fills
  ! unused space with 0. If not, it results in a
  ! segmentation fault. Thus the index of the strided elements
  ! that are used to form T_ijOpt and tau_ijOpt arrays should be
  ! mapped with _boxCenter indices to give a contiguous array (17,17,17)
  !  
  ! $
  ! PRINT H_IJ
  ! if (check_h_ij) then
  !    print*, 'size(h_ij):', size(h_ij)
  !    do i = 1, size(h_ij,dim=1)
  !       print*, i, h_ij(i,:)
  !    end do
  !    check_h_ij = 0
  ! end if
  !
  !  $$
  !  FIND TRAINING ERROR FOR EACH BOX [THERE WILL BE TOO MANY BOXES; USE IF CONDITION]
  !      call trainingError(T_ijOpt,   T_ij,    error_cross_T_ij,   'plot', cross_csv_T_ij   )
  !      call trainingError(tau_ijOpt, tau_ij,  error_cross_tau_ij, 'plot', cross_csv_tau_ij )
  !
  !
  !   ^^ CHECK V,T with results from MATLAB code [BOULDER]
  !            if (stress.eq.'dev') then
  !              if (withPressure.eqv..false.) then          
  !                 if (dataset.eq.'jhu256') then
  !                    if (V(1500,2000).ne.2.0009431419772586d-2) then
  !                       print*, "Error! Check sorting order in  V matrix!"
  !                       print*, 'V(1500,2000)', V(1500,2000)
  !
  !                       if (T(3,1).ne.8.7759832493259110d-2)then
  !                          print*, "Error! Check sorting order in  T vector!"
  !                          print*, T(3,1)
  !                       end if
  !                    end if
  !                 end if
  !              end if
  !           end if
  !
  !
  !----------------------------------------------------------------
  
  subroutine autonomicClosure(u_f, u_t, tau_ij, T_ij, h_ij, tau_ijOpt, T_ijOpt)
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(1:,extLower:,extLower:,z_extLower:), intent(in) :: u_f
    real(8), dimension(1:,extLower:,extLower:,z_extLower:), intent(in) :: u_t
    real(8), dimension(1:,extLower:,extLower:,z_extLower:), intent(in) :: T_ij
    real(8), dimension(:,:,:,:), intent(in) :: tau_ij

    real(8), dimension(:,:),     intent(out):: h_ij
    real(8), dimension(1:,1:,1:,zLower:), intent(out):: tau_ijOpt
    real(8), dimension(1:,1:,1:,zLower:), intent(out):: T_ijOpt
    !
    !    ..LOCAL ARRAYS..
    real(8), dimension(:,:),     allocatable :: V  
    real(8), dimension(:,:),     allocatable :: T  
    !
    !    .. NON-COLOCATED..
    real(8), dimension(:), allocatable :: u_n
    integer :: non_col_1, non_col_2
    !
    !    ..RANDOM TRAINING POINTS (STENCIL-CENTERS)..
    integer :: rand_count
    integer,dimension(:), allocatable :: randMask
    !
    !    ..INDICES..
    integer :: i_train,     j_train,     k_train  
    integer :: i_stencil,   j_stencil,   k_stencil 
    integer :: row_index, col_index, row, col 
    integer :: u_comp, uu_comp 

    real(8) :: error_cross_T_ij, error_cross_tau_ij
    !
    !    ..DEBUG..
    logical :: printval
    logical :: check_h_ij = 1

! ##
    allocate (V (M, N) )
    allocate (T (M, P) )
! ###

    ! COLOCATED FORMULATION WITH RANDOMLY SELECTED TRANING POINTS:
    if (formulation.eq.'colocated') then

!call cpu_time(tic)

       ! COMPUTE VELOCITY PRODUCTS: TAKE THIS STEP OUT - DON'T NEED TO SAVE EVERYTHING IN AN ARRAY
       ! BUT IF I DO IT THEN I WILL HAVE TO USE AN "IF" STATEMENT INSIDE THE LOOP. SO THIS MIGHT
       ! TAKE UP EXTRA MEMORY BUT CUTS COMPUTATION TIME.
       if (order == 2) then
          if (allocated(uu_t).eqv..false.)allocate(uu_t(n_uu, extLower:extUpper,extLower:extUpper,z_extLower:z_extUpper))
          if (allocated(uu_f).eqv..false.)allocate(uu_f(n_uu, extLower:extUpper,extLower:extUpper,z_extLower:z_extUpper))
          call secondOrderProducts(uu_t, u_t)
          call secondOrderProducts(uu_f, u_f)
       end if


       ! Z-MIDPLANE COMPUTATION: 
       do k_boxCenter = zLower, zUpper
       do j_boxCenter = boxFirst, boxLast, boxCenterSkip
       do i_boxCenter = boxFirst, boxLast, boxCenterSkip

!print*, i_boxCenter, j_boxCenter, k_boxCenter
          if (trainingPoints.eq.'random') then
             call randTrainingSet(randMask)
             rand_count = 0
          end if
          row_index  = 0 

             ! VISIT M TRANING POINTS:
             do k_train = k_boxCenter-boxLower, k_boxCenter+boxUpper, trainingPointSkip
             do j_train = j_boxCenter-boxLower, j_boxCenter+boxUpper, trainingPointSkip
             do i_train = i_boxCenter-boxLower, i_boxCenter+boxUpper, trainingPointSkip
             
                if (trainingPoints.eq.'random') then
                   rand_count = rand_count + 1
                   if (any(randMask.eq.rand_count)) cycle
                end if

!print*, row_index, i_train, j_train, k_train

                col_index = 0 
                row_index = row_index + 1


                ! ZERO ORDER TERMS:
                col_index = col_index + 1
                V(row_index, col_index) = 1.d0

                ! BUILD 3x3x3 STENCIL AT Delta_test SCALE:
                do k_stencil = k_train-Delta_test, k_train+Delta_test, Delta_test ! Vectorize using (PACK/UNPACK)
                do j_stencil = j_train-Delta_test, j_train+Delta_test, Delta_test
                do i_stencil = i_train-Delta_test, i_train+Delta_test, Delta_test
              
                   ! FIRST ORDER TERMS:
                   do u_comp = 1, n_u ! 1 to 3 -> 3x(3x3x3) = 81
                      col_index = col_index+1
                      V(row_index,col_index) = u_t(u_comp,i_stencil,j_stencil,k_stencil)
                   end do

                   ! SECOND ORDER TERMS: 6x(3x3x3) = 162 (GIVES A TOTAL OF 243 TERMS)
                   if (order == 2) then
                      do uu_comp = 1, n_uu 
                         col_index = col_index+1
                         V(row_index,col_index) = uu_t(uu_comp,i_stencil,j_stencil,k_stencil)
                      end do
                   end if

                end do 
                end do 
                end do ! STENCIL
       
                T(row_index,:) = T_ij(:,i_train,j_train,k_train) 

             end do
             call progressBar(j_boxCenter, boxLast)
             end do
             end do ! DONE VISITING ALL RANDOM TRANING POINTS IN A BOUNDING BOX

!call cpu_time(toc)
!print*, 'Time to build V matrix, t1 = ', toc-tic          

             !
             ! BEGIN SUPERVISED TRAINING: FEATURE SELECTION 
!             do iter = 1, n_lambda
!                lambda = lambda_0(1) * 10**(iter-1)
!call cpu_time (tic)
                ! CALL SOLVER:
                if (solutionMethod.eq.'LU') then
                   call LU(V, T, h_ij)                       ! DAMPED LEAST SQUARES 
                elseif(solutionMethod.eq.'SVD') then
                   call SVD(V, T, h_ij, printval)             ! TSVD
                else
                   print*, 'Choose correct solver: LU, SVD'
                   stop
                end if

!call cpu_time(toc)
!print*, 'Time to invert, t2 = ', toc-tic

!$
                call computedStress (u_f, u_t, h_ij, T_ijOpt, tau_ijOpt)
! $$

!             end do ! DONE COMPUTING OPTIMIZED STRESS. MOVE ON TO THE NEXT BOUNDING BOX       
! &&&     

       end do
       end do
       end do ! BOX. DONE COMPUTING OPTIMIZED STRESSES IN ALL BOUNDING BOXES. 

       
       ! COMPUTE OPTIMIZED STRESS USING h_ij AT A GIVEN lambda
        if (plot_Stress)                                        call plotComputedStress(lambda,'All')     
        if (production_Term) then
           call productionTerm(Pij_fOpt, tau_ijOpt, Sij_f)
           call productionTerm(Pij_tOpt, T_ijOpt,   Sij_t)
           if (save_ProductionTerm)                             call plotProductionTerm(lambda)
        end if

       
    
    else

       ! NON-COLOCATED FORMULATION WITH ORDERED TRAINING POINTS:

       allocate (u_n(stencil_size))

       ! WHOLE DOMAIN COMPUTATION: 
!call cpu_time(tic)
!       do k_boxCenter = boxFirst, boxLast, boxCenterSkip
       do k_boxCenter = z_plane, z_plane                    
       do j_boxCenter = boxFirst, boxLast, boxCenterSkip
       do i_boxCenter = boxFirst, boxLast, boxCenterSkip

          if (trainingPoints.eq.'random') then
             call randTrainingSet(randMask)
             rand_count = 0
          end if
          row_index  = 0 

           ! VISIT TRAINING POINT: C-ORDER
          do k_train = k_boxCenter-boxLower, k_boxCenter+boxUpper, trainingPointSkip       
          do j_train = j_boxCenter-boxLower, j_boxCenter+boxUpper, trainingPointSkip
          do i_train = i_boxCenter-boxLower, i_boxCenter+boxUpper, trainingPointSkip       

              
              if (trainingPoints.eq.'random') then
                 rand_count = rand_count + 1
                 if (any(randMask.eq.rand_count)) cycle
              end if

             ! Replace this loop with subroutine build_V()
             col_index = 0 
             row_index = row_index + 1
!print*, row_index, i_train, j_train, k_train


             ! ZERO ORDER TERMS: 80 C 0
             col_index = col_index + 1
             V(row_index, col_index) = 1.d0


             ! ENTER 3x3x3 STENCIL: C-ORDER
             u_n = reshape(u_t (:, i_train-Delta_test : i_train+Delta_test : Delta_test, & 
                                   j_train-Delta_test : j_train+Delta_test : Delta_test, &
                                   k_train-Delta_test : k_train+Delta_test : Delta_test), [stencil_size])

             ! FIRST ORDER TERMS: 81 C 1      
             do non_col_1 = 1,stencil_size 
                col_index = col_index + 1
                V(row_index,col_index) = u_n(non_col_1)
             end do

             ! SECOND ORDER TERMS: 82 C 2
             if(order == 2) then
             do non_col_1 = 1, stencil_size
             do non_col_2 = non_col_1, stencil_size
                col_index = col_index + 1
                V(row_index,col_index) = u_n(non_col_1) * u_n(non_col_2)
             end do
             end do
             end if

             T(row_index,:) = T_ij(:,i_train,j_train,k_train) !Change 1 to (1-6) here. !THIS ONE IS CORRECT; KEEP IT.

          end do
          end do
          end do ! DONE VISITING ALL TRAINING POINTS IN A BOUNDING BOX
          
!call cpu_time(toc)
!print*, 'Time to build V matrix, t1 = ', toc-tic          

          ! CHECK V,T: ^^
!print*, 'V(1500,2000) ',V(1500,2000)
!print*, 'T(3,1) ', T(3,1)

          !
          ! BEGIN SUPERVISED TRAINING: FEATURE SELECTION 
!          do iter = 1, n_lambda
!             lambda = lambda_0(1) * 10**(iter-1)
!             lambda = lambda_0(iter)
!             print('(a8,ES10.2)'), 'lambda ',lambda
!call cpu_time (tic)
             !
             ! CALL SOLVER: LATER SPLIT IT INTO 1) FACTORIZATION STEP (OUT OF LOOP)
             ! AND 2] SOLUTION USING LAMBDA (WITHIN LOOP)
             if (solutionMethod.eq.'LU') then
                call LU(V, T, h_ij)                       ! DAMPED LEAST SQUARES 
             elseif(solutionMethod.eq.'SVD') then
                call SVD(V, T, h_ij, printval)             ! TSVD
             else
                print*, 'Choose correct solver: LU, SVD'
                stop
             end if
!call cpu_time(toc)
!print*, 'Time to invert, t2 = ', toc-tic
             ! CHECK h_ij:
             if (dataset.eq.'jhu256'.and.solutionMethod.eq.'SVD'.and.stress.eq.'dev') then
                if (h_ij(350,1).ne.-4.5121154730201521d-2)then
                   print*, "Error! Check lambda, method or sorting order for h_ij computation:"
                   print*, h_ij(350,1)
                   !stop
                else 
                   print*,'SVD check ... Passed'
                end if
             end if

!              ! COMPUTE OPTIMIZED STRESS USING h_ij AT A GIVEN lambda
              call computedStress (u_f, u_t, h_ij, T_ijOpt, tau_ijOpt)

              ! KEEP IT FOR THE LAMBDA LOOP:
             
!              if (production_Term) then
!                 call productionTerm(Pij_fOpt, tau_ijOpt, Sij_f)
!                 call productionTerm(Pij_tOpt, T_ijOpt,   Sij_t)
!                 if (save_ProductionTerm)                             call plotProductionTerm(lambda)
!              end if

             !  FIND TRAINING ERROR FOR EACH BOX [THERE WILL BE TOO MANY BOXES; USE IF CONDITION]
!             call trainingerror(T_ijOpt,   T_ij,    error_cross_T_ij,   'plot', cross_csv_T_ij   )
!             call trainingError(tau_ijOpt, tau_ij,  error_cross_tau_ij, 'plot', cross_csv_tau_ij )

!          end do ! DONE COMPUTING OPTIMIZED STRESS. MOVE ON TO THE NEXT BOUNDING BOX

       end do
       end do
       end do ! BOX. DONE COMPUTING OPTIMIZED STRESSES IN ALL BOUNDING BOXES. 



       
       ! COMPUTE OPTIMIZED STRESS USING h_ij AT A GIVEN lambda
        if (plot_Stress)                                        call plotComputedStress(lambda,'All')     
        if (production_Term) then
           call productionTerm(Pij_fOpt, tau_ijOpt, Sij_f)
           call productionTerm(Pij_tOpt, T_ijOpt,   Sij_t)
           if (save_ProductionTerm)                             call plotProductionTerm(lambda)
        end if


    end if

  end subroutine autonomicClosure
  
  
  !****************************************************************
  !                        COMPUTED STRESS                        !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Computes tau_ij and T_ij from calculated h_ij
  !      
  !      
  
  !      
  ! FORM: subroutine computedStress (u_f, u_t, h_ij, T_ijOpt, tau_ijOpt)
  !       
  !      
  ! BEHAVIOR: T_ijOpt and tau_ijOpt should be allocated before. 
  !           ZERO ORDER terms are always 1.
  !           Indices are in C-order for direct comparison with MATLAB.
  !
  !      *** INDEX Nomenclature:
  !           _boxCenter: for a point in the TEST-scale field
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
  
  
  subroutine computedStress(u_f, u_t, h_ij, T_ijOpt, tau_ijOpt)
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(1:,extLower:,extLower:,z_extLower:), intent(in) :: u_f
    real(8), dimension(1:,extLower:,extLower:,z_extLower:), intent(in) :: u_t
    real(8), dimension(:,:),     intent(in) :: h_ij
    real(8), dimension(1:,1:,1:,zLower:), intent(out):: T_ijOpt
    real(8), dimension(1:,1:,1:,zLower:), intent(out):: tau_ijOpt 
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
    integer :: i_opt,     j_opt,     k_opt  
    integer :: i_stencil, j_stencil, k_stencil
    integer :: non_col_1, non_col_2

    if (formulation.eq.'colocated') then
! Whole domain:
!       do k_opt = k_boxCenter-optLower, k_boxCenter+optUpper
! Cross-validation points [Preferred]:
       do k_opt = k_boxCenter-3*smallHalf(N_cr), k_boxCenter+3*smallHalf(N_cr), 3
       do j_opt = j_boxCenter-optLower, j_boxCenter+optUpper
       do i_opt = i_boxCenter-optLower, i_boxCenter+optUpper
         
!call cpu_time (tic)

          ! T_ij^F: Test scale
          col_index = 0 
          
          ! ZERO ORDER TERMS: [1]
          col_index = col_index + 1
          T_ijOpt   (:,i_opt, j_opt, k_opt) = h_ij(col_index,:) 
          
          do k_stencil = k_opt-Delta_test, k_opt+Delta_test, Delta_test
          do j_stencil = j_opt-Delta_test, j_opt+Delta_test, Delta_test
          do i_stencil = i_opt-Delta_test, i_opt+Delta_test, Delta_test

             ! FIRST ORDER TERMS:             
             do u_comp = 1, n_u
                col_index = col_index + 1

                T_ijOpt (:,i_opt,j_opt,k_opt) = T_ijOpt (:,i_opt,j_opt,k_opt)             &
                                        +                                          &
                     (u_t(u_comp,i_stencil,j_stencil,k_stencil)) * h_ij(col_index,:)
             end do

             
             ! SECOND ORDER TERMS: 
             if (order == 2) then
                do uu_comp = 1, n_uu
                   col_index = col_index + 1

                   T_ijOpt (:,i_opt,j_opt,k_opt) = T_ijOpt(:,i_opt,j_opt,k_opt)               &
                        +                                            &
                        (uu_t(uu_comp,i_stencil,j_stencil,k_stencil)) * h_ij(col_index,:)
                end do
             end if

          end do
          end do
          end do

          ! tau_ij^F: LES scale
          col_index = 0
          
          ! ZERO ORDER TERMS: [1]
          col_index = col_index + 1
          tau_ijOpt (:,i_opt, j_opt, k_opt) = h_ij(col_index,:) 
          
          do k_stencil = k_opt-Delta_LES, k_opt+Delta_LES, Delta_LES
          do j_stencil = j_opt-Delta_LES, j_opt+Delta_LES, Delta_LES
          do i_stencil = i_opt-Delta_LES, i_opt+Delta_LES, Delta_LES

             ! FIRST ORDER TERMS:             
             do u_comp = 1, n_u
                col_index = col_index + 1
                tau_ijOpt(:,i_opt,j_opt,k_opt) = tau_ijOpt(:,i_opt,j_opt,k_opt)         &
                     +                                          &
                     (u_f(u_comp,i_stencil,j_stencil,k_stencil)) * h_ij(col_index,:)
             end do
             
             ! SECOND ORDER TERMS: 
             if (order == 2) then
                do uu_comp = 1, n_uu
                   col_index = col_index + 1
                   tau_ijOpt(:,i_opt,j_opt,k_opt) = tau_ijOpt(:,i_opt,j_opt,k_opt) &
                        +                                &
                        (uu_f(uu_comp,i_stencil,j_stencil,k_stencil)) * h_ij(col_index,:)
                end do
             end if
!call cpu_time(toc)
!print*, 'Time to compute one cell, t3 = ', toc-tic
!stop

          end do
          end do
          end do

       end do
       end do
       end do

    else
       
       ! NON-COLOCATED

       ! ENTER STENCIL-CENTER POINTS: C-ORDER
! Whole domain:
          ! do k_opt = k_boxCenter-optLower, k_boxCenter+optUpper 
! Cross-validation points [Preferred]:
       do k_opt = k_boxCenter-3*smallHalf(N_cr), k_boxCenter+3*smallHalf(N_cr), 3
       do j_opt = j_boxCenter-optLower, j_boxCenter+optUpper
       do i_opt = i_boxCenter-optLower, i_boxCenter+optUpper

!call cpu_time (tic)
          col_index = 0 
         
          ! ENTER 3x3x3 STENCIL: C-ORDER                                      
          u_t_stencil = reshape(u_t (:, i_opt-Delta_test : i_opt+Delta_test : Delta_test, &
                                        j_opt-Delta_test : j_opt+Delta_test : Delta_test, &
                                        k_opt-Delta_test : k_opt+Delta_test : Delta_test),  [stencil_size]) 

          u_f_stencil = reshape(u_f (:, i_opt-Delta_LES : i_opt+Delta_LES : Delta_LES, &
                                        j_opt-Delta_LES : j_opt+Delta_LES : Delta_LES, &
                                        k_opt-Delta_LES : k_opt+Delta_LES : Delta_LES),  [stencil_size]) 

          ! ZERO ORDER TERMS: WILL BE 1
          col_index = col_index + 1

          T_ijOpt   (:,i_opt, j_opt, k_opt) = h_ij(col_index,:) 
          tau_ijOpt (:,i_opt, j_opt, k_opt) = h_ij(col_index,:) 

          ! FIRST ORDER TERMS:             
          do non_col_1 = 1, stencil_size 
             col_index = col_index + 1

             T_ijOpt(:,i_opt,j_opt,k_opt) = T_ijOpt (:,i_opt,j_opt,k_opt)             &
                                         +                                          &
                                         (u_t_stencil(non_col_1) * h_ij(col_index,:))

             tau_ijOpt(:,i_opt,j_opt,k_opt) = tau_ijOpt(:,i_opt,j_opt,k_opt)         &
                                          +                                          &
                                          (u_f_stencil(non_col_1) * h_ij(col_index,:))
          end do

          ! SECOND ORDER TERMS:
          if (order == 2) then
          do non_col_1 = 1, stencil_size
          do non_col_2 = non_col_1, stencil_size
             col_index = col_index+1

             T_ijOpt(:,i_opt,j_opt,k_opt) = T_ijOpt(:,i_opt,j_opt,k_opt) &
                                         +                             &
                                         (u_t_stencil(non_col_1) * u_t_stencil(non_col_2) * h_ij(col_index,:))

             tau_ijOpt(:,i_opt,j_opt,k_opt) = tau_ijOpt(:,i_opt,j_opt,k_opt) &
                                            +                                &
                                            (u_f_stencil(non_col_1) * u_f_stencil(non_col_2) * h_ij(col_index,:))
          end do
          end do
          end if
!call cpu_time(toc)
!print*, 'Time to compute one cell, t3 = ', toc-tic
!stop
       end do
       end do
       end do ! DONE COMPUTING STRESSES USING NON-COLOCATED FORM

    end if

  end subroutine computedStress


  
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
  ! FORM: subroutine LU (V, T_ij, h_ij)
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
  
  subroutine LU (V,  T_ij, h_ij)
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:), intent(in)  :: V 
    real(8), dimension(:,:), intent(in)  :: T_ij
    real(8), dimension(:,:), intent(out) :: h_ij
    !
    !    ..SCALAR ARGUMENTS..
!    real(8),optional :: lambda
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
    integer :: LWMAX
    integer :: NRHS 
    integer :: INFO

    LWMAX = M * N
    NRHS = P
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
  ! FORM: subroutine SVD (A, T, h_ij, [printval])
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
  
  subroutine SVD (A, T, h_ij, printval)
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:), intent(in)  :: A
    real(8), dimension(:,:), intent(in)  :: T
    real(8), dimension(:,:), intent(out) :: h_ij
    !
    !    ..SCALAR ARGUMENTS..
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
    integer :: LWMAX 
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

 
    LWMAX =3 * M * N
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
       call DGESDD('S', M, N, A, M, S, U, M, VT, N, WORK, LWORK, IWORK, info) ! TAKES ABOUT 45 SECS
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
!    lambda = 1.d-1
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

