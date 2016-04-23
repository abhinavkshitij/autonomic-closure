!****************************************************************
!                              SOLVER
!****************************************************************

!----------------------------------------------------------------
! USE: Solver
!      
!      
!      
! FORM: module solver
!          contains
!       
!
! BEHAVIOR: 
!           
!
! STATUS : Refactoring this unit.
! 
!----------------------------------------------------------------

!! Linsolve is a collection of subroutines that perform linear algebra operations
!! and related operations to solve the inverse problem.

module solver

 ! use fileio
  use global

  ! Stencil parameters:
  integer,parameter :: stride = 1 ! Is the ratio between LES(taken as 1) and test scale
  integer,parameter :: skip = 10
  integer,parameter :: X = 1     ! Number of realizations
  integer,parameter :: n_DAMP = 1  ! Number of lambda's
  
 
  ! Bounding Box parameters:  
  integer,parameter :: box       = 252
  integer,parameter :: boxSize   = box**3
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


subroutine cutout(array,n_comp)
  implicit none

  integer, intent(in) :: n_comp
  real(8), allocatable, dimension(:,:,:,:),intent(inout) :: array
  real(8), allocatable, dimension(:,:,:,:):: temp

  allocate (temp(n_comp,testcutSize,testcutSize,testcutSize))
  temp = array(:,lBound:uBound,lBound:uBound,lBound:uBound)
  deallocate(array)

  allocate(array(n_comp,testcutSize,testcutSize,testcutSize))
  array = temp   
  deallocate(temp)

  return
end subroutine cutout


! RANDOM NUMBER GENERATION
subroutine init_random_seed()

  integer, dimension(:), allocatable :: seed
  integer :: i, n_seed, clock

  call RANDOM_seed(size = n_seed) !Generate random seeds
  allocate(seed(n_seed))
  call SYSTEM_clock(COUNT=clock)
  seed = clock + 37 * (/ (i-1, i=1,n_seed) /)
  call RANDOM_seed(PUT = seed)
  deallocate(seed)
  return
end subroutine init_random_seed

subroutine randTrainingSet(randMask)
! Status: Passed testing 
! Notes: Creates a mask to skip cells randomly.
!      : Always retains the center point. So M-1 points are randomly selected
!      : whereas the center point is ALWAYS retained.
implicit none

real(8)                 :: random_0_to_1 !stores real value from 0 to 1
integer                 :: randomInt     !stores integer value  from 1 to 512
integer,parameter       :: maskSize = boxSize - M ! 512-Training points(243) = 269
integer,intent(out)     :: randMask(maskSize) 

! DEBUG VARIABLES:
logical,parameter       :: debugRandom = .FALSE.
integer                 :: i, j, k
integer                 :: count
integer                 :: randomMatrix(box,box,box)

randMask=0; i=1          ! i starts at 1 because (243 - 1) points are randomly selected

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

return       
end subroutine randTrainingSet

!  COMPUTE SGS STRESS:
! VERSION : Performs entire domain computation with just a single damping value and
!           realization. Suppress all printing and outputs the data files for images.
! Enter only one DAMPING value.

subroutine synStress(u_f, u_t, tau_ij, T_ij, h_ij)
implicit none


! ARGUMENTS:
real(8),dimension(:,:,:,:),intent(in) :: u_f,    u_t
real(8),dimension(:,:,:,:),intent(in) :: tau_ij, T_ij
real(8),dimension(:,:),intent(out) :: h_ij

! PARAMETERS:
integer,parameter :: n_u=3, n_uu=6
real(8) :: lambda = 1.d-1

! LOCAL VARIABLES:
real(8),allocatable,dimension(:,:,:,:) :: uu_t, uu_f
real(8),allocatable,dimension(:,:,:,:) :: TijOpt, tau_ijOpt !Computed stresses

! FOR LEAST SQUARES PROBLEM:
!real(8),allocatable,dimension(:,:) :: V  ! Non-linear combination of u_i, u_iu_j
!real(8),allocatable,dimension(:,:)   :: T  ! Input training data in LAPACK routines.
real(8), dimension(M,N) :: V = 0.
real(8), dimension(M,P) :: T = 0.
integer                :: row_index, col_index, row, col ! Indices to build V,h_ij arrays.
integer                :: u_comp, uu_comp ! Indices to select u_i,u_iu_j components

! FOR NON-COLOCATED FORMULATION:
logical :: coloc = 0            !0 means use non-coloc
integer,parameter :: stencil_size = 3*(3*3*3) !3 * 27 = 81
real(8) :: u_n(stencil_size)
integer :: non_col_1, non_col_2                   ! To run non coloc combinantions

! FOR RANDOM TRAINING POINTS (STENCIL-CENTERS):
integer :: rand_count
integer :: randMask(boxSize - M) !512-243=269

! LIST OF INDICES FOR LOOPS:
integer :: i_test,    j_test,    k_test ! for a point in the TEST-scale field
integer :: i_box,     j_box,     k_box  ! for a point in bounding BOX (stencil-center)
integer :: i_stencil, j_stencil, k_stencil ! for a point in the 3x3x3 STENCIL
integer :: i_opt,     j_opt,     k_opt     ! for computed(OPTIMIZED)values at LES and Test scales.
integer :: i_proj,    j_proj,    k_proj    ! to PROJECT final computed data for postprocessing


! EXTERNAL FILES:
integer,parameter      ::  n_files = 1     !CHANGE HERE
character(50)::  PATH=trim(RUN_DIR)//'run17/' !CHANGE HERE FOR DIFFERENT EXPERIMENTS
character(10)::  l_val = "l_0_10/1/" ! l_var stands for variable lambda.

logical      ::  writeStress = 0    
character(4) ::  z_plane,idx
integer      ::  vizPlaneC, vizPlaneF

! DEBUG SWITCHES:
logical :: debug(4) = [1,1,0,0]
! 1 - Takes just one (tau_ij_11) for testing purpose. 
! 2 - Prints cutout shapes to check if cutout ops have gone well.
! 3 - Performs all computations on just the first point.
! 4 - To check if values from the stencil and bounding box are working well.
integer :: lim

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

!allocate (V(M,N),T(M,P))
print*, shape(V)
TijOpt=0.
tau_ijOpt=0.

if(debug(2)) then
   print*,'shape tau_ij cutout:    ',shape(tau_ij)
   print*,'shape uu_t cutout:      ',shape(uu_t)
   print*,'shape tau_ijOpt cutout: ',shape(tau_ijOpt)
end if


! COLOCATED FORMULATION:
if (coloc) then

! ! Compute velocity products:(Lower triangle order for ij)
! k = 0
! do j=1,n_u
! do i=1,n_u
!    if (i.ge.j) then          
!       k=k+1
!       uu_t(k,:,:,:) = u_t(i,:,:,:) * u_t(j,:,:,:)
!       uu_f(k,:,:,:) = u_f(i,:,:,:) * u_f(j,:,:,:)
!    end if
! end do
! end do


! ! WHOLE DOMAIN COMPUTATION: 
! do k_test = testLower, testUpper, stride 
! do j_test = testLower, testUpper, stride
! do i_test = testLower, testUpper, stride ! i_test = 11,43,2

!       row_index  = 0 

!       ! ENTER STENCIL-CENTER POINTS:
!       do k_box = k_test-boxLower, k_test+boxUpper, skip
!       do j_box = j_test-boxLower, j_test+boxUpper, skip
!       do i_box = i_test-boxLower, i_test+boxUpper, skip  ! i_box = 3,49,2

!          col_index = 0 
!          row_index = row_index + 1
         
!          ! ENTER 3x3x3 STENCIL:
!          do k_stencil = k_box-stride, k_box+stride, stride ! Vectorize using (PACK/UNPACK)
!          do j_stencil = j_box-stride, j_box+stride, stride
!          do i_stencil = i_box-stride, i_box+stride, stride
         
!             ! ZERO ORDER TERMS:
!             V(row_index, col_index) = 1.d0

!             ! FIRST ORDER TERMS:
!             do u_comp = 1,n_u ! 1 to 3 -> 3x(3x3x3) = 81
!                col_index = col_index+1
!                V(row_index,col_index) = u_t(u_comp,i_stencil,j_stencil,k_stencil)
!             end do

!             ! SECOND ORDER TERMS: 6x(3x3x3) = 162 (GIVES A TOTAL OF 243 TERMS)
!             do uu_comp = 1,n_uu ! 1 to 6
!                col_index = col_index+1
!                V(row_index,col_index) = uu_t(uu_comp,i_stencil,j_stencil,k_stencil)
!             end do

!          end do 
!          end do 
!          end do ! STENCIL

!          T(row_index) = T_ij(n,i_box,j_box,k_box) !Change 1 to (1-6) here.

!       end do
!       end do
!       end do ! BOUNDING BOX


!       ! 1) Until here, the mechanism to compute the coefficient matrix
!       ! is tested and verified. There are no issues with the basic mechanism.
!       ! Only variable strides, for more than 2 is to be adjusted in the code later.

!       ! 2) For testing, only the '_11' component is evaluated. Later extend to all 6 ij's
!       ! V matrix generated can be used to determine Volterra coefficients for all six ij components.


!       !call LU(V,T,h_ij,lambda)     ! Least squares by LU decomposition


!       ! Test scale computations will involve creation
!       ! of a strided matrix. Since an array stores contiguous
!       ! set of memory, two things will happen if the array index
!       ! is not carefully set - if there is a regular chuck of
!       ! heap memory, then it creates a (51,51,51) array and fills
!       ! unused space with 0. If not, it results in a
!       ! segmentation fault. Thus the index of the strided elements
!       ! that are used to form TijOpt and tau_ijOpt arrays should be
!       ! mapped with _test indices to give a contiguous array (17,17,17)

!       i_proj = (i_test-testLower)/stride + 1
!       j_proj = (j_test-testLower)/stride + 1
!       k_proj = (k_test-testLower)/stride + 1

!       ! COARSE STENCIL: Project T_ij back to itself to compare with the original T_ij field
!       ! VECTORIZE THIS PART FOR SPEED
!       col = 0
!       do k_opt = k_test-stride, k_test+stride, stride
!       do j_opt = j_test-stride, j_test+stride, stride
!       do i_opt = i_test-stride, i_test+stride, stride

!          do u_comp = 1,n_u ! 1 to 3
!             col = col+1
!             TijOpt(1,i_proj,j_proj,k_proj) = TijOpt(1,i_proj,j_proj,k_proj)&
!                  + u_t(u_comp,i_opt,j_opt,k_opt) * h_ij(col)
!          end do
!          do uu_comp = 1,n_uu ! 1 to 6
!             col = col+1
!             TijOpt(1,i_proj,j_proj,k_proj) = TijOpt(1,i_proj,j_proj,k_proj)&
!                  + uu_t(uu_comp,i_opt,j_opt,k_opt) * h_ij(col)
!          end do

!       end do
!       end do
!       end do ! coarse

!       ! FINE STENCIL : Calculate the SGS stress with h_ij;compare with the original tau_ij     
!       col = 0 
!       do k_opt = k_test-1, k_test+1
!       do j_opt = j_test-1, j_test+1
!       do i_opt = i_test-1, i_test+1

!          do u_comp = 1,n_u ! 1 to 3
!             col = col+1
!             tau_ijOpt(1,i_proj,j_proj,k_proj) = tau_ijOpt(1,i_proj,j_proj,k_proj)&
!                  + u_f(u_comp,i_opt,j_opt,k_opt) * h_ij(col)
!          end do
!          do uu_comp = 1,n_uu ! 1 to 6
!             col = col+1
!             tau_ijOpt(1,i_proj,j_proj,k_proj) = tau_ijOpt(1,i_proj,j_proj,k_proj)&
!                  + uu_f(uu_comp,i_opt,j_opt,k_opt) * h_ij(col)
!          end do

!       end do
!       end do
!       end do  ! fine                

! end do
! end do
! end do ! test

else


    ! WHOLE DOMAIN COMPUTATION: 
    do i_test = 129, 129, 2 
    do j_test = 129, 129, 2
    do k_test = 129, 129, 2 ! i_test = 11,43,2

       row_index  = 0 

       ! ENTER STENCIL-CENTER POINTS: C-ORDER
       do i_box = i_test-126, i_test+125, 10
       do j_box = j_test-126, j_test+125, 10
       do k_box = k_test-126, k_test+125, 10 ! i_box = 3,49,2
          col_index = 0 
          row_index = row_index + 1
         
         ! ENTER 3x3x3 STENCIL: C-ORDER
                                                                            
          u_n = reshape(u_t (:,i_box-2:i_box+2:2, & 
                               j_box-2:j_box+2:2, &
                               k_box-2:k_box+2:2),&
                               [stencil_size]   )

          ! ZERO ORDER TERMS:
          col_index = col_index+1
          V(row_index, col_index) = 1.d0

          
          ! FIRST ORDER TERMS:          
          do non_col_1 = 1,stencil_size 
             col_index = col_index+1
             V(row_index,col_index) = u_n(non_col_1)
          end do

          ! SECOND ORDER TERMS: 6x(3x3x3) = 162 (GIVES A TOTAL OF 243 TERMS)
          do non_col_1 = 1, stencil_size
             do non_col_2 = non_col_1, stencil_size
                col_index = col_index+1
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
         print*,T(3,1)
         stop
      else 
         print*,'T vector check ... Passed'
      end if
      !print*, 'T(1,2)', T(2)
      !print*, 'T(1,3)', T(3)
      print*, 'T(1,4)', T(1,4)
      !print*, 'T(1,5)', T(5)
      
      
!      if (computeVolterra) then

         !         call SVD(V,T,h_ij,lambda,printval) ! TSVD
         call LU(V,T,h_ij,lambda)           !Damped least squares
         ! PRINT h_ij:
         print*, 'h_ij(1,1):',h_ij(1,1)
         print*, 'h_ij(20,1):', h_ij(20,1)
         print*, 'h_ij(350,1):', h_ij(350,1)


         ! CHECK h_ij:
         if (h_ij(350,1).ne.-4.5121154730201521d-2)then
            print*, "Error! Check lambda, method or sorting order for h_ij computation:"
            print*,h_ij(350,1)
            stop
         else 
            print*,'SVD check ... Passed'
         end if
!      end if

      ! SAVE/READ h_ij:
      open(1,file='./run/apples/h_ij.dat')
      read(1,*) h_ij
      print*,'h_ij(350,1) from file: ',h_ij(350,1)
      close(1)
   end do
   end do
   end do ! test


end if


! !  POST-PROCESSING:
! !  Write 3-plane slices of all 4 stresses for MATLAB plots.
! !  Convert to subroutine later; combine with matrixview()
! if(writeStress) then  
!    do k = 1,3
!       vizPlaneC = nint(0.25 * k * testSize) ! Mark three slices at 1/4, 1/2, 3/4 z-planes.
!       vizPlaneF = (vizPlaneC - 1)*2 + testLower
!       write(z_plane, '(i1)'), k

!       open(1,file= trim(PATH)//trim(l_val)//'T_ij_z'//trim(z_plane)//'.dat')
!       open(2,file= trim(PATH)//trim(l_val)//'TijOpt_z'//trim(z_plane)//'.dat')
!       open(3,file= trim(PATH)//trim(l_val)//'tau_ij_z'//trim(z_plane)//'.dat')
!       open(4,file= trim(PATH)//trim(l_val)//'tau_ijOpt_z'//trim(z_plane)//'.dat')

!       write(1,*), (T_ij      (n_files, i, testLower:testUpper:stride,vizPlaneF), i=testLower,testUpper,stride) !Change 1 here
!       write(2,*), (TijOpt    (1, i, :, vizPlaneC), i=1,testSize)
!       write(3,*), (tau_ij    (n_files, i, testLower:testUpper:stride,vizPlaneF), i=testLower,testUpper,stride) !Change 1 here
!       write(4,*), (tau_ijOpt (1, i, :, vizPlaneC), i=1,testSize)

!       close(1)
!       close(2)
!       close(3)
!       close(4)
!    end do
! end if


! deallocate(uu_t,uu_f,TijOpt,tau_ijOpt)
! return

end subroutine synStress



 subroutine optimizedTij(u_f,u_t,h_ij,TijOpt,tau_ijOpt)

   ! ARGUMENTS:
   real(8), dimension(:,:,:,:), intent(in) :: u_f, u_t
   real(8), dimension(:,:),intent(in) :: h_ij

   ! LOCAL VARIABLES:
   real(8),allocatable,dimension(:,:,:,:),intent(out) :: TijOpt, tau_ijOpt !Computed stresses

   integer                :: col_index
   integer                :: u_comp, uu_comp ! Indices to select u_i,u_ij components


   ! FOR NON-COLOCATED FORMULATION:
   integer,parameter :: stencil_size = 3*(3*3*3) !3 * 27 = 81
   real(8),dimension(stencil_size) :: u_f_stencil, u_t_stencil
   integer :: non_col_1, non_col_2  ! To run non coloc combinantions   

   ! LIST OF INDICES FOR LOOPS:
   integer :: i_test,    j_test,    k_test ! for a point in the TEST-scale field
   integer :: i_opt,     j_opt,     k_opt  ! for a point in bounding BOX (stencil-center)
   integer :: i_stencil, j_stencil, k_stencil


   allocate (TijOpt(1,i_GRID,j_GRID,k_GRID))
   allocate (tau_ijOpt(1,i_GRID,j_GRID,k_GRID))

   ! WHOLE DOMAIN COMPUTATION: 
   do i_test = 129, 129
   do j_test = 129, 129
   do k_test = 129, 129

      ! ENTER STENCIL-CENTER POINTS: C-ORDER
      do i_opt = i_test-127, i_test+126
      do j_opt = j_test-127, j_test+126
      do k_opt = k_test-127, k_test+126

         col_index = 0 
         
         ! ENTER 3x3x3 STENCIL: C-ORDER                                       
          u_t_stencil = reshape(u_t(:,i_opt-2:i_opt+2:2, j_opt-2:j_opt+2:2,k_opt-2:k_opt+2:2),(/stencil_size/)) 
          u_f_stencil = reshape(u_f(:,i_opt-1:i_opt+1:1, j_opt-1:j_opt+1:1,k_opt-1:k_opt+1:1),(/stencil_size/)) 

          ! ZERO ORDER TERMS:
          col_index = col_index+1
          TijOpt(:,i_opt,j_opt,k_opt) = h_ij(col_index,:) ! = 1
          tau_ijOpt(:,i_opt,j_opt,k_opt) = h_ij(col_index,:) ! = 1

          ! FIRST ORDER TERMS:             
          do non_col_1 = 1,stencil_size 
             col_index = col_index+1
             TijOpt(:,i_opt,j_opt,k_opt) = TijOpt(:,i_opt,j_opt,k_opt) + (u_t_stencil(non_col_1)*h_ij(col_index,:))
             tau_ijOpt(:,i_opt,j_opt,k_opt) = tau_ijOpt(:,i_opt,j_opt,k_opt) + (u_f_stencil(non_col_1)*h_ij(col_index,:))
          end do

          ! SECOND ORDER TERMS: 6x(3x3x3) = 162 (GIVES A TOTAL OF 243 TERMS)
          do non_col_1 = 1, stencil_size
          do non_col_2 = non_col_1, stencil_size
             col_index = col_index+1
             TijOpt(:,i_opt,j_opt,k_opt) = TijOpt(:,i_opt,j_opt,k_opt) &
                  + (u_t_stencil(non_col_1) * u_t_stencil(non_col_2)*h_ij(col_index,:))
             tau_ijOpt(:,i_opt,j_opt,k_opt) = tau_ijOpt(:,i_opt,j_opt,k_opt) &
                  + (u_f_stencil(non_col_1) * u_f_stencil(non_col_2)*h_ij(col_index,:))
          end do
          end do
       end do
       end do
       end do ! BOUNDING BOX
    end do
    end do
    end do ! test

    ! WRITE COMPUTED STRESS:
    open(1,file='./run/apples/T_ijOpt.dat')
    open(2,file='./run/apples/tau_ijOpt.dat')
    write(1,*) TijOpt(1,:,:,129)
    write(2,*) tau_ijOpt(1,:,:,129)
    close(1)
    close(2)
    print*,'TijOpt(1,15,24,10):',TijOpt(1,15,24,10)
    print*,'tau_ijOpt(1,15,24,10):',tau_ijOpt(1,15,24,10)
    return
  end subroutine optimizedTij


  
  !****************************************************************
  !                                LU                             !
  !****************************************************************
  subroutine LU(V, T_ij, h_ij, lambda)
    implicit none

    ! dgesv destroys the original matrix. So V is copied into 'a' matrix. This is the LU product matrix.
    ! dgesv returns x vector. First it copies the values from 'T_ij' vector and then computes 'h_ij'.
    ! 'h_ij' is the h_ij vector.

    ! ARGUMENTS:
    real(8), dimension(:,:), intent(in) :: V 
    real(8), dimension(:,:), intent(in) :: T_ij
    real(8), dimension(:,:), intent(out) :: h_ij
    real(8), intent(in) :: lambda

    ! DGESV ARGUMENTS:
    integer, parameter        :: LDA = N, LDB = N, nrhs = P
    real(8), dimension(:,:),allocatable :: A,VT ! EYE - IDENTITY MATRIX
    real(8), dimension(:,:),allocatable :: b
    integer, dimension(:), allocatable  :: ipiv
    integer                   :: info
    real(8) :: tic, toc


    !DGEMM ARGUMENTS:
    real(8):: alpha=1.d0, beta=0.d0

    allocate (A(N,N), b(N,P), ipiv(N))

    ! VT = transpose(V)
    call cpu_time(tic)
    call dgemm('T','N',N,N,M, alpha, V,M, V,M, beta, A,N)
    call cpu_time(toc)
    print*,'Elapsed time', toc-tic


    !A
    forall(i=1:N) A(i,i) = A(i,i) + lambda 
    print*,'A:'

    call printplane(A,frameLim=4)
    !b
    call dgemm('T','N',N,P,M, alpha, V,M, T_ij,M, beta, b,N)
    print*,'b:'

    call printplane(b,frameLim=4)
    call cpu_time(tic)
    !$omp parallel 
    call DGESV(N, nrhs, A, LDA, ipiv, b, LDB, info) !A is a symmetric matrix

    !$omp end parallel 
    call cpu_time(toc)
    print*,'Elapsed time:', toc-tic

    h_ij = b

    deallocate (A,b,ipiv)
    return
  end subroutine LU

  !****************************************************************
  !                             SVD                               !
  !****************************************************************

  subroutine SVD(A,T,h_ij,lambda,printval)

    ! ARGUMENTS:
    real(8),dimension(:,:),intent(in)  :: A
    real(8),dimension(:,:),intent(in)  :: T
    real(8),intent(in)                 ::lambda
    real(8),dimension(:,:),intent(out) :: h_ij


    ! DGESVD ARGUMENTS:                                                                                                
    integer, parameter :: LDA = M
    integer, parameter :: LDU = M
    integer, parameter :: LDVT = N
    integer :: LWMAX
    integer :: info, LWORK

    real(8),dimension(:,:),allocatable :: U
    real(8),dimension(:,:),allocatable :: VT 
    real(8),dimension(:,:),allocatable :: D, Vinv
    real(8),dimension(:), allocatable :: S, work


    logical :: printval


    LWMAX = M*N
    allocate (U(LDU,M), VT(N,N), S(N), D(N,M), Vinv(N,M), work(LWMAX))


    if (printval) then
       print*,'A'                      ! A MATRIX
       call printplane(A,frameLim=4)
    end if


    ! CALL SVD ROUTINE:
    LWORK = -1
    call dgesvd('All','All', M,N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, info)
    LWORK = min( LWMAX, int(WORK(1)) ) 
    call dgesvd('All','All',M,N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, info) !BOTTLENECK = 850 SECS

    ! Convergence check:
    if (info.gt.0) then
       print*, 'Failed to converge'
       stop
    end if


    if (printval) then
       print*,'A'                      ! A MATRIX
       call printplane(A,size(A,dim=1))
       print*,'U'                      ! U MATRIX
       call printplane(U,size(U,dim=1))
       print*,'S'                      ! S VECTOR
       print*, S(1:3)
       print*,'VT'                     ! VT MATRIX
       call printplane(VT,size(VT,dim=1))
    end if


    ! COMPUTE PSEUDOINVERSE:
    D = 0.d0
    forall(i=1:N) D(i,i) = S(i) / (S(i)**2 + lambda**2)
    Vinv = matmul(matmul(transpose(VT),D),transpose(U)) !BOTTLENECK = 1300 SECS

    if(printval) then
       print*,'D'                     ! D MATRIX - DIAG MATRIX
       call printplane(D,size(D,dim=1))

       print*,'Vinv'                     ! Vinv MATRIX - DIAG MATRIX
       call printplane(Vinv,size(Vinv,dim=1))
    end if

    h_ij = matmul(Vinv,T)
    return
  end subroutine SVD



end module solver

