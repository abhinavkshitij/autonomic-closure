module linsolve

  !! Linsolve is a collection of subroutines that perform linear algebra operations
  !! and related operations to solve the inverse problem.
 
  use fileio

  integer,parameter :: stride = 2 ! Is the ratio between LES(taken as 1) and test scale
  
  !!$ Stencil parameters:
  integer,parameter :: M = 300   !!$ Number of training points 3x3x3x9
  integer,parameter :: N = 27 * 9

  !!$ Bounding Box parameters:
  real,   parameter :: eps = 1e-3 !Do not remove eps from here. Ensures correct integer values since
  ! both ceiling and floor functions can give wrong integer values due to conversion of
  ! an integer value into its machine representation.
  integer,parameter :: box = 8
  integer,parameter :: boxSize = box**3
  integer,parameter :: bigHalf   = ceiling(0.5*real(box) + eps) !!$ for 8->5
  integer,parameter :: smallHalf = floor(0.5*real(box) + eps)   !!$ for 8->4
  integer,parameter :: boxCenter = smallHalf * box*(box + 1) + bigHalf
  integer,parameter :: boxLower  = stride * (bigHalf - 1)
  integer,parameter :: boxUpper  = stride * (box - bigHalf)

  !!$ Test field parameters: 
  integer,parameter :: testSize = 17
  integer,parameter :: testcutSize = stride * (testSize + box) + 1
  integer,parameter :: testLower = stride * bigHalf + 1
  integer,parameter :: testUpper = stride * (bigHalf - 1 + testSize) + 1

  !!$ Cutout parameters:
  integer,parameter :: lBound = 0.5*(GRID - testcutSize)
  integer,parameter :: uBound = 0.5*(GRID + testcutSize) - 1

contains

  subroutine printParams()
    implicit none
    print*, 'Stencil parameters:'
    print*, 'Training points :    ',M
    print*, 'Features :           ',N
    print*, ''
    print*, 'Box parameters:      '
    print*, 'Bounding box:        ',box
    print*, 'bigHalf:             ',bigHalf
    print*, 'smallHalf:           ',smallHalf
    print*, 'boxCenter:           ',boxCenter
    print*, 'boxLower:            ',boxLower
    print*, 'boxUpper:            ',boxUpper
    print*, ''
    print*, 'Test field parameters:'
    print*, 'testcutSize:         ',testcutSize
    print*, 'testLower:           ',testLower
    print*, 'testUpper:           ',testUpper
    print*, ''
    print*, 'Cutout parameters:'
    print*, 'lower bound:         ',lBound
    print*, 'upper bound:         ',uBound
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

  
  !!$ RANDOM NUMBER GENERATION
  subroutine init_random_seed()
          
    integer, dimension(:), allocatable :: seed
    integer :: i, n_seed, clock

    call RANDOM_seed(size = n_seed) !!$Generate random seeds
    allocate(seed(n_seed))
    call SYSTEM_clock(COUNT=clock)
    seed = clock + 37 * (/ (i-1, i=1,n_seed) /)
    call RANDOM_seed(PUT = seed)
    deallocate(seed)
    return
  end subroutine init_random_seed

  subroutine randTrainingSet(randMask)
    !!$ Status: Passed testing 
    !!$ Notes: Creates a mask to skip cells randomly.
    !!$      : Always retains the center point. So M-1 points are randomly selected
    !!$      : whereas the center point is ALWAYS retained.
    implicit none

    real(8)                 :: random_0_to_1 !!$stores real value from 0 to 1
    integer                 :: randomInt     !!$stores integer value  from 1 to 512
    integer,parameter       :: maskSize = boxSize - M !!$ 512-Training points(243) = 269
    integer,intent(out)     :: randMask(maskSize) 

    !!$ Debug variables:
    logical,parameter       :: debugRandom = .FALSE.
    integer                 :: i, j, k, count
    integer                 :: randomMatrix(box,box,box)

    randMask=0; i=1          !!$ i starts at 1 because (243 - 1) points are randomly selected

    call init_random_seed()
    do
       call random_number(random_0_to_1) !!$returns a random value from 0 to 1
       randomInt = nint( (boxSize-1) * random_0_to_1 ) + 1 !!$ (511 * x) + 1
       if(any(randMask.eq.randomInt).or.(randomInt.eq.boxCenter)) cycle !!$ Enforce element uniqueness
       randMask(i) = randomInt       
       i = i+1   
       if (i.gt.maskSize) exit   
    end do

 
    !!$ DEBUG:
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
 
       !!$ Print random training points indices:
       print*,'randomMatrix'
       do k=1,box
       do i=1,box
          print*,randomMatrix(i,:,k)
       end do
          print*,''
       end do
    end if

  
  return       
  end subroutine randTrainingSet


!!$  COMPUTE SGS STRESS:
  
  subroutine synStress(u_f, u_t, tau_ij, T_ij, n_u, n_uu)
    implicit none

    !!$ ARGUMENT LIST:
    real(8),dimension(:,:,:,:),intent(in) :: u_f,    u_t
    real(8),dimension(:,:,:,:),intent(in) :: tau_ij, T_ij
    integer,intent(in) :: n_u, n_uu

    !!$ LOCAL VARAIBLES:
    real(8),allocatable,dimension(:,:,:,:) :: uu_t, uu_f
    real(8),allocatable,dimension(:,:,:,:) :: TijOpt, tau_ijOpt !Computed stresses

    !!$ FOR THE LEAST SQUARES PROBLEM:
    real(8),dimension(M,N) :: V = 0. !!$ Non-linear combination of u_i, u_iu_j
    real(8),dimension(N)   :: h_ij = 0. !!$ T takes the T_ij values and copies to h_ij in inverse().
    real(8),dimension(M)   :: T = 0.
    real(8)                :: lambda(5) = (/ 1.d2, 1.1d2, 1.2d2, 1.3d2, 1.4d2 /) 
    integer                :: row_index, col_index, row, col !!$ Indices to build V,h_ij arrays.
    integer                :: u_comp, uu_comp !!$ Indices to select u_i,u_iu_j components

    !!$ FOR RANDOM TRAINING POINTS (STENCIL-CENTERS):
    integer :: rand_count
    integer :: randMask(boxSize-M) !!$512-243=269

    !!$ LIST OF INDICES FOR LOOPS:
    integer :: i_test,    j_test,    k_test !for a point in the TEST-scale field
    integer :: i_box,     j_box,     k_box  !for a point in bounding BOX (stencil-center)
    integer :: i_stencil, j_stencil, k_stencil !for a point in the 3x3x3 STENCIL
    integer :: i_opt,     j_opt,     k_opt     !for computed(OPTIMIZED)values at LES and Test scales.
    integer :: i_proj,    j_proj,    k_proj    !to PROJECT final computed data for postprocessing
    integer :: i,         j,         k     !General indices

    !!$ TO SAVE RESULTS:
    character(50)::  PATH="./testResults/dampedLeast/bin4020/"
    character(10)::  l_val = "l_4000/" !!$ l_var stands for variable lambda.
    logical      ::  writeStress = .FALSE.    
    character(4) ::  z_plane
    integer      ::  vizPlaneC, vizPlaneF
       
    !!$ DEBUG SWITCHES:
    integer,dimension(4) :: debug=(/1,0,1,0/)
    integer :: lim, p


   print*, 'Computing SGS stress...'

   !!$ To determine stresses at coarse and fine stencil,
   !!$ velocities and their products must be known at LES
   !!$ scale. Though test scale will skip every other point,
   !!$ the fine stencil wiil need LES scale values to compute tau_ijOpt

   allocate(uu_t(n_uu,testcutSize,testcutSize,testcutSize)) !(6,51,51,51)
   allocate(uu_f(n_uu,testcutSize,testcutSize,testcutSize))
   
   if (debug(1).eq.1) then
      allocate(TijOpt(1,testSize,testSize,testSize)) !(1,17,17,17)
      allocate(tau_ijOpt(1,testSize,testSize,testSize))
   else
      allocate(TijOpt(n_uu,testSize,testSize,testSize)) !(6,17,17,17)
      allocate(tau_ijOpt(n_uu,testSize,testSize,testSize))
   end if

   if(debug(2).eq.1) then
      print *,'shape tau_ij cutout: ',shape(tau_ij)
      print*,'shape uu_t cutout:',shape(uu_t)
      print*,'shape tau_ijOpt cutout',shape(tau_ijOpt)
   end if

   !!$ Compute velocity products:
   k=0
   do j=1,n_u
   do i=1,n_u
      if (i.ge.j) then          !!$ Lower triangle order for ij
         k=k+1
         uu_t(k,:,:,:) = u_t(i,:,:,:) * u_t(j,:,:,:)
         uu_f(k,:,:,:) = u_f(i,:,:,:) * u_f(j,:,:,:)
      end if
   end do
   end do
    

!! SINGLE POINT COMPUTATION: COMMENT THIS PART OUT FOR A COMPLETE RUN:
   if(debug(3).eq.1)then
      lim = testLower
      print*, 'Check for the first element...(11,11,11)'
      print('(2(a15,f20.15))'), 'T_ij:',T_ij(1,lim,lim,lim),'tau_ij:',tau_ij(1,lim,lim,lim)
   end if   
   do k_test = lim,lim,stride
   do j_test = lim,lim,stride
   do i_test = lim,lim,stride 


!! WHOLE DOMAIN COMPUTATION: COMMENT THIS PART OUT FOR SINGLE POINT COMPUTATION.   
!!$     do k_test = testLower, testUpper, stride 
!!$     do j_test = testLower, testUpper, stride
!!$     do i_test = testLower, testUpper, stride !!$ i_test = 11,43,2


     !!$ GENERATE RANDOM STENCIL CENTER POINTS:
     call randTrainingSet(randMask)   
     rand_count = 0 
     row_index  = 0 
     
     do k_box = k_test-boxLower, k_test+boxUpper, stride
     do j_box = j_test-boxLower, j_test+boxUpper, stride
     do i_box = i_test-boxLower, i_test+boxUpper, stride !!$ i_box = 3,49,2
     
        rand_count = rand_count + 1
        if (any(randMask.eq.rand_count)) cycle         !  Skip if the point is listed on randMask

        col_index = 0 
        row_index = row_index + 1

        do k_stencil = k_box-stride, k_box+stride, stride
        do j_stencil = j_box-stride, j_box+stride, stride
        do i_stencil = i_box-stride, i_box+stride, stride

           do u_comp = 1,n_u ! 1 to 3
              col_index = col_index+1
              V(row_index,col_index) = u_t(u_comp,i_stencil,j_stencil,k_stencil)
           end do
    
           do uu_comp = 1,n_uu ! 1 to 6
              col_index = col_index+1
              V(row_index,col_index) = uu_t(uu_comp,i_stencil,j_stencil,k_stencil)
           end do

        end do
        end do
        end do ! stencil        

        T(row_index) = T_ij(1,i_box,j_box,k_box)

     end do
     end do
     end do ! box
     
    
     !!$ 1) Until here, the mechanism to compute the coefficient matrix
     !!$ is tested and verified. There are no issues with the basic mechanism.
     !!$ Only variable strides, for more than 2 is to be adjusted in the code later.

     !!$ 2) For testing, only the '_11' component is evaluated. Later extend to all 6 ij's
     !!$ V matrix generated can be used to determine Volterra coefficients for all six ij components.


     !!$ inverse takes in T as T_ij for each of the point in the bounding box.
     lam:do i=1,5
        TijOpt = 0.d0
        tau_ijOpt = 0.d0
        
        call inverse(V,T,h_ij,lambda(i))     !!$ SOLVE THE INVERSE PROBLEM:



     !!$ The test scale computations will involve creation
     !!$ of a strided matrix. Since an array stores contiguous
     !!$ set of memory, two things will happen if the array index
     !!$ is not carefully set - if there is a regular chuck of
     !!$ heap memory, then it creates a (51,51,51) array and fills
     !!$ unused space with 0. If not, it results in a
     !!$ segmentation fault. Thus the index of the strided elements
     !!$ that are used to form TijOpt and tau_ijOpt arrays should be
     !!$ mapped with _test indices to give a contiguous array (17,17,17)

     i_proj = (i_test-testLower)/stride + 1
     j_proj = (j_test-testLower)/stride + 1
     k_proj = (k_test-testLower)/stride + 1
 
     !!$ COARSE STENCIL: Project T_ij back to itself to compare with the original T_ij field
     col = 0
     do k_opt = k_test-stride, k_test+stride, stride
     do j_opt = j_test-stride, j_test+stride, stride
     do i_opt = i_test-stride, i_test+stride, stride

        do u_comp = 1,n_u !!$ 1 to 3
           col = col+1
           TijOpt(1,i_proj,j_proj,k_proj) = TijOpt(1,i_proj,j_proj,k_proj)&
                + u_t(u_comp,i_opt,j_opt,k_opt) * h_ij(col)
        end do
        do uu_comp = 1,n_uu !!$ 1 to 6
           col = col+1
           TijOpt(1,i_proj,j_proj,k_proj) = TijOpt(1,i_proj,j_proj,k_proj)&
                + uu_t(uu_comp,i_opt,j_opt,k_opt) * h_ij(col)
        end do

     end do
     end do
     end do 

     !!$ 1) Calculate bias and varience and repeat until both are minimized.
     !!$ Establish metrics to compute the figure of merit. Compute R-value, and change lambda
     !!$ to see when it apporaches the "right" value for each of the point.
     !!$ 
     !!$ 2) Need an error measuring parameter. From Hastie, it is the EPE (Expected
     !!$ Prediction Error). EPE is related to Bayesian probabilistic distributions.


     !!$ FINE STENCIL : Calculate the SGS stress with h_ij;compared with the original tau_ij     
     col = 0 
     do k_opt = k_test-1, k_test+1
     do j_opt = j_test-1, j_test+1
     do i_opt = i_test-1, i_test+1

        do u_comp = 1,n_u !!$ 1 to 3
           col = col+1
           tau_ijOpt(1,i_proj,j_proj,k_proj) = tau_ijOpt(1,i_proj,j_proj,k_proj)&
                + u_f(u_comp,i_opt,j_opt,k_opt) * h_ij(col)
        end do
        do uu_comp = 1,n_uu !!$ 1 to 6
           col = col+1
           tau_ijOpt(1,i_proj,j_proj,k_proj) = tau_ijOpt(1,i_proj,j_proj,k_proj)&
                + uu_f(uu_comp,i_opt,j_opt,k_opt) * h_ij(col)
        end do
     end do
     end do
     end do 
     
     print('(2(a15,f20.15),a10,f9.1)'),'TijOpt:',TijOpt(1,i_proj,j_proj,k_proj),&
          'tau_ijOpt:',tau_ijOpt(1,i_proj,j_proj,k_proj),'lambda:',lambda(i)

     end do lam

  end do
  end do
  end do !!$ test


!!$ FOR POST-PROCESSING:
!!$ Files to save results in:
!!$ Convert to a subroutine later 
  if(writeStress) then
     do k = 1,3
        vizPlaneC = nint(0.25 * k * testSize) ! Mark three slices at 1/4, 1/2, 3/4 z-planes.
        vizPlaneF = (vizPlaneC - 1)*2 + testLower
        write(z_plane, '(i1)'), k

        open(1,file= trim(PATH)//trim(l_val)//'T_ij_z'//trim(z_plane)//'.dat',status='replace')
        open(2,file= trim(PATH)//trim(l_val)//'TijOpt_z'//trim(z_plane)//'.dat',status='replace')
        open(3,file= trim(PATH)//trim(l_val)//'tau_ij_z'//trim(z_plane)//'.dat',status='replace')
        open(4,file= trim(PATH)//trim(l_val)//'tau_ijOpt_z'//trim(z_plane)//'.dat',status='replace')

        write(1,*), (T_ij(1,i,testLower:testUpper:stride,vizPlaneF), i=testLower,testUpper,stride)
        write(2,*), (TijOpt(1,i,:,vizPlaneC), i=1,testSize)
        write(3,*), (tau_ij(1,i,testLower:testUpper:stride,vizPlaneF), i=testLower,testUpper,stride)
        write(4,*), (tau_ijOpt(1,i,:,vizPlaneC), i=1,testSize)
    
        close(1)
        close(2)
        close(3)
        close(4)
     end do
  end if




!!$ ***************DEBUG********************
  if(debug(4).eq.1) then

  if (any(randMask.eq.rand_count)) then
  if ((i_stencil-stride).eq.testCutSize)then
     print*, "Error : Flag in last stencil cell!!$"
  end if
  end if

  print*,'Testing for V:'
  i_stencil=0;j_stencil=0;k_stencil=0
  k=0 ; p=0
  do k_stencil = 3-stride,3+stride,stride
  do j_stencil = 3-stride,3+stride,stride
  do i_stencil = 3-stride,3+stride,stride
     p=p+1
     print*,p
     do i=1,n_u
        k=k+1
        print*, 'u_t',i,i_stencil,j_stencil,k_stencil,u_t(i,i_stencil,j_stencil,k_stencil),'V:',k,V(1,k)
     end do
     do i=1,n_uu
        k=k+1
        print*, 'uu_t',i,i_stencil,j_stencil,k_stencil,uu_t(i,i_stencil,j_stencil,k_stencil),'V:',k,V(2,k)
     end do
     print*,''
  end do
  end do
  end do

  end if
  !!$ ****************************************


  deallocate(uu_t,uu_f,TijOpt)



 
  return
  end subroutine synStress




  
  subroutine inverse(V,T_ij,h_ij,lambda)
    implicit none

  !  integer, parameter :: n = M !!$ M is calculated right in the beginning
    real(8), dimension(M,N),intent(in) :: V 
    real(8), dimension(M),intent(in) :: T_ij
    real(8), intent(in) :: lambda
    real(8), dimension(N),intent(out) :: h_ij
    real(8), dimension(N,N) :: A,eye
    real(8), dimension(N)   :: b,work

    real(8) :: errnorm, xnorm, rcond, anorm, colsum
    integer :: i,j,info, LDA = N, LDB = N, nrhs = 1
    integer, dimension(N) :: ipiv
    !!$integer, dimension(:),allocatable :: iwork
    integer :: d
    character, dimension(1) :: norm
    logical :: debug = .false., print_h = .false., cond=.false.

    !!$ Create matrix for computing the direct inverse:
    forall(d = 1:N) eye(d,d) = 1.d0 !!$ create Identity matrix
    !!$ Use the SAVE attribute or something to avoid repeated construction.
    
    A = matmul(transpose(V),V) + lambda*eye !!$ V'V + lambda*I
    b = matmul(transpose(V),T_ij) !!$ Pass initial value from T_ij to h_ij for LU decomposition

    if (debug) then
       do j=1,18
          write(*,"(4(F16.2))") V(j,1:3),T_ij(j)
       end do
    end if

    
    !!$ compute 1-norm needed for condition number
   if(cond)then
      anorm = 0.d0
      do j=1,N
         colsum = 0.d0
         do i=1,N
            colsum = colsum + abs(A(i,j))
         enddo
         anorm = max(anorm, colsum)
      enddo
   end if

       
    !!$ dgesv destroys the orginal matrix. So V is copied into 'a' matrix. This is the LU product matrix.
    !!$ dgesv returns x vector. First it copies the values from 'T_ij' vector and then computes 'h_ij'.
    !!$ 'h_ij' is the h_ij vector.
    
    call dgesv(N, nrhs, A, LDA, ipiv, b, LDB, info)
    h_ij = b

    if(print_h) then
       print*,'lambda',lambda
       do j=1,N
          write(*,"(i3,',',f22.15)") j,h_ij(j)
       end do
       print*,''
    end if

    if (debug) then
       print*,"After dgesv operation:"
       do j=1,18
          write(*,"(5(F16.4) )") a(j,1:3),T_ij(j),h_ij(j)
       end do
    end if



 
   
!!$$    ! compute 1-norm of error
!!$    errnorm = 0.d0
!!$    xnorm = 0.d0
!!$    do i=1,M
!!$        errnorm = errnorm + abs(h_ij(i)-T_ij(i))
!!$        xnorm = xnorm + abs(h_ij(i))
!!$    enddo
!!$
!!$    ! relative error in 1-norm:
!!$    errnorm = errnorm / xnorm


    !!$ compute condition number of matrix:
    !!$ note: uses V returned from dgesv with L,U factors:

!!$     allocate(work(4*M))
!!$     allocate(iwork(M))
!!$     norm = '1'  !!$ use 1-norm
!!$     call dgecon(norm,M,a,lda,anorm,rcond,work,iwork,info)

!!$     if (info /= 0) then
!!$        print *, "*** Error in dgecon: info = ",info
!!$     endif

!!$     if(debug) then
!!$     print 201, M, 1.d0/rcond
!!$     201 format("For M = ",i4," the approx. condition number is ",e10.3,/, &
!!$            " and the relative error in 1-norm is ",e10.3)    
!!$  end if
 


    !!$deallocate(work,iwork)
        
    return
  end subroutine inverse
  



  
  
  subroutine bubblesort(array)
    implicit none
    integer, dimension(:)::array
    integer :: i,j,n_sort,temp

    n_sort = size(array,dim=1)

     bubble:do i=1,n_sort
            do j=1,n_sort-i
                if( array(j+1).lt.array(j))then
                   temp = array(j)
                   array(j) = array(j+1)
                   array(j+1) = temp
                end if
             end do
             end do bubble
     return
   end subroutine bubblesort





  
  
end module linsolve

