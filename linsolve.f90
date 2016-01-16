module linsolve

  !! Linsolve is a collection of subroutines that perform linear algebra operations
  !! and related operations to solve the inverse problem.
 
  use fileio

  integer,parameter :: stride = 2
  
  ! Stencil parameters:
  integer,parameter :: coloc2 = 27 * 9   ! For a 3x3x3 stencil
  
  ! Bounding Box parameters:
  real,   parameter :: eps = 1e-3
  integer,parameter :: box = 8
  integer,parameter :: boxSize = box**3
  integer,parameter :: bigHalf   = ceiling(0.5*real(box)+eps) ! for 8->5
  integer,parameter :: smallHalf = floor(0.5*real(box)+eps)   ! for 8->4
  integer,parameter :: boxCenter = smallHalf * box*(box+1) + bigHalf
  integer,parameter :: boxLower  = stride * (bigHalf-1)
  integer,parameter :: boxUpper  = stride * (box-bigHalf)

  ! Test field parameters: 
  integer,parameter :: testSize = 17
  integer,parameter :: testcutSize = stride * (testSize+box) + 1
  integer,parameter :: testLower = stride * bigHalf + 1
  integer,parameter :: testUpper = stride * (bigHalf-1+testSize) + 1

  ! Cutout parameters:
  integer,parameter :: lBound = 0.5*(GRID - testcutSize)
  integer,parameter :: uBound = 0.5*(GRID + testcutSize) - 1

contains

  subroutine printParams()
    implicit none
    print*, 'Stencil parameters:'
    print*, 'colocated variables :',coloc2
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
    print*, ''
    print*, 'lower bound:         ',lBound
    print*, 'upper bound:         ',uBound
    return
  end subroutine printParams
  
  subroutine cutout(array,n)
    implicit none

    integer, intent(in) :: n
    real(8), allocatable, dimension(:,:,:,:),intent(inout) :: array
    real(8), allocatable, dimension(:,:,:,:):: temp
   
    allocate (temp(n,testcutSize,testcutSize,testcutSize))
    temp = array(:,lBound:uBound,lBound:uBound,lBound:uBound)
    deallocate(array)
    
    allocate(array(n,testcutSize,testcutSize,testcutSize))
    array = temp   
    deallocate(temp)
    
    return
  end subroutine cutout

  
  ! RANDOM NUMBER GENERATION
  subroutine init_random_seed()
      
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed
  
    call RANDOM_seed(size = n)
    allocate(seed(n))
    call SYSTEM_clock(COUNT=clock)
    seed = clock + 37 * (/ (i-1, i=1,n) /)
    call RANDOM_seed(PUT = seed)
    deallocate(seed)
    return
  end subroutine init_random_seed

  subroutine randAlloc(num)
    ! Status: Integrating with bounding box.
    ! Notes: Creates a mask to skip cells randomly.
    !      : Always retains the center point.
    implicit none

    real(8)                 :: x 
    integer, parameter      :: n = boxSize - coloc2 ! 512-378 = 134
    integer                 :: randomNum, num(n)
    integer,dimension(box,box,box):: randomMatrix
    integer                 :: i, j, k, index, c
    logical,parameter       :: debug = .false.

    num = 0 ; i=1

    call init_random_seed()
    do
       call random_number(x)
       randomNum = nint( (boxSize-1) * x ) + 1
       
       ! Enforce element uniqueness and keep boxCenter -> see Notes
       if(any(num.eq.randomNum).or.(randomNum.eq.boxCenter)) cycle
       
       num(i) = randomNum       
       i = i+1   
       
       if (i.gt.n) exit   
    end do
    call bubblesort(num) ! Remove in real case.Keep for testing
    print*,num    
   !! ******************DEBUG**********************
    if (debug) then
    c=0;index=1;randomMatrix=0
    do k=1,box
    do j=1,box
    do i=1,box
       c = c+1
       if (c.eq.num(index))then
          randomMatrix(i,j,k) = num(index)
          index = index+1
       end if
    enddo
    enddo
    enddo
 
    print*,'randomMatrix'
    do k=1,box
    do i=1,box
       print*,randomMatrix(i,:,k)
    end do
    print*,''
    end do
    print*,'boxCenter',boxCenter
    end if
  !! ***********************************************
  
  return       
  end subroutine randAlloc


!  SGS STRESS
  
  subroutine synStress(u_f,u_t,tau_ij,T_ij,n_u,n_uu)
    implicit none

    real(8),dimension(:,:,:,:),intent(in)  :: u_f, u_t ,tau_ij,T_ij
    real(8),allocatable,dimension(:,:,:,:) :: uu_t,uu_f,TijOpt,tau_ijOpt
    real(8),dimension(coloc2,coloc2) :: V = 0. ! Is the non-linear combination of the velocities.
    real(8),dimension(coloc2) :: T = 0. , h_ij = 0. ! T takes the T_ij values and copies to h_ij in inverse().
    real(8) :: lambda(5) = (/ 1.d4, 1.1d2, 1.2d2, 1.3d2, 1.4d2 /) 

    character(50)::  PATH="./testResults/dampedLeast/bin4020/"
    character(10)::  l_val = "l_4000/" ! l_var stands for variable lambda.
    character(4) ::  z_plane
    integer :: n_u
    integer :: n_uu

    integer :: randMask(boxSize-coloc2) !512-378=134; 512-270=240

    integer :: i_test,    j_test,    k_test 
    integer :: i_box,     j_box,     k_box  
    integer :: i_stencil, j_stencil, k_stencil
    integer :: i_opt,     j_opt,     k_opt
    integer :: i_proj,    j_proj,    k_proj
    
    
    integer :: rand_count,lim
    integer :: row_index, col_index , row, col
    integer :: u_comp,    uu_comp

    integer,dimension(4) :: debug=(/0,1,0,0/)

    integer :: i,j,k,p=0,z,vizPlaneC, vizPlaneF
    logical :: writeStress=.false.
 
   print*, 'Computing SGS stress...'
   print *,'shape tau_ij cutout: ',shape(tau_ij)

   ! To determine stresses at corase and fine stencil,
   ! velocities and their products must tbe known at LES
   ! scale. Though test scale will skip every other point,
   ! the fine stencil needed to compute synthetic SGS scale
   ! will be required. 
   allocate(uu_t(n_uu,testcutSize,testcutSize,testcutSize))
   allocate(uu_f(n_uu,testcutSize,testcutSize,testcutSize))
   
   ! TijOpt is synthetic Test scale stress, while tau_ijOpt
   ! is the computed stress at LES scale. 
   allocate(TijOpt(n_uu,testSize,testSize,testSize)) 
   allocate(tau_ijOpt(n_uu,testSize,testSize,testSize)) 
   print*,'shape uu_t cutout:',shape(uu_t)
   print*,'shape tau_ijOpt cutout',shape(tau_ijOpt)
   ! I could have used a 5-D matrix and then the FORALL statement
   ! But this makes it more complicated to manage. Better to stick
   ! with 4-D matrix.
   print*, 'Compute velocity products:'
   k=0
   do j=1,n_u
   do i=1,n_u
      if (i.ge.j) then          ! Lower triangle 
         k=k+1
         uu_t(k,:,:,:) = u_t(i,:,:,:) * u_t(j,:,:,:)
         uu_f(k,:,:,:) = u_f(i,:,:,:) * u_f(j,:,:,:)
      end if
   end do
   end do
    

  
! 2*(5+1)-1 for 8 = 11 ; 5 for box,1 for stencil
   ! As a modification, it should be independent of stride=2 , as well

   if(debug(3).eq.1)then
      lim=testUpper
!      print*, 'Check for the last element..(43,43,43)'
!      print*, u_t(1,uBound,uBound,uBound)
   else
      lim=testLower
!      print*, 'Check for the first element...(11,11,11)'
!      print*, T_ij(1,testLower,testLower,testLower),         &
!            tau_ij(1,testLower,testLower,testLower)
      print*,''
   end if

   lim = 11
   print*, T_ij(1,lim,lim,lim), tau_ij(1,lim,lim,lim)
   print*, ''
   

   
!     do k_test = testLower, testUpper, stride 
!     do j_test = testLower, testUpper, stride
!     do i_test = testLower, testUpper, stride ! i_test = 11,43,2


  do k_test = lim,lim,stride
  do j_test = lim,lim,stride
  do i_test = lim,lim,stride 

     call randAlloc(randMask)   !Each point will have different stencil-center points.
     rand_count = 0 
     row_index  = 0 
     
     open(111,file='./testResults/h_ij/randomPoints.csv',status='replace') !Save random point co-od on CSV file to plot in Paraview

     do k_box = k_test-boxLower, k_test+boxUpper, stride
     do j_box = j_test-boxLower, j_test+boxUpper, stride
     do i_box = i_test-boxLower, i_test+boxUpper, stride ! i_box = 3,49,2
        rand_count = rand_count + 1
        ! Skip if the point is listed in randMask
        if (any(randMask.eq.rand_count)) cycle
        write(111,*) rand_count
        ! Reset pointer to the first column after each stencil operation
        col_index = 0 
        row_index = row_index + 1

        do k_stencil = k_box-stride, k_box+stride, stride
        do j_stencil = j_box-stride, j_box+stride, stride
        do i_stencil = i_box-stride, i_box+stride, stride

           ! Compute non-linear combinations to generate V matrix:
           ! Get all three velocity components for the V matrix
           do u_comp = 1,n_u ! 1 to 3
              col_index = col_index+1
              V(row_index,col_index) = u_t(u_comp,i_stencil,j_stencil,k_stencil)
           end do
           ! Get all six velocity products for the V matrix:                  
           do uu_comp = 1,n_uu ! 1 to 6
              col_index = col_index+1
              V(row_index,col_index) = uu_t(uu_comp,i_stencil,j_stencil,k_stencil)
           end do

        end do
        end do
        end do !stencil

        ! At this point, the V matrix generated can be used to
        ! determine Volterra coefficients for all six ij components.
        ! For testing, only 11 component evaluated. Later extend to all 6 ij's
        T(row_index) = T_ij(1,i_box,j_box,k_box)  ! Testing for T_11

     end do
     end do
     end do !box
     
     close(111)                 !111 is the code for CSV files

     ! Until here, the mechanism to compute the coefficient matrix
     ! is tested and verified. There are no issues with the basic mechanism.
     ! Only variable strides, for more than 2 is to be adjusted in the code later.

     ! SOLVE THE INVERSE PROBLEM:
     ! inverse takes in T as T_ij for each of the point in the bounding box.
     lam:do i=1,1
        TijOpt = 0.d0
        tau_ijOpt = 0.d0
        
        call inverse(V,T,h_ij,lambda(i))

 !!$         p = p+1 !should be --> 4913
 !!$         if(p.eq.1.or.p.eq.10.or.p.eq.20.or.p.eq.40)


     ! Check if h_ij computed gives close value for T_ij:
     ! This is done at test scale
    
     ! The test scale computations will involved creation
     ! of a strided matrix. Since an array stores contiguous
     ! set of memory, two things will happen if the array index
     ! is not carefully set - if there is a regular chuck of
     ! heap memory, then it creates a (51,51,51) array and fills
     ! unsused space with 0. If not, then it results into a
     ! segmentation fault. Thus the index of the strided elements
     ! that are used to form TijOpt and tau_ijOpt arrays should be
     ! mapped with _test indices to give a continous array (17,17,17)
     col = 0
     i_proj = (i_test-testLower)/stride + 1
     j_proj = (j_test-testLower)/stride + 1
     k_proj = (k_test-testLower)/stride + 1
 
     ! COARSE STENCIL : Project T_ij back to itself
     ! to compare with the original T_ij field
     do k_opt = k_test-stride, k_test+stride, stride
     do j_opt = j_test-stride, j_test+stride, stride
     do i_opt = i_test-stride, i_test+stride, stride

        do u_comp = 1,n_u ! 1 to 3
           col = col+1
           TijOpt(1,i_proj,j_proj,k_proj) = TijOpt(1,i_proj,j_proj,k_proj)&
                + u_t(u_comp,i_opt,j_opt,k_opt) * h_ij(col)
        end do
        do uu_comp = 1,n_uu ! 1 to 6
           col = col+1
           TijOpt(1,i_proj,j_proj,k_proj) = TijOpt(1,i_proj,j_proj,k_proj)&
                + uu_t(uu_comp,i_opt,j_opt,k_opt) * h_ij(col)
        end do

     end do
     end do
     end do 

     ! FINE STENCIL : Calculate the SGS stress with determined features
     ! and this will be compared with the original tau_ij

     ! Begin auto-tune : Compute tau_ijOpt and determine if it is good enough for the point
     ! or not. Calculate bias and varience and repeat until both are minimized.
     ! Establish metrics to compute the figure of merit. Compute R-value, and change lambda
     ! to see when it apporaches the "right" value for each of the point.
     ! 
     ! Need a error measuring parameter. From Hastie, it is the EPE (Expected
     ! Prediction Error). Now EPE is related with Bayesian probabilistic distributions.
     col = 0 
     do k_opt = k_test-1, k_test+1
     do j_opt = j_test-1, j_test+1
     do i_opt = i_test-1, i_test+1

        do u_comp = 1,n_u ! 1 to 3
           col = col+1
           tau_ijOpt(1,i_proj,j_proj,k_proj) = tau_ijOpt(1,i_proj,j_proj,k_proj)&
                + u_f(u_comp,i_opt,j_opt,k_opt) * h_ij(col)
        end do
        do uu_comp = 1,n_uu ! 1 to 6
           col = col+1
           tau_ijOpt(1,i_proj,j_proj,k_proj) = tau_ijOpt(1,i_proj,j_proj,k_proj)&
                + uu_f(uu_comp,i_opt,j_opt,k_opt) * h_ij(col)
        end do

     end do
     end do
     end do 

!      print'(9i5,4f12.4)',i_test,j_test,k_test,&
!           i_box,j_box,k_box,&
!           i_proj,j_proj,k_proj,&
!           maxval(h_ij),minval(h_ij),maxval(V),minval(V)

     
     print*, 'TijOpt:',TijOpt(1,i_proj,j_proj,k_proj),&
          'tau_ij:',tau_ijOpt(1,i_proj,j_proj,k_proj)
     print*,'M:', row_Index
  end do lam
     ! print*,'Tij',      T_ij     (1,i_test,j_test,k_test)
!      print*,'TijOpt',   TijOpt   (1,i_proj,j_proj,k_proj)
!      print*,'tau_ij',   tau_ij   (1,i_test,j_test,k_test)
!      print*,'tau_ijOpt',tau_ijOpt(1,i_proj,j_proj,k_proj)
!      print*,i_proj,j_proj,k_proj

 !!$       end if

  end do
  end do
  end do ! test

 ! Files to save results in:

  if(writeStress) then
  do z=1,3
   vizPlaneC = nint(0.25*z*testSize)  
   vizPlaneF = (vizPlaneC-1)*2+testLower
   write(z_plane, '(i1)'), z

   open(1,file= trim(PATH)//trim(l_val)//'T_ij_z'//trim(z_plane)//'.dat',status='new',action='write')
   open(2,file= trim(PATH)//trim(l_val)//'TijOpt_z'//trim(z_plane)//'.dat',status='new',action='write')
   open(3,file= trim(PATH)//trim(l_val)//'tau_ij_z'//trim(z_plane)//'.dat',status='new',action='write')
   open(4,file= trim(PATH)//trim(l_val)//'tau_ijOpt_z'//trim(z_plane)//'.dat',status='new',action='write')

  ! Print results:
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




!! ***************DEBUG********************
  if(debug(4).eq.1) then
  print*,''
  print*, 'Exclusion list:'
  print*, randMask
  print*,''
  print *, 'Number of random cells sampled:',rand_count

  if (any(randMask.eq.rand_count)) then
  if ((i_stencil-stride).eq.testCutSize)then
     print*, "Error : Flag in last stencil cell!"
  end if
  end if

  print*, randMask   
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
  !! ****************************************


  deallocate(uu_t,uu_f,TijOpt)

!!$open(1,file="eig.dat")
!!$do j=1,243
!!$   write(1,202) V(j,:)
!!$end do     
!!$202 format(243f16.4)                                          
!!$close(1)                       
 
  return
  end subroutine synStress




  
  subroutine inverse(V,T_ij,h_ij,lambda)
    implicit none

    integer, parameter :: n = coloc2 ! coloc2 is calculated right in the beginning
    real(8), dimension(n,n),intent(in) :: V 
    real(8), dimension(n),intent(in) :: T_ij
    real(8), intent(in) :: lambda
    real(8), dimension(n),intent(out) :: h_ij
    real(8), dimension(n,n) :: a,eye
    real(8), dimension(n)   :: b,work

    real(8) :: errnorm, xnorm, rcond, anorm, colsum
    integer :: i,info, lda, ldb, nrhs, j
    integer, dimension(n) :: ipiv
    !integer, dimension(:),allocatable :: iwork
    integer :: d
    character, dimension(1) :: norm
    logical :: debug = .false., print_h = .false., cond=.false.

    ! Create matrix for computing the direct inverse:
    forall(d = 1:n) eye(d,d) = 1.d0 ! create Identity matrix
    ! Use the SAVE attribute or something to avoid repeated construction.
    
    a = matmul(transpose(V),V) + lambda*eye ! V'V + lambda*I
    b = matmul(transpose(V),T_ij) ! Pass initial value from T_ij to h_ij for LU decomposition

    if (debug) then
       do j=1,18
          write(*,"(4(F16.2))") V(j,1:3),T_ij(j)
       end do
    end if

    
    ! compute 1-norm needed for condition number
   if(cond)then
      anorm = 0.d0
      do j=1,n
         colsum = 0.d0
         do i=1,n
            colsum = colsum + abs(a(i,j))
         enddo
         anorm = max(anorm, colsum)
      enddo
   end if

       
    nrhs = 1 ! number of right hand sides in T_ij  ! extend upto 6 for 6 tau_ijs
    lda = n  ! leading dimension of a
    ldb = n  ! leading dimension of T_ij

    ! dgesv destroys the orginal matrix. So V is copied into 'a' matrix. This is the LU product matrix.
    ! dgesv returns x vector. First it copies the values from 'T_ij' vector and then computes 'h_ij'.
    ! 'h_ij' is the h_ij vector.
    
    call dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
    h_ij = b

    if(print_h) then
       print*,'lambda',lambda
       do j=1,n
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



 
   
!!$    ! compute 1-norm of error
!!$    errnorm = 0.d0
!!$    xnorm = 0.d0
!!$    do i=1,n
!!$        errnorm = errnorm + abs(h_ij(i)-T_ij(i))
!!$        xnorm = xnorm + abs(h_ij(i))
!!$    enddo
!!$
!!$    ! relative error in 1-norm:
!!$    errnorm = errnorm / xnorm


    ! compute condition number of matrix:
    ! note: uses V returned from dgesv with L,U factors:

!     allocate(work(4*n))
!     allocate(iwork(n))
!     norm = '1'  ! use 1-norm
!     call dgecon(norm,n,a,lda,anorm,rcond,work,iwork,info)

!     if (info /= 0) then
!        print *, "*** Error in dgecon: info = ",info
!     endif

!     if(debug) then
!     print 201, n, 1.d0/rcond
!     201 format("For n = ",i4," the approx. condition number is ",e10.3,/, &
!            " and the relative error in 1-norm is ",e10.3)    
!  end if
 


    !deallocate(work,iwork)
        
    return
  end subroutine inverse
  

  
  
  subroutine bubblesort(array)
    implicit none
    integer, dimension(:)::array
    integer :: i,j,n,temp

    n = size(array,dim=1)

     bubble:do i=1,n
             do j=1,n-i
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

