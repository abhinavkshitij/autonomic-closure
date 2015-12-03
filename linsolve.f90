module linsolve

  !! Linsolve is a collection of subroutines that perform linear algebra operations
  !! and related operations to solve the inverse problem.
 
  use fileio

  integer,parameter :: stride = 2
  
  ! Stencil parameters:
  integer,parameter :: coloc2 = 27*9   ! For a 3x3x3 stencil
  
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
    real(kind=8), allocatable, dimension(:,:,:,:),intent(inout) :: array
    real(kind=8), allocatable, dimension(:,:,:,:):: temp
   
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

    real(kind=8)            :: x 
    integer, parameter      :: n = boxSize - coloc2 ! 512-378 = 134
    integer                 :: randomNum, num(n)
    integer,dimension(box,box,box):: randomMatrix
    integer                 :: i, j, k, index, c
    logical,parameter       :: debug = .false.

    num = 0 ; i=0

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
    !call bubblesort(num) ! Remove in real case.Keep for testing
        
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


  ! SGS STRESS
  
  subroutine synStress(u_f,u_t,tau_ij,T_ij,n_u,n_uu)
    implicit none

    real(kind=8),dimension(:,:,:,:),intent(in)  :: u_f, u_t ,tau_ij,T_ij
    real(kind=8),allocatable,dimension(:,:,:,:) :: uu_t,uu_f,T_SGS
    real(kind=8),dimension(coloc2,coloc2) :: A = 0.
    real(kind=8),dimension(coloc2) :: b = 0.
    real(kind=8),dimension(33,33) :: T_res
    
    integer :: n_u
    integer :: n_uu

    integer :: randMask(boxSize-coloc2) !512-378=134; 512-270=240

    integer :: i_test,    j_test,    k_test 
    integer :: i_box,     j_box,     k_box  
    integer :: i_stencil, j_stencil, k_stencil
    integer :: i_sgs,     j_sgs,     k_sgs
    
    integer :: rand_count,lim
    integer :: row_index, col_index , row, col
    integer :: u_comp,    uu_comp

    integer,dimension(4) :: debug=(/0,1,0,0/)

    integer :: i,j,k,p=0
  
 
   print*, 'Computing SGS stress...'
   print *,'shape u cutout: ',shape(u_t)

   allocate(uu_t(n_uu,testcutSize,testcutSize,testcutSize))
   allocate(uu_f(n_uu,testcutSize,testcutSize,testcutSize))
   allocate(T_SGS(n_uu,testcutSize,testcutSize,testcutSize))
   print*,'shape uu cutout:',shape(uu_t)
   T_SGS = 0.

   print*, 'Compute velocity products:'
   k=0
   do j=1,n_u
      do i=1,n_u
         if (i.ge.j) then
            k=k+1
            uu_t(k,:,:,:) = u_t(i,:,:,:)*u_t(j,:,:,:)
            uu_f(k,:,:,:) = u_f(i,:,:,:)*u_f(j,:,:,:)
         end if       
      end do
   end do

   

! 2*(5+1)-1 for 8 = 11 ; 5 for box,1 for stencil
   ! As a modification, it should be independent of stride=2 , as well

!!$   if(debug(3).eq.1)then
!!$      lim=testUpper
!!$      print*, 'Check for the last element..(43,43,43)'
!!$      print*, u_t(1,uBound,uBound,uBound)
!!$   else
!!$      lim=testLower
!!$      print*, 'Check for the first element...(11,11,11)'
!!$      print*, u_t(1,testLower,testLower,testLower)
!!$   end if

do k_test = testLower, testUpper, stride 
   do j_test = testLower, testUpper, stride
      do i_test = testLower, testUpper, stride ! i_test = 11,43,2


!!$ do k_test = lim,lim,stride
!!$   do j_test = lim,lim,stride
!!$      do i_test = lim,lim,stride         
        
         rand_count = 0 
         row_index  = 0 

         call randAlloc(randMask)
          do k_box = k_test-boxLower, k_test+boxUpper,stride
            do j_box = j_test-boxLower, j_test+boxUpper,stride
               do i_box = i_test-boxLower,i_test+boxUpper,stride ! i_box = 3,49,2
                  rand_count = rand_count+1
                  ! Skip if the point is listed in randMask
                  if (any(randMask.eq.rand_count)) cycle
                  
                  ! Reset pointer to the first column after each stencil operation
                  col_index = 0 
                  row_index = row_index + 1

                   do k_stencil = k_box-stride,k_box+stride,stride
                     do j_stencil = j_box-stride,j_box+stride,stride
                       do i_stencil = i_box-stride,i_box+stride,stride

                          ! Compute non-linear combinations:
                           do u_comp = 1,n_u ! 1 to 3
                              col_index = col_index+1
                              A(row_index,col_index) = u_t(u_comp,i_stencil,j_stencil,k_stencil)
                           end do
                           do uu_comp = 1,n_uu ! 1 to 6
                              col_index = col_index+1
                              A(row_index,col_index) = uu_t(uu_comp,i_stencil,j_stencil,k_stencil)
                           end do                     
                        
                        end do
                     end do
                  end do !stencil
                  
                  b(row_index) = T_ij(1,i_box,j_box,k_box)
                  
               end do
            end do
         end do !box
         
         ! Solve the inverse problem
         call inverse(A,b)
         
!!$         p = p+1 !should be --> 4913
!!$         if(p.eq.1.or.p.eq.10.or.p.eq.20.or.p.eq.40)


         ! Compute SGS stress: NEEDS REVISION and TESTS
         col = 0
         do k_sgs = k_test-1,k_test+1
            do j_sgs = j_test-1,j_test+1
               do i_sgs = i_test-1,i_test+1 
                  
                  do u_comp = 1,n_u ! 1 to 3
                     col = col+1
                      T_SGS(1,i_test,j_test,k_test) = T_SGS(1,i_test,j_test,k_test)&
                           + u_f(u_comp,i_sgs,j_sgs,k_sgs)*b(col)
                   end do
                   do uu_comp = 1,n_uu ! 1 to 6
                      col = col+1
                      T_SGS(1,i_test,j_test,k_test) = T_SGS(1,i_test,j_test,k_test)&
                           + uu_f(uu_comp,i_sgs,j_sgs,k_sgs)*b(col)
                   end do
                end do
             end do
          end do ! SGS

!!$          print*,'T_SGS', T_SGS(1,i_test,j_test,k_test)
!!$          print*,'tau_ij',tau_ij(1,i_test,j_test,k_test)
!!$       end if
 
      end do
   end do
end do ! test


!!$open(1,file="T_SGS.dat")
!!$do j=1,33
!!$   do i=1,33
!!$      if (T_SGS(1,i,j,testLower).ne.0) then
!!$         T_res(i,j) = T_SGS(1,i,j,testLower)
!!$        write(1,*) T_res(i,j)
!!$      end if
!!$   end do
!!$   
!!$end do
!!$close(1)
   
! Comprison of results:
!!$open(1,file="T_SGS.dat")
!!$open(2,file="tau_ij.dat")                                       
!!$      
!!$! T_SGS FIELD:
!!$
!!$do j=testLower,testUpper
!!$   do i=testLower, testUpper
!!$      write(1,*) T_SGS(1,i,j,testLower)
!!$      write(2,*) tau_ij(1,i,j,testLower)
!!$   end do
!!$end do
!!$
!!$
!!$close(2)



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
print*,'Testing for A:'
i_stencil=0;j_stencil=0;k_stencil=0
k=0 ; p=0
do k_stencil = 3-stride,3+stride,stride
   do j_stencil = 3-stride,3+stride,stride
      do i_stencil = 5-stride,5+stride,stride
         p=p+1
         print*,p

do i=1,n_u
      k=k+1
      print*, 'u_t',i,i_stencil,j_stencil,k_stencil,u_t(i,i_stencil,j_stencil,k_stencil),'A:',k,A(1,k)
   end do

   do i=1,n_uu
      k=k+1
      print*, 'uu_t',i,i_stencil,j_stencil,k_stencil,uu_t(i,i_stencil,j_stencil,k_stencil),'A:',k,A(2,k)
   end do
   print*,''
         
         
      end do
   end do
end do

end if
!! ****************************************


deallocate(uu_t,uu_f,T_SGS)

!!$open(1,file="eig.dat")
!!$do j=1,243
!!$   write(1,202) A(j,:)
!!$end do     
!!$202 format(243f16.4)                                          
!!$close(1)                       
 
    return
  end subroutine synStress




  
  subroutine inverse(A_arr,b)
    implicit none

    integer, parameter :: n = coloc2
    real(kind=8), dimension(n,n),intent(in) :: A_arr
    real(kind=8), dimension(n),intent(inout) :: b
     real(kind=8), dimension(:),allocatable :: x
    real(kind=8), dimension(:,:),allocatable :: a
    real(kind=8), dimension(:),allocatable :: work
    
    real(kind=8) :: errnorm, xnorm, rcond, anorm, colsum
    integer :: i, info, lda, ldb, nrhs, j
    integer, dimension(:),allocatable :: ipiv
    integer, dimension(:),allocatable :: iwork
    character, dimension(1) :: norm
    logical :: debug = .false.

    allocate(a(n,n))
    allocate(x(n),ipiv(n))
    
    a = A_arr
    x = 0.

    if (debug) then
    do j=1,18
      write(*,"(4(F16.2) )") A_arr(j,1:3),b(j)
   end do
end if

    
    ! compute 1-norm needed for condition number

    anorm = 0.d0
    do j=1,n
        colsum = 0.d0
        do i=1,n
            colsum = colsum + abs(a(i,j))
            enddo
        anorm = max(anorm, colsum)
        enddo

        ! extend upto 6 for 6 tau_ijs
    nrhs = 1 ! number of right hand sides in b 
    lda = n  ! leading dimension of a
    ldb = n  ! leading dimension of b

    call dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)

    if (debug) then
   Print*,"LU:"
   do j=1,18
      write(*,"(4(F16.4) )") a(j,1:3),b(j)
   end do
end if




 
   
   
!!$    ! compute 1-norm of error
!!$    errnorm = 0.d0
!!$    xnorm = 0.d0
!!$    do i=1,n
!!$        errnorm = errnorm + abs(x(i)-b(i))
!!$        xnorm = xnorm + abs(x(i))
!!$    enddo
!!$
!!$    ! relative error in 1-norm:
!!$    errnorm = errnorm / xnorm


    ! compute condition number of matrix:
    ! note: uses A returned from dgesv with L,U factors:

    allocate(work(4*n))
    allocate(iwork(n))
    norm = '1'  ! use 1-norm
    call dgecon(norm,n,a,lda,anorm,rcond,work,iwork,info)

    if (info /= 0) then
       print *, "*** Error in dgecon: info = ",info
    endif

    if(debug) then
    print 201, n, 1.d0/rcond
    201 format("For n = ",i4," the approx. condition number is ",e10.3,/, &
           " and the relative error in 1-norm is ",e10.3)    
 end if
 

    deallocate(a,x,ipiv,work,iwork)
        

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

