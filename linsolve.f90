module linsolve
  use fileio

  integer,parameter :: stride = 2
  
  ! Stencil attributes:
  integer, parameter :: coloc2 = 27*9   ! For a 3x3x3 stencil
  
  ! Box attributes:
  real , parameter  :: eps = 1e-3
  integer,parameter :: box = 8
  integer,parameter :: boxSize = box**3
  integer,parameter :: bigHalf   = ceiling(0.5*real(box)+eps) ! for 8->5
  integer,parameter :: smallHalf = floor(0.5*real(box)+eps)   ! for 8->4
  integer,parameter :: boxCenter = smallHalf * box*(box+1) + bigHalf
  integer,parameter :: boxLower =  stride * (bigHalf-1)
  integer,parameter :: boxUpper =  stride * (box-bigHalf)

  ! Test field attributes: 
  integer,parameter :: testSize = 17
  integer,parameter :: testcutSize = stride*(testSize+box) + 1
  integer,parameter :: testLower = stride*bigHalf + 1
  integer,parameter :: testUpper = stride*(bigHalf-1+testSize) + 1
 
contains

subroutine cutout(array,n_u)
    implicit none
   
    real(kind=8), allocatable, dimension(:,:,:,:),intent(inout) :: array
    integer, intent(in) :: n_u
    
    real(kind=8), allocatable, dimension(:,:,:,:):: temp
    integer :: lBound,uBound

    lBound = 0.5*(GRID - testcutSize)
    uBound = 0.5*(GRID + testcutSize) - 1
    
    print*, array(n_u,lBound,lBound,lBound)
    print*, 'lBound',lBound , 'uBound',uBound
    
    allocate (temp(n_u,testcutSize,testcutSize,testcutSize))
    temp = array(:,lBound:uBound,lBound:uBound,lBound:uBound)
    deallocate(array)
    allocate(array(n_u,testcutSize,testcutSize,testcutSize))
    array = temp
    
    print *, array(n_u,1,1,1)
    return
  end subroutine cutout
  
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
    integer, parameter      :: n = boxSize - coloc2 ! n=(512-378) = 134
    integer                 :: randomNum, num(n)
    integer,dimension(box,box,box):: randomMatrix
    integer                 :: i, j, k, index, c
    logical,parameter       :: debug = .false.

    num = 0 ; i=0
    
    call init_random_seed()
    do
       call random_number(x)
       randomNum = nint( (boxSize-1) * x ) + 1    
      
       if(any(num.eq.randomNum).or.(randomNum.eq.boxCenter)) cycle
          
       num(i) = randomNum
       i = i+1
       if (i.gt.n) exit   
    end do

    
     call bubblesort(num) ! Remove in real case.Keep for testing

    !!$    ! Activate for debugging:
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
     
  return       
  end subroutine randAlloc
  

  subroutine condition(A)
    implicit none

    integer, parameter :: n = coloc2
    real(kind=8), dimension(n,n) :: A
    real(kind=8), dimension(n) :: x,b
    real(kind=8), dimension(:),allocatable ::work
    
    real(kind=8) :: errnorm, xnorm, rcond, anorm, colsum
    integer :: i, info, lda, ldb, nrhs, j
    integer, dimension(:),allocatable :: ipiv
    integer, dimension(:),allocatable :: iwork
    character, dimension(1) :: norm

   
    allocate(ipiv(n))
    b = 1.d0
    
    ! compute 1-norm needed for condition number

    anorm = 0.d0
    do j=1,n
        colsum = 0.d0
        do i=1,n
            colsum = colsum + abs(a(i,j))
            enddo
        anorm = max(anorm, colsum)
        enddo
    
    nrhs = 1 ! number of right hand sides in b
    lda = n  ! leading dimension of a
    ldb = n  ! leading dimension of b

    call dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)

    ! compute 1-norm of error
    errnorm = 0.d0
    xnorm = 0.d0
    do i=1,n
        errnorm = errnorm + abs(x(i)-b(i))
        xnorm = xnorm + abs(x(i))
    enddo

    ! relative error in 1-norm:
    errnorm = errnorm / xnorm


    ! compute condition number of matrix:
    ! note: uses A returned from dgesv with L,U factors:

    allocate(work(4*n))
    allocate(iwork(n))
    norm = '1'  ! use 1-norm
    call dgecon(norm,n,a,lda,anorm,rcond,work,iwork,info)

    if (info /= 0) then
        print *, "*** Error in dgecon: info = ",info
        endif

    print 201, n, 1.d0/rcond, errnorm
201 format("For n = ",i4," the approx. condition number is ",e10.3,/, &
           " and the relative error in 1-norm is ",e10.3)    

    deallocate(ipiv)
    deallocate(work,iwork)

  end subroutine condition
  




  
  
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

