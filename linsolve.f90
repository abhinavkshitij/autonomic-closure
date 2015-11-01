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

  ! Cutout atributes:
  integer,parameter :: lBound = 0.5*(GRID - testcutSize)
  integer,parameter :: uBound = 0.5*(GRID + testcutSize) - 1

 
contains

subroutine cutout(array,dim)
    implicit none
   
    real(kind=8), allocatable, dimension(:,:,:,:),intent(inout) :: array
    integer, intent(in) :: dim
    
    real(kind=8), allocatable, dimension(:,:,:,:):: temp
    
    print*, 'lBound',lBound , 'uBound',uBound
    print*, 'u_t(1,1,1)', array(1,lBound,lBound,lBound)
    print*, 'Shape(u_t):', shape(array)
    
    print*, "Begin cutout ..."
    
    allocate (temp(dim,testcutSize,testcutSize,testcutSize))
    temp = array(:,lBound:uBound,lBound:uBound,lBound:uBound)
    deallocate(array)
    allocate(array(dim,testcutSize,testcutSize,testcutSize))
    array = temp

    print*, ''
    print*, "Done cutout ..."
    print*, 'u_t(1,1,1)',array(1,1,1,1)
    print*, 'Shape(u_t):',shape(array)
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





  subroutine synStress(u,n_u,n_uu)
    implicit none

    real(kind=8),allocatable,dimension(:,:,:,:),intent(inout) :: u
    real(kind=8),allocatable,dimension(:,:,:,:) :: uu
    real(kind=8),dimension(coloc2,coloc2) :: A = 0

    integer :: n_u
    integer :: n_uu

    integer :: randMask(boxSize-coloc2) !512-378=134; 512-270=240
    integer :: cell_flag(testcutSize,testcutSize,testcutSize)

    integer :: i_test,    j_test,    k_test
    integer :: i_box,     j_box,     k_box
    integer :: i_stencil, j_stencil, k_stencil
    
    integer :: rand_count,lim
    integer :: row_index, col_index
    integer :: u_comp,    uu_comp

    integer,dimension(4) :: debug=(/0,1,0,1/)

    integer :: i,j,k,p,DIM
  
  ! call cutout(u,n_u)
   print*, 'Transfer to synStress ... check first element'
   print *, 'u(1,1,1) ',u(1,1,1,1)
   print *, 'shape u cutout: ',shape(u)

   allocate(uu(n_uu,testcutSize,testcutSize,testcutSize))
   print*,'shape uu cutout:',shape(uu)

   print*, 'Compute velocity products:'
   k=0
   do j=1,n_u
      do i=1,n_u
         if (i.ge.j) then
            k=k+1
            uu(k,:,:,:) = u(i,:,:,:)*u(j,:,:,:)
         end if       
      end do
   end do

   print*,''
   print*,'boxCenter index:',boxCenter
   print*,'smallHalf,bighalf:',smallHalf,bigHalf



! 2*(5+1)-1 for 8 = 11 ; 5 for box,1 for stencil
   ! As a modification, it should be independent of stride=2 , as well

   if(debug(3).eq.1)then
      lim=testUpper
      print*, 'Check for the last element..(43,43,43)'
      print*, u(1,uBound,uBound,uBound)
   else
      lim=testLower
      print*, 'Check for the first element... (11,11,11)'
      print*, u(1,1,1,1)
   end if

!!$do k_test = testLower, testUpper, stride 
!!$   do j_test = testLower, testUpper, stride
!!$      do i_test = testLower, testUpper, stride ! i_test = 11,43,2


 do k_test = lim,lim,stride
   do j_test = lim,lim,stride
      do i_test = lim,lim,stride
         
         cell_flag  = 0   
         rand_count = 0 
         row_index  = 0 ! Reset pointer to the first position after moving
                        ! to the next cell.

         call randAlloc(randMask)
          do k_box = k_test-boxLower, k_test+boxUpper,stride
            do j_box = j_test-boxLower, j_test+boxUpper,stride
               do i_box = i_test-boxLower,i_test+boxUpper,stride ! i_box = 3,49,2
                  
                  rand_count = rand_count+1
                  if (any(randMask.eq.rand_count)) cycle   ! Skip if the point is listed in randMask
                    
                  ! cell_flag(i_box,j_box,k_box) = 1
                  
                  col_index = 0 ! Reset pointer to the first column after each stencil operation
                  row_index = row_index + 1

                   do k_stencil = k_box-stride,k_box+stride,stride
                     do j_stencil = j_box-stride,j_box+stride,stride
                        do i_stencil = i_box-stride,i_box+stride,stride
                           
                          cell_flag(i_stencil,j_stencil,k_stencil) = 2 ! Flag in stencil-->check:first cell
                                                    
                           do u_comp = 1,n_u ! 1 to 3
                              col_index = col_index+1
                              A(row_index,col_index)=u(u_comp,i_stencil,j_stencil,k_stencil)                           
                           end do
                           do uu_comp = 1,n_uu ! 1 to 6
                              col_index = col_index+1
                              A(row_index,col_index)=uu(uu_comp,i_stencil,j_stencil,k_stencil)
                           end do

                          
                        
                        end do
                     end do
                  end do !stencil
                  
                 ! B            
                  
               end do
            end do
         end do !box
         
         p=p+1
         print*,'Begin Inverse operation:'
         if(p.eq.1.or.p.eq.100.or.p.eq.200) then
            print*,'Check A(5,5):',A(5,5)
            call inverse(A)
            print*,'Check A(5,5):',A(5,5)
         end if
         
         ! Compute SGS stress:
         

         
      end do
   end do
end do ! test

if(debug(4).eq.1) then
   
print*,''
print*, 'Exclusion list:'
print*, randMask
print*,''
print *, 'Number of random cells sampled:',rand_count

print*, 'Check for first cell flag:',cell_flag(1,1,1)
print*, 'Check for last cell:',cell_flag(testcutSize,testcutSize,testcutSize)
        
print*,''

print*,'Last limit of test array:',i_test-stride 
print*,'Last limit of bounding box:',i_box-stride ! This one will always count till the end.
print*,'Last limit of stencil:',i_stencil-stride
print*,''


!print*,uu(6,testcutSize,testcutSize,testcutSize),A(coloc2,coloc2) !--> check for the last cell 
!print*,uu(6,testcutSize,testcutSize-stride,testcutSize),A(coloc2,coloc2) !--> check for the last cell

end if

deallocate(uu)

!!$open(1,file="eig.dat")
!!$do j=1,243
!!$   write(1,202) A(j,:)
!!$end do     
!!$202 format(243f16.4)                                          
!!$close(1)                       

  
    return
  end subroutine synStress




  
  subroutine inverse(A_arr)
    implicit none

    integer, parameter :: n = coloc2
    real(kind=8), dimension(n,n),intent(in) :: A_arr
    
    real(kind=8), dimension(:,:),allocatable :: a
    real(kind=8), dimension(:),allocatable :: x,b
    real(kind=8), dimension(:),allocatable :: work
    
    real(kind=8) :: errnorm, xnorm, rcond, anorm, colsum
    integer :: i, info, lda, ldb, nrhs, j
    integer, dimension(:),allocatable :: ipiv
    integer, dimension(:),allocatable :: iwork
    character, dimension(1) :: norm


    allocate(a(n,n))
    allocate(b(n),x(n),ipiv(n))
    
    a = A_arr
    b = 1.d0 ; x = 1.d0  !! change this for real b,x
    
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

    print*,'a(5,5)',a(5,5)
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


    deallocate(a,x,b,ipiv,work,iwork)
        

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

