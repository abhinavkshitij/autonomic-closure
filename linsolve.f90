module linsolve
  integer,parameter :: box = 8
  integer,parameter :: boxCenter = (0.5*box+1)*(box**2) + (0.5*box)*(box)+ (0.5*box+1)
  contains
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


  subroutine randAlloc()

    ! Status:Begin testing for inverse random
    ! Notes: Need 378 random numbers from 1 to 512
    !      : Check for repetition
    implicit none
    integer, parameter     :: nmax = 1000
    integer, parameter     :: rand_n = box**3-378
    real(kind=8) :: x    
    integer      :: num(rand_n)
    integer,dimension(box,box,box)::randomMatrix
    integer      :: i,j,k,index,randomNum,c

    c=0;num = 0 ; i=0
    call init_random_seed()
    do
       call random_number(x)
       randomNum = nint(511*x)+1     !! 511*x + 1
      
       if(any(num.eq.randomNum).or.(randomNum.eq.boxCenter)) then
          c = c+1
          cycle
       end if
       

     
       num(i) = randomNum
       i = i+1
       if (i.gt.rand_n)exit
       
    end do

     
    

    print *, "         x   -initial      "
    print *, num
    
    call bubblesort(num)
     print *, "         x -after unique sort        "
     print *, num

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
 
    print*,count(randomMatrix.gt.0.)
    
  return       
  end subroutine randAlloc
  

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

