module linsolve
  contains
subroutine init_random_seed()
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed
  
    call RANDOM_seed(size = n)
    allocate(seed(n))
  
    call SYSTEM_clock(COUNT=clock)
  
    seed = clock + 37 * (/ (i-1, i=1,n) /)
    call RANDOM_seed(PUT = seed)
  
    print *, "Using random seed = ", seed
    print *, " "

    deallocate(seed)
    return
  end subroutine init_random_seed


  subroutine randAlloc()

    ! Status:Begin testing for inverse random
    ! Notes: Need 378 random numbers from 1 to 512
    !      : Check for repetition
    implicit none
    integer, parameter     :: nmax = 1000
    integer, parameter     :: rand_n = 378
    real(kind=8) :: x(rand_n)    
    integer      :: num(rand_n)
    integer,dimension(8,8,8)::randomMatrix
    integer      :: i,j,k,count,index

    
    call init_random_seed()   
    call random_number(x(1:rand_n))
    num = nint(511*x)+1
    print *, "         x   -initial      "
    print *, num
    

     ! Sort random numbers in ascending order
     call unique_sort(num,rand_n)
     print *, "         x -after unique sort        "
     print *, num

    count=0;index=1;randomMatrix=0
    do k=1,8
       do j=1,8
          do i=1,8
             count = count+1
             if (count.eq.num(index))then
                randomMatrix(i,j,k) = num(index)
                index = index+1
             end if  
          enddo
       enddo
    enddo
    print*,'randomMatrix'
    do k=1,8
    do i=1,8
       print*,randomMatrix(i,:,k)
    end do
    print*,''
    end do
  return       
  end subroutine randAlloc
  

  subroutine unique_sort(array,n)
    !! This subroutine will sort and ensure that all elements in the random array are unique.
    implicit none

    
    integer,intent(in) :: n 
    integer, dimension(1:n),intent(inout) :: array
    integer, allocatable,dimension(:) :: bank,list
    integer :: i,j,k,p,n_list,bankIndex
    allocate (bank(n))
    bank = 0
   

    print *, "         x-after sort         "
    call bubblesort(array,n)
    print *, array

   
   ! Create a list of non-selected numbers -- bank
    bankIndex = 0
    do j = 1,n   
       if ( any(array.eq.j)) cycle
       bankIndex = bankIndex+1
       bank(bankIndex) = j             
    end do

    n_list = count(bank.gt.0)
    
       print*, "length of the non-selected elements in the bank:",bankIndex,n_list
       print*, "Bank elements:"
       print *, bank
      

       allocate(list(1:n_list))
       do i =1,n_list
          list(i) = bank(i)
       end do
       deallocate(bank)
       
       print*, list
       
    k=1    
    do p=2,n 
       if ((array(p).eq. array(p+1)).or.(array(p).eq.array(p-1)))then
          
          replace:do i=2,n-1
             if ((array(i).eq. array(i+1)).or.(array(i).eq.array(i-1)))then
                array(i) = list(n_list-k)
                k=k+1
             end if
          end do replace          
       end if
    end do
    
    call bubblesort(array,n)

   
    return
  end subroutine unique_sort

  subroutine bubblesort(array,n)
    implicit none
    integer, dimension(:)::array
    integer :: i,j,n,temp


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

