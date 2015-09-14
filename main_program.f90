module simple_ops
contains
  subroutine sincos(a)
    implicit none
    real(kind=8),intent(in), dimension(:,:)::a
    real(kind=8),allocatable, dimension (:,:)::b
    
    b = sin(a)+cos(a)
    call printdata(b)
    deallocate(b)
    return
  end subroutine sincos

  subroutine printdata(a)
    ! Prints a matrix in a neat way.
    ! CAUTION: Change the entry in format statement for different matrix sizes.
    implicit none
    real(kind=8),intent(inout),dimension(:,:) ::a
    integer :: dims(2)
    integer :: i

    dims = shape(a)
    do i=1,dims(1)
       write(*,20) a(i,:) ! Print the ith row
    end do
20  format (4(f10.4))
    
  return  
  end subroutine printdata
end module simple_ops

  
program array_module
  use simple_ops

  implicit none
  integer, parameter :: N=4
  real(kind=8), dimension(N,N) :: x,y
  integer ::i,j

  ! Initialize x:
  do j=1,N
     do i=1,N
        x(i,j) = 2*i+j
     end do
  end do

  print*, "x="
  print*, ""
  call printdata(x)
  print*, "sin(x)+cos(x) ="
  print*, ""
  call sincos(x)
  
  end program array_module
  

