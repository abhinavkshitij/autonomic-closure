module actools

use fourier
use fileio

contains

!****************************************************************
!                             S_ij                              !
!****************************************************************
  
  subroutine computeSij(u,Sij)
    implicit none

    real(8), dimension(:,:,:,:),  intent(in) :: u
    real(8), dimension(:,:,:,:,:),intent(out) :: Sij 
    real(8), allocatable, dimension(:,:,:,:,:) :: A
    real(8), allocatable, dimension(:,:,:,:) :: grad_u

    integer :: i,j !DO NOT REMOVE THIS LINE. i IS PASSED INTO gradient() AND ITS VALUE IS MODIFIED.

    allocate (A(3,3,i_GRID,j_GRID,k_GRID))
    allocate (grad_u(3,i_GRID,j_GRID,k_GRID))

    
    write(*,'(a32)',ADVANCE='NO'), 'Compute velocity gradients:'
    call cpu_time(tic)
    do i=1,3
       call gradient(u(i,:,:,:),grad_u)
       A(i,:,:,:,:) = grad_u        ! A(1),A(2),A(3) -> grad_u_x, grad_u_y, grad_u_z
    end do
    call cpu_time(toc)
    print*, toc-tic
    deallocate(grad_u)

    write(*,'(a16)',ADVANCE='NO'), "Compute Sij:"
    call cpu_time(tic)
    ! COMPUTE S_ij:
    do j=1,3
       do i=1,3
          Sij(i,j,:,:,:) = 0.5d0 * (A(i,j,:,:,:) + A(j,i,:,:,:))
       end do
    end do
    call cpu_time(toc)
    print*, toc-tic

    deallocate (A)

    return
  end subroutine computeSij
  
  
  subroutine gradient(f, grad_f)
    implicit none

    real(8), dimension (:,:,:), intent(in) ::  f
    real(8), dimension (:,:,:,:),intent(out) :: grad_f

    real(8), allocatable,dimension (:,:,:) :: u
    real(8) :: dx = 2.d0*pi/dble(i_GRID) !Only for JHU data. Change for others.
    real(8) :: dx_inv

    dx_inv = 1.d0/(2.d0*dx)

    allocate(u (0:i_GRID+1, 0:j_GRID+1, 0:k_GRID+1))
    u(1:i_GRID,1:j_GRID,1:k_GRID) = f

    ! CREATE GHOST CELLS: (BASED ON PERIODIC BC)
    ! NEED AT 6 FACES:
    u(0,:,:) = u(i_GRID,:,:) ; u(i_GRID+1,:,:) = u(1,:,:)
    u(:,0,:) = u(:,j_GRID,:) ; u(:,j_GRID+1,:) = u(:,1,:)
    u(:,:,0) = u(:,:,k_GRID) ; u(:,:,k_GRID+1) = u(:,:,1)


    ! APPLY 3D CENTRAL DIFFERENCING SCHEME AT INTERIOR POINTS:
    do k=1,k_GRID
    do j=1,j_GRID
    do i=1,i_GRID
       grad_f(1,i,j,k) = dx_inv * (  u(i+1,j,k) - u(i-1,j,k) )
       grad_f(2,i,j,k) = dx_inv * (  u(i,j+1,k) - u(i,j-1,k) )
       grad_f(3,i,j,k) = dx_inv * (  u(i,j,k+1) - u(i,j,k-1) )
    end do
    end do
    end do

    return
  end subroutine gradient

end module actools
