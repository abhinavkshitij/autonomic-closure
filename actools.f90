!****************************************************************
!                              ACTOOLS
!****************************************************************

!----------------------------------------------------------------
! USE:  1) Compute strain field - S_ij                  : FILTER
!       2) Compute gradient with central differencing   : FILTER
!       3) [SAVE] final results in RUN dir              : SINK 
!
! FORM: module actools
!          contains
!       subroutine energySpectra       [FILTER]
!       subroutine computeS_ij         [FILTER]
!       subroutine gradient            [FILTER]    
!       subroutine productionTerm      [FILTER]
!       subroutine velocityProducts    [FILTER]
!       subroutine createPDF           [FILTER]
! BEHAVIOR: 
!           
!           
!          
!
! STATUS : Reverse indices to gain performance. 
!          Latest test had errors in S_ij. 
! 
!----------------------------------------------------------------

module actools

  use fourier
  use fileio

contains
  

  !****************************************************************
  !                        ENERGY SPECTRA
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE:  Calculates and plots energy spectra 
  !       
  !       
  !  
  !
  ! FORM: subroutine energySpectra(u)
  !       
  !
  ! BEHAVIOR: 
  !           
  !           
  !           
  !          
  !
  ! STATUS : 
  ! 
  !----------------------------------------------------------------
  
  subroutine energySpectra(u)
    implicit none
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:,:,:), intent(in) :: u
    !
    !    ..LOCAL VARS..
    real(8), dimension(:), allocatable :: Ek

    !print*, center
    
  end subroutine energySpectra

  !****************************************************************
  !                        COMPUTE S_IJ                           !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Calculates strain field S_ij 
  !       - calls gradient() to compute grad_x, grad_y, grad_z for 
  !       each u_i.
  !
  ! FORM: subroutine computeSij (u, Sij)
  !       
  !
  ! BEHAVIOR: Each velocity component u_i has three gradient components:
  !           grad_u_i_x, grad_u_i_y, grad_u_i_z, returned by gradient() 
  !           For each i, the three grad components are loaded to A(i)
  !           So A(1) has all three grad components of u1
  !              A(2) has all three grad components of u2
  !              A(3) has all three grad components of u3
  !          
  !
  ! STATUS : Reverse indices to improve speed. 
  !          Need to time this section to track  performance.
  !----------------------------------------------------------------
  
  subroutine computeSij(u, Sij)
    implicit none
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:,:,:),  intent(in) :: u
    real(8), dimension(:,:,:,:,:),intent(out) :: Sij 
    !
    !    ..WORK ARRAYS..
    real(8), allocatable, dimension(:,:,:,:,:)  :: A
    real(8), allocatable, dimension(:,:,:,:)    :: grad_u
    !
    !    ..INDICES..
    integer :: i,j !DO NOT REMOVE THIS LINE. i IS PASSED INTO gradient() AND ITS VALUE IS MODIFIED.

    allocate (A(3,3,i_GRID,j_GRID,k_GRID))
    allocate (grad_u(3,i_GRID,j_GRID,k_GRID))

    !
    ! COMPUTE VELOCITY GRADIENTS:
    write(*,'(a32)',ADVANCE='NO'), adjustl('        Compute velocity gradients:')
    call cpu_time(tic)
    do i=1,3
       call gradient(u(i,:,:,:),grad_u)
       A(i,:,:,:,:) = grad_u        ! A(1),A(2),A(3) -> grad_u_x, grad_u_y, grad_u_z
    end do
    call cpu_time(toc)
    write(*,*), toc-tic
    deallocate(grad_u)

    !
    ! COMPUTE S_ij:
    write(*,'(a16)',ADVANCE='NO'), adjustl('      Compute Sij:')
    call cpu_time(tic)
    do j=1,3
       do i=1,3
          Sij(i,j,:,:,:) = 0.5d0 * (A(i,j,:,:,:) + A(j,i,:,:,:))
       end do
    end do
    call cpu_time(toc)
    write(*,*), toc-tic

    deallocate (A)

    return
  end subroutine computeSij


  
  !****************************************************************
  !                              GRADIENT
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE:  Calculates gradients using central order differencing
  !       Returns three gradient vector components: 
  !                 grad_f(1), grad_f(2), grad_f(3)
  !  
  !
  ! FORM: subroutine gradient(f, grad_f)
  !       
  !
  ! BEHAVIOR: Uses ghost cells
  !           
  !           
  !           
  !          
  !
  ! STATUS : Extend this capability to all possible BC
  ! 
  !----------------------------------------------------------------
  
  subroutine gradient(f, grad_f)
    implicit none
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension (:,:,:), intent(in) ::  f
    real(8), dimension (:,:,:,:),intent(out) :: grad_f
    !
    !    ..LOCAL ARRAY..
    real(8), allocatable,dimension (:,:,:) :: u
    !
    !    ..LOCAL SCALARS..
    real(8) :: dx_inv

    !
    ! FIND 1/dx FOR CENTRAL DIFFERENCING
    dx_inv = 1.d0/(2.d0*dx)

    !
    ! LOAD f ON WORK ARRAY u WITH EXTENDED DOMAIN TO CONTAIN GHOST CELLS
    allocate(u (0:i_GRID+1, 0:j_GRID+1, 0:k_GRID+1))
    u(1:i_GRID,1:j_GRID,1:k_GRID) = f

    ! CREATE GHOST CELLS: (BASED ON PERIODIC BC)
    ! NEED AT 6 FACES:
    u(0,:,:) = u(i_GRID,:,:) ; u(i_GRID+1,:,:) = u(1,:,:)
    u(:,0,:) = u(:,j_GRID,:) ; u(:,j_GRID+1,:) = u(:,1,:)
    u(:,:,0) = u(:,:,k_GRID) ; u(:,:,k_GRID+1) = u(:,:,1)

    !
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
  end subroutine gradient


  
  !****************************************************************
  !                        PRODUCTION TERM                        !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Calculates production term P = tau_ij * S_ij
  !       
  !       
  !
  ! FORM: subroutine productionTerm(P, tau_ij, S_ij)
  !       
  !
  ! BEHAVIOR: 
  !           
  !           
  !           
  !           
  !          
  !          
  !
  ! STATUS : 
  !          
  !----------------------------------------------------------------

  subroutine productionTerm(P, tau_ij, S_ij)
    implicit none
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:,:,:), intent(in) :: tau_ij
    real(8), dimension(:,:,:,:,:), intent(in) :: S_ij
    real(8), dimension(:,:,:), intent(out) :: P
    !
    !   .. LOCAL VARS..
    integer :: count 
    
    count = 0
    P = 0
    do j = 1,3
       do i = 1,3
          if(i.ge.j) then
             count = count + 1
             P(:,:,:) = P(:,:,:) + tau_ij(count,:,:,:) * S_ij (i,j,:,:,:) 
          end if
       end do
    end do
  end subroutine productionTerm

  
  !****************************************************************
  !                        VELOCITY PRODUCTS                      !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Calculates velocity products uiuj
  !       
  !       
  !
  ! FORM: subroutine productionTerm(P, tau_ij, S_ij)
  !       
  !
  ! BEHAVIOR: 
  !           
  !           
  !           
  !           
  !          
  !          
  !
  ! STATUS : 
  !          
  !----------------------------------------------------------------

  subroutine velocityProducts(uiuj,ui)
    implicit none
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:,:,:), intent(in) :: ui
    real(8), dimension(:,:,:,:), intent(out) :: uiuj
    !
    !   .. LOCAL VARS..
    integer :: count = 0

    do j = 1,3
       do i = 1,3
          if(i.ge.j) then
             count = count + 1
             uiuj(count,:,:,:) = ui(i,:,:,:) * ui(j,:,:,:)
          end if
       end do
    end do
  end subroutine velocityProducts

  
  !****************************************************************
  !                        CREATE PDF                             !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Samples a field and computes the probability density function PDF
  !       
  !       
  !
  ! FORM: subroutine createPDF(field, [plot], [fieldname])
  !       
  !
  ! BEHAVIOR: 
  !           
  !           
  !           
  !           
  !          
  !          
  !
  ! STATUS :  
  !
  ! CHECK GAUSSIAN BY MEASURING THE THIRD ORDER MOMENT:
  !            print*, sum((field-avg)**3)/size(field)
  !
  !          
  !----------------------------------------------------------------

  subroutine createPDF(field, plotOption, fieldname)
    implicit none
    !
    !    ..ARRAY ARGUMENT..
    real(8), dimension(:,:,:), intent(in) :: field
    !
    !    ..PLOT ARGUMENT..
    character(*), optional, intent(in) :: plotOption
    character(*), optional, intent(in) :: fieldname

    !
    !    ..LOCAL ARRAYS..
    real(8), dimension(:), allocatable :: x
    real(8), dimension(:), allocatable :: pdf
    !
    !    ..LOCAL VARS..
    real(8) :: delta, avg, SD


    allocate (x(samples))
    allocate (pdf(samples))
    
    !
    ! LINSPACE-LIKE FUNCTION - IND. DOMAIN
    delta = (maxval(field) - minval(field)) / (samples-1)
    forall(i=1:samples)  x(i) = minval(field) + (delta * (i-1))

    !
    ! CALCULATE FIELD MEAN AND STD DEV
    avg = mean(field)
    SD = stdev(field)

    !
    ! PDF
    pdf = exp(-((x-avg)**2 / (2*SD**2))) / (sqrt (2*pi)*SD) 

    !
    ! PLOT:
    if(plotOption.eq.'plot'.and.present(fieldname)) then
       call plotPDF(x,pdf,fieldname,avg,SD)
    else
       print*, 'Need fieldname argument to name the file'
    end if
  end subroutine createPDF
  

  !****************************************************************
  !                        CROSS VALIDATION                       !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Compute cross validation error at test scale
  !       
  !       
  !
  ! FORM: subroutine createPDF(field, [plot], [fieldname])
  !       
  !
  ! BEHAVIOR: 
  !           
  !           
  !           
  !           
  !          
  !          
  !
  ! STATUS :  
  !
  ! CHECK GAUSSIAN BY MEASURING THE THIRD ORDER MOMENT:
  !            print*, sum((field-avg)**3)/size(field)
  !
  !          
  !----------------------------------------------------------------
  
  subroutine trainingError(T_ijOpt, T_ij, error, plotOption,fID)
    implicit none
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:,:,:), intent(in) :: T_ijOpt
    real(8), dimension(:,:,:,:), intent(in) :: T_ij
    !
    !    ..SCALAR ARGUMENTS..
    real(8), intent(out) :: error
    character(*), optional, intent(in) :: plotOption
    integer,optional, intent(in) :: fID   
    !
    !   ..LOCAL VARS..
    integer :: i
    real(8) :: error_i(6)


    do i = 1,6
       error_i(i) =  norm  (T_ijOpt(i, bigHalf(i_GRID)-3 : bigHalf(i_GRID)+3 : 3, &
                                       bigHalf(j_GRID)-3 : bigHalf(j_GRID)+3 : 3, &
                                       bigHalf(k_GRID)-3 : bigHalf(k_GRID)+3 : 3) &

                          - T_ij   (i, bigHalf(i_GRID)-3 : bigHalf(i_GRID)+3 : 3, &
                                       bigHalf(j_GRID)-3 : bigHalf(j_GRID)+3 : 3, &
                                       bigHalf(k_GRID)-3 : bigHalf(k_GRID)+3 : 3) )  
    end do

    error_i = error_i / 27
    error = maxval(error_i)

    !
    !    PLOT RESULTS:
    if (present(plotOption).and.plotOption.eq.'plot') then
       write(fID,"( 8(ES16.7,',') )"), lambda, error_i, error    
    end if

  end subroutine trainingError

end module actools
