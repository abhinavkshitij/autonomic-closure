!****************************************************************
!                              ACTOOLS
!****************************************************************

!----------------------------------------------------------------
! USE:  1) Compute strain field - S_ij                  : FILTER
!       2) Compute gradient with central differencing   : FILTER
!       3) [SAVE] final results in RESULTS dir          : SINK 
!
! FORM: module actools
!          contains
!       subroutine energySpectra       [FILTER]
!       subroutine computeS_ij         [FILTER]
!       subroutine gradient            [FILTER]    
!       subroutine productionTerm      [FILTER]
!       subroutine dissipationRate     
!       function   turbulentKE         
!       subroutine secondOrderProducts [FILTER]
!       subroutine extendDomain        [FILTER]
!       subroutine createPDF           [FILTER]
!       subroutine trainingError
!       
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
  implicit none

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
  ! STASH:
  !  ***
  !      write(*,'(a32)',ADVANCE='NO'), adjustl('        Compute velocity gradients:')
  !      call cpu_time(tic)
  !      ...
  !      call cpu_time(toc)
  !      write(*,*), toc-tic
  !      
  !      write(*,'(a16)',ADVANCE='NO'), adjustl('      Compute Sij:')
  !      call cpu_time(tic)
  !      ...
  !      call cpu_time(toc)
  !      write(*,*), toc-tic
  !     
  !
  !----------------------------------------------------------------
  
  subroutine computeSij(u, Sij)
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:,:,:),  intent(in) :: u
    real(8), dimension(:,:,:,:),intent(out) :: Sij 
    !
    !    ..WORK ARRAYS..
    real(8), allocatable, dimension(:,:,:,:,:)  :: A
    real(8), allocatable, dimension(:,:,:,:)    :: grad_u
    !
    !    ..INDICES..
    integer :: i,j,k !DO NOT REMOVE THIS LINE. i IS PASSED INTO gradient() AND ITS VALUE IS MODIFIED.

    allocate (A(3,3,i_GRID,j_GRID,k_GRID))
    allocate (grad_u(3,i_GRID,j_GRID,k_GRID))

    !
    ! COMPUTE VELOCITY GRADIENTS:
    do i=1,3
       call gradient(u(i,:,:,:),grad_u)
       A(i,:,:,:,:) = grad_u        ! A(1),A(2),A(3) -> grad_u_x, grad_u_y, grad_u_z
    end do
    deallocate(grad_u)

    !
    ! COMPUTE S_ij:
    k = 0
    do j=1,3
       do i=1,3

          if (i.ge.j) then
             k = k + 1
             Sij(k,:,:,:) = 0.5d0 * (A(i,j,:,:,:) + A(j,i,:,:,:))
          end if

       end do
    end do
    deallocate (A)

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
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:,:,:), intent(in) :: tau_ij
    real(8), dimension(:,:,:,:), intent(in) :: S_ij
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
             if (i.eq.j) then
                P(:,:,:) = P(:,:,:) + tau_ij(count,:,:,:) * S_ij (count,:,:,:) 
             else
                P(:,:,:) = P(:,:,:) + 2.d0 * (tau_ij(count,:,:,:) * S_ij(count,:,:,:))
             end if
          end if
       end do
    end do

  end subroutine productionTerm

  
  !****************************************************************
  !                        DISSIPATION RATE                        !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Calculates dissipation rate: epsilon = nu * S_ij * S_ij
  !       
  !       
  !
  ! FORM: subroutine dissipationRate(S_ij, epsilon)
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
  !       ##     print*, 'SijSij: ',mean(temp), temp (10,25,24), temp(3,4,5)
  !          
  !----------------------------------------------------------------


  subroutine dissipationRate(S_ij, epsilon)
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:,:,:), intent(in) :: S_ij
    real(8), intent(out) :: epsilon
    !
    !   .. LOCAL VARS..
    real(8), dimension(:,:,:), allocatable :: temp
    integer :: count 
    
    allocate (temp(i_GRID, j_GRID, k_GRID))

    
    count = 0
    do j = 1,3
       do i = 1,3
          if(i.ge.j) then
             count = count + 1
             if (i.eq.j) then
                temp(:,:,:) = temp(:,:,:) + S_ij(count,:,:,:) * S_ij (count,:,:,:) 
             else
                temp(:,:,:) = temp(:,:,:) + 2.d0 * (S_ij(count,:,:,:) * S_ij(count,:,:,:))
             end if
          end if
       end do
    end do

!## Print pointwise values
    epsilon = 2.d0 * nu * mean(temp)

  end subroutine dissipationRate


  
  !****************************************************************
  !                   TURBULENT KINETIC ENERGY                    !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Calculates turbulent KE = 0.5 * (u_1'^2 + u_2'^2 + u_3'^3)
  !       
  !       
  !
  ! FORM: subroutine turbulentKE(u, TKE)
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



  function turbulentKE(u)
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:,:,:), intent(in) :: u
    real(8)  :: turbulentKE
           
    turbulentKE = 0.5d0 * mean( (u(1,:,:,:) - mean(u(1,:,:,:)))**2 &
                               +(u(2,:,:,:) - mean(u(2,:,:,:)))**2 &
                               +(u(3,:,:,:) - mean(u(3,:,:,:)))**2 ) 

  end function turbulentKE



  
  !****************************************************************
  !                        VELOCITY PRODUCTS                      !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Calculates velocity products uiuj
  !       
  !       
  !
  ! FORM: subroutine secondOrderProducts(uiuj, ui)
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

  subroutine secondOrderProducts(uiuj,ui)
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(1:,extLower:,extLower:,extLower:), intent(in) :: ui
    real(8), dimension(1:,extLower:,extLower:,extLower:), intent(out) :: uiuj
    !
    !   .. LOCAL VARS..
    integer :: count 

    count = 0
    do j = 1, n_u
       do i = 1, n_u
          if(i.ge.j) then
             count = count + 1
             uiuj(count,:,:,:) = ui(i,:,:,:) * ui(j,:,:,:)
          end if
       end do
    end do
  end subroutine secondOrderProducts


  
  !****************************************************************
  !                        EXTEND DOMAIN
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Use ghost cells to extend domain 
  !       
  !       
  !
  ! FORM: subroutine extendDomain()
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

  subroutine extendDomain(array, varName)
    !
    !   .. ARRAY ARGUMENTS ..
    real(8), dimension (:,:,:,:), allocatable, intent(inout) :: array
    !
    !   .. SCALAR ARGUMENTS..
    character(*), intent(in), optional :: varName
    !
    !   .. LOCAL VARS..
    real(8), dimension(:,:,:,:), allocatable :: temp
    integer :: n_extendedLayers
    integer :: tempDim(4)
    integer :: i,j,k

    tempDim = shape(array)
    allocate (temp (tempDim(1), tempDim(2), tempDim(3), tempDim(4)))

    ! REALLOCATE ARRAY SIZE
    temp = array
    deallocate(array)
    allocate(array(tempDim(1), extLower:extUpper, extLower:extUpper, extLower:extUpper))
    array(:, 1:i_GRID, 1:j_GRID, 1:k_GRID) = temp(:, 1:i_GRID, 1:j_GRID, 1:k_GRID)
    
    ! USE PERIODIC CONDITIONS ON GHOST CELLS:
    if (scheme.eq.'global') then
       n_extendedLayers = 2
    else
        n_extendedLayers = boxLower + 2
    end if
    
    do k = 1, n_extendedLayers
       do j = 1, n_extendedLayers
          do i = 1, n_extendedLayers
             array(:, 1-i, :, :) = array(:, i_GRID-i+1, :, :) ! x-plane 
             array(:, :, 1-j, :) = array(:, :, j_GRID-j+1, :) ! y-plane 
             array(:, :, :, 1-k) = array(:, :, :, k_GRID-k+1) ! z-plane 

             array(:, i_GRID+i, :, :) = array(:, i, :, :)     ! x-plane 
             array(:, :, j_GRID+j, :) = array(:, :, j, :)     ! y-plane 
             array(:, :, :, k_GRID+k) = array(:, :, :, k)     ! z-plane 
          end do
       end do
    end do

  end subroutine extendDomain

  
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


    allocate (x(n_bins))
    allocate (pdf(n_bins))
    
    !
    ! LINSPACE-LIKE FUNCTION - IND. DOMAIN
    delta = (maxval(field) - minval(field)) / (n_bins-1)
    forall(i=1:n_bins)  x(i) = minval(field) + (delta * (i-1))

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
  ! FORM: subroutine trainingError(T_ijOpt, T_ij, error, plotOption,fID)
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
  ! STASH:
  !
  ! Ensemble mean of L2 norm error < ||.||2 > N_cr
  !  error_ij = error_ij / N_cr**3
  !
  ! Non-dimensionalize error_ij by ||T_ij||_2/(nx*ny) -> Mean value at z-midplane [using numerical integration]
  !  do ij = 1,6
  !      error_ij(ij) = error_ij(ij) / sqrt(sum(T_ij(ij,:,:,bigHalf(k_GRID))**2)/(i_GRID*j_GRID))         
  !  end do
  !
  !   
  ! USING NORM: [OLD METHOD]
  !      error_ij(ij) = norm((T_ij_F_cr(ij,:,:,:) - T_ij_cr(ij,:,:,:))) / norm(T_ij_cr(ij,:,:,:))
  !      error_ij = error_ij/(N_cr)
  !
  ! PRINT T_ij_cr and T_ij_F_cr
  ! !      print*, T_ij_F_cr(ij,6,6,6), T_ij_cr(ij,6,6,6)      
  !----------------------------------------------------------------
  
  subroutine trainingError(T_ijOpt, T_ij, error, plotOption,fID)
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:,:,:), intent(in) :: T_ijOpt
    real(8), dimension(:,:,:,:), intent(in) :: T_ij
    !
    !    ..SCALAR ARGUMENTS..
    real(8), intent(out) :: error
    real(8), dimension (:,:,:,:), allocatable :: T_ij_F_cr, T_ij_cr
    character(*), optional, intent(in) :: plotOption
    integer, optional, intent(in) :: fID   
    !
    !   ..LOCAL VARS..
    integer :: ij, i,j,k
    real(8) :: error_ij(6)

    allocate (T_ij_F_cr (6, N_cr, N_cr, N_cr))
    allocate (T_ij_cr   (6, N_cr, N_cr, N_cr))

    !
    ! MARK CROSS-VALIDATION POINTS:
     T_ij_F_cr = T_ijOpt(:,  bigHalf(i_GRID)-3*smallHalf(N_cr) : bigHalf(i_GRID)+3*smallHalf(N_cr) : 3, &
                             bigHalf(j_GRID)-3*smallHalf(N_cr) : bigHalf(j_GRID)+3*smallHalf(N_cr) : 3, &
                             bigHalf(k_GRID)-3*smallHalf(N_cr) : bigHalf(k_GRID)+3*smallHalf(N_cr) : 3) 

     T_ij_cr = T_ij(:,       bigHalf(i_GRID)-3*smallHalf(N_cr) : bigHalf(i_GRID)+3*smallHalf(N_cr) : 3, &
                             bigHalf(j_GRID)-3*smallHalf(N_cr) : bigHalf(j_GRID)+3*smallHalf(N_cr) : 3, &
                             bigHalf(k_GRID)-3*smallHalf(N_cr) : bigHalf(k_GRID)+3*smallHalf(N_cr) : 3) 

     !
     ! TRAINING ERROR (_ij):
     do ij = 1,6
        error_ij(ij) = sqrt(sum((T_ij_F_cr(ij,:,:,:) - T_ij_cr(ij,:,:,:)) ** 2) / sum(T_ij_cr(ij,:,:,:) ** 2))
     end do

     ! ERROR (_rms):
     error = sqrt(sum(error_ij**2)/6.d0)
     
     open(77, file = trim(RES_PATH)//'TrainingError27.csv', status = 'replace')
     do j=1,N_cr
     do i=1,N_cr
!     write(77,"( 2(ES16.7,','),ES16.7 )"), T_ij_F_cr(1,i,j,6), T_ij_cr(1,i,j,6), error_ij(1)
     end do
     end do
     close(77)

     ! PLOT RESULTS:
     if (present(plotOption).and.plotOption.eq.'plot') then
        write(fID,"( ES8.1,',',6(ES16.7,','),ES16.7 )"), lambda, error_ij, error    
     end if

  end subroutine trainingError


end module actools
