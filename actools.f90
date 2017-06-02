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
!       subroutine vorticity           [FILTER]
!       subroutine dissipationRate     
!       function   turbulentKE         
!       subroutine secondOrderProducts [FILTER]
!       subroutine extendDomain        [FILTER]
!       subroutine createPDF           [FILTER]
!       subroutine trainingError
!       subroutine bandPassFilter      [FILTER]
!       subroutine rotateX             [FILTER]
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
    real(8), dimension(1:,1:,1:,zLower:),intent(out) :: Sij 
    !
    !    ..WORK ARRAYS..
    real(8), allocatable, dimension(:,:,:,:,:)  :: A
    real(8), allocatable, dimension(:,:,:,:)    :: grad_u
    !
    !    ..INDICES..
    integer :: i,j,k !DO NOT REMOVE THIS LINE. i IS PASSED INTO gradient() AND ITS VALUE IS MODIFIED.

    allocate (A(3,3,i_GRID,j_GRID,zLower:zUpper))
    allocate (grad_u(3,i_GRID,j_GRID,zLower:zUpper))

    !
    ! COMPUTE VELOCITY GRADIENTS:
    do i=1,3
       call gradient(u(i,:,:,:),grad_u) ! Returns dudx, dudy, dudz for each ui
       A(i,:,:,:,zLower:zUpper) = grad_u(:,:,:,zLower:zUpper)  ! A(1),A(2),A(3) -> grad_u_x, grad_u_y, grad_u_z
    end do
    deallocate(grad_u)

!    print*, A(:,:,15,24,zLower)
 
    !
    ! COMPUTE S_ij:
    k = 0
    do j=1,3
       do i=1,3

          if (i.ge.j) then
             k = k + 1
             Sij(k,:,:,zLower:zUpper) = 0.5d0 * (A(i,j,:,:,zLower:zUpper) + A(j,i,:,:,zLower:zUpper))
          end if

       end do
    end do
    deallocate (A)
    
!    print*, 'Sij(2,15,24,zLower), Sij(2,15,24,z_plane), Sij(2,15,24,zUpper)'
!    print*, Sij(2,15,24,zLower), Sij(2,15,24,z_plane), Sij(2,15,24,zUpper)
!    stop

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
    real(8), dimension (1:,1:,1:,zLower:),intent(out) :: grad_f
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
    allocate(u (0:i_GRID+1, 0:j_GRID+1, zLower-1:zUpper+1))

    if (compDomain.eq.'all') then
       ! COPY INTERIOR CELLS IN X-, Y-, Z-DIRS:
       u(1:i_GRID, 1:j_GRID, 1:k_GRID) = f(1:i_GRID, 1:j_GRID, 1:k_GRID)

       ! CREATE GHOST CELLS: (BASED ON PERIODIC BC) AT 6 FACES:
       u(0,:,:) = u(i_GRID,:,:) ; u(i_GRID+1,:,:) = u(1,:,:)
       u(:,0,:) = u(:,j_GRID,:) ; u(:,j_GRID+1,:) = u(:,1,:)
       u(:,:,0) = u(:,:,k_GRID) ; u(:,:,k_GRID+1) = u(:,:,1) ! [ALL] Apply Periodic BC in z-dir

    else if (compDomain.eq.'plane') then
       ! COPY INTERIOR CELLS IN X- AND Y-DIRS: Retain ONE plane above and below z_plane
       u(1:i_GRID, 1:j_GRID, zLower-1:zUpper+1) = f(1:i_GRID, 1:j_GRID, zLower-1:zUpper+1)

       ! CREATE GHOST CELLS: (BASED ON PERIODIC BC) AT 4 FACES:
       u(0,:,:) = u(i_GRID,:,:) ; u(i_GRID+1,:,:) = u(1,:,:)
       u(:,0,:) = u(:,j_GRID,:) ; u(:,j_GRID+1,:) = u(:,1,:)

    end if

    !
    ! APPLY 3D CENTRAL DIFFERENCING SCHEME AT INTERIOR POINTS:
    do k=zLower, zUpper
    do j=1,j_GRID
    do i=1,i_GRID
       grad_f(1,i,j,k) = dx_inv * (  u(i+1,j,k) - u(i-1,j,k) )
       grad_f(2,i,j,k) = dx_inv * (  u(i,j+1,k) - u(i,j-1,k) )
       grad_f(3,i,j,k) = dx_inv * (  u(i,j,k+1) - u(i,j,k-1) )
    end do
    end do
    end do

!    print*, grad_f(:,15,24,zLower)

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
  ! BEHAVIOR: While passing an array section, ensure that 
  !           zLower maps to zLower. 
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
    real(8), dimension(1:,1:,1:,zLower:), intent(in) :: tau_ij
    real(8), dimension(1:,1:,1:,zLower:), intent(in) :: S_ij
    real(8), dimension(1:,1:,zLower:), intent(out) :: P
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
                P(:,:,zLower:zUpper) = P(:,:,zLower:zUpper) + &
                     tau_ij(count,:,:,zLower:zUpper) * S_ij (count,:,:,zLower:zUpper) 
             else
                P(:,:,zLower:zUpper) = P(:,:,zLower:zUpper) + &
                     2.d0 * (tau_ij(count,:,:,zLower:zUpper) * S_ij(count,:,:,zLower:zUpper))
             end if
          end if
       end do
    end do

  end subroutine productionTerm

  
  !****************************************************************
  !                              VORTICITY
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE:  Calculates CURL using central order differencing
  !       Returns vorticity, omega as curl of velocity
  !                 
  !  
  !
  ! FORM: subroutine vorticity(f, curl_f)
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

  subroutine vorticity(f, curl_f)
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension (:,:,:,:), intent(in) ::  f
    real(8), dimension (1:,1:,1:,zLower:), intent(out) :: curl_f
    !
    !    ..LOCAL ARRAY..
    real(8), allocatable,dimension (:,:,:,:) :: u
    !
    !    ..LOCAL SCALARS..
    real(8) :: dx_inv

    !
    ! FIND 1/dx FOR CENTRAL DIFFERENCING
    dx_inv = 1.d0/(2.d0*dx)

    !
    ! LOAD f ON WORK ARRAY u WITH EXTENDED DOMAIN TO CONTAIN GHOST CELLS
    allocate(u (1:3,0:i_GRID+1, 0:j_GRID+1, zLower-1:zUpper+1))

    if (compDomain.eq.'all') then
       ! COPY INTERIOR CELLS IN X-, Y-, Z-DIRS:
       u(:,1:i_GRID, 1:j_GRID, 1:k_GRID) = f(:,1:i_GRID, 1:j_GRID, 1:k_GRID)

       ! CREATE GHOST CELLS: (BASED ON PERIODIC BC) AT 6 FACES:
       u(:,0,:,:) = u(:,i_GRID,:,:) ; u(:,i_GRID+1,:,:) = u(:,1,:,:)
       u(:,:,0,:) = u(:,:,j_GRID,:) ; u(:,:,j_GRID+1,:) = u(:,:,1,:)
       u(:,:,:,0) = u(:,:,:,k_GRID) ; u(:,:,:,k_GRID+1) = u(:,:,:,1) ! [ALL] Apply Periodic BC in z-dir

    else if (compDomain.eq.'plane') then
       ! COPY INTERIOR CELLS IN X- AND Y-DIRS: Retain ONE plane above and below z_plane
       u(:,1:i_GRID, 1:j_GRID, zLower-1:zUpper+1) = f(:,1:i_GRID, 1:j_GRID, zLower-1:zUpper+1)

       ! CREATE GHOST CELLS: (BASED ON PERIODIC BC) AT 4 FACES:
       u(:,0,:,:) = u(:,i_GRID,:,:) ; u(:,i_GRID+1,:,:) = u(:,1,:,:)
       u(:,:,0,:) = u(:,:,j_GRID,:) ; u(:,:,j_GRID+1,:) = u(:,:,1,:)

    end if

    !
    ! APPLY 3D CENTRAL DIFFERENCING SCHEME AT INTERIOR POINTS:
    do k=zLower, zUpper
    do j=1,j_GRID
    do i=1,i_GRID
       curl_f(1,i,j,k) = dx_inv * ( ( u(3,i,j+1,k) - u(3,i,j-1,k) ) - ( u(2,i,j,k+1) - u(2,i,j,k-1) ) )
       curl_f(2,i,j,k) = dx_inv * ( ( u(1,i,j,k+1) - u(1,i,j,k-1) ) - ( u(3,i+1,j,k) - u(3,i-1,j,k) ) )
       curl_f(3,i,j,k) = dx_inv * ( ( u(2,i+1,j,k) - u(2,i-1,j,k) ) - ( u(1,i,j+1,k) - u(1,i,j-1,k) ) )
    end do
    end do
    end do

!    print*, grad_f(:,15,24,zLower)
 
  end subroutine vorticity
  



  
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
  !           This is whole domain computation. Applies only when 
  !           compDomain is [ALL]
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
  ! STATUS :
  ! 
  !          
  !  STASH:
  !  uiuj(count,:,:,z_extLower:z_extUpper) = ui(i,:,:,z_extLower:z_extUpper)
  !                                        * ui(j,:,:,z_extLower:z_extUpper)
  !----------------------------------------------------------------

  subroutine secondOrderProducts(uiuj,ui)
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(1:,extLower:,extLower:,z_extLower:), intent(in) :: ui
    real(8), dimension(1:,extLower:,extLower:,z_extLower:), intent(out) :: uiuj
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
  ! USE: Use ghost cells to extend domain. Implements cyclic  
  !       padding for periodic BC.
  !       
  !
  ! FORM: subroutine extendDomain()
  !       
  !
  ! BEHAVIOR: 
  !           1. Copies interior cells to a temp array
  !           2. If ALL, then pads in x,y,z direction
  !           3. If PLANE, then pads in x,y direction
  !           or z-direction if needed.
  !           4. n_extendedLayers is determined right 
  !              in the beginning [global.f90] 
  !          
  !
  ! STATUS : 1. Adding capability for other planes.
  !          2. Generalizing ALL/PLANE case.
  !
  ! STASH :
  ! *
  !     allocate (temp (tempDim(1), extLower:extUpper, &
  !    extLower:extupper, &
  !    min(1,z_extLower):max(k_GRID,z_extUpper)))
  ! **
  !   temp(:, extLower:0,:,:) = array(:, i_GRID+extLower:i_GRID,:,:)
  !   temp(:, i_GRID+1:extUpper,:,:) = array(:,1:extUpper-i_GRID,:,:)
  !   temp(:, :, extLower:0,:) = array(:, :, j_GRID+extLower:j_GRID,:)
  !   temp(:, :,j_GRID+1:extUpper,:) = array(:,:,1:extUpper-j_GRID,:)
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
    integer :: tempDim(4)
    integer :: i,j

    tempDim = shape(array) !*
    allocate (temp (tempDim(1), extLower:extUpper, &
                                extLower:extUpper, &
                                min(1,z_extLower):max(k_GRID,z_extUpper)))

    temp(:, 1:i_GRID, 1:j_GRID, 1:k_GRID) = array(:, 1:i_GRID, 1:j_GRID, 1:k_GRID)
    deallocate(array)

    ! TEMP: CYCLIC PADDING
    ! 1) X,Y DIRECTION **
    ! X:
    temp (:,extLower:0,:,:) = temp(:,i_GRID+extLower:i_GRID,:,:)
    temp (:,i_GRID+1:extUpper,:,:) = temp(:,1:extUpper-i_GRID,:,:)

    ! Y:
    temp (:,:,extLower:0,:) = temp(:,:,j_GRID+extLower:j_GRID,:)
    temp (:,:,j_GRID+1:extUpper,:) = temp(:,:,1:extUpper-j_GRID,:)
    
    ! 2) Z-DIRECTION
    if (z_extLower < 1)        temp (:,:,:,z_extLower:0) = temp(:,:,:,k_GRID+z_extLower:k_GRID)
    if (z_extUpper > k_GRID)   temp (:,:,:,k_GRID+1:z_extUpper) = temp(:,:,:,1:z_extUpper-k_GRID)

    ! REALLOCATE ARRAY WITH PADDED ELEMENTS:
    allocate(array(tempDim(1), extLower:extUpper, extLower:extUpper, z_extLower:z_extUpper))
    array(:,:,:,z_extLower:z_extUpper) = temp(:,:,:, z_extLower:z_extUpper)

  end subroutine extendDomain


  
  !****************************************************************
  !                      PRECOMPUTE S_f_Sij_f_t                    !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Precompute S_f_Sij_f_t term 
  !       
  !       
  !
  ! FORM: subroutine computeS_f_Sij_f_t(Sij_f, test, S_f_Sij_f_t)
  !       
  !
  ! BEHAVIOR: 1] Performs whole domain tau_DS and saves it
  !          
  ! 
  ! 1) Multiply S_f(scalar) and Sij_f (vector)
  ! 2) Needs whole domain computations for Sij_f - use computeSij_all ! NEW SUBROUTINE
  ! 3) Filter at test scale -- requires a 3D sharp filter at test scale
  ! 
  !          
  ! STATUS :  
  !
  !          
  !----------------------------------------------------------------

  subroutine computeS_f_Sij_f_t(Sij_f, test, S_f_Sij_f_t)
    implicit none
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(1:,1:,1:,zLower:), intent(in) :: Sij_f
    real(8), dimension(   :,:,:)        , intent(in) :: test
    real(8), dimension(1:,1:,1:,zLower:), intent(out) :: S_f_Sij_f_t
    !
    !    ..WORK ARRAY - SCALARS..
    real(8), dimension(:,:,:), allocatable :: S_f
    !
    !   .. LOCAL VARS..
    integer :: ij 
     
    ! Find S_f = sqrt(sum (Sij_f))
    allocate (S_f(i_GRID,j_GRID,zLower:zUpper))
    S_f = sqrt(2.d0 * sum_ij(Sij_f,Sij_f))

    ! Compute S_f_Sij_f_t
    do ij = 1,6
       S_f_Sij_f_t(ij,:,:,:) = sharpFilter(S_f * Sij_f(ij,:,:,:),test)
    end do
 
  end subroutine computeS_f_Sij_f_t
  
  

  
  !****************************************************************
  !                      DYNAMIC SMAGORINSKY                      !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Calculate Dynamic Smagorinsky stress (tau_DS)
  !       
  !       
  !
  ! FORM: subroutine dyn_Smag(tau_DS)
  !       
  !
  ! BEHAVIOR: 1] Performs whole domain tau_DS and saves it
  !          
  ! 
  !
  ! STATUS :  
  !          
  !----------------------------------------------------------------

  subroutine dyn_Smag(Sij_f, Sij_t, S_f_Sij_f_t, T_ij, tau_DS)
    implicit none
    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(1:,1:,1:,zLower:), intent(in) :: T_ij
    real(8), dimension(1:,1:,1:,zLower:), intent(in) :: Sij_f
    real(8), dimension(1:,1:,1:,zLower:), intent(in) :: Sij_t
    real(8), dimension(1:,1:,1:,zLower:), intent(in) :: S_f_Sij_f_t
    real(8), dimension(1:,1:,1:,zLower:), intent(out):: tau_DS
    !
    !    ..WORK ARRAY - SCALARS..
    real(8), dimension(:,:,:), allocatable :: S_f
    real(8), dimension(:,:,:), allocatable :: S_t
    real(8), dimension(:,:,:,:), allocatable :: M_ij
    real(8), dimension(:,:,:), allocatable :: Cs_2        ! DS CONSTANT

    !
    !   .. LOCAL VARS..
    real(8) :: Delta_tilde ! LES GRID SPACING 
    real(8) :: Delta_hat   ! TEST GRID SPACING
    !
    integer :: count 
    
    ! Find Delta_tilde, Delta_hat
    ! dx is known in the beginning - global.f90
!    Delta_tilde = 2.d0*PI/dble(LES_scale) * dx ![Delta_LES = 3]
!    Delta_hat = 2.d0*PI/dble(test_scale) * dx ! [Delta_test = 6]

    Delta_tilde = 2.d0*PI/dble(LES_scale)  ![Delta_LES = 3]
    Delta_hat = 2.d0*PI/dble(test_scale)   ! [Delta_test = 6]

    ! Find Sij_t:
    ! Filter Sij_f at test scale
    !
    ! 1) Sij_t = filter(Sij_f, test_scale)
    ! 
    ! (Not required since the low pass sharp filter
    ! at the test scale is smaller than the LES_scale) 


    !
    ! Find S_f = sqrt(sum (Sij_f))
    !      S_t = sqrt(sum (Sij_t))            [S_f, S_t are scalars]
    allocate (S_f(i_GRID,j_GRID,zLower:zUpper))
    allocate (S_t(i_GRID,j_GRID,zLower:zUpper))

    S_f = sqrt(2.d0 * sum_ij(Sij_f,Sij_f))
    S_t = sqrt(2.d0 * sum_ij(Sij_t,Sij_t))

!    print*, 'S_f(15,24,129)', S_f(15,24,129)
!    print*, 'S_t(15,24,129)', S_t(15,24,129)

    !
    ! Find M_ij = 2*[ Delta_hat**2 * S_t * Sij_t - Delta_tilde**2 * S_f_Sij_f_t]
    allocate (M_ij (6,i_GRID,j_GRID,zLower:zUpper))
    do count = 1,6
       M_ij(count,:,:,:) =  2.d0 * ( Delta_hat**2 * S_t * Sij_t(count,:,:,:) -  &
                                   Delta_tilde**2 * S_f_Sij_f_t(count,:,:,:) )
    end do
!    print*, 'M_ij', M_ij(2,15,24,129)
    
    !
    ! Find Cs_2(:,:,:) = sum(T_ij * M_ij) / sum(M_ij *  M_ij)
    allocate (Cs_2 (i_GRID,j_GRID,zLower:zUpper))
    Cs_2 =  sum_ij(M_ij,T_ij) / sum_ij(M_ij,M_ij)
!    print*, 'Cs_2', Cs_2(15,24,129)

    !
    ! Find tau_DS = 2 * (Cs_2 * Delta_tilde)**2 * S_f * Sij_f
    !      T_DS   = 2 * (Cs_2 * Delta_hat)**2   * S_t * Sij_t
    do count = 1,6
       tau_DS(count,:,:,:) = 2.d0 * & 
                      Cs_2 * Delta_tilde ** 2 * S_f * Sij_f(count,:,:,:)
     end do
     print*, 'tau_DS(2,15,24,129) ',tau_DS(2,15,24,129)
       

  end subroutine dyn_Smag


  !****************************************************************
  !                        SUM IJ                                 !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Adds up all the elements of a symmetric matrix 
  !       
  !       
  ! FORM: real(8) function sum_ij(matrix)
  !       
  !
  ! BEHAVIOR: Multiplies off-diagonal elements by 2 due to symmetry 
  !           
  !          
  !
  ! STATUS :  
  !
  ! 
  !          
  !----------------------------------------------------------------

  
  function sum_ij(mat1, mat2)
    !
    !    ..ARRAY ARGUMENT..
    real(8), dimension(1:,1:,1:,zLower:), intent(in) :: mat1
    real(8), dimension(1:,1:,1:,zLower:), intent(in) :: mat2
    real(8), dimension(i_GRID,j_GRID,zLower:zUpper) :: sum_ij

                  sum_ij = sum_ij           &
                         +                  &
         ( mat1(1,:,:,:) * mat2(1,:,:,:) +  &
   2.d0 *  mat1(2,:,:,:) * mat2(2,:,:,:) +  &
   2.d0 *  mat1(3,:,:,:) * mat2(3,:,:,:) +  &
           mat1(4,:,:,:) * mat2(4,:,:,:) +  &
   2.d0 *  mat1(5,:,:,:) * mat2(5,:,:,:) +  &
           mat1(6,:,:,:) * mat2(6,:,:,:) )

  end function sum_ij

  
  
    
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


  !****************************************************************
  !                        ROTATE X                               !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Performs +90 deg rotation on a rank-4 array along x-axis (i-dir)       
  !       
  !
  ! FORM: subroutine rotateX(a)
  !       
  !
  ! BEHAVIOR: The original y-axis now becomes z-axis running from
  !           right to left. Old indices (2,15,24,129) -> New (2,15,128,24) 
  !           
  !           
  !           
  ! STATUS :  
  !
  !
  !----------------------------------------------------------------
  
  subroutine rotateX(a)

    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:,:,:), intent(inout):: a
    !
    !   ..LOCAL VARIABLES..
    real(8), dimension(:,:,:,:), allocatable :: temp
    integer :: j,k
    integer :: dims(4)

    dims = shape(a)
    allocate(temp(dims(1),dims(2), dims(3), dims(4)))
    
    do k = 1, size(a,dim=4)
       do j = 1, size(a,dim=3)
          temp (:,:,j_GRID-j+1,k) = a(:,:,k,j)
       end do
    end do
    a = temp
  end subroutine rotateX

  
  !****************************************************************
  !                        ROTATE Y                               !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Performs +90 deg rotation on a rank-4 array along y-axis (j-dir)       
  !       
  !
  ! FORM: subroutine rotateY(a)
  !       
  !
  ! BEHAVIOR: The original y-axis now becomes z-axis running from
  !           right to left. Old indices (2,15,24,129) -> New (2,15,128,24) 
  !           
  !           
  !           
  ! STATUS :  
  !
  !
  !----------------------------------------------------------------
  
  subroutine rotateY(a)

    !
    !    ..ARRAY ARGUMENTS..
    real(8), dimension(:,:,:,:), intent(inout):: a
    !
    !   ..LOCAL VARIABLES..
    real(8), dimension(:,:,:,:), allocatable :: temp
    integer :: i,k
    integer :: dims(4)

    dims = shape(a)
    allocate(temp(dims(1),dims(2), dims(3), dims(4)))

    do k = 1, size(a,dim=4)
       do i = 1, size(a,dim=2)
          temp (:,i,:,k_GRID-k+1) = a(:,k,:,i)
       end do
    end do

    a = temp
  end subroutine rotateY


end module actools
