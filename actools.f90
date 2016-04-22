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

    allocate (A(3,3,grid,grid,grid))
    allocate (grad_u(3,grid,grid,grid))

    
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
    real(8) :: dx = 2.d0*pi/dble(grid) !Only for JHU data. Change for others.
    real(8) :: dx_inv

    dx_inv = 1.d0/(2.d0*dx)

    allocate(u (0:grid+1, 0:grid+1, 0:grid+1))
    u(1:grid,1:grid,1:grid) = f

    ! CREATE GHOST CELLS: (BASED ON PERIODIC BC)
    ! NEED AT 6 FACES:
    u(0,:,:) = u(grid,:,:) ; u(grid+1,:,:) = u(1,:,:)
    u(:,0,:) = u(:,grid,:) ; u(:,grid+1,:) = u(:,1,:)
    u(:,:,0) = u(:,:,grid) ; u(:,:,grid+1) = u(:,:,1)


    ! APPLY 3D CENTRAL DIFFERENCING SCHEME AT INTERIOR POINTS:
    do k=1,grid
    do j=1,grid
    do i=1,grid
       grad_f(1,i,j,k) = dx_inv * (  u(i+1,j,k) - u(i-1,j,k) )
       grad_f(2,i,j,k) = dx_inv * (  u(i,j+1,k) - u(i,j-1,k) )
       grad_f(3,i,j,k) = dx_inv * (  u(i,j,k+1) - u(i,j,k-1) )
    end do
    end do
    end do

    return
  end subroutine gradient
  
!   subroutine testNonCo(u_f,u_t,tau_ij,T_ij,h_ij)

!     integer,parameter :: n_u=3, n_uu=6

!     ! PARAMETERS:
!     real(8) :: lambda = 1.d-1

!     ! ARGUMENTS:
!     real(8), dimension(:,:,:,:), intent(in) :: u_f, u_t
!     real(8), dimension(:,:,:,:), intent(in) :: tau_ij, T_ij
!     real(8), dimension(:,:),intent(out) :: h_ij

!     ! LOCAL VARIABLES:
!     real(8),allocatable,dimension(:,:,:,:) :: TijOpt, tau_ijOpt !Computed stresses

!     ! FOR LEAST SQUARES PROBLEM:
!     real(8),dimension(M,N)  :: V = 0. ! Non-linear combination of u_i, u_iu_j
!     real(8),dimension(M,P)   :: T = 0. ! Input training data in LAPACK routines.

!     integer                :: row_index, col_index, row, col ! Indices to build V,h_ij arrays.
!     integer                :: u_comp, uu_comp ! Indices to select u_i,u_ij components


!     ! FOR NON-COLOCATED FORMULATION:                                                                                        
!     logical :: coloc = 0            !0 means use non-coloc                                                                  
!     integer,parameter :: stencil_size = 3*(3*3*3) !3 * 27 = 81
!     real(8) :: u_n(stencil_size)
!     integer :: non_col_1, non_col_2  ! To run non coloc combinantions   

!     ! LIST OF INDICES FOR LOOPS:
!     integer :: i_test,    j_test,    k_test ! for a point in the TEST-scale field
!     integer :: i_box,     j_box,     k_box  ! for a point in bounding BOX (stencil-center)
!     integer :: i_stencil, j_stencil, k_stencil


!     ! DEBUG SWITCHES:
!     logical :: debug(4)=[1,0,0,0]
!     ! 1 - Takes just one (tau_ij_11) for testing purpose. 
!     ! 2 - Prints cutout shapes to check if cutout ops have gone well.
!     ! 3 - Performs all computations on just the first point.
!     ! 4 - To check if values from the stencil and bounding box are working well.
!     integer :: lim
!     logical :: printval



!     ! WHOLE DOMAIN COMPUTATION: 
!     do i_test = 129, 129, 2 
!     do j_test = 129, 129, 2
!     do k_test = 129, 129, 2 ! i_test = 11,43,2

!        row_index  = 0 

!        ! ENTER STENCIL-CENTER POINTS: C-ORDER
!        do i_box = i_test-126, i_test+125, 10
!        do j_box = j_test-126, j_test+125, 10
!        do k_box = k_test-126, k_test+125, 10 ! i_box = 3,49,2

!           col_index = 0 
!           row_index = row_index + 1
         
!          ! ENTER 3x3x3 STENCIL: C-ORDER
                                                                            
!           u_n = reshape(u_t(:,i_box-2:i_box+2:2, j_box-2:j_box+2:2,k_box-2:k_box+2:2),(/stencil_size/)) 

!           ! ZERO ORDER TERMS:
!           col_index = col_index+1
!           V(row_index, col_index) = 1.d0

          
!           ! FIRST ORDER TERMS:          
!           do non_col_1 = 1,stencil_size 
!              col_index = col_index+1
!              V(row_index,col_index) = u_n(non_col_1)
!           end do

!           ! SECOND ORDER TERMS: 6x(3x3x3) = 162 (GIVES A TOTAL OF 243 TERMS)
!           do non_col_1 = 1, stencil_size
!              do non_col_2 = non_col_1, stencil_size
!                 col_index = col_index+1
!                 V(row_index,col_index) = u_n(non_col_1) * u_n(non_col_2)
!              end do
!           end do

!           T(row_index,:) = T_ij(:,i_box,j_box,k_box) !Change 1 to (1-6) here. !THIS ONE IS CORRECT; KEEP IT.

!       end do
!       end do
!       end do ! BOUNDING BOX

!       print*,row_index,col_index
!       ! CHECK V:
!       if (V(1500,2000).ne.2.0009431419772586d-2) then
!          print*, "Error! Check sorting order in  V matrix!"
!          print*, 'V(1500,2000)', V(1500,2000)
!          stop
!       else 
!          print*,'V matrix check ... Passed'
!       end if


!       !CHECK T:
!       if (T(3,1).ne.8.7759832493259110d-2)then
!          print*, "Error! Check sorting order in  T vector!"
!          print*,T(3,1)
!          stop
!       else 
!          print*,'T vector check ... Passed'
!       end if
!       !print*, 'T(1,2)', T(2)
!       !print*, 'T(1,3)', T(3)
!       print*, 'T(1,4)', T(1,4)
!       !print*, 'T(1,5)', T(5)
      
      
! !      if (computeVolterra) then

!          !         call SVD(V,T,h_ij,lambda,printval) ! TSVD
!          call LU(V,T,h_ij,lambda)           !Damped least squares
!          ! PRINT h_ij:
!          print*, 'h_ij(1,1):',h_ij(1,1)
!          print*, 'h_ij(20,1):', h_ij(20,1)
!          print*, 'h_ij(350,1):', h_ij(350,1)

!          stop
!          ! CHECK h_ij:
!          if (h_ij(350,1).ne.-4.5121154730201521d-2)then
!             print*, "Error! Check lambda, method or sorting order for h_ij computation:"
!             print*,h_ij(350,1)
!             stop
!          else 
!             print*,'SVD check ... Passed'
!          end if
! !      end if

!       ! SAVE/READ h_ij:
!       open(1,file='./run/apples/h_ij.dat')
!       read(1,*) h_ij
!       print*,'h_ij(350,1) from file: ',h_ij(350,1)
!       close(1)
!    end do
!    end do
!    end do ! test

!    return
!  end subroutine testNonCO


 
!  subroutine optimizedTij(u_f,u_t,h_ij,TijOpt,tau_ijOpt)

!    ! ARGUMENTS:
!    real(8), dimension(:,:,:,:), intent(in) :: u_f, u_t
!    real(8), dimension(:,:),intent(in) :: h_ij

!    ! LOCAL VARIABLES:
!    real(8),allocatable,dimension(:,:,:,:),intent(out) :: TijOpt, tau_ijOpt !Computed stresses

!    integer                :: col_index
!    integer                :: u_comp, uu_comp ! Indices to select u_i,u_ij components


!    ! FOR NON-COLOCATED FORMULATION:
!    integer,parameter :: stencil_size = 3*(3*3*3) !3 * 27 = 81
!    real(8),dimension(stencil_size) :: u_f_stencil, u_t_stencil
!    integer :: non_col_1, non_col_2  ! To run non coloc combinantions   

!    ! LIST OF INDICES FOR LOOPS:
!    integer :: i_test,    j_test,    k_test ! for a point in the TEST-scale field
!    integer :: i_opt,     j_opt,     k_opt  ! for a point in bounding BOX (stencil-center)
!    integer :: i_stencil, j_stencil, k_stencil


!    allocate (TijOpt(1,grid,grid,grid))
!    allocate (tau_ijOpt(1,grid,grid,grid))
!    ! WHOLE DOMAIN COMPUTATION: 
!    do i_test = 129, 129
!    do j_test = 129, 129
!    do k_test = 129, 129

!       ! ENTER STENCIL-CENTER POINTS: C-ORDER
!       do i_opt = i_test-127, i_test+126
!       do j_opt = j_test-127, j_test+126
!       do k_opt = k_test-127, k_test+126

!          col_index = 0 
         
!          ! ENTER 3x3x3 STENCIL: C-ORDER                                       
!           u_t_stencil = reshape(u_t(:,i_opt-2:i_opt+2:2, j_opt-2:j_opt+2:2,k_opt-2:k_opt+2:2),(/stencil_size/)) 
!           u_f_stencil = reshape(u_f(:,i_opt-1:i_opt+1:1, j_opt-1:j_opt+1:1,k_opt-1:k_opt+1:1),(/stencil_size/)) 

!           ! ZERO ORDER TERMS:
!           col_index = col_index+1
!           TijOpt(:,i_opt,j_opt,k_opt) = h_ij(col_index,:) ! = 1
!           tau_ijOpt(:,i_opt,j_opt,k_opt) = h_ij(col_index,:) ! = 1

!           ! FIRST ORDER TERMS:             
!           do non_col_1 = 1,stencil_size 
!              col_index = col_index+1
!              TijOpt(:,i_opt,j_opt,k_opt) = TijOpt(:,i_opt,j_opt,k_opt) + (u_t_stencil(non_col_1)*h_ij(col_index,:))
!              tau_ijOpt(:,i_opt,j_opt,k_opt) = tau_ijOpt(:,i_opt,j_opt,k_opt) + (u_f_stencil(non_col_1)*h_ij(col_index,:))
!           end do

!           ! SECOND ORDER TERMS: 6x(3x3x3) = 162 (GIVES A TOTAL OF 243 TERMS)
!           do non_col_1 = 1, stencil_size
!           do non_col_2 = non_col_1, stencil_size
!              col_index = col_index+1
!              TijOpt(:,i_opt,j_opt,k_opt) = TijOpt(:,i_opt,j_opt,k_opt) &
!                   + (u_t_stencil(non_col_1) * u_t_stencil(non_col_2)*h_ij(col_index,:))
!              tau_ijOpt(:,i_opt,j_opt,k_opt) = tau_ijOpt(:,i_opt,j_opt,k_opt) &
!                   + (u_f_stencil(non_col_1) * u_f_stencil(non_col_2)*h_ij(col_index,:))
!           end do
!           end do
!        end do
!        end do
!        end do ! BOUNDING BOX
!     end do
!     end do
!     end do ! test

!     ! WRITE COMPUTED STRESS:
!     open(1,file='./run/apples/T_ijOpt.dat')
!     open(2,file='./run/apples/tau_ijOpt.dat')
!     write(1,*) TijOpt(1,:,:,129)
!     write(2,*) tau_ijOpt(1,:,:,129)
!     close(1)
!     close(2)
!     print*,'TijOpt(1,15,24,10):',TijOpt(1,15,24,10)
!     print*,'tau_ijOpt(1,15,24,10):',tau_ijOpt(1,15,24,10)
!     return
!   end subroutine optimizedTij


  
!   !****************************************************************
!   !                                LU                             !
!   !****************************************************************
!   subroutine LU(V, T_ij, h_ij, lambda)
!     implicit none

!     ! dgesv destroys the original matrix. So V is copied into 'a' matrix. This is the LU product matrix.
!     ! dgesv returns x vector. First it copies the values from 'T_ij' vector and then computes 'h_ij'.
!     ! 'h_ij' is the h_ij vector.

!     ! ARGUMENTS:
!     real(8), dimension(:,:), intent(in) :: V 
!     real(8), dimension(:,:), intent(in) :: T_ij
!     real(8), dimension(:,:), intent(out) :: h_ij
!     real(8), intent(in) :: lambda

!     ! DGESV ARGUMENTS:
!     integer, parameter        :: LDA = N, LDB = N, nrhs = P
!     real(8), dimension(:,:),allocatable :: A,VT ! EYE - IDENTITY MATRIX
!     real(8), dimension(:,:),allocatable :: b
!     integer, dimension(:), allocatable  :: ipiv
!     integer                   :: info
!     real(8) :: tic, toc


!     !DGEMM ARGUMENTS:
!     real(8):: alpha=1.d0, beta=0.d0

!     allocate (A(N,N), b(N,P), ipiv(N))

!     ! VT = transpose(V)
!     call cpu_time(tic)
!     call dgemm('T','N',N,N,M, alpha, V,M, V,M, beta, A,N)
!     call cpu_time(toc)
!     print*,'Elapsed time', toc-tic


!     !A
!     forall(i=1:N) A(i,i) = A(i,i) + lambda 
!     print*,'A:'

!     call printplane(A,frameLim=4)
!     !b
!     call dgemm('T','N',N,P,M, alpha, V,M, T_ij,M, beta, b,N)
!     print*,'b:'

!     call printplane(b,frameLim=4)
!     call cpu_time(tic)
!     !$omp parallel 
!     call DGESV(N, nrhs, A, LDA, ipiv, b, LDB, info) !A is a symmetric matrix

!     !$omp end parallel 
!     call cpu_time(toc)
!     print*,'Elapsed time:', toc-tic

!     h_ij = b

!     deallocate (A,b,ipiv)
!     return
!   end subroutine LU

!   !****************************************************************
!   !                             SVD                               !
!   !****************************************************************

!   subroutine SVD(A,T,h_ij,lambda,printval)

!     ! ARGUMENTS:
!     real(8),dimension(:,:),intent(in)  :: A
!     real(8),dimension(:,:),intent(in)  :: T
!     real(8),intent(in)                 ::lambda
!     real(8),dimension(:,:),intent(out) :: h_ij


!     ! DGESVD ARGUMENTS:                                                                                                
!     integer, parameter :: LDA = M
!     integer, parameter :: LDU = M
!     integer, parameter :: LDVT = N
!     integer :: LWMAX
!     integer :: info, LWORK

!     real(8),dimension(:,:),allocatable :: U
!     real(8),dimension(:,:),allocatable :: VT 
!     real(8),dimension(:,:),allocatable :: D, Vinv
!     real(8),dimension(:), allocatable :: S, work


!     logical :: printval


!     LWMAX = M*N
!     allocate (U(LDU,M), VT(N,N), S(N), D(N,M), Vinv(N,M), work(LWMAX))


!     if (printval) then
!        print*,'A'                      ! A MATRIX
!        call printplane(A,frameLim=4)
!     end if


!     ! CALL SVD ROUTINE:
!     LWORK = -1
!     call dgesvd('All','All', M,N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, info)
!     LWORK = min( LWMAX, int(WORK(1)) ) 
!     call dgesvd('All','All',M,N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, info) !BOTTLENECK = 850 SECS

!     ! Convergence check:
!     if (info.gt.0) then
!        print*, 'Failed to converge'
!        stop
!     end if


!     if (printval) then
!        print*,'A'                      ! A MATRIX
!        call printplane(A,size(A,dim=1))
!        print*,'U'                      ! U MATRIX
!        call printplane(U,size(U,dim=1))
!        print*,'S'                      ! S VECTOR
!        print*, S(1:3)
!        print*,'VT'                     ! VT MATRIX
!        call printplane(VT,size(VT,dim=1))
!     end if


!     ! COMPUTE PSEUDOINVERSE:
!     D = 0.d0
!     forall(i=1:N) D(i,i) = S(i) / (S(i)**2 + lambda**2)
!     Vinv = matmul(matmul(transpose(VT),D),transpose(U)) !BOTTLENECK = 1300 SECS

!     if(printval) then
!        print*,'D'                     ! D MATRIX - DIAG MATRIX
!        call printplane(D,size(D,dim=1))

!        print*,'Vinv'                     ! Vinv MATRIX - DIAG MATRIX
!        call printplane(Vinv,size(Vinv,dim=1))
!     end if

!     h_ij = matmul(Vinv,T)
!     return
!   end subroutine SVD

end module actools
