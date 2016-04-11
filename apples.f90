program apples

use global
use fourier
use fileio

implicit none
integer,parameter           :: LES_scale=40, test_scale=20
character(3) :: stress = 'dev'

!integer, parameter :: M=26**3, N=3403, P=6  ! Define array size here.
real(8) :: lambda = 0.1d0         ! lambda, damping factor


! DEFINE VELOCITIES:
real(8),allocatable,dimension(:,:,:,:) :: u, u_f, u_t
real(8),allocatable,dimension(:,:,:,:) :: tau_ij,T_ij
real(8),allocatable,dimension(:,:,:,:) :: TijOpt, tau_ijOpt
real(8),allocatable,dimension(:,:,:) :: LES,test

! DEFINE STRAIN RATES:
real(8),allocatable,dimension(:,:,:,:,:) :: Sij_f, Sij_t
real(8),allocatable,dimension(:,:):: h_ij

integer :: n_u, n_uu

character(50):: CUT_DATA = '../derived_data/cutout-valid/jhu/' !Change bin4020 by 'sed' in shell script
character(10):: f_CUT 
character(3) :: d_set = 'jhu'   ! for binRead()


! DEBUG FLAGS:
logical :: debug(2) = (/0,2/)
! 1- To select only one velocity component and 3 velocity products.
! 2- Choose between NRL/JHU256 database[1] or JHU(HDF5) database[0]

logical :: computeStresses = 0
logical :: computeStrain = 0
logical :: filterVelocities = 0
logical :: readFile = 0
logical :: computeVolterra = 1

real(4) :: tic, toc
real(8) :: pre_cut, post_cut

integer :: i,j,k,d=0


! FORMAT:
3015 format(a30,f22.15)
307 format(a30,f22.7)
507 format(a50,f22.7)
!call system('clear')
!call printParams()

!! Set debug flags for velocity components:
n_u=3; n_uu = 6
if (debug(1)) then
   n_u = 1
   n_uu = 3
   print*, 'Debug mode for velocity components...'
end if


!! Select file to read:
if(readFile)then
allocate(u(n_u,GRID,GRID,GRID))

fileSelect:if (debug(2)) then
   call binRead(u,  d_set,  DIM=n_u)
   write(f_CUT,'(a3, 2(i2), a1)') 'bin',LES_scale,test_scale,'/' !Write dirname
   !print*, trim(CUT_DATA)//trim(f_CUT) !Test pathname
else
   call hdf5Read() 
end if fileSelect
end if


allocate(u_f(n_u,GRID,GRID,GRID))
allocate(u_t(n_u,GRID,GRID,GRID))
if (filterVelocities) then
   !! Create filters:
   allocate(LES(GRID,GRID,GRID))
   allocate(test(GRID,GRID,GRID))
   call createFilter(LES,LES_scale)
   call createFilter(test,test_scale)
   call fftshift(LES)
   call fftshift(test)

   ! FILTER VELOCITIES:

   print*,'Filter velocities ... '

   filter:do i=1,n_u
      u_f(i,:,:,:) = sharpFilter(u(i,:,:,:),LES) ! Speed up this part -- Bottleneck
      u_t(i,:,:,:) = sharpFilter(u_f(i,:,:,:),test) ! Speed up this part -- Bottleneck
   end do filter


   ! CHECK FFT:
   if (u_t(1,15,24,10).ne.-0.48241021987284982d0) then
      print*, 'Precision test for FFT: Failed'
      print*, 'Check precision or data for testing is not JHU 256'
      print*, u_t(1,15,24,10)
      stop
   else
      print*, 'Precision test for FFT: Passed'
   end if
end if

! ASSERT 6 COMPONENTS FOR ij TO COMPUTE STRESS:
if (n_uu.ne.6) then
   print*,"Need all 6 ij components to compute stress... Aborting"
   stop
end if

   allocate(tau_ij(n_uu,GRID,GRID,GRID))
   allocate(T_ij(n_uu,GRID,GRID,GRID))

! COMPUTE STRESS:
if(computeStresses)then

   call cpu_time(tic)
   print*,'Compute stress:',stress
   call computeStress(u,u_f,u_t,tau_ij,T_ij,n_u,n_uu,LES,test,stress)
   call cpu_time(toc)
   write(*,307),'computeStress - time elapsed:',toc-tic
   deallocate(LES,test)

else

! READ STRESS files:
   open(1,file="../derived_data/cutout-valid/jhu/bin4020/u_f.bin",form="unformatted")
   open(2,file="../derived_data/cutout-valid/jhu/bin4020/u_t.bin",form="unformatted")
   open(3,file="../derived_data/cutout-valid/jhu/bin4020/tau_ij.bin",form="unformatted")
   open(4,file="../derived_data/cutout-valid/jhu/bin4020/T_ij.bin",form="unformatted")

   read(1) u_f
   read(2) u_t
   read(3) tau_ij
   read(4) T_ij

   close(1)
   close(2)
   close(3)
   close(4)

   ! CHECK INPUT DATA:
if (u_t(1,15,24,10).ne.-0.48241021987284982d0) then
         print*, 'ERROR READING DATA'
         print*, u_t(1,15,24,10)
         stop
      elseif(T_ij(1,15,24,10).ne.-5.2544371578038401d-3) then
            print*, 'ERROR COMPUTING T_ij'
            print*,'T_ij(1,15,24,10)',T_ij(1,15,24,10)
            stop
         else
         print*, 'Read data saved from main.f90: Passed'
         end if
      end if


!open(1,file='./run/apples/T_ij.dat')
!open(2,file='./run/apples/tau_ij.dat')
!write(1,*) T_ij(1,:,:,129)
!write(2,*) tau_ij(1,:,:,129)
!close(1)
!close(2)



! ! COMPUTE STRAIN RATE:
if(computeStrain)then
   allocate (Sij_f(3,3,grid,grid,grid))
   allocate (Sij_t(3,3,grid,grid,grid))
   print*, 'Compute strain rate'
!   call computeSij(u_f, Sij_f)
!   call computeSij(u_t, Sij_t)

   print*,'Sij_fl(1,2,15,24,10)',Sij_f(1,2,15,24,10)
   print*,'Sij_ft(1,2,15,24,10)',Sij_t(1,2,15,24,10)
   print*,'Sij_ft(1,2,1,3,5)',Sij_t(1,2,1,3,5)

   deallocate (Sij_f, Sij_t)
end if


if (allocated(h_ij).neqv..true.) allocate(h_ij(N,P))
! COMPUTE h_ij:
call testNonCo(u_f,u_t,tau_ij,T_ij,h_ij)

!open(1,file='./run/apples/h_ij.dat')
!read(1,*) h_ij

!CHECK h_ij:
if (h_ij(350,1).ne.-4.5121154730201521d-2 ) then
   print*,'Error reading h_ij'
   print*,h_ij(350,1)
   stop
else
   print*,'Reading h_ij(1) ... Success'
end if

!close(1)
if (allocated(TijOpt).neqv..true.) allocate(TijOpt(n_uu,grid,grid,grid))
if (allocated(tau_ijOpt).neqv..true.) allocate(tau_ijOpt(n_uu,grid,grid,grid))

stop

call optimizedTij(u_f,u_t,h_ij,TijOpt,tau_ijOpt)
open(11,file='./run/apples/TijOpt.dat')
open(12,file='./run/apples/tau_ijOpt.dat')
read(11,*)TijOpt(1,:,:,129)
read(12,*)tau_ijOpt(1,:,:,129)

stop
!PRODUCTION FIELD
open(13,file='./run/apples/P11_t.dat')
open(14,file='./run/apples/P11_f.dat')

write(13,*) TijOpt(1,:,:,129)*Sij_t(1,1,:,:,129)
write(14,*) tau_ijOpt(1,:,:,129)*Sij_f(1,1,:,:,129)

close(11)
close(12)
close(13)
close(14)

stop
contains

!****************************************************************
!                             S_ij                              !
!****************************************************************

! subroutine computeSij(u,Sij)
! implicit none

! real(8), dimension(:,:,:,:),  intent(in) :: u
! real(8), dimension(:,:,:,:,:),intent(out) :: Sij

! real(8), dimension(3,3,grid,grid,grid) :: A
! real(8), dimension(:,:,:,:),allocatable :: grad_u
! integer :: i,j



! print*, 'Compute velocity gradients:'
! allocate (grad_u(3,grid,grid,grid))
! do i=1,3
!    call gradient(u(i,:,:,:),grad_u)
!    A(i,:,:,:,:) = grad_u        ! A(1),A(2),A(3) -> grad_u_x, grad_u_y, grad_u_z
! end do
!    deallocate(grad_u)

! ! COMPUTE S_ij:
! do i=1,3
! do j=1,3
!     Sij(i,j,:,:,:) = 0.5d0 * (A(i,j,:,:,:) + A(j,i,:,:,:))
! end do
! end do

! return
! end subroutine computeSij
  



subroutine gradient(f, grad_f)
implicit none

real(8), dimension (grid,grid,grid), intent(in) ::  f
real(8), dimension (3,grid, grid, grid),intent(out) :: grad_f
real(8), dimension (0:grid+1,0:grid+1,0:grid+1) :: u
real(8), parameter :: pi = 4.d0 * atan(1.d0)
real(8) :: dx = 2.d0*pi/dble(grid) !Only for JHU data. Change for others.
real(8) :: dx_inv
logical :: printval
integer :: i,j,k

dx_inv = 1.d0/(2.d0*dx)

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

subroutine testNonCo(u_f,u_t,tau_ij,T_ij,h_ij)

integer :: i,j,k     ! General indices
integer,parameter :: n_u=3, n_uu=6

! PARAMETERS:
real(8) :: lambda = 1.d-1

! ARGUMENTS:
real(8), dimension(:,:,:,:), intent(in) :: u_f, u_t
real(8), dimension(:,:,:,:), intent(in) :: tau_ij, T_ij
real(8), dimension(:,:),intent(out) :: h_ij

! LOCAL VARIABLES:
real(8),allocatable,dimension(:,:,:,:) :: TijOpt, tau_ijOpt !Computed stresses

! FOR LEAST SQUARES PROBLEM:
real(8),dimension(M,N)  :: V = 0. ! Non-linear combination of u_i, u_iu_j
real(8),dimension(M,P)   :: T = 0. ! Input training data in LAPACK routines.

integer                :: row_index, col_index, row, col ! Indices to build V,h_ij arrays.
integer                :: u_comp, uu_comp ! Indices to select u_i,u_ij components


! FOR NON-COLOCATED FORMULATION:                                                                                        
logical :: coloc = 0            !0 means use non-coloc                                                                  
integer,parameter :: stencil_size = 3*(3*3*3) !3 * 27 = 81
real(8) :: u_n(stencil_size)
integer :: non_col_1, non_col_2  ! To run non coloc combinantions   

! LIST OF INDICES FOR LOOPS:
integer :: i_test,    j_test,    k_test ! for a point in the TEST-scale field
integer :: i_box,     j_box,     k_box  ! for a point in bounding BOX (stencil-center)
integer :: i_stencil, j_stencil, k_stencil


! DEBUG SWITCHES:
logical :: debug(4)=(/1,0,0,0/)
! 1 - Takes just one (tau_ij_11) for testing purpose. 
! 2 - Prints cutout shapes to check if cutout ops have gone well.
! 3 - Performs all computations on just the first point.
! 4 - To check if values from the stencil and bounding box are working well.
integer :: lim, p
logical :: printval




! WHOLE DOMAIN COMPUTATION: 
do i_test = 129, 129, 2 
do j_test = 129, 129, 2
do k_test = 129, 129, 2 ! i_test = 11,43,2

      row_index  = 0 

      ! ENTER STENCIL-CENTER POINTS: C-ORDER
      do i_box = i_test-126, i_test+125, 10
      do j_box = j_test-126, j_test+125, 10
      do k_box = k_test-126, k_test+125, 10 ! i_box = 3,49,2

         col_index = 0 
         row_index = row_index + 1
         
         ! ENTER 3x3x3 STENCIL: C-ORDER
                                                                            
          u_n = reshape(u_t(:,i_box-2:i_box+2:2, j_box-2:j_box+2:2,k_box-2:k_box+2:2),(/stencil_size/)) 

            ! ZERO ORDER TERMS:
             col_index = col_index+1
             V(row_index, col_index) = 1.d0


            ! FIRST ORDER TERMS:
             
            do non_col_1 = 1,stencil_size 
               col_index = col_index+1
               V(row_index,col_index) = u_n(non_col_1)
            end do

            ! SECOND ORDER TERMS: 6x(3x3x3) = 162 (GIVES A TOTAL OF 243 TERMS)

            do non_col_1 = 1, stencil_size
            do non_col_2 = non_col_1, stencil_size
               col_index = col_index+1
               V(row_index,col_index) = u_n(non_col_1) * u_n(non_col_2)
            end do
            end do

         T(row_index,:) = T_ij(:,i_box,j_box,k_box) !Change 1 to (1-6) here. !THIS ONE IS CORRECT; KEEP IT.


      end do
      end do
      end do ! BOUNDING BOX

      print*,row_index,col_index
      ! CHECK V:
      if (V(1500,2000).ne.2.0009431419772586d-2) then
         print*, "Error! Check sorting order in  V matrix!"
         print*, 'V(1500,2000)', V(1500,2000)
         stop
      else 
         print*,'V matrix check ... Passed'

      end if


      !CHECK T:
      if (T(3,1).ne.8.7759832493259110d-2)then
         print*, "Error! Check sorting order in  T vector!"
         print*,T(3,1)
         stop
      else 
         print*,'T vector check ... Passed'
      end if
      !print*, 'T(1,2)', T(2)
      !print*, 'T(1,3)', T(3)
      print*, 'T(1,4)', T(1,4)
      !print*, 'T(1,5)', T(5)
      

      if (computeVolterra) then
        
!         call SVD(V,T,h_ij,lambda,printval) ! TSVD
         call LU(V,T,h_ij,lambda)           !Damped least squares
         ! PRINT h_ij:                                                                                                                
         print*, 'h_ij(1,1):',h_ij(1,1)
         print*, 'h_ij(20,1):', h_ij(20,1)
         print*, 'h_ij(350,1):', h_ij(350,1)

stop
         ! CHECK h_ij:
         if (h_ij(350,1).ne.-4.5121154730201521d-2)then
            print*, "Error! Check lambda, method or sorting order for h_ij computation:"
            print*,h_ij(350,1)
            stop
         else 
            print*,'SVD check ... Passed'
         end if
      end if

      ! SAVE/READ h_ij:
      open(1,file='./run/apples/h_ij.dat')
      read(1,*) h_ij
      print*,'h_ij(350,1) from file: ',h_ij(350,1)
      close(1)
end do
end do
end do ! test


return
end subroutine testNonCO



subroutine optimizedTij(u_f,u_t,h_ij,TijOpt,tau_ijOpt)

! ARGUMENTS:
real(8), dimension(:,:,:,:), intent(in) :: u_f, u_t
real(8), dimension(:,:),intent(in) :: h_ij

! LOCAL VARIABLES:
real(8),allocatable,dimension(:,:,:,:),intent(out) :: TijOpt, tau_ijOpt !Computed stresses

integer                :: col_index
integer                :: u_comp, uu_comp ! Indices to select u_i,u_ij components


! FOR NON-COLOCATED FORMULATION:
integer,parameter :: stencil_size = 3*(3*3*3) !3 * 27 = 81
real(8),dimension(stencil_size) :: u_f_stencil, u_t_stencil
integer :: non_col_1, non_col_2  ! To run non coloc combinantions   

! LIST OF INDICES FOR LOOPS:
integer :: i_test,    j_test,    k_test ! for a point in the TEST-scale field
integer :: i_opt,     j_opt,     k_opt  ! for a point in bounding BOX (stencil-center)
integer :: i_stencil, j_stencil, k_stencil


allocate (TijOpt(1,grid,grid,grid))
allocate (tau_ijOpt(1,grid,grid,grid))
! WHOLE DOMAIN COMPUTATION: 
do i_test = 129, 129
do j_test = 129, 129
do k_test = 129, 129

      ! ENTER STENCIL-CENTER POINTS: C-ORDER
      do i_opt = i_test-127, i_test+126
      do j_opt = j_test-127, j_test+126
      do k_opt = k_test-127, k_test+126

         col_index = 0 
         
         ! ENTER 3x3x3 STENCIL: C-ORDER                                       
          u_t_stencil = reshape(u_t(:,i_opt-2:i_opt+2:2, j_opt-2:j_opt+2:2,k_opt-2:k_opt+2:2),(/stencil_size/)) 
          u_f_stencil = reshape(u_f(:,i_opt-1:i_opt+1:1, j_opt-1:j_opt+1:1,k_opt-1:k_opt+1:1),(/stencil_size/)) 

            ! ZERO ORDER TERMS:
             col_index = col_index+1
             TijOpt(:,i_opt,j_opt,k_opt) = h_ij(col_index,:) ! = 1
             tau_ijOpt(:,i_opt,j_opt,k_opt) = h_ij(col_index,:) ! = 1

            ! FIRST ORDER TERMS:             
            do non_col_1 = 1,stencil_size 
               col_index = col_index+1
                TijOpt(:,i_opt,j_opt,k_opt) = TijOpt(:,i_opt,j_opt,k_opt) + (u_t_stencil(non_col_1)*h_ij(col_index,:))
                tau_ijOpt(:,i_opt,j_opt,k_opt) = tau_ijOpt(:,i_opt,j_opt,k_opt) + (u_f_stencil(non_col_1)*h_ij(col_index,:))
            end do

            ! SECOND ORDER TERMS: 6x(3x3x3) = 162 (GIVES A TOTAL OF 243 TERMS)
            do non_col_1 = 1, stencil_size
            do non_col_2 = non_col_1, stencil_size
               col_index = col_index+1
               TijOpt(:,i_opt,j_opt,k_opt) = TijOpt(:,i_opt,j_opt,k_opt) &
                                           + (u_t_stencil(non_col_1) * u_t_stencil(non_col_2)*h_ij(col_index,:))
               tau_ijOpt(:,i_opt,j_opt,k_opt) = tau_ijOpt(:,i_opt,j_opt,k_opt) &
                                           + (u_f_stencil(non_col_1) * u_f_stencil(non_col_2)*h_ij(col_index,:))
            end do
            end do
      end do
      end do
      end do ! BOUNDING BOX
end do
end do
end do ! test

! WRITE COMPUTED STRESS:
open(1,file='./run/apples/T_ijOpt.dat')
open(2,file='./run/apples/tau_ijOpt.dat')
write(1,*) TijOpt(1,:,:,129)
write(2,*) tau_ijOpt(1,:,:,129)
close(1)
close(2)
print*,'TijOpt(1,15,24,10):',TijOpt(1,15,24,10)
print*,'tau_ijOpt(1,15,24,10):',tau_ijOpt(1,15,24,10)
return
end subroutine optimizedTij



!****************************************************************
!                                LU                             !
!****************************************************************
subroutine LU(V, T_ij, h_ij, lambda)
  implicit none
 
 ! dgesv destroys the original matrix. So V is copied into 'a' matrix. This is the LU product matrix.
 ! dgesv returns x vector. First it copies the values from 'T_ij' vector and then computes 'h_ij'.
 ! 'h_ij' is the h_ij vector.
 
  ! ARGUMENTS:
  real(8), dimension(:,:), intent(in) :: V 
  real(8), dimension(:,:), intent(in) :: T_ij
  real(8), dimension(:,:), intent(out) :: h_ij
  real(8), intent(in) :: lambda
 
 ! DGESV ARGUMENTS:
  integer, parameter        :: LDA = N, LDB = N, nrhs = P
  real(8), dimension(:,:),allocatable :: A,VT ! EYE - IDENTITY MATRIX
  real(8), dimension(:,:),allocatable :: b
  integer, dimension(:), allocatable  :: ipiv
  integer                   :: info
  real(8) :: tic, toc
  integer :: i
 
  !DGEMM ARGUMENTS:
  real(8):: alpha=1.d0, beta=0.d0

 allocate (A(N,N), b(N,P), ipiv(N))

! VT = transpose(V)
 call cpu_time(tic)
 call dgemm('T','N',N,N,M, alpha, V,M, V,M, beta, A,N)
 call cpu_time(toc)
 print*,'Elapsed time', toc-tic


!A
 forall(i=1:N) A(i,i) = A(i,i) + lambda 
 print*,'A:'
 call printmatrix(A,size(A,dim=1)) 
!b
 call dgemm('T','N',N,P,M, alpha, V,M, T_ij,M, beta, b,N)
 print*,'b:'
 call printmatrix(b,size(b,dim=1))
 
call cpu_time(tic)
!$omp parallel 
call DGESV(N, nrhs, A, LDA, ipiv, b, LDB, info) !A is a symmetric matrix

!$omp end parallel 
call cpu_time(toc)
print*,'Elapsed time:', toc-tic

h_ij = b

deallocate (A,b,ipiv)
return
end subroutine LU

!****************************************************************
!                             SVD                               !
!****************************************************************

subroutine SVD(A,T,h_ij,lambda,printval)

! ARGUMENTS:
real(8),dimension(:,:),intent(in)  :: A
real(8),dimension(:,:),intent(in)  :: T
real(8),intent(in)                 ::lambda
real(8),dimension(:,:),intent(out) :: h_ij


! DGESVD ARGUMENTS:                                                                                                
integer, parameter :: LDA = M
integer, parameter :: LDU = M
integer, parameter :: LDVT = N
integer :: LWMAX
integer :: i,j, info, LWORK

real(8),dimension(:,:),allocatable :: U
real(8),dimension(:,:),allocatable :: VT 
real(8),dimension(:,:),allocatable :: D, Vinv
real(8),dimension(:), allocatable :: S, work


logical :: printval


LWMAX = M*N
allocate (U(LDU,M), VT(N,N), S(N), D(N,M), Vinv(N,M), work(LWMAX))


if (printval) then
print*,'A'                      ! A MATRIX
call printmatrix(A,size(A,dim=1))
end if


! CALL SVD ROUTINE:
LWORK = -1
call dgesvd('All','All', M,N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, info)
LWORK = min( LWMAX, int(WORK(1)) ) 
call dgesvd('All','All',M,N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, info) !BOTTLENECK = 850 SECS

! Convergence check:
if (info.gt.0) then
   print*, 'Failed to converge'
   stop
end if


if (printval) then
print*,'A'                      ! A MATRIX
call printmatrix(A,size(A,dim=1))
print*,'U'                      ! U MATRIX
call printmatrix(U,size(U,dim=1))
print*,'S'                      ! S VECTOR
print*, S(1:3)
print*,'VT'                     ! VT MATRIX
call printmatrix(VT,size(VT,dim=1))
end if


! COMPUTE PSEUDOINVERSE:
D = 0.d0
forall(i=1:N) D(i,i) = S(i) / (S(i)**2 + lambda**2)
Vinv = matmul(matmul(transpose(VT),D),transpose(U)) !BOTTLENECK = 1300 SECS

if(printval) then
print*,'D'                     ! D MATRIX - DIAG MATRIX
call printmatrix(D,size(D,dim=1))


print*,'Vinv'                     ! Vinv MATRIX - DIAG MATRIX
call printmatrix(Vinv,size(Vinv,dim=1))
end if


h_ij = matmul(Vinv,T)
return
end subroutine SVD


!****************************************************************
!                     SUBROUTINE PRINTMATRIX                    !            
!****************************************************************
subroutine printmatrix(A,LDA)
implicit none

! ARGUMENTS:
real(8), dimension(:,:) :: A
integer :: LDA,width
integer :: i,j


if (LDA.gt.8) then
   LDA=8
   width=8
end if

if (size(A,dim=2).lt.8) width = size(A,dim=2)

do i=1,LDA
   print*, A(i,1:width)
end do

return
end subroutine printmatrix



end program apples




