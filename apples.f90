program apples

use actools
use linsolve

implicit none


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

character(50):: CUT_DATA = '../temp/cutout-valid/jhu/' !Change bin4020 by 'sed' in shell script
character(64) :: filename
character(10):: f_CUT 


! DEBUG FLAGS:
! 1- To select only one velocity component and 3 velocity products.
! 2- Choose between NRL/JHU256 database[1] or JHU(HDF5) database[0]
logical :: debug(2) = [0,2]

logical :: computeStresses = 0
logical :: computeStrain = 0
logical :: filterVelocities = 0
logical :: readFile = 0
logical :: computeVolterra = 1

real(8) :: pre_cut, post_cut


write(f_CUT,'(a3, 2(i2), a1)') 'bin',LES_scale,test_scale,'/' !Write dirname


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
   call readData(u,  DIM=n_u)
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
   call computeStress(u,u_f,u_t,tau_ij,T_ij,n_u,n_uu,LES,test)
   call cpu_time(toc)
   write(*,307),'computeStress - time elapsed:',toc-tic
   deallocate(LES,test)

else

! READ STRESS files:
   do i = 1,size(var_FFT)
      filename = trim(CUT_DATA)//trim(f_CUT)//trim(var_FFT(i)%name)//'.bin'
      print*, filename
      open(i,file = filename,form='unformatted')

      if (i.eq.1) read(i) u_f
      if (i.eq.2) read(i) u_t
      if (i.eq.3) read(i) tau_ij
      if (i.eq.4) read(i) T_ij
      close(i)
   end do
   
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
   call computeSij(u_f, Sij_f)
   call computeSij(u_t, Sij_t)

   print*,'Sij_fl(1,2,15,24,10)',Sij_f(1,2,15,24,10)
   print*,'Sij_ft(1,2,15,24,10)',Sij_t(1,2,15,24,10)
   print*,'Sij_ft(1,2,1,3,5)',Sij_t(1,2,1,3,5)

   deallocate (Sij_f, Sij_t)
end if


if (allocated(h_ij).neqv..true.) allocate(h_ij(N,P))
print*,shape(h_ij)
! COMPUTE h_ij:
!call testNonCo(u_f,u_t,tau_ij,T_ij,h_ij)
call synStress(u_f,u_t,tau_ij,T_ij,h_ij)


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
stop
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

end program apples



