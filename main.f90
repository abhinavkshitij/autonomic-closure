program main
  
! STATUS : Integrated with readfile.f90 as a subroutine. 9/11/2015 (Passed)
!          Integrated with two .mod files. 9/13/2015.(Passed)
!          Incorporated HDF5 read capability. 10/22/2015(Passed)
! Result : Passed 
! Notes  : 

use fourier
use fileio
use linsolve

implicit none
integer,parameter           :: LES_scale=40, test_scale=20
character(3) :: stress = 'dev'


! Define velocities:
real(8),allocatable,dimension(:,:,:,:) :: u, u_f, u_t
real(8),allocatable,dimension(:,:,:,:) :: tau_ij,T_ij
real(8),allocatable,dimension(:,:,:) :: LES,test

integer :: n_u, n_uu

character(50):: CUT_DATA = '../derived_data/cutout64/jhu/' !Change bin4020 by 'sed' in shell script
character(10):: f_CUT 
character(3) :: d_set = 'jhu'   ! for binRead()


! DEBUG FLAGS:
logical :: debug(2) = (/1,1/)
! 1- To select only one velocity component and 3 velocity products.
! 2- Choose between NRL/JHU256 database[1] or JHU(HDF5) database[0]


real :: tic, toc
real(8):: pre_cut, post_cut

integer                     :: i,j,k,d=0

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
allocate(u(n_u,GRID,GRID,GRID))

fileSelect:if (debug(2)) then
   call binRead(u,  d_set,  DIM=n_u)
   write(f_CUT,'(a3,2(i2),a1)') 'bin',LES_scale,test_scale,'/' !Write dirname
   !print*, trim(CUT_DATA)//trim(f_CUT) !Test pathname
else
   call hdf5Read() 
end if fileSelect



!! Create filters:
allocate(LES(GRID,GRID,GRID))
allocate(test(GRID,GRID,GRID))
call createFilter(LES,LES_scale)
call createFilter(test,test_scale)
call fftshift(LES)
call fftshift(test)


! FILTER VELOCITIES:
print*,'Filter velocities ... '
allocate(u_f(n_u,GRID,GRID,GRID))
allocate(u_t(n_u,GRID,GRID,GRID))

filter:do i=1,n_u
   u_f(i,:,:,:) = sharpFilter(u(i,:,:,:),LES) ! Speed up this part -- Bottleneck
   u_t(i,:,:,:) = sharpFilter(u_f(i,:,:,:),test) ! Speed up this part -- Bottleneck
end do filter

print*,precision(-0.48241021987284982d0),precision(u_t(1,15,24,10))
! CHECK FFT:(based on u_t)
if (u_t(1,15,24,10).ne.dble(-0.48241021987284982d0)) then
         print*, 'Precision test for FFT: Failed'
         print*, 'Check precision or data for testing is not JHU 256'
         print*, u_t(1,15,24,10)
         stop
      else
         print*, 'Precision test for FFT: Passed'
      end if


! ASSERT 6 COMPONENTS FOR ij TO COMPUTE STRESS:
if (n_uu.ne.6) then
   print*,"Need all 6 ij components to compute stress... Aborting"
   stop
end if


! COMPUTE STRESS:
allocate(tau_ij(n_uu,GRID,GRID,GRID))
allocate(T_ij(n_uu,GRID,GRID,GRID))
call cpu_time(tic)
print*,'Compute stress:',stress
call computeStress(u,u_f,u_t,tau_ij,T_ij,n_u,n_uu,LES,test,stress)
call cpu_time(toc)
write(*,307),'computeStress - time elapsed:',toc-tic
deallocate(LES,test)


! CHECK STRESS: (based on T_ij)

if (T_ij(1,15,24,10).ne.-5.2544371578038401E-003) then
         print*, 'Precision test for stress: Failed'
         print*, 'Check precision or data for testing is not JHU 256'
         print*, T_ij(1,15,24,10)
         stop
      else
         print*, 'Precision test for FFT: Passed'
   end if



! CUTOUT:
pre_cut = u_t(1,testLower+lBound-1,testLower+lBound-1,testLower+lBound-1)
call cutout(u_f,   n_u )
call cutout(u_t,   n_u )
call cutout(tau_ij,n_uu)
call cutout(T_ij,  n_uu)
post_cut = u_t(1,testLower,testLower,testLower)
if (pre_cut.ne.post_cut) then
   print*,"Cutout: Error in cutout!"
   stop
   else
      print*,"Cutout: Passed"
   end if

stop
! WRITE VELOCITIES & STRESSES TO FILES:
print*,'Write cutout to file ... '
open(1,file=trim(CUT_DATA)//trim(f_CUT)//'u_f.dat')
open(2,file=trim(CUT_DATA)//trim(f_CUT)//'u_t.dat')
open(3,file=trim(CUT_DATA)//trim(f_CUT)//'tau_ij.dat')
open(4,file=trim(CUT_DATA)//trim(f_CUT)//'T_ij.dat')

 write(1,*) u_f
 write(2,*) u_t
 write(3,*) tau_ij
 write(4,*) T_ij

close(1)
close(2)
close(3)
close(4)

stop

! COMPUTE STRESS:
call synStress(u_f,u_t,tau_ij,T_ij,n_u,n_uu)
!print*,'tau_ij',tau_ij(1,testLower,testLower,testLower) !! Check if input arg get altered 


end program main
