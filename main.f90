program main
  
! STATUS : Integrated with readfile.f90 as a subroutine. 9/11/2015 (Passed)
!          Integrated with two .mod files. 9/13/2015.(Passed)
!          Incorporated HDF5 read capability. 10/22/2015(Passed)
! Result : Passed 
! Notes  : 

use actools
use linsolve

implicit none


! Define velocities:
real(8),allocatable,dimension(:,:,:,:) :: u, u_f, u_t
real(8),allocatable,dimension(:,:,:,:) :: tau_ij,T_ij
real(8),allocatable,dimension(:,:,:) :: LES,test

real(8) :: pre_cut, post_cut

integer :: n_u, n_uu

character(50):: CUT_DATA = '../temp/cutout-valid/jhu/' !Change bin4020 by 'sed' in shell script
character(10):: f_CUT 

! DEBUG FLAGS:
! 1- To select only one velocity component and 3 velocity products.
! 2- Choose between NRL/JHU256 database[1] or JHU(HDF5) database[0]
! 3- Write as binary file [0]
logical :: debug(3) = [0,1,0]



! FORMAT:
3015 format(a30,f22.15)
307 format(a30,f22.7)
507 format(a50,f22.7)
!call system('clear')
call printParams()


!! Set debug flags for velocity components:
n_u=3; n_uu = 6
if (debug(1)) then
   n_u = 1
   n_uu = 3
   print*, 'Debug mode for velocity components...'
end if


!! Select file to read:
allocate(u(n_u,GRID,GRID,GRID))

call readData(u, DIM = n_u)
write(f_CUT,'(a3,2(i2),a1)') 'bin',LES_scale,test_scale,'/' !Write dirname
call printplane(u(1,:,:,:),frameLim=4)


!! Create filters:
allocate(LES(GRID,GRID,GRID))
allocate(test(GRID,GRID,GRID))
call createFilter(LES,LES_scale)
call createFilter(test,test_scale)
call fftshift(LES)
call fftshift(test)


! FILTER VELOCITIES AND CHECK:
print*,'Filter velocities ... '
allocate(u_f(n_u,GRID,GRID,GRID))
allocate(u_t(n_u,GRID,GRID,GRID))

filter:do i=1,n_u
   u_f(i,:,:,:) = sharpFilter(u(i,:,:,:),LES) ! Speed up this part -- Bottleneck
   u_t(i,:,:,:) = sharpFilter(u_f(i,:,:,:),test) ! Speed up this part -- Bottleneck
end do filter
call check_FFT(u_t(1,15,24,10))


! ASSERT 6 COMPONENTS FOR ij TO COMPUTE STRESS:
if (n_uu.ne.6) then
   print*,"Need all 6 ij components to compute stress... Aborting"
   stop
end if


! COMPUTE STRESS AND CHECK:
allocate(tau_ij(n_uu,GRID,GRID,GRID))
allocate(T_ij(n_uu,GRID,GRID,GRID))
call cpu_time(tic)
print*,'Compute stress:',stress
call computeStress(u,u_f,u_t,tau_ij,T_ij,n_u,n_uu,LES,test)
call cpu_time(toc)
write(*,307),'computeStress - time elapsed:',toc-tic
deallocate(LES,test)




 ! COMPUTE STRAIN RATE:
 call computeSij(u_f, Sij_f)
 call computeSij(u_t, Sij_t)


! CUTOUT:
! pre_cut = u_t(1,testLower+lBound-1,testLower+lBound-1,testLower+lBound-1)
! call cutout(u_f,   n_u )
! call cutout(u_t,   n_u )
! call cutout(tau_ij,n_uu)
! call cutout(T_ij,  n_uu)
! post_cut = u_t(1,testLower,testLower,testLower)
! if (pre_cut.ne.post_cut) then
!    print*,"Cutout: Error in cutout!"
!    stop
!    else
!       print*,"Cutout: Passed"
!    end if

stop
! WRITE TO FILES:
if(debug(3)) then
   print*,'Write cutout to .dat file ... '
   do i = 1,4
      open(i,file = trim(CUT_DATA)//trim(f_CUT)//trim(var_FFT(i)%name)//'.dat')
      if (i.eq.1) write(i,*) u_f
      if (i.eq.2) write(i,*) u_t
      if (i.eq.3) write(i,*) tau_ij
      if (i.eq.4) write(i,*) T_ij
      close(i)
   end do
else
   print*,'Write u_f, u_t, tau_ij, T_ij to  .bin files ... '
   do i = 1,4
      open(i,file = trim(CUT_DATA)//trim(f_CUT)//trim(var_FFT(i)%name)//'.bin',form='unformatted')
      if (i.eq.1) write(i) u_f
      if (i.eq.2) write(i) u_t
      if (i.eq.3) write(i) tau_ij
      if (i.eq.4) write(i) T_ij
      close(i)
   end do
end if

! COMPUTE STRESS:
!call synStress(u_f,u_t,tau_ij,T_ij,n_u,n_uu)
!print*,'tau_ij',tau_ij(1,testLower,testLower,testLower) !! Check if input arg get altered 

stop

end program main
