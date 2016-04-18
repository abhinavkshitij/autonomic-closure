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


! Define velocities:
real(8),allocatable,dimension(:,:,:,:) :: u, u_f, u_t
real(8),allocatable,dimension(:,:,:,:) :: tau_ij,T_ij
real(8),allocatable,dimension(:,:,:) :: LES,test

real(8) :: pre_cut, post_cut

integer :: n_u, n_uu

character(50):: CUT_DATA = '../derived_data/cutout-valid/jhu/' !Change bin4020 by 'sed' in shell script
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


! FILTER VELOCITIES:
print*,'Filter velocities ... '
allocate(u_f(n_u,GRID,GRID,GRID))
allocate(u_t(n_u,GRID,GRID,GRID))

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
call computeStress(u,u_f,u_t,tau_ij,T_ij,n_u,n_uu,LES,test)
call cpu_time(toc)
write(*,307),'computeStress - time elapsed:',toc-tic
deallocate(LES,test)


! CHECK STRESS: (based on T_ij)
if (T_ij(1,15,24,10).ne.-5.2544371578038401d-3) then
         print*, 'Precision test for stress: Failed'
         print*, 'Check precision or data for testing is not JHU 256'
         print*, T_ij(1,15,24,10)
         stop
      else
         print*, 'Precision test for FFT: Passed'
      end if

! ! COMPUTE STRAIN RATE:
! call computeSij(u_f, Sij_f)
! call computeSij(u_t, Sij_t)


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


! WRITE TO FILES:
if(debug(3)) then
   print*,'Write cutout to .dat file ... '
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

else

   print*,'Write u_f, u_t, tau_ij, T_ij to  .bin files ... '
   open(1,file=trim(CUT_DATA)//trim(f_CUT)//'u_f.bin',form='unformatted',action='write') 
   open(2,file=trim(CUT_DATA)//trim(f_CUT)//'u_t.bin',form='unformatted',action='write') 
   open(3,file=trim(CUT_DATA)//trim(f_CUT)//'tau_ij.bin',form='unformatted',action='write') 
   open(4,file=trim(CUT_DATA)//trim(f_CUT)//'T_ij.bin',form='unformatted',action='write') 

   write(1) u_f
   write(2) u_t
   write(3) tau_ij
   write(4) T_ij
   
   close(1)
   close(2)
   close(3)
   close(4)

end if

! COMPUTE STRESS:
!call synStress(u_f,u_t,tau_ij,T_ij,n_u,n_uu)
!print*,'tau_ij',tau_ij(1,testLower,testLower,testLower) !! Check if input arg get altered 

stop

end program main
