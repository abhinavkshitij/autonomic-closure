!****************************************************************
!                              MAIN
!****************************************************************

!----------------------------------------------------------------
! USE: 1) Main program for autonomic closure.
!      2) Reads datasets and saves filtered datasets in temp/ 
!      
!      
!
! FORM: program main
!       
!       
!       
!       
!       
!
! BEHAVIOR: Depends on actools and solver
!           
!
! STATUS : 
! 
!----------------------------------------------------------------

program main
  
  use actools
  use solver

  implicit none
  

  real(8) :: pre_cut, post_cut
  character(50):: CUT_DATA = '../temp/cutout-valid/jhu/' !Change bin4020 by 'sed' in shell script
  character(10):: f_CUT 

  ! DEBUG FLAGS:
  ! 1- To select only one velocity component and 3 velocity products.
  ! 2- Write as binary file [0]

  logical :: debug(2) = [1,0]

  ! FORMAT:
3015 format(a30,f22.15)
307 format(a30,f22.7)
507 format(a50,f22.7)
  !call system('clear')
  call printParams()


  !! Set debug flags for velocity components:
  if (debug(1)) then
     n_u = 1
     n_uu = 3
     print*, 'Debug mode for velocity components...'
  end if


  !! Select file to read:
  allocate(u(n_u,i_GRID,j_GRID,k_GRID))

  call readData(u, DIM = n_u)
  write(f_CUT,'(a3,2(i2),a1)') 'bin',LES_scale,test_scale,'/' !Write dirname
  call printplane(u(1,:,:,:),frameLim=4)


  !! Create filters:
  allocate(LES(f_GRID,f_GRID,f_GRID))
  allocate(test(f_GRID,f_GRID,f_GRID))
  call createFilter(LES,LES_scale)
  call createFilter(test,test_scale)
  call fftshift(LES)
  call fftshift(test)


  ! FILTER VELOCITIES AND CHECK:
  print*,'Filter velocities ... '
  allocate(u_f(n_u,i_GRID,j_GRID,k_GRID))
  allocate(u_t(n_u,i_GRID,j_GRID,k_GRID))

  call cpu_time(tic)
  filter:do i=1,n_u
     u_f(i,:,:,:) = sharpFilter(u(i,:,:,:),LES) ! Speed up this part -- Bottleneck
     u_t(i,:,:,:) = sharpFilter(u_f(i,:,:,:),test) ! Speed up this part -- Bottleneck
  end do filter
  call check_FFT(u_t(1,15,24,10))
  call cpu_time(toc)
  print*, 'Elapsed time FFT:' ,toc-tic

  ! ASSERT 6 COMPONENTS FOR ij TO COMPUTE STRESS:
  if (n_uu.ne.6) then
     print*,"Need all 6 ij components to compute stress... Aborting"
     stop
  end if


  ! COMPUTE STRESS AND CHECK:
  allocate(tau_ij(n_uu,i_GRID,j_GRID,k_GRID))
  allocate(T_ij(n_uu,i_GRID,j_GRID,k_GRID))
  call cpu_time(tic)
  print*,'Compute stress:',stress
  call computeStress(u,u_f,u_t,tau_ij,T_ij,n_u,n_uu,LES,test)
  call cpu_time(toc)
  write(*,307),'computeStress - time elapsed:',toc-tic
  deallocate(LES,test)




  ! COMPUTE STRAIN RATE:
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

  stop
  ! WRITE TO FILES:
  if(debug(2)) then
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
  !call autonomicClosure(u_f,u_t,tau_ij,T_ij,n_u,n_uu)
  !print*,'tau_ij',tau_ij(1,testLower,testLower,testLower) !! Check if input arg get altered 

  stop

end program
