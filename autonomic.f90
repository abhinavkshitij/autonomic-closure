!****************************************************************
!                            AUTONOMIC
!****************************************************************

!----------------------------------------------------------------
! USE: 1)
!      
!      
!      
!
! FORM: program autonomic
!       
!       
!       
!       
!       
!
! BEHAVIOR: 
!           
! DEBUG FLAGS:
!     useTestData:  Use only one velocity component and 3 velocity products
!     writeStressBin:  Write stress binary file [0]
!
!           
!
! STATUS : 
!    +  allocate(u (n_u, i_GRID,j_GRID,k_GRID))
!    +  call readData(u, DIM=n_u)
!   ++  call printplane(u_f(1,:,:,:),frameLim=4)
!----------------------------------------------------------------

program autonomic
  
  use fourier
  use actools
  use solver

  implicit none


  character(64) :: filename
  character(10) :: scale
  !
  !    ..CONTROL SWITCHES..
  logical :: useTestData      =  0
  logical :: readFile         =  0
  logical :: filterVelocities =  0
  logical :: computeOrigStress   =  0
  logical :: save_FFT_data    =  0
  logical :: plot_FFT_data    =  1
  logical :: computeStrain    =  0
  logical :: computeVolterra  =  0


  ! FORMAT:
3015 format(a30,f22.15)
307 format(a30,f22.7)
507 format(a50,f22.7)
  !call system('clear')
  !call printParams()

  !! Set debug flags for velocity components:
  if (useTestData) then
     n_u = 1
     n_uu = 3
     print*, 'Debug mode for velocity components...'
  end if


  ! LOAD DATASET: ALTERNATE METHOD STASHED ABOVE. +
  if(readFile) call readData(DIM=n_u)

  
  ! FILTER VELOCITIES:

  ! SET READ/WRITE PATH FOR FFT_DATA:
  write(scale,'(2(i2))') LES_scale, test_scale 
  TEMP_PATH = trim(TEMP_DIR)//trim(d_set)//'/'//'bin'//trim(scale)//'/'
  RES_PATH =  trim(RES_DIR)//trim(d_set)//'/'//'dat'//trim(scale)//'/'


  allocate(u_f (n_u, i_GRID,j_GRID,k_GRID))
  allocate(u_t (n_u, i_GRID,j_GRID,k_GRID))
  if (filterVelocities) then
     print*
     print*,'Filter velocities ... '

     !! Create filters:
     allocate(LES(f_GRID,f_GRID,f_GRID))
     allocate(test(f_GRID,f_GRID,f_GRID))
     call createFilter(LES,LES_scale)
     call createFilter(test,test_scale)
     call fftshift(LES)
     call fftshift(test)

     filter:do i=1,n_u
        u_f(i,:,:,:) = sharpFilter(u(i,:,:,:),LES) ! Speed up this part -- Bottleneck
        u_t(i,:,:,:) = sharpFilter(u_f(i,:,:,:),test) ! Speed up this part -- Bottleneck
     end do filter
     call check_FFT(u_t(1,15,24,10))
  end if

  ! ASSERT 6 COMPONENTS FOR ij TO COMPUTE STRESS:
  if (n_uu.ne.6) then
     print*,"Need all 6 ij components to compute stress... Aborting"
     stop
  end if

  allocate(tau_ij(n_uu,i_GRID,j_GRID,k_GRID))
  allocate(T_ij(n_uu,i_GRID,j_GRID,k_GRID))

  ! COMPUTE STRESS:
  if(computeOrigStress)then

     call cpu_time(tic)
     print*
     print*,'Compute stress:',stress
     call computeStress(u, u_f, u_t, tau_ij, T_ij, LES, test)
     call cpu_time(toc)
     write(*,307),'computeStress - time elapsed:',toc-tic
     deallocate(LES,test)
     if (save_FFT_DATA) then 
        call saveFFT_data()
     end if

  else

     ! LOAD SAVED FFT_DATA : Filtered velocities and stress files: 
     call loadFFT_data()

     ! CHECK INPUT DATA:
     ! ++
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
  
  if (plot_FFT_data) call plotFFT_data() 
  


stop

  ! COMPUTE STRAIN RATE:
  if(computeStrain)then
     allocate (Sij_f(n_u,n_u,i_GRID,j_GRID,k_GRID))
     allocate (Sij_t(n_u,n_u,i_GRID,j_GRID,k_GRID))
     print*
     print*, 'Compute strain rate'
     call computeSij(u_f, Sij_f)
     call computeSij(u_t, Sij_t)

     print*,'Sij_fl(1,2,15,24,10)',Sij_f(1,2,15,24,10)
     print*,'Sij_ft(1,2,15,24,10)',Sij_t(1,2,15,24,10)

     deallocate (Sij_f, Sij_t)
  end if
  
  ! COMPUTE h_ij using autonomic closure:
  if (allocated(h_ij).neqv..true.) allocate(h_ij(N,P))
  call autonomicClosure(u_f, u_t, tau_ij, T_ij, h_ij)

  
  ! PRINT h_ij:
  print*, 'h_ij(1,1):',h_ij(1,1)
  print*, 'h_ij(20,1):', h_ij(20,1)
  print*, 'h_ij(350,1):', h_ij(350,1)


  ! SAVE/READ h_ij: in '../run/apples/h_ij.dat'


  ! CHECK h_ij:
  if (h_ij(350,1).ne.-4.5121154730201521d-2)then
     print*, "Error! Check lambda, method or sorting order for h_ij computation:"
     print*,h_ij(350,1)
     stop
  else 
     print*,'SVD check ... Passed'
  end if


  if (allocated(TijOpt).neqv..true.) allocate(TijOpt(n_uu,i_GRID,j_GRID,k_GRID))
  if (allocated(tau_ijOpt).neqv..true.) allocate(tau_ijOpt(n_uu,i_GRID,j_GRID,k_GRID))

  
  call cpu_time(tic)
  call optimizedTij(u_f, u_t, h_ij, TijOpt, tau_ijOpt)
  call cpu_time(toc)
  print*,'Elapsed time', toc-tic


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

end program autonomic



