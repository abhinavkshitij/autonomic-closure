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
  logical :: useTestData          =  0
  logical :: readFile             =  0
  logical :: filterVelocities     =  0
  logical :: plot_Velocities      =  0
  logical :: computeFFT_data      =  0
  logical :: save_FFT_data        =  1
  logical :: plot_Stresses        =  1
  logical :: production_Term      =  1
  logical :: save_ProductionTerm  =  0
  logical :: computeVolterra      =  1

  !
  !    ..TEST VARAIBLES FOR DEBUGGING..
  real(8) :: error
  real(8) :: test_lam = 1d-11
  character(64) :: test_ch

  
  write(test_ch, '(ES6.0E2)') test_lam
  write(*, '(ES6.0E2)') test_lam
  print*, trim(test_ch(4:6))


  call setEnv()
  !
  !    ..INIT POSTPROCESSING..
  open(22, file = trim(RES_DIR)//'path.txt')


 
  ! FORMAT:
3015 format(a30,f22.15)
307 format(a30,f22.7)
507 format(a50,f22.7)
  !call system('clear')
  call printParams()


  !! Set debug flags for velocity components:
  if (useTestData) then
     n_u = 1
     n_uu = 3
     print*, 'Debug mode for velocity components... \n'
  end if


  ! 1] LOAD DATASET: ALTERNATE METHOD STASHED ABOVE. +
  if(readFile) call readData(DIM = n_u)
  if(allocated(u).eqv..false.) allocate(u(n_u,i_GRID,j_GRID,k_GRID))

  ! GET STATISTICS OF INITIAL VELOCITY
  !  call createPDF(u(1,:,:,:),plot=.true.,fieldname='Velocities')


  ! INITIALIZE PATH: LEVEL 1 (DATASET)
  write(22,*) trim(DATA_PATH)
  TEMP_PATH = trim(TEMP_DIR)//trim(dataset)//'/'
  RES_PATH =  trim(RES_DIR)//trim(dataset)//'/'
  if (dataset.eq.'hst') then
     TEMP_PATH = trim(TEMP_PATH)//trim(hst_set)//'/'
     RES_PATH = trim(RES_PATH)//trim(hst_set)//'/'
  end if
  write(22,*) RES_PATH
  call system ('mkdir -p '//trim(TEMP_PATH))
  call system ('mkdir -p '//trim(RES_PATH))

  ! PLOT Energy spectra:
  call energySpectra(u(1:3,:,:,:))


  ! ************** LEVEL 2 ****************!
  !
  ! ADD PATH DEPTH : (SCALE)
  write(scale,'(2(i0))') LES_scale, test_scale 
  TEMP_PATH = trim(TEMP_PATH)//'bin'//trim(scale)//'/'
  RES_PATH =  trim(RES_PATH)//'dat'//trim(scale)//'/'
  write(22,*) RES_PATH
  call system ('mkdir -p '//trim(TEMP_PATH))
  call system ('mkdir -p '//trim(RES_PATH))



  ! 2] FILTER VELOCITIES [AND PRESSURE]:
  allocate(u_f (n_u, i_GRID,j_GRID,k_GRID))
  allocate(u_t (n_u, i_GRID,j_GRID,k_GRID))


  if (filterVelocities) then
     print*
     write(*, '(a32)', ADVANCE='NO'), 'Filter velocities ... '

     !! Create filters:
     allocate(LES (f_GRID,f_GRID,f_GRID))
     allocate(test(f_GRID,f_GRID,f_GRID))
     call createFilter(LES,LES_scale)
     call createFilter(test,test_scale)
     call fftshift(LES)
     call fftshift(test)

     do i=1,n_u
        u_f(i,:,:,:) = sharpFilter(u(i,:,:,:),LES) ! Speed up this part -- Bottleneck
        u_t(i,:,:,:) = sharpFilter(u_f(i,:,:,:),test) ! Speed up this part -- Bottleneck
     end do 
     call check_FFT(u_t(1,15,24,10))  
     print*, 'OK'
     if (plot_Velocities) call plotVelocities()
  end if

  ! CHECK ALL 6ij:
  if (n_uu.ne.6) then
     print*, 'Cannot proceed further without 6 ij'
     stop
  end if


  ! 3] GET FFT_DATA:
  allocate (tau_ij (n_uu,i_GRID,j_GRID,k_GRID))
  allocate (T_ij   (n_uu,i_GRID,j_GRID,k_GRID))
  allocate (Sij_f  (3,3, i_GRID,j_GRID,k_GRID))
  allocate (Sij_t  (3,3, i_GRID,j_GRID,k_GRID))


  if(computeFFT_data)then

     ! STRESS:
     print*,'Compute stress:',stress
     call computeStress(u, u_f, u_t, tau_ij, T_ij, LES, test)
     deallocate(LES,test)

     if (save_FFT_DATA) call saveFFT_data()

  else

     ! LOAD SAVED FFT_DATA: Filtered velocities, stress and strain rates 
     call loadFFT_data()
!     call checkFFT_data()

     ! CHECK INPUT DATA:
     ! ++
     print*, 'check input'
     if (withPressure) then
     if (u_f(4,15,24,10).ne.0.12164158031779296d0) then
        print*, 'ERROR READING DATA'
        print*, u_f(4,15,24,10)
        stop
     end if
     end if

     if (u_t(1,15,24,10).ne.-0.48241021987284982d0) then
        print*, 'ERROR READING DATA'
        print*, u_t(1,15,24,10)
        stop
     elseif(T_ij(1,15,24,10).ne.-5.2544371578038401d-3) then
        print*, 'ERROR COMPUTING T_ij'
        print*,'T_ij(1,15,24,10)',T_ij(1,15,24,10)
        stop
     else
        print*, 'Read data saved from main.f90: Passed \n'
     end if

  end if

  ! PLOT STRESSES - ORIGINAL ON LEVEL 2
  if (plot_Stresses) call plotOriginalStress()

  
  ! STRAIN:
  print*, 'Compute strain rate \n'
  call computeSij(u_f, Sij_f)
  call computeSij(u_t, Sij_t)

  print*,'Sij_f(1,2,15,24,10)',Sij_f(1,2,15,24,10)
  print*,'Sij_f(1,2,15,24,10)',Sij_t(1,2,15,24,10),'\n'     



!  5] PRODUCTION FIELD - ORIGINAL 
  if(production_Term) then
     allocate (Pij_f (i_GRID, j_GRID, k_GRID))
     allocate (Pij_t (i_GRID, j_GRID, k_GRID))
     call productionTerm(Pij_f, tau_ij, Sij_f)
     call productionTerm(Pij_t, T_ij,   Sij_t)
     print*,'Pij_f(15,24,10)',Pij_f(15,24,10)
     print*,'Pij_t(15,24,10)',Pij_t(15,24,10),'\n'
     if (save_ProductionTerm) then
        call plotProductionTerm()
     end if
  end if

stop

  ! ************** LEVEL 3 ****************!
  !
  ! ADD PATH DEPTH : (METHOD) - LU or SVD
  TEMP_PATH = trim(TEMP_PATH)//trim(solutionMethod)//'/'
  RES_PATH =  trim(RES_PATH)//trim(solutionMethod)//'/'
  call system ('mkdir -p '//trim(TEMP_PATH))
  call system ('mkdir -p '//trim(RES_PATH))
  write(22,*) RES_PATH


  ! 6] AUTONOMICALLY TUNED LAMBDA
  allocate(h_ij(N,P))
  if (allocated(Pij_f).eqv..false.)  allocate (Pij_f (i_GRID, j_GRID, k_GRID))
  if (allocated(Pij_t).eqv..false.)  allocate (Pij_t (i_GRID, j_GRID, k_GRID))

  fileID = 81
  open(fileID,file=trim(RES_PATH)//'crossValidationError.csv')
  do iter = 1, n_DAMP
     lambda = 1.d-10 * 10**(iter-1)

     call autonomicClosure(u_f, u_t, tau_ij, T_ij, h_ij)
     ! OPTIMIZED STRESS:
     if (allocated(T_ijOpt).eqv..false.)    allocate(T_ijOpt   (n_uu,i_GRID,j_GRID,k_GRID))
     if (allocated(tau_ijOpt).eqv..false.)  allocate(tau_ijOpt (n_uu,i_GRID,j_GRID,k_GRID))
     call computedStress(u_f, u_t, h_ij, T_ijOpt, tau_ijOpt)
     ! PLOT COMPUTED STRESSES 
     if (plot_Stresses) call plotComputedStress(lambda)
     call trainingerror(T_ijOpt, T_ij, error,'plot',fileID)


 ! 7] PRODUCTION FIELD FROM COMPUTED STRESSES:
  if(production_Term) then
     call productionTerm(Pij_f, tau_ijOpt, Sij_f)
     call productionTerm(Pij_t, T_ijOpt,   Sij_t)
     if (save_ProductionTerm) call plotProductionTerm(lambda)
  end if


     ! GET STATISTICS:
     !  print*, 'testError'
     !  call createPDF(,plot=.true.,fieldname='errorTest')
     !  call createPDF(,plot=.true.,fieldname='errorSGS')
  end do
  close(fileID)

  if (allocated(Sij_f) ) deallocate (Sij_f, Sij_t)


end program autonomic
