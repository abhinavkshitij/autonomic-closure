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
!   +   call createPDF(u(1,:,:,:),plot=.true.,fieldname='Velocities')
!   ##  call energySpectra(u(1:3,:,:,:))
!   ++  
!     print*,'Sij_f(1,2,15,24,10)',Sij_f(1,2,15,24,10)
!     print*,'Sij_t(1,2,15,24,10)',Sij_t(1,2,15,24,10),'\n'     
!  +++ 
!     print*,'Pij_f(15,24,10)',Pij_f(15,24,10)
!     print*,'Pij_t(15,24,10)',Pij_t(15,24,10),'\n'
!     
! ### Speed up this part -- Bottleneck  
! 
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
  logical :: readFile             =  1
  logical :: filterVelocities     =  1
  logical :: plot_Velocities      =  0
  logical :: computeFFT_data      =  1 ! **** ALWAYS CHECK THIS ONE BEFORE A RUN**** !
  logical :: save_FFT_data        =  0

  logical :: plot_Stresses        =  0
  logical :: production_Term      =  0
  logical :: save_ProductionTerm  =  0

  integer :: time_index
  real(8) :: error_cross


  if (computeFFT_data.eqv..false.) then
     useTestData      = 0
     readfile         = 0
     filterVelocities = 0
     plot_Velocities  = 0
     save_FFT_data    = 0
  end if

  !    ..INIT POSTPROCESSING..
  open(path_txt, file = trim(RES_DIR)//'path.txt')

  call setEnv()
  call printParams()
  print*, 'Dataset: ', dataset, '\n'


  ! TEST DATA:
  if (useTestData) then
     n_u = 1
     n_uu = 3
     print*, 'Debug mode for velocity components... \n'
  end if


  time_loop: do time_index = time_init, time_final, time_incr


     ! SET TIME PARAMS:
     write(time,'(i0)') time_index !num2str
     if ((len(trim(time))-1).lt.2) time = trim('0')//trim(time)

     ! 1] LOAD DATASET:
     if(allocated(u).eqv..false.)        allocate(u(n_u,i_GRID,j_GRID,k_GRID))
     if(readFile)                        call readData(DIM = n_u)
     if (dataset.eq.'hst')               u(:,:,256:130:-1,:) = u(:,:,2:128,:)


     ! + GET STATISTICS OF INITIAL VELOCITY:

     
     ! ************** LEVEL 1 ****************!
     !
     ! ADD PATH DEPTH: DATASET
     
     TEMP_PATH = trim(TEMP_DIR)//trim(dataset)//'/'
     RES_PATH  = trim(RES_DIR)//trim(dataset)//'/'
     if (dataset.eq.'hst') then
        TEMP_PATH = trim(TEMP_PATH)//trim(hst_set)//'/'
        RES_PATH  = trim(RES_PATH)//trim(hst_set)//'/'
     end if
     write(path_txt,*) trim(DATA_PATH)
     write(path_txt,*) RES_PATH
     call system ('mkdir -p '//trim(TEMP_PATH))
     call system ('mkdir -p '//trim(RES_PATH))

     ! ## PLOT Energy spectra:


     ! ************** LEVEL 2 ****************!
     !
     ! ADD PATH DEPTH : SCALE
     write(scale,'(2(i0))') LES_scale, test_scale 
     TEMP_PATH = trim(TEMP_PATH)//'bin'//trim(scale)//'/'
     RES_PATH =  trim(RES_PATH)//'dat'//trim(scale)//'/'
     write(path_txt,*) RES_PATH
     call system ('mkdir -p '//trim(TEMP_PATH))
     call system ('mkdir -p '//trim(RES_PATH))



     ! 2] FILTER VELOCITIES [AND PRESSURE]:
     if(allocated(u_f).eqv..false.)        allocate(u_f (n_u, i_GRID,j_GRID,k_GRID))
     if(allocated(u_t).eqv..false.)        allocate(u_t (n_u, i_GRID,j_GRID,k_GRID))


     if (filterVelocities) then
        write(*, '(a32)', ADVANCE='NO'), adjustl('        Filter velocities ... ')

        ! CREATE FILTERS:
        allocate(LES (f_GRID,f_GRID,f_GRID))
        allocate(test(f_GRID,f_GRID,f_GRID))
        call createFilter(LES,LES_scale)
        call createFilter(test,test_scale)
        call fftshift(LES)
        call fftshift(test)

        do i=1,n_u
           u_f(i,:,:,:) = sharpFilter(u  (i,:,:,:),LES)
           u_t(i,:,:,:) = sharpFilter(u_f(i,:,:,:),test) 
        end do
        ! ###
        print*, 'Completed'
        call check_FFT(u_t(1,15,24,10))  
        if (plot_Velocities)                                        call plotVelocities()
     end if



     ! BREAKPOINT 1:
     if (useTestData) stop


     ! 3] GET FFT_DATA:
     if(allocated(tau_ij).eqv..false.)     allocate (tau_ij (n_uu,i_GRID,j_GRID,k_GRID))
     if(allocated(T_ij).eqv..false.)       allocate (T_ij   (n_uu,i_GRID,j_GRID,k_GRID))


     ! COMPUTE ORIGINAL STRESS [SAVE]
     if(computeFFT_data) then
        print*,'Compute original stress:',stress
        call computeStress(u, u_f, u_t, tau_ij, T_ij, LES, test)
        deallocate(LES,test)
        if (save_FFT_DATA) call saveFFT_data()
     else
        ! LOAD SAVED FFT_DATA ../temp/ [CHECK]
        call loadFFT_data()
        call checkFFT_data()
     end if
     if (plot_Stresses)                                            call plotOriginalStress()



     if(allocated(Sij_f).eqv..false.)     allocate (Sij_f  (3,3, i_GRID,j_GRID,k_GRID))
     if(allocated(Sij_t).eqv..false.)     allocate (Sij_t  (3,3, i_GRID,j_GRID,k_GRID))
     ! 5] ORIGINAL PRODUCTION FIELD 
     if(production_Term) then

        ! STRAIN RATE:
        print*, 'Compute strain rate \n'
        call computeSij(u_f, Sij_f)
        call computeSij(u_t, Sij_t)
        ! ++ CHECK S_ij

        if(allocated(Pij_f).eqv..false.)          allocate (Pij_f (i_GRID, j_GRID, k_GRID))
        if(allocated(Pij_t).eqv..false.)          allocate (Pij_t (i_GRID, j_GRID, k_GRID))
        call productionTerm(Pij_f, tau_ij, Sij_f)
        call productionTerm(Pij_t, T_ij,   Sij_t)
        ! +++  CHECK P_ij
        if (save_ProductionTerm)                                   call plotProductionTerm()     
        deallocate (Pij_f, Pij_t)
     end if


     ! ************** LEVEL 3 ****************!
     !
     ! ADD PATH DEPTH : (METHOD) - LU or SVD
     TEMP_PATH = trim(TEMP_PATH)//trim(solutionMethod)//'/'
     RES_PATH =  trim(RES_PATH)//trim(solutionMethod)//'/'
     write(path_txt,*) RES_PATH
     call system ('mkdir -p '//trim(TEMP_PATH))
     call system ('mkdir -p '//trim(RES_PATH))



     ! 6] AUTONOMICALLY TUNED LAMBDA
     if(allocated(T_ijOpt).eqv..false.)            allocate (T_ijOpt   (n_uu,i_GRID,j_GRID,k_GRID))
     if(allocated(tau_ijOpt).eqv..false.)          allocate (tau_ijOpt (n_uu,i_GRID,j_GRID,k_GRID))
     if(allocated(h_ij).eqv..false.)               allocate (h_ij      (N,P))

     if (production_Term) then
        if(allocated(Pij_fOpt).eqv..false.)        allocate (Pij_fOpt  (i_GRID, j_GRID, k_GRID))
        if(allocated(Pij_tOpt).eqv..false.)        allocate (Pij_tOpt  (i_GRID, j_GRID, k_GRID))
     end if


     print*, 'Autonomic closure ... '
     open(cross_csv, file=trim(RES_PATH)//trim('crossValidationError')//trim(time)//trim('.csv'))
     do iter = 1, n_lambda
        if (mod(n_lambda,2).eq.1)                   lambda = lambda_0(1) * 10**(iter-1)
        if (mod(n_lambda,2).eq.0)                   lambda = lambda_0(2) * 10**(iter-1)

        call autonomicClosure (u_f, u_t, tau_ij, T_ij, h_ij)
        call computedStress   (u_f, u_t, h_ij, T_ijOpt, tau_ijOpt)
        if (plot_Stresses)                                      call plotComputedStress(lambda)
        call trainingerror(T_ijOpt, T_ij, error_cross,'plot',cross_csv)

        ! 7] PRODUCTION FIELD - COMPUTED 
        if(production_Term) then
           call productionTerm(Pij_fOpt, tau_ijOpt, Sij_f)
           call productionTerm(Pij_tOpt, T_ijOpt,   Sij_t)
           if (save_ProductionTerm)                             call plotProductionTerm(lambda)
        end if


        !     GET STATISTICS:
        !      print*, 'testError'
        !      call createPDF(,plot=.true.,fieldname='errorTest')
        !      call createPDF(,plot=.true.,fieldname='errorSGS')
     end do
     close(cross_csv)

  end do time_loop
end program
