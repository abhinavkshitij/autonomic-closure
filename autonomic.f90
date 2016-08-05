!****************************************************************
!                            AUTONOMIC
!****************************************************************

!----------------------------------------------------------------
! USE: This is the MAIN() equivalent of the code. Loads library 
!      procedures that are defined in MODULES. 
!      
!      
!
! FORM: program autonomic
!       
!
! BEHAVIOR: Currently has no capability to read command line inputs.
!
!
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
! $$$ lambda (coarse points vs fine points)
!         if (mod(iter,2).eq.1) then
!            lambda = lambda_0(1) * 10**((iter-1)/2)
!         else
!            lambda = lambda_0(2) * 10**((iter-1)/2)
!         end if
!
! %%%
!   GET STATISTICS:
!       print*, 'testError'
!       call createPDF(var,plot=.true.,fieldname='errorTest')
!       call createPDF(var,plot=.true.,fieldname='errorSGS')
!
!----------------------------------------------------------------

program autonomic

  use fourier
  use actools
  use solver

  implicit none

  character(64) :: filename
  character(10) :: scale

  integer :: time_index
  real(8) :: u_rms, epsilon, TKE


  if (computeFFT_data.eqv..false.) then
     useTestData      = 0
     readfile         = 0
     filterVelocities = 0
     plot_Velocities  = 0
     save_FFT_data    = 0
  end if

  !    ..INIT POSTPROCESSING..
  open(path_txt, file = trim(RES_DIR)//'path.txt', status = 'replace', action = 'write')

  call setEnv()
  call printParams('display')
  print*, 'Dataset: ', dataset, '\n'

stop
!print*, n_u, n_uu

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

     ! TURBULENT FIELD STATISTICS:
     if (turbulentStats) then
     ! ## PLOT Energy spectra
         TKE  = turbulentKE(u)
         u_rms = sqrt( (mean(u(1,:,:,:)**2) + mean(u(2,:,:,:)**2) + mean(u(3,:,:,:)**2))/3.d0 )
         print*, 'turbulentKE (k) = ',TKE
         print*, 'rms velocity (u_rms) = ', u_rms 
         if(allocated(Sij).eqv..false.)     allocate (Sij  (n_uu, i_GRID,j_GRID,k_GRID))
         print*, 'nu = ', nu
         call computeSij(u, Sij)
         call dissipationRate(Sij, epsilon)
         print*, 'Dissipation rate (epsilon) = ', epsilon
     end if

!stop
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
        if (plot_Velocities) then
                                                                       call plotVelocities('All')
        if (withPressure)                                              call plotPressure()
        end if
     end if

!stop

     ! BREAKPOINT 1:
     if (useTestData) stop


     ! 3] GET FFT_DATA:
     if(allocated(tau_ij).eqv..false.)     allocate (tau_ij (6, i_GRID, j_GRID, k_GRID))
     if(allocated(T_ij).eqv..false.)       allocate (T_ij   (6, i_GRID, j_GRID, k_GRID))


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
     if (plot_Stress)                                            call plotOriginalStress('All')


     if(allocated(Sij_f).eqv..false.)     allocate (Sij_f  (6, i_GRID,j_GRID,k_GRID))
     if(allocated(Sij_t).eqv..false.)     allocate (Sij_t  (6, i_GRID,j_GRID,k_GRID))


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
!     TEMP_PATH = trim(TEMP_PATH)//trim(solutionMethod)//'/' 
!     RES_PATH =  trim(RES_PATH)//trim(solutionMethod)//'/'  
     TEMP_PATH = trim(TEMP_PATH)//trim(CASE_NAME)//'/'
     RES_PATH =  trim(RES_PATH)//trim(CASE_NAME)//'/'


     write(path_txt,*) RES_PATH
!     call system ('mkdir -p '//trim(TEMP_PATH))
     call system ('mkdir -p '//trim(RES_PATH))


     ! 6] AUTONOMICALLY TUNED LAMBDA
     if(allocated(T_ijOpt).eqv..false.)            allocate (T_ijOpt   (6,i_GRID,j_GRID,k_GRID))
     if(allocated(tau_ijOpt).eqv..false.)          allocate (tau_ijOpt (6,i_GRID,j_GRID,k_GRID))
     if(allocated(h_ij).eqv..false.)               allocate (h_ij      (N,P))

     if (production_Term) then
        if(allocated(Pij_fOpt).eqv..false.)        allocate (Pij_fOpt  (i_GRID, j_GRID, k_GRID))
        if(allocated(Pij_tOpt).eqv..false.)        allocate (Pij_tOpt  (i_GRID, j_GRID, k_GRID))
     end if


     ! EXTEND DOMAIN:
     call extendDomain(u_f)
     call extendDomain(u_t)
     call extendDomain(T_ij) 

     print*, 'Autonomic closure ... '
     open(cross_csv_T_ij,   file=trim(RES_PATH)//trim('crossValidationError_T_ij')//trim(time)//trim('.csv'))
     open(cross_csv_tau_ij, file=trim(RES_PATH)//trim('crossValidationError_tau_ij')//trim(time)//trim('.csv'))
! $$$
     call autonomicClosure (u_f, u_t, tau_ij, T_ij, h_ij, tau_ijOpt, T_ijOpt)

     close(cross_csv_T_ij)
     close(cross_csv_tau_ij)
! %%
  end do time_loop

  close (path_txt)
end program
