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
!  $$ Create files for cross-validation
!     open(cross_csv_T_ij,   file=trim(RES_PATH)//trim('crossValidationError_T_ij')//trim(time)//trim('.csv'))
!     open(cross_csv_tau_ij, file=trim(RES_PATH)//trim('crossValidationError_tau_ij')//trim(time)//trim('.csv'))
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

  integer :: time_index, plane_idx
  integer :: p_levels(3) = [43, 129, 212]
  real(8) :: u_rms, epsilon, TKE
  character(1) :: idx

  logical :: debug_PrintFilters = 1

  real(8),allocatable,dimension(:,:,:)   :: dev_t


  if (computeFFT_data.eqv..false.) then
     useTestData      = 0
     readfile         = 0
     filterVelocities = 0
     plot_Velocities  = 0
     save_FFT_data    = 0
  end if
 if (computeFFT_data) multiFilter  = 0

  !    ..INIT POSTPROCESSING..

  open(path_txt, file = trim(RES_DIR)//'path.txt', status = 'replace', action = 'write')

  call setEnv()
  call printParams('display')
  print*, "CASE_NAME:", CASE_NAME
!  call memRequirement()
  print*, boxSize, maskSize
  print*, zLower, zUpper
! stop 


  ! TEST DATA:
  if (useTestData) then
     n_u = 1;   n_uu = 3
     print*, 'Debug mode for velocity components... \n'
  end if

  time_loop: do time_index = time_init, time_final, time_incr

     ! SET TIME PARAMS:
     write(time,'(i0)') time_index !num2str
     if ((len(trim(time))-1).lt.2) time = trim('0')//trim(time)

     ! 1] LOAD DATASET:
!     if(allocated(u).eqv..false.)        allocate(u(n_u,i_GRID,j_GRID,k_GRID))
     if(readFile)                        call readData(DIM = n_u)
     if (dataset.eq.'hst')               u(:,:,256:130:-1,:) = u(:,:,2:128,:) ! CHANGE THIS PART

    !print*, 'u(:,1,1,1): ', u(:,1,1,1)
    ! print*,'u(:,3,15,17): ', u(:,3,15,17)
     !print*, u(:,3,4,23)
     !print*, bigHalf(i_GRID)
     !stop

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
     write(path_txt,*) trim(RES_PATH)
     call system ('mkdir -p '//trim(TEMP_PATH))
     call system ('mkdir -p '//trim(RES_PATH))

  !    ! PLOT VELOCITY Z-MIDPLANE:
!      do i = 1,3
!         write(idx, '(i0)') i
!         open(20,file = trim(RES_PATH)//'u_'//trim(idx)//'.dat')
!         write(20,*) u(i,:,:,z_plane)
!         close(20)
!      end do
     


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

     ! VORTICITY:
     if(compute_vorticity) then
        allocate (omega (n_u, i_GRID, j_GRID, zLower:zUpper))
        call vorticity(u, omega)
        call plotVorticity()
        print*, omega(1,15,24,zLower)
     end if

     ! ************** LEVEL 2 ****************!
     !
     ! ADD PATH DEPTH : SCALE
     write(scale,'(2(i0))') LES_scale, test_scale 
     TEMP_PATH = trim(TEMP_PATH)//'bin'//trim(scale)//'/'
     RES_PATH =  trim(RES_PATH)//'dat'//trim(scale)//'/'//&
                 trim(LESfilterType)//'/'//&
                 trim(TestfilterType)//'/'//&
                 trim(rotationPlane)//&
                 trim(z_plane_name)// '/'
     write(path_txt,*) RES_PATH
     call system ('mkdir -p '//trim(TEMP_PATH))
     call system ('mkdir -p '//trim(RES_PATH))

     ! 2] FILTER VELOCITIES [AND PRESSURE]:
     if(allocated(u_f).eqv..false.)        allocate(u_f (n_u, i_GRID,j_GRID,k_GRID))
     if(allocated(u_t).eqv..false.)        allocate(u_t (n_u, i_GRID,j_GRID,k_GRID))
     if (multiFilter) then
        if(allocated(u_tB).eqv..false.)        allocate(u_tB (n_u, i_GRID,j_GRID,k_GRID))
       ! if(allocated(u_tG).eqv..false.)        allocate(u_tG (n_u, i_GRID,j_GRID,k_GRID))
     end if


     if (filterVelocities) then
        write(*, '(a32)', ADVANCE='NO'), adjustl('        Filter velocities ... ')

        ! CREATE FILTERS:
        allocate(LES (f_GRID,f_GRID,f_GRID))
        allocate(test(f_GRID,f_GRID,f_GRID))

        call createFilter(LES,LES_scale,LESfilterType)
        call createFilter(test,test_scale,TestfilterType)


        !DEBUG : Print filters 
        if (debug_PrintFilters) call printFilters()
        !stop
           
        LES = fftshift(LES)
        test = fftshift(test)

 
        ! Apply filter in Fourier domain (DEFAULT:sharp)
        do i=1,n_u
           u_f(i,:,:,:) = sharpFilter(u  (i,:,:,:),LES)
           u_t(i,:,:,:) = sharpFilter(u_f(i,:,:,:),test) 
        end do
        ! ###
        print*, 'Completed'
!        call check_FFT(u_t(1,15,24,10))  
        if (plot_Velocities) then
                                                                       call plotVelocities('All')
        if (withPressure)                                              call plotPressure()
        end if
     end if

     ! print*, 'P (:,1,1,1)' 
     ! print*, u (:,1,1,1)
     ! print*, u_f (:,1,1,1)
     ! print*, u_t (:,1,1,1)

     ! print*, 'P (:,3,15,17)' 
     ! print*, u (:,3,15,17)
     ! print*, u_f (:,3,15,17)
     ! print*, u_t (:,3,15,17)


!stop
     ! BREAKPOINT 1:
     if (useTestData) stop


     ! 3] GET FFT_DATA:
     if(allocated(tau_ij).eqv..false.)     allocate (tau_ij (6, i_GRID, j_GRID, k_GRID))
     if(allocated(T_ij).eqv..false.)       allocate (T_ij   (6, i_GRID, j_GRID, k_GRID))
     if (multiFilter) then
        if(allocated(T_ijB).eqv..false.)        allocate(T_ijB (6, i_GRID,j_GRID,k_GRID))
!        if(allocated(T_ijG).eqv..false.)        allocate(T_ijG (6, i_GRID,j_GRID,k_GRID))
     end if


     ! COMPUTE ORIGINAL STRESS [ROTATE][SAVE]
     if(computeFFT_data) then
        print*,'Compute original stress:',stress
        call computeStress(u, u_f, u_t, tau_ij, T_ij, LES, test)
        deallocate(LES,test)
        !-     
     else
        ! LOAD SAVED FFT_DATA ../temp/ [CHECK]
        if (multiFilter) then 
          call loadFFT_data3(LESFilterType)
        else
          call loadFFT_data()
        end if
!        call checkFFT_data()
     end if

     !-
     if (rotationAxis == 'X') then
        print*, 'Rotate array along x-axis'
        call rotateX(u_f) 
        call rotateX(u_t) 
        call rotateX(tau_ij) 
        call rotateX(T_ij)
        if (multiFilter) then
            call rotateX(u_tB) 
 !           call rotateX(u_tG) 
            call rotateX(T_ijB)
 !           call rotateX(T_ijG)
        endif 
     else if (rotationAxis == 'Y') then
        print*, 'Rotate array along y-axis'
        call rotateY(u_f) 
        call rotateY(u_t) 
        call rotateY(tau_ij) 
        call rotateY(T_ij)
        if (multiFilter) then
            call rotateY(u_tB) 
   !         call rotateY(u_tG) 
            call rotateY(T_ijB)
  !          call rotateY(T_ijG)
        endif 
     end if


     ! SAVE FFT DATA. THEN LOAD AND ROTATE IT. 
     if (save_FFT_DATA) call saveFFT_data()

     print*, 'tau_ij(1,1,1,1):', tau_ij(1,1,1,1)
     print*, 'tau_ij(2,3,15,17):', tau_ij(2,3,15,17)
     print*, 'T_ij(1,1,1,1):', T_ij(1,1,1,1)
     print*, 'T_ij(2,3,15,17):', T_ij(2,3,15,17)

     !-
     !if (stress.eq.'dev') then
      if (make_Deviatoric) then
        print*, 'Convert to deviatoric stress'
        call makeDeviatoric (tau_ij)
        call makeDeviatoric (T_ij)
        if (multiFilter) call makeDeviatoric (T_ijB)
      end if
      print*, 'tau_ij_dev(1,1,1,1):', tau_ij(1,1,1,1)
      print*, 'T_ij_dev(2,3,15,17):', T_ij(2,3,15,17)
      !print*, 'T_ij_devB(1,15,24,129):', T_ijB(1,15,24,129), '\n'
      
!    stop

    if (plot_Stress)                    call plotOriginalStress('All')
      

     ! if(allocated(Sij_f).eqv..false.)     allocate (Sij_f  (6, i_GRID,j_GRID,zLower:zUpper))
     ! if(allocated(Sij_t).eqv..false.)     allocate (Sij_t  (6, i_GRID,j_GRID,zLower:zUpper))

     
     ! 5] ORIGINAL PRODUCTION FIELD 
     if(production_Term) then

        ! STRAIN RATE:
        print*, 'Compute strain rate \n'
        call computeSij(u_f, Sij_f)
        call computeSij(u_t, Sij_t)
        
        ! ++ 
        if(allocated(Pij_f).eqv..false.)          allocate (Pij_f (i_GRID, j_GRID, zLower:zUpper))
        if(allocated(Pij_t).eqv..false.)          allocate (Pij_t (i_GRID, j_GRID, zLower:zUpper))

        ! tau_ij and T_ij are declared with their limits from 1 to i_GRID. 
        ! If passed without specifying indices, it will map 1 to 129 instead
        ! 129 to 129.
        call productionTerm(Pij_f, tau_ij(:,:,:,zLower:zUpper), Sij_f)
        call productionTerm(Pij_t, T_ij  (:,:,:,zLower:zUpper), Sij_t)

        ! +++  CHECK P_ij
!        call check_beforeExtension()

        if (save_ProductionTerm)                                   call plotProductionTerm()     
        deallocate (Pij_f, Pij_t)
     end if
     

     !  DYNAMIC SMAGORINSKY : PERFORM PLANEWISE COMPUTATION
     ! if (computeDS) then
     !    allocate(S_f_Sij_f_t (6, i_GRID, j_GRID, zLower:zUpper))
     !    call loadPrecomputedStress(S_f_Sij_f_t)  ! LOAD PRECOMPUTED S_f_Sij_f_t: PLANE SLICE
     !    call dyn_Smag(Sij_f, Sij_t, S_f_Sij_f_t, T_ij, tau_DS)
     ! end if


     ! ************** Level 3 ****************!
     !
     ! ADD PATH DEPTH : CASE-NAME
     TEMP_PATH = trim(TEMP_PATH)//trim(CASE_NAME)//'/'
     RES_PATH =  trim(RES_PATH)//trim(CASE_NAME)//'/'

     write(path_txt,*) RES_PATH
     call system ('mkdir -p '//trim(RES_PATH))

!stop
     ! 6] AUTONOMICALLY TUNED LAMBDA
     if(allocated(T_ijOpt).eqv..false.)            allocate (T_ijOpt   (6,i_GRID,j_GRID,zLower:zUpper))
     if(allocated(tau_ijOpt).eqv..false.)          allocate (tau_ijOpt (6,i_GRID,j_GRID,zLower:zUpper))
     if(allocated(h_ij).eqv..false.)               allocate (h_ij      (N,1))

     if (production_Term) then
        if(allocated(Pij_fOpt).eqv..false.)        allocate (Pij_fOpt  (i_GRID, j_GRID, zLower:zUpper))
        if(allocated(Pij_tOpt).eqv..false.)        allocate (Pij_tOpt  (i_GRID, j_GRID, zLower:zUpper))
     end if

     ! EXTEND DOMAIN:
     print*,'Extend domain:'
     call extendDomain(u_f)
     call extendDomain(u_t)
     call extendDomain(T_ij)   
     if (multiFilter) then
        call extendDomain(u_tB)
  !      call extendDomain(u_tG)
        call extendDomain(T_ijB) 
 !       call extendDomain(T_ijG) 
    endif 
 
     call cpu_time(tic)

     ! AUTONOMIC CLOSURE:
     lambda_loop: do iter = 1, size(lambda_0)
        lambda = lambda_0(iter)
        print('(a32,ES8.1)'), 'Autonomic closure, lambda = ', lambda, '\n'
        ! $$
        ! $$$ 
!        allocate (Sij_t  (6, i_GRID,j_GRID,z_extLower:z_extUpper))
!        call computeSij(u_t, Sij_t)
        call autonomicClosure (u_f, u_t, tau_ij, T_ij, T_ijB, h_ij, tau_ijOpt, T_ijOpt)
 !       call check_afterExtension()
     end do lambda_loop


     call cpu_time(toc)
     print*,'Total Elapsed time', toc-tic
     
     ! %%
  end do time_loop

  close (path_txt)

  ! Ring an alarm:
   print*, '\a'
!   print*, 'tau_ijOpt(:,1,1,1):', tau_ijOpt(:,1,1,1)
   print*, 'tau_ijOpt(:,23,23,23):', tau_ijOpt(:,23,23,23)
!   print*, 'tau_ijOpt(:,3,15,17):', tau_ijOpt(:,3,15,17)
contains 


  !****************************************************************
  !                        DEBUG: PRINT FILTERS                   !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Prints out a midplane slice (z=129) of the LES and test 
  !     filter
  !      
  !      
  ! FORM:   subroutine debug_PrintFilters()
  !      
  ! BEHAVIOR: 
  !          
  ! STATUS : Temporary utility procedure for debugging
  !          
  !
  !----------------------------------------------------------------
  
  subroutine printFilters()
    implicit none

    open(20,file=trim(RES_PATH)//'../../test_filter.dat',status='replace')
    open(21,file=trim(RES_PATH)//'../../LES_filter.dat',status='replace')
       
    write(20,*) test (:,:,z_plane)
    write(21,*) LES  (:,:,z_plane)
       
    close(20)
    close(21)

  end subroutine printFilters
      
  !****************************************************************
  !                        CHECK: BEFORE EXTENSION                !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Dumps values at point P(15,24,129) before domain extension.
  !      
  !      
  ! FORM:   subroutine check_beforeExtension()
  !      
  ! BEHAVIOR: 
  !          
  ! STATUS : Temporary utility procedure for debugging
  !          
  !
  !----------------------------------------------------------------
  
  subroutine check_beforeExtension()
    implicit none

    print*, 'Sij_f(2,15,24,129)', Sij_f(2,15,24,z_plane)!, Sij_f(2,15,24,zLower), Sij_t(2,15,24,zUpper)
    print*, 'Sij_t(2,15,24,129)', Sij_t(2,15,24,z_plane), '\n'
    print*, 'tau_ij(2,15,24,129)', tau_ij(2,15,24,z_plane)
    print*, 'T_ij(2,15,24,129)', T_ij(2,15,24,z_plane) 
    !    print*, 'T_ij(2,15,24,z_plane-boxLower-Delta_test)',T_ij(2,15,24,z_plane-boxLower-Delta_test)
    !    print*, 'T_ij(2,15,24,z_plane+boxUpper+Delta_test)',T_ij(2,15,24,z_plane+boxUpper+Delta_test)
    print*, 'u_f(2,15,24,129)', u_f(2,15,24,z_plane)
    !    print*, 'u_f(2,15,24,z_plane-boxLower-Delta_test)',u_f(2,15,24,z_plane-boxLower-Delta_test)
    !    print*, 'u_f(2,15,24,z_plane+boxUpper+Delta_test)',u_f(2,15,24,z_plane+boxUpper+Delta_test)
    print*, 'u_t(2,15,24,129)', u_t(2,15,24,z_plane)
    !    print*, 'u_t(2,15,24,z_plane-boxLower-Delta_test)',u_t(2,15,24,z_plane-boxLower-Delta_test)
    !    print*, 'u_t(2,15,24,z_plane+boxUpper+Delta_test)',u_t(2,15,24,z_plane+boxUpper+Delta_test)

    print*, 'Pij_f(15,24,129)', Pij_f(15,24,z_plane)
    print*, 'Pij_t(15,24,129)', Pij_t(15,24,z_plane)

    !print*, count(u_t == 0), count(u_f == 0), count(T_ij == 0)

  end subroutine check_beforeExtension


  !****************************************************************
  !                        CHECK: AFTER EXTENSION                !
  !****************************************************************
  
  !----------------------------------------------------------------
  ! USE: Dumps values at point P(15,24,129) after domain extension.
  !      
  !      
  ! FORM:   subroutine check_afterExtension()
  !      
  ! BEHAVIOR: 
  !          
  ! STATUS : Temporary utility procedure for debugging
  !          
  !
  !----------------------------------------------------------------

  subroutine check_afterExtension()
    implicit none

    print*, 'T_ij(2,15,24,129)', T_ij(2,15,24,zLower) 
    !    print*, 'T_ij(2,15,24,z_plane-boxLower-Delta_test)',T_ij(2,15,24,z_plane-boxLower-Delta_test)
    !    print*, 'T_ij(2,15,24,z_plane+boxUpper+Delta_test)',T_ij(2,15,24,z_plane+boxUpper+Delta_test)
    print*, 'u_f(2,15,24,129)', u_f(2,15,24,zLower)
    !    print*, 'u_f(2,15,24,z_plane-boxLower-Delta_test)',u_f(2,15,24,z_plane-boxLower-Delta_test)
    !    print*, 'u_f(2,15,24,z_plane+boxUpper+Delta_test)',u_f(2,15,24,z_plane+boxUpper+Delta_test)
    print*, 'u_t(2,15,24,129)', u_t(2,15,24,zLower)
    !    print*, 'u_t(2,15,24,z_plane-boxLower-Delta_test)',u_t(2,15,24,z_plane-boxLower-Delta_test)
    !    print*, 'u_t(2,15,24,z_plane+boxUpper+Delta_test)',u_t(2,15,24,z_plane+boxUpper+Delta_test)

    print*, 'tau_ijOpt(2,15,24,129)', tau_ijOpt(2,15,24,z_plane)
    print*, 'T_ijOpt(2,15,24,129)', T_ijOpt(2,15,24,z_plane)


    print*, 'Pij_fOpt(15,24,129)', Pij_fOpt(15,24,z_plane)
    print*, 'Pij_tOpt(15,24,129)', Pij_tOpt(15,24,z_plane)

  end subroutine check_afterExtension
  
end program



