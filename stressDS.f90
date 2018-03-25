!****************************************************************
!                            STRESS DS
!****************************************************************

!----------------------------------------------------------------
! USE: Precomputes terms for DS stresses. Saves tau_DS
!      
!      
!
! FORM: program stressDS
!       
!
! BEHAVIOR: Load and filter velocities.
!           Computes Sij_f, Sij_t. 
!           Saves tau_DS for the whole domain to a binary file
!
!           
! DEBUG FLAGS:
!
!           
!
! STATUS : + Change suffix from _dev to _GS, etc in line 169.
!
! STASH :
! **
!        allocate(LES (f_GRID,f_GRID,f_GRID))
!        call createFilter(LES,LES_scale)
!        call fftshift(LES)
!----------------------------------------------------------------

program stressDS

  use fourier
  use actools
  use solver

  implicit none

  character(64) :: filename
  character(10) :: scale

  integer :: time_index
  integer :: ij
  integer :: ferr
  character(1) :: idx

  logical :: save_tau_DS = 1
  !    ..INIT POSTPROCESSING..
  call setEnv()
  call printParams('display')

  !    ..CHECK WHOLE DOMAIN CONDITION..
  if (compDomain.ne.'all') then
     print*, 'CompDomain is not set to ALL'
     stop
  end if

  time_loop: do time_index = time_init, time_final, time_incr

     ! SET TIME PARAMS:
     write(time,'(i0)') time_index !num2str
     if ((len(trim(time))-1).lt.2) time = trim('0')//trim(time)

     ! 1] LOAD DATASET:
!     if(allocated(u).eqv..false.)        allocate(u(n_u,i_GRID,j_GRID,k_GRID))
!     call readData(DIM = n_u)
!     if (dataset.eq.'hst')               u(:,:,256:130:-1,:) = u(:,:,2:128,:) ! CHANGE THIS PART


     ! ************** LEVEL 1 ****************!
     !
     ! ADD PATH DEPTH: DATASET
     
     TEMP_PATH = trim(TEMP_DIR)//trim(dataset)//'/'
     RES_PATH  = trim(RES_DIR)//trim(dataset)//'/'
     if (dataset.eq.'hst') then
        TEMP_PATH = trim(TEMP_PATH)//trim(hst_set)//'/'
        RES_PATH  = trim(RES_PATH)//trim(hst_set)//'/'
     end if

     call system ('mkdir -p '//trim(TEMP_PATH))
     call system ('mkdir -p '//trim(RES_PATH))


     ! ************** LEVEL 2 ****************!
     !
     ! ADD PATH DEPTH : SCALE
     write(scale,'(2(i0))') LES_scale, test_scale 
     TEMP_PATH = trim(TEMP_PATH)//'bin'//trim(scale)//'/'
     RES_PATH =  trim(RES_PATH)//'dat'//trim(scale)//'/'//&
                 trim(LESfilterType)//'/'//&
                 trim(TestfilterType)//'/'//&
                 trim(rotationPlane) !// &
                 !trim(z_plane_name)// '/'

     call system ('mkdir -p '//trim(TEMP_PATH))
     call system ('mkdir -p '//trim(RES_PATH))


     ! 3] GET FFT_DATA:
     if(allocated(u_f).eqv..false.)        allocate (u_f    (n_u,i_GRID, j_GRID, k_GRID))
     if(allocated(u_t).eqv..false.)        allocate (u_t    (n_u,i_GRID, j_GRID, k_GRID))
     if(allocated(tau_ij).eqv..false.)     allocate (tau_ij (6,  i_GRID, j_GRID, k_GRID))
     if(allocated(T_ij).eqv..false.)       allocate (T_ij   (6,  i_GRID, j_GRID, k_GRID))

     ! LOAD SAVED FFT_DATA ../temp/ [CHECK]
     call loadFFT_data()
!     call checkFFT_data()        
     deallocate(tau_ij)     ! Don't need tau_ij
     
     !->>
     if (rotationAxis == 'X') then
        print*, 'Rotate array along x-axis'
        call rotateX(u_f) 
        call rotateX(u_t) 
        call rotateX(T_ij)
     else if (rotationAxis == 'Y') then
        print*, 'Rotate array along y-axis'
        call rotateY(u_f) 
        call rotateY(u_t) 
        call rotateY(T_ij)
     end if

     if (stress.eq.'dev') then
        print*, 'Convert to deviatoric stress'
        call makeDeviatoric (T_ij)
      end if

          
     ! STRAIN RATE [ALL]: zLower=1, zUpper=256
     print*, 'computeSij \n'
     if(allocated(Sij_f).eqv..false.)     allocate (Sij_f  (6, i_GRID,j_GRID,zLower:zUpper))
     if(allocated(Sij_t).eqv..false.)     allocate (Sij_t  (6, i_GRID,j_GRID,zLower:zUpper))

     call computeSij(u_f, Sij_f) 
     call computeSij(u_t, Sij_t) 

!     print*, 'Sij_f(2,15,24,129)', Sij_f(2,15,24,129)
!     print*, 'Sij_t(2,15,24,129)', Sij_t(2,15,24,129)

     
     ! CREATE TEST FILTER: **
     allocate(test(f_GRID,f_GRID,f_GRID))
     call createFilter(test, test_scale, TestfilterType)        
     test = fftshift(test)
 
     ! PRECOMPUTE S_f_Sij_f_t [SAVE]: 
     print*
     print*, 'Precompute S_f_Sij_f_t'
     print*

     allocate(S_f_Sij_f_t (6, i_GRID, j_GRID, zLower:zUpper)) 
!call cpu_time(tic)
     call computeS_f_Sij_f_t(Sij_f, test, S_f_Sij_f_t) ! SAVE S_f_Sij_f_t WITHIN THE SUBROUTINE
     deallocate(test)

     !  DYNAMIC SMAGORINSKY:
     print*, 'Compute Dynamic Smagorinsky stress \n'
     allocate (tau_DS (6, i_GRID, j_GRID, zLower:zUpper))
     call dyn_Smag(Sij_f, Sij_t, S_f_Sij_f_t, T_ij, tau_DS)
!call cpu_time(toc)
!print*, 'Elapsed time', toc-tic
!stop
    
     ! COMPUTE Pij_DS:
     allocate (Pij_DS (i_GRID, j_GRID, zLower:zUpper))
     call productionTerm(Pij_DS, tau_DS, Sij_f)

      ! SAVE DYN SMAG:
     if(save_tau_DS) call plotDynSmag()

      print*,'Pij_DS(15,24,129)', Pij_DS(15,24,129)
     ! print*,'Pij_DS(max)', maxval(Pij_DS(:,:,129)), 'at', maxloc(Pij_DS(:,:,129))
     ! print*,'Sij_f(max)', maxval(Sij_f(:,:,:,129)), 'at', maxloc(Sij_f(:,:,:,129))
!stop
    
  end do time_loop

  
end program stressDS




