!****************************************************************
!                            STRESS BD
!****************************************************************

!----------------------------------------------------------------
! USE: Precomputes terms for BD stresses. Saves tau_BD
!      
!      
!
! FORM: program stressBD
!       
!
! BEHAVIOR: Load and filter velocities.
!           Computes Sij_f, Sij_t. 
!           Saves tau_BD for the whole domain to a binary file
!
!           
! DEBUG FLAGS:
!
!           
!
! STATUS :
!
! STASH :
! **
!        allocate(LES (f_GRID,f_GRID,f_GRID))
!        call createFilter(LES,LES_scale)
!        call fftshift(LES)
!----------------------------------------------------------------

program stressBD

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
  character(*), parameter :: BD_CASE = 'a1' ! a1B[ox], a2G[auss]

  real(8), parameter :: C_Bardina = 0.9d0

  logical :: save_BD = 1

  
  !    ..INIT POSTPROCESSING..
  call setEnv()
  call printParams('display')


  time_loop: do time_index = time_init, time_final, time_incr

     ! SET TIME PARAMS:
     write(time,'(i0)') time_index !num2str
     if ((len(trim(time))-1).lt.2) time = trim('0')//trim(time)


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
     !!RES_PATH =  trim(RES_PATH)//'dat'//trim(scale)//'/'
     RES_PATH =  trim(RES_PATH)//'dat'//trim(scale)//'/'//&
                 trim(LESfilterType)//&
                 trim(TestfilterType)//&
                 trim(rotationPlane)//&
                 trim(z_plane_name)// '/'

     call system ('mkdir -p '//trim(TEMP_PATH))
     call system ('mkdir -p '//trim(RES_PATH))

     ! 3] GET FFT_DATA:
     if(allocated(u_f).eqv..false.)        allocate (u_f    (n_u,i_GRID, j_GRID, k_GRID))
     if(allocated(u_t).eqv..false.)        allocate (u_t    (n_u,i_GRID, j_GRID, k_GRID))


     ! Can load and filter velocities here or do the same with autonomic.f90. 

     ! if (filterVelocities) then
     !    write(*, '(a32)', ADVANCE='NO'), adjustl('        Filter velocities ... ')

     !    ! CREATE FILTERS:
     !    allocate(LES (f_GRID,f_GRID,f_GRID))
     !    allocate(test(f_GRID,f_GRID,f_GRID))
     !    call createFilter(LES,LES_scale)
     !    call createFilter(test,LES_scale,'Gauss')
     !    call fftshift(LES)
     !    call fftshift(test)

     !    ! Apply filter in Fourier domain 
     !    do i=1,n_u
     !     u_f(i,:,:,:) = sharpFilter(u  (i,:,:,:),LES)
     !     u_t(i,:,:,:) = sharpFilter(u_f(i,:,:,:),test) 
     !   end do
     !   ! ###
     !   print*, 'Completed'
     !   !        call check_FFT(u_t(1,15,24,10))  
     !   if (plot_Velocities) then
     !                                                                  call plotVelocities('All')
     !     if (withPressure)                                            call plotPressure()
     !   end if
     !  end if

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

          
     ! STRAIN RATE [ALL]: zLower=1, zUpper=256
     print*, 'computeSij \n'
     if(allocated(Sij_f).eqv..false.)     allocate (Sij_f  (6, i_GRID,j_GRID,zLower:zUpper))
!     if(allocated(Sij_t).eqv..false.)     allocate (Sij_t  (6, i_GRID,j_GRID,zLower:zUpper))

 
     call computeSij(u_f, Sij_f) 
!     call computeSij(u_t, Sij_t) 

!     print*, 'Sij_f(2,15,24,129)', Sij_f(2,15,24,129)
          
     !  FIND BARDINA STRESS:
     print*, 'Compute SGS stress using Bardina model \n'
    if(allocated(tau_BD).eqv..false.) allocate (tau_BD (6, i_GRID, j_GRID, zLower:zUpper))
     tau_BD = C_Bardina * T_ij

    
     ! COMPUTE Pij_BD:
      if(allocated(Pij_BD).eqv..false.) allocate (Pij_BD (i_GRID, j_GRID, zLower:zUpper))
     call productionTerm(Pij_BD, tau_BD, Sij_f)
     print*,'tau_a2_BD(15,24,129)', tau_BD(2,15,24,z_plane)
     print*,'Pij_BD(15,24,129)', Pij_BD(15,24,z_plane)


     ! SAVE BARDINA STRESS:
     if(save_BD) call plotBardina(BD_CASE)

     ! ! SAVE Pij_BD:
     !  print*,'Saving BD production field in', RES_PATH
     !  open(53, file = trim(RES_PATH)//'Pij_'//trim(BD_CASE)//'_BD.dat', iostat=ferr)
     !  write(53,*) Pij_BD(:,:,z_plane)
     !  close(53)

  end do time_loop

end program stressBD




