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

program jhu_dealias

  use fileio
  use fourier

  implicit none


  character(64) :: filename
  character(10) :: scale

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
  call printParams('display')
  print*, 'Dataset: ', dataset, '\n'

    
  ! SET TIME PARAMS:
  write(time,'(i0)') time_index !num2str
  if ((len(trim(time))-1).lt.2) time = trim('0')//trim(time)

  ! 1] LOAD DATASET:
  if(allocated(u).eqv..false.)        allocate(u(n_u,i_GRID,j_GRID,k_GRID))
  if(readFile)                        call readData(DIM = n_u)
  print*, u(1,15,24,10)
 print*,'u(1,9,9,129)', u(1,9,9,129) ! corresponds to u(1,2,2,32)


  ! 2] FILTER VELOCITIES [AND PRESSURE]:
  if(allocated(u_f).eqv..false.)        allocate(u_f (n_u, i_GRID,j_GRID,k_GRID))
  write(*, '(a32)', ADVANCE='NO'), adjustl('        Filter velocities ... ')


  ! CREATE AND APPLY DE-ALIASING FILTER:
  print*,'allocate'
  allocate(LES (f_GRID,f_GRID,f_GRID))
  print*, 'create filter'
  call createFilter(LES,LES_scale)
  print*, 'fftshift'
  call fftshift(LES)


  do i=1,n_u
     u_f(i,:,:,:) = sharpFilter(u  (i,:,:,:),LES)
  end do
  print*,'u_f(1,9,9,129)',u_f(1,9,9,129)

  print*, 'Completed'
  call check_FFT(u_t(1,15,24,10))  
  if (plot_Velocities)                                        call plotVelocities()

  ! 3] SAVE TO BIN FILE:


end program jhu_dealias


