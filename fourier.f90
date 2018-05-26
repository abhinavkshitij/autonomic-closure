!****************************************************************
!                              FOURIER
!****************************************************************

!----------------------------------------------------------------
! USE: Signal processing toolbox
!       1) Generate filter (spectrally-sharp)           : SOURCE  
!       2) Apply filter (coarse-graining)               : FILTER
!       3) [SAVE] final results in RUN dir              : SINK 
!
! FORM: module fourier
!          contains
!       subroutine createFilter       [SOURCE]
!       subroutine fftshift           [FILTER]
!       function sharpFilter          [FILTER]
!          - subroutine check_FFT      [CHECK] 
!       subroutine computeStress      [FILTER]
!          - subroutine check_Stress   [CHECK]
!
! BEHAVIOR: Needs C-binding for FFTW libraries. 
!           
!
! STATUS : Check for non-uniform domains
! 
!----------------------------------------------------------------

module fourier
  use global

  use, intrinsic :: iso_c_binding
  include 'fftw3.f03' 
 
contains

  !****************************************************************
  !                         CREATEFILTER
  !****************************************************************

  !----------------------------------------------------------------
  ! USE: Create an explicit filter.
  !      
  ! FORM: subroutine createFilter(array, scale,[filterOption-default:Sharp])
  !       
  ! BEHAVIOR: GRID must be cubic for spherical symmetry.
  !          
  ! Layout: The center should lie at an unit offset from the geometric center.
  !         This is because the DC component occupies the first index. 
  ! 
  !    1             129          256 
  !  1  |-------------|------------|
  !     |             |            |
  !     |             |            |
  !     |      A      |     B      |
  !  128|             |            |
  !  129--------------o------------|
  !     |             |            |
  !     |      D      |     C      |
  !     |             |            |
  !  256--------------|------------|
  !
  !  STATUS : > Passed test for FFTW layout.
  !          > Passed test for center in spectral space.
  !          > Check for the axes -- should be in alignment with FFTW output.
  !
  ! Notes  : 1) Create and define variables in a module. 
  !          2) The filter creation can be parallelized (not needed right now).
  !          3) Include options for other filters.   
  !----------------------------------------------------------------
  
  subroutine createFilter(filter,scale,filterOption)    
    implicit none
    !
    !    ..ARGUMENTS..
    real(8),dimension(:,:,:),intent(inout) :: filter
    integer, intent(in)                    :: scale
    character(*), optional, intent(in)     :: filterOption
    !
    !    ..LOCAL VARS.. 
    real(8) :: distance
    real(8) :: slope  
    integer :: last, center

    last = size(filter,dim=1)
    center = bigHalf(last)

    if (present(filterOption).and.filterOption.eq.'Gauss') then
      !  Create Gauss filter:
      do k = 1,last
         do j = 1,last
            do i = 1,last

               distance = sqrt( dble((i - center)**2) &
                    +           dble((j - center)**2) &
                    +           dble((k - center)**2) )

               filter(i,j,k) = exp(-5.d-1 * (distance / scale)**2)

            end do
         end do
      end do
    elseif (present(filterOption).and.filterOption.eq.'Box') then
      !  Create Box filter:
      do k = 1,last
         do j = 1,last
            do i = 1,last

               distance = sqrt( dble((i - center)**2) &
                    +           dble((j - center)**2) &
                    +           dble((k - center)**2) )

               !  filter(i,j,k) = sin( (distance+eps) / (scale+eps) ) / &
               ! &                   ( (distance+eps) / (scale+eps) )

               filter(i,j,k) = sin( (PI*(distance+eps)) / (scale+eps) ) / &
               &                  ( (PI*(distance+eps)) / (scale+eps) )

            end do
         end do
      end do
      elseif (present(filterOption).and.filterOption.eq.'Tri') then
      !  Create Tri filter:
      do k = 1,last
         do j = 1,last
            do i = 1,last

               distance = sqrt( dble((i - center)**2) &
                    +           dble((j - center)**2) &
                    +           dble((k - center)**2) )

                filter(i,j,k) = sin( (PI*(distance+eps)) / (scale+eps) * 0.5d0) ** 2 / &
               &                   ( (PI*(distance+eps)) / (scale+eps) * 0.5d0) ** 2

            end do
         end do
      end do
      elseif (present(filterOption).and.filterOption.eq.'Custom') then
      !  Create Custom filter:
      do k = 1,last
         do j = 1,last
            do i = 1,last

               distance = sqrt( dble((i - center)**2) &
                    +           dble((j - center)**2) &
                    +           dble((k - center)**2) )

    
               if (distance.ge.0.d0.and.distance.lt.test_scale) then
                  slope = (0.75d0 - 1.d0) / real (test_scale)
                  filter(i,j,k) = slope * distance + 1.d0
               elseif (distance.ge.test_scale.and.distance.lt.LES_scale) then
                  filter(i,j,k) = 0.75d0  
               elseif (distance.ge.LES_scale.and.distance.lt.3*LES_scale)  then
                  filter(i,j,k) = 0.33d0
               endif        

            end do
         end do
      end do  
    else
      !  Create spectrally sharp filter:
      do k = 1,last
         do j = 1,last
            do i = 1,last

               distance = sqrt( dble((i - center)**2) &
                    +           dble((j - center)**2) &
                    +           dble((k - center)**2) )
               if (distance.le.scale) filter(i,j,k) = 1.d0

            end do
         end do
      end do

      ! ! Check spherical symmetry:
      ! if ( (filter (center+scale+1,center,center) .eq. 0)      .and.  &
      !      (filter (center,center-scale,center)   .eq. 1)      .and.  &
      !      (filter (center+scale,center,center)   .eq. 1)      .and.  &
      !      (filter (center,center,center+scale)   .eq. 1)      .and.  &
      !      (filter (center,center-scale-1,center) .eq. 0) )    then
      ! else
      !    print*, 'Failed creating spherically symmetric filter ... Aborting'
      !    stop
      ! end if
    endif

    return
  end subroutine createFilter

  !****************************************************************
  !                            FFTSHIFT
  !****************************************************************

  !----------------------------------------------------------------
  ! USE: Peforms fftshift operation.
  !      
  ! FORM: subroutine fftshift   [FILTER]
  !       
  ! BEHAVIOR: Brings DC component (k=0) at the center. 
  !           Nyquist freqeuncy at the top left corner.
  !         
  ! 
  ! Layout:                     k 
  !             E       F       |
  !           A       B         |
  !                             /----> j
  !             H       G      /
  !           D       C       i
  !  
  !
  ! In clockwise direction from the topleft corner, 
  ! Front half - (i = 129,256) :  A, B, C, D
  ! Rear half  - (i = 1,128)   :  E, F, G, H
  !
  ! SWAP RULES : A <-> G, C <-> E, D <-> F, B <-> H
  !----------------------------------------------------------------
  
  function fftshift(array)
    !
    !    ..ARRAY ARGUMENTS..
    real(8),dimension(:,:,:):: array
    real(8),dimension(:,:,:), allocatable:: fftshift
    !
    !    ..WORK ARRAYS..
    real(8),allocatable,dimension(:,:,:) :: temp
    real(8),allocatable,dimension(:,:,:) :: A,B,C,D,E,F,G,H 
    !
    !   ..LOCAL VARS..
    integer :: last,center
    
    last = size(array,dim=1)
    center = bigHalf(last)
    
    allocate (fftshift(last,last,last))

    ! CREATE SUBARRAYS:
    allocate(temp(1:(center-1), 1:(center-1) , 1:(center-1) ) )
    allocate(A ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
    allocate(B ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
    allocate(C ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
    allocate(D ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
    allocate(E ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
    allocate(F ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
    allocate(G ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
    allocate(H ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )


    A = array( center:last  ,  1:(center-1) ,  center:last  )
    B = array( center:last  ,  center:last  ,  center:last  )
    C = array( center:last  ,  center:last  ,  1:(center-1) )
    D = array( center:last  ,  1:(center-1) ,  1:(center-1) )
    E = array( 1:(center-1) ,  1:(center-1) ,  center:last  )
    F = array( 1:(center-1) ,  center:last  ,  center:last  )
    G = array( 1:(center-1) ,  center:last  ,  1:(center-1) )
    H = array( 1:(center-1) ,  1:(center-1) ,  1:(center-1) )

    ! SWAP SUBARRAYS: (_SHIFT)
    temp = A ; A = G ; G = temp
    temp = C ; C = E ; E = temp
    temp = D ; D = F ; F = temp
    temp = B ; B = H ; H = temp

    ! RECREATE array WITH FFTSHIFT:
    fftshift( center:last  , 1:(center-1) , center:last  ) = A
    fftshift( center:last  , center:last  , center:last  ) = B
    fftshift( center:last  , center:last  , 1:(center-1) ) = C
    fftshift( center:last  , 1:(center-1) , 1:(center-1) ) = D
    fftshift( 1:(center-1) , 1:(center-1) , center:last  ) = E
    fftshift( 1:(center-1) , center:last  , center:last  ) = F
    fftshift( 1:(center-1) , center:last  , 1:(center-1) ) = G
    fftshift( 1:(center-1) , 1:(center-1) , 1:(center-1) ) = H

    deallocate (temp,A,B,C,D,E,F,G,H)
    return
  end function fftshift

  !****************************************************************
  !                            FFTSHIFT_C
  !****************************************************************

  !----------------------------------------------------------------
  ! USE: Peforms fftshift operation on a complex array 
  !      
  ! FORM: subroutine fftshift_c   [FILTER]
  !       
  ! BEHAVIOR: Brings DC component (k=0) at the center. 
  !           Nyquist freqeuncy at the top left corner.
  !         
  ! 
  ! Layout:                     k 
  !             E       F       |
  !           A       B         |
  !                             /----> j
  !             H       G      /
  !           D       C       i
  !  
  !
  ! In clockwise direction from the topleft corner, 
  ! Front half - (i = 129,256) :  A, B, C, D
  ! Rear half  - (i = 1,128)   :  E, F, G, H
  !
  ! SWAP RULES : A <-> G, C <-> E, D <-> F, B <-> H
  !----------------------------------------------------------------

  function fftshift_c(array)
    !
    !    ..ARRAY ARGUMENTS..
    complex(C_DOUBLE_COMPLEX),dimension(:,:,:):: array
    complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:,:) :: fftshift_c
    !
    !    ..WORK ARRAYS..
    complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:,:) :: temp
    complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:,:) :: A,B,C,D,E,F,G,H 
    !
    !   ..LOCAL VARS..
    integer :: last,center
    
    last = size(array,dim=1)
    center = bigHalf(last)
    allocate (fftshift_c(last,last,last))

    ! CREATE SUBARRAYS:
    allocate(temp(1:(center-1), 1:(center-1) , 1:(center-1) ) )
    allocate(A ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
    allocate(B ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
    allocate(C ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
    allocate(D ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
    allocate(E ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
    allocate(F ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
    allocate(G ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
    allocate(H ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )


    A = array( center:last  ,  1:(center-1) ,  center:last  )
    B = array( center:last  ,  center:last  ,  center:last  )
    C = array( center:last  ,  center:last  ,  1:(center-1) )
    D = array( center:last  ,  1:(center-1) ,  1:(center-1) )
    E = array( 1:(center-1) ,  1:(center-1) ,  center:last  )
    F = array( 1:(center-1) ,  center:last  ,  center:last  )
    G = array( 1:(center-1) ,  center:last  ,  1:(center-1) )
    H = array( 1:(center-1) ,  1:(center-1) ,  1:(center-1) )

    ! SWAP SUBARRAYS: (_SHIFT)
    temp = A ; A = G ; G = temp
    temp = C ; C = E ; E = temp
    temp = D ; D = F ; F = temp
    temp = B ; B = H ; H = temp

    ! RECREATE array WITH FFTSHIFT:
    fftshift_c( center:last  , 1:(center-1) , center:last  ) = A
    fftshift_c( center:last  , center:last  , center:last  ) = B
    fftshift_c( center:last  , center:last  , 1:(center-1) ) = C
    fftshift_c( center:last  , 1:(center-1) , 1:(center-1) ) = D
    fftshift_c( 1:(center-1) , 1:(center-1) , center:last  ) = E
    fftshift_c( 1:(center-1) , center:last  , center:last  ) = F
    fftshift_c( 1:(center-1) , center:last  , 1:(center-1) ) = G
    fftshift_c( 1:(center-1) , 1:(center-1) , 1:(center-1) ) = H

    deallocate (temp,A,B,C,D,E,F,G,H)
    return
  end function fftshift_c


  !****************************************************************
  !                         FORWARD FT
  !****************************************************************

  !----------------------------------------------------------------
  ! USE: Takes forward FFT
  !      
  ! FORM:  function forwardFFT(array)  
  !       
  ! BEHAVIOR: array_work must be of size(f_GRID,f_GRID,f_GRID)
  !         
  ! STATUS : 
  ! Notes  : 
  !         
  !----------------------------------------------------------------

  function forwardFT(array)
    implicit none
    !
    !    ..ARRAY ARGUMENTS..
    complex(C_DOUBLE_COMPLEX),dimension(:,:,:) :: array
    complex(C_DOUBLE_COMPLEX),dimension(f_GRID,f_GRID,f_GRID) :: forwardFT
    !
    !    ..LOCAL VARS..
    type(C_PTR) :: plan

    ! FFT:
    call dfftw_plan_dft_3d(plan,f_GRID,f_GRID,f_GRID,array,forwardFT,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute(plan)    
    call dfftw_destroy_plan(plan)
    
  end function forwardFT

  !****************************************************************
  !                         INVERSE FFT
  !****************************************************************

  !----------------------------------------------------------------
  ! USE: Takes inverse FFT
  !      
  ! FORM:  function inverseFFT(array)  
  !       
  ! BEHAVIOR: array_work must be of size(f_GRID,f_GRID,f_GRID)
  !         
  ! STATUS : 
  ! Notes  : 
  !         
  !----------------------------------------------------------------

  function inverseFT(array)
    implicit none
    !
    !    ..ARRAY ARGUMENTS..
    complex(C_DOUBLE_COMPLEX),dimension(:,:,:) :: array
    complex(C_DOUBLE_COMPLEX),dimension(f_GRID,f_GRID,f_GRID) :: inverseFT
    !
    !    ..LOCAL VARS..
    type(C_PTR) :: plan

    ! IFT:
    call dfftw_plan_dft_3d(plan,f_GRID,f_GRID,f_GRID,array,inverseFT,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_execute(plan)    
    call dfftw_destroy_plan(plan)
    
  end function inverseFT

  !****************************************************************
  !                         SHARPFILTER
  !****************************************************************

  !----------------------------------------------------------------
  ! USE: Convolves physical field with a filter 
  !      
  ! FORM:  function sharpFilter(array_work,filter)  
  !       
  ! BEHAVIOR: array_work must be of size(f_GRID,f_GRID,f_GRID)
  !         
  ! STATUS : > Passed test with MATLAB results (10/21/2015).    
  ! Notes: 
  !         1) dfftw_execute(plan) takes 2 secs
  !         2) Uses FFTW libraries.
  !         3) Use C_DOUBLE, C_INT, etc for C-compatibility on any platform.
  !         4) Can use ABS(in_cmplx) but residual imaginary values 
  !            are not exactly 0.d0.
  !         5) Normalization depends on total number of points
  !            including the padded section. So f_GRID^3 is the
  !            correct normalization factor. 
  !
  ! %%%
  !     print*,'Write FFT files from sharpFilter()'
  !     open(1, file= trim(RES_PATH)//'out_cmplx1.dat')
  !     open(11,file= trim(RES_PATH)//'out_cmplx2.dat')
  !     open(2, file= trim(RES_PATH)//'in_cmplx.dat')
  !     write(1,*) abs(out_cmplx(:,1,:))
  !     write(11,*) abs(out_cmplx(:,:,1))
  !     write(2,*) real(in_cmplx(:,:,129))
  !     close(1)
  !     close(11)
  !     close(2)
  !----------------------------------------------------------------

  function sharpFilter(u_real,filter)
    implicit none
    !
    !    ..ARRAY ARGUMENTS..
    real(C_DOUBLE),dimension(:,:,:), intent(in) :: u_real, filter
    real(C_DOUBLE),dimension(f_GRID,f_GRID,f_GRID) :: sharpFilter 
    !
    !    ..WORK ARRAYS..
    complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:,:):: u_cmplx
    !
    !    ..LOCAL VARS..
    type(C_PTR) :: plan
   
    allocate(u_cmplx (f_GRID,f_GRID,f_GRID))

    u_cmplx = dcmplx (u_real) / (dble(f_GRID**3)) 
    sharpFilter = real(inverseFT(forwardFT(u_cmplx) * filter)) 
    
  end function sharpFilter

  !****************************************************************
  !                         CHECK FFT
  !****************************************************************

  !----------------------------------------------------------------
  ! USE: Check for FFT opertion.
  !
  ! FORM: subroutine check_FFT(value)
  !       
  ! BEHAVIOR: Optionally filters the resulting product. 
  !         
  ! STATUS: Redundant. Replace with a checksum.  
  ! 
  !----------------------------------------------------------------

  subroutine check_FFT(value)
    implicit none
    real(8) :: value
    if (dataset.eq.'jhu256') then
       if (abs(value - (-0.48241021987284982d0)).gt.tiny(0.d0)) then
          print*, 'Precision test for FFT: Failed'
          print*, 'Check precision or data for testing is not JHU 256'
          print*, value
!          stop
       else
          print*, 'Precision test for FFT: Passed'
       end if
    end if
  end subroutine check_FFT

  !****************************************************************
  !                         DEALIASED PRODUCTS
  !****************************************************************

  !----------------------------------------------------------------
  ! USE: Convolves TWO products at 2x grid resolution and then 
  !      projects resulting product at the original (1x) grid 
  !      resolution.
  !      
  ! FORM: function dealiasedProducts(u1,u2,[filter])  
  !       
  ! BEHAVIOR: Optionally filters the resulting product. 
  !         
  ! STATUS: > Unit testing with full run.   
  ! 
  ! Notes: The 3 DP complex arrays take 6GB memory
  !         
  !----------------------------------------------------------------

  function dealiasedProducts(u1_real, u2_real, filter)
    implicit none
    !
    !    ..ARRAY ARGUMENTS..
    real(C_DOUBLE),dimension(:,:,:), intent(in) :: u1_real, u2_real
    real(C_DOUBLE),dimension(:,:,:), intent(in), optional :: filter
    real(C_DOUBLE),dimension(f_GRID,f_GRID,f_GRID) :: dealiasedProducts
    !
    !    ..WORK ARRAYS..
    complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:,:):: u1_cmplx, u2_cmplx
    complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:,:):: temp
    !
    !    ..LOCAL VARS..
    type(C_PTR) :: plan
   
    allocate(u1_cmplx (f_GRID,f_GRID,f_GRID))
    allocate(u2_cmplx (f_GRID,f_GRID,f_GRID))


    u1_cmplx = dcmplx (u1_real) / (dble(f_GRID**3)) 
    u2_cmplx = dcmplx (u2_real) / (dble(f_GRID**3)) 

    ! FT: [FORTRAN is CASE-INSENSITIVE]
    !U1_cmplx = fftshift_c(forwardFT(u1_cmplx))
    !U2_cmplx = fftshift_c(forwardFT(u2_cmplx))

    ! print*,'Write FFT files from dealiasedProducts()'
    !   open(1, file= trim(RES_PATH)//'U1_cmplx_256.dat')
    !   open(2, file= trim(RES_PATH)//'U2_cmplx_256.dat')
    !   write(1,*) abs(U1_cmplx(:,:,129))
    !   write(2,*) abs(U2_cmplx(:,:,129))
    !   close(1)
    !   close(2)

    !stop

    print*, 2*f_GRID

    ! ZERO-PAD U(256) from k = 129 to 256 -> U(512)
    allocate (temp (2*f_GRID, 2*f_GRID, 2*f_GRID))
    temp = U1_cmplx
    deallocate (U1_cmplx)
    allocate (U1_cmplx(2*f_GRID, 2*f_GRID, 2*f_GRID))
    !U1_cmplx(129:384, 129:384, 129:384) = temp(1:256, 1:126, 1:256)

    temp = U2_cmplx
    deallocate(U2_cmplx)
    allocate(U2_cmplx(2*f_GRID, 2*f_GRID, 2*f_GRID))
    !U2_cmplx(129:384, 129:384, 129:384) = temp(1:256, 1:256, 1:256)

    ! IFFT:
    
 
    !dealiasedProducts = real(in_cmplx) 
    
  end function dealiasedProducts




  !****************************************************************
  !                         COMPUTESTRESS
  !****************************************************************

  !----------------------------------------------------------------
  ! USE: Computes tau_ij and T_ij for all ij 
  !      
  ! FORM: subroutine computeStress(u, u_f, u_t, tau_ij, T_ij,  LES, test)
  !       
  ! BEHAVIOR: 
  !         
  ! STATUS : > Passed test with MATLAB results (10/21/2015).    
  !         
  ! Boulder CODE: [for code validation]
  !  ! DEVIATORIC STRESS:
  !     elseif (stress.eq.'dev')then
  !        k = 1
  !        do j=1,3
  !           do i=1,3
  !              if (i.le.3.and.j.le.3.and.i.ge.j) then
  !                 print *, 'tau(', i, ',', j, ')'
  !                 tau_ij(k,:,:,:) = sharpFilter(u(i,:,:,:) * u(j,:,:,:),LES)     &
  !                                           - u_f(i,:,:,:) * u_f(j,:,:,:)
  !                 print *, 'T(', i, ',', j, ')'
  !                    T_ij(k,:,:,:) = sharpFilter(u(i,:,:,:) * u(j,:,:,:),test)   &  
  !               k = k + 1
  !              end if
  !           end do
  !        end do
  !       
  !----------------------------------------------------------------


  subroutine computeStress(u, u_f, u_t, tau_ij, T_ij, LES, test)
    implicit none
    !
    !    ..ARRAY ARGUMENTS..
    real(8),dimension(:,:,:,:),intent(in)  :: u, u_f, u_t
    real(8),dimension(  :,:,:),intent(in)  :: LES, test
    real(8),dimension(:,:,:,:),intent(out) :: tau_ij, T_ij
    !
    !    ..WORK ARRAY..
    real(8),allocatable,dimension(:,:,:)   :: dev_t

    ! ABSOLUTE STRESS:
    k = 1
     do j=1,3
        do i=1,3
           if (i.le.3.and.j.le.3.and.i.ge.j) then
              print *, 'tau(', i, ',', j, ')' !<-- CHECK ORDER OF i,j,k...affects performance!!!
              tau_ij(k,:,:,:) = sharpFilter(u(i,:,:,:) * u(j,:,:,:), LES)       &
                                        - u_f(i,:,:,:) * u_f(j,:,:,:)
              print *, 'T(', i, ',', j, ')'
              T_ij(k,:,:,:) = sharpFilter(u_f(i,:,:,:) * u_f(j,:,:,:), test)    &
                                        - u_t(i,:,:,:) * u_t(j,:,:,:)
              k = k + 1
           end if
        end do
     end do

    ! DEVIATORIC STRESS:
    if (stress.eq.'dev') then
     allocate(dev_t(i_GRID,j_GRID,k_GRID))

     dev_t = (tau_ij(1,:,:,:) + tau_ij(4,:,:,:) + tau_ij(6,:,:,:)) / 3.d0
     tau_ij(1,:,:,:) = tau_ij(1,:,:,:) - dev_t
     tau_ij(4,:,:,:) = tau_ij(4,:,:,:) - dev_t
     tau_ij(6,:,:,:) = tau_ij(6,:,:,:) - dev_t

     dev_t = (T_ij(1,:,:,:) + T_ij(4,:,:,:) + T_ij(6,:,:,:)) / 3.d0
     T_ij(1,:,:,:) = T_ij(1,:,:,:) - dev_t
     T_ij(4,:,:,:) = T_ij(4,:,:,:) - dev_t
     T_ij(6,:,:,:) = T_ij(6,:,:,:) - dev_t
    end if

  contains

    subroutine check_Stress(value)
      implicit none
      real(8) :: value
      if (dataset.eq.'jhu256') then
         if( ( (stress.eq.'dev') .and. (abs(value - (-5.2544371578038401d-3)).gt.tiny(0.d0)) )         &
              .or.( (stress.eq.'abs') .and. (abs(value - (-4.0152351790891661d-3)).gt.tiny(0.d0)) ) )  then
            print*, 'Precision test for stress: Failed.'
            print*, 'Check precision or dataset is not JHU 256'
            print*, value
!            stop
         else
            print*, 'Precision test for stress: Passed'
         end if
      end if
    end subroutine check_Stress

  end subroutine computeStress

end module fourier
