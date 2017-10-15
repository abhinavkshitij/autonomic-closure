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
  ! FORM: subroutine createFilter(filter, scale,[filterOption-default:Sharp])
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

    if (present(filterOption).and.filterOption.eq.'Gauss') then
      !  Create Gauss filter:
      do k = 1,f_GRID
         do j = 1,f_GRID
            do i = 1,f_GRID

               distance = sqrt( dble((i - center)**2) &
                    +           dble((j - center)**2) &
                    +           dble((k - center)**2) )

               filter(i,j,k) = exp(-5.d-1 * (distance / scale)**2)

            end do
         end do
      end do
    elseif (present(filterOption).and.filterOption.eq.'Box') then
      !  Create Box filter:
      do k = 1,f_GRID
         do j = 1,f_GRID
            do i = 1,f_GRID

               distance = sqrt( dble((i - center)**2) &
                    +           dble((j - center)**2) &
                    +           dble((k - center)**2) )

                filter(i,j,k) = sin( (distance+eps) / (scale+eps) ) / &
               &                   ( (distance+eps) / (scale+eps) )

            end do
         end do
      end do
      elseif (present(filterOption).and.filterOption.eq.'Tri') then
      !  Create Tri filter:
      do k = 1,f_GRID
         do j = 1,f_GRID
            do i = 1,f_GRID

               distance = sqrt( dble((i - center)**2) &
                    +           dble((j - center)**2) &
                    +           dble((k - center)**2) )

                filter(i,j,k) = sin( (distance+eps) / (scale+eps) * 0.5d0) ** 2 / &
               &                   ( (distance+eps) / (scale+eps) * 0.5d0) ** 2

            end do
         end do
      end do
      elseif (present(filterOption).and.filterOption.eq.'Custom') then
      !  Create Custom filter:
      do k = 1,f_GRID
         do j = 1,f_GRID
            do i = 1,f_GRID

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
      do k = 1,f_GRID
         do j = 1,f_GRID
            do i = 1,f_GRID

               distance = sqrt( dble((i - center)**2) &
                    +           dble((j - center)**2) &
                    +           dble((k - center)**2) )
               if (distance.le.scale) filter(i,j,k) = 1.d0

            end do
         end do
      end do

      ! Check spherical symmetry:
      if ( (filter (center+scale+1,center,center) .eq. 0)      .and.  &
           (filter (center,center-scale,center)   .eq. 1)      .and.  &
           (filter (center+scale,center,center)   .eq. 1)      .and.  &
           (filter (center,center,center+scale)   .eq. 1)      .and.  &
           (filter (center,center-scale-1,center) .eq. 0) )    then
      else
         print*, 'Failed creating spherically symmetric filter ... Aborting'
         stop
      end if
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
  
  subroutine fftshift(filter)
    !
    !    ..ARRAY ARGUMENTS..
    real(8),dimension(:,:,:),intent(inout):: filter
    !
    !    ..WORK ARRAYS..
    real(8),allocatable,dimension(:,:,:) :: temp
    real(8),allocatable,dimension(:,:,:) :: A,B,C,D,E,F,G,H 

    

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


    A = filter( center:f_GRID  ,  1:(center-1) ,  center:f_GRID  )
    B = filter( center:f_GRID  ,  center:f_GRID  ,  center:f_GRID  )
    C = filter( center:f_GRID  ,  center:f_GRID  ,  1:(center-1) )
    D = filter( center:f_GRID  ,  1:(center-1) ,  1:(center-1) )
    E = filter( 1:(center-1) ,  1:(center-1) ,  center:f_GRID  )
    F = filter( 1:(center-1) ,  center:f_GRID  ,  center:f_GRID  )
    G = filter( 1:(center-1) ,  center:f_GRID  ,  1:(center-1) )
    H = filter( 1:(center-1) ,  1:(center-1) ,  1:(center-1) )

    ! SWAP SUBARRAYS: (_SHIFT)
    temp = A ; A = G ; G = temp
    temp = C ; C = E ; E = temp
    temp = D ; D = F ; F = temp
    temp = B ; B = H ; H = temp

    ! RECREATE FILTER WITH FFTSHIFT:
    filter( center:f_GRID  , 1:(center-1) , center:f_GRID  ) = A
    filter( center:f_GRID  , center:f_GRID  , center:f_GRID  ) = B
    filter( center:f_GRID  , center:f_GRID  , 1:(center-1) ) = C
    filter( center:f_GRID  , 1:(center-1) , 1:(center-1) ) = D
    filter( 1:(center-1) , 1:(center-1) , center:f_GRID  ) = E
    filter( 1:(center-1) , center:f_GRID  , center:f_GRID  ) = F
    filter( 1:(center-1) , center:f_GRID  , 1:(center-1) ) = G
    filter( 1:(center-1) , 1:(center-1) , 1:(center-1) ) = H

    deallocate (temp,A,B,C,D,E,F,G,H)
    return
  end subroutine fftshift


  !****************************************************************
  !                         FORWARD FFT
  !****************************************************************

  !----------------------------------------------------------------
  ! USE: Take forward FFT
  !      
  ! FORM:  function forwardFFT(array)  
  !       
  ! BEHAVIOR: array_work must be of size(f_GRID,f_GRID,f_GRID)
  !         
  ! STATUS : 
  ! Notes  : 
  !         
  !----------------------------------------------------------------

!   function forwardFFT(array)
!     implicit none
!     !
!     !    ..ARRAY ARGUMENTS..
!     real(C_DOUBLE),dimension(:,:,:) :: array, forwardFFT
!     !
!     !    ..WORK ARRAYS..
!     complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:,:):: in_cmplx, out_cmplx
!     !
!     !    ..LOCAL VARS..
!     type(C_PTR) :: plan
   
!     allocate(in_cmplx(f_GRID,f_GRID,f_GRID))
!     allocate(out_cmplx(f_GRID,f_GRID,f_GRID))

!     in_cmplx = dcmplx (array) / (dble(f_GRID**3)) 

!     ! FFT:
!     call dfftw_plan_dft_3d(plan,f_GRID,f_GRID,f_GRID,in_cmplx,out_cmplx,FFTW_FORWARD,FFTW_ESTIMATE)
!     call dfftw_execute(plan)    
!     call dfftw_destroy_plan(plan)

!     forwardFFT = abs(out_cmplx) 
    
!   end function forwardFFT

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
  ! Notes  : 
  !         1) dfftw_execute(plan) takes 2 secs
  !         2) Uses FFTW libraries.
  !         3) Use C_DOUBLE, C_INT, etc for C-compatibility on any platform.
  !         4) Can use ABS(in_cmplx) but residual imaginary values 
  !            are not exactly 0.d0.
  !         5) Normalization depends on total number of points
  !            including the padded section. So f_GRID^3 is the
  !            correct normalization factor.
  !----------------------------------------------------------------

  function sharpFilter(array_work,filter)
    implicit none
    !
    !    ..ARRAY ARGUMENTS..
    real(C_DOUBLE),dimension(:,:,:) :: array_work, filter
    real(C_DOUBLE),dimension(f_GRID,f_GRID,f_GRID) :: sharpFilter 
    !
    !    ..WORK ARRAYS..
    complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:,:):: in_cmplx, out_cmplx
    !
    !    ..LOCAL VARS..
    type(C_PTR) :: plan
   
    allocate(in_cmplx (f_GRID,f_GRID,f_GRID))
    allocate(out_cmplx(f_GRID,f_GRID,f_GRID))


    in_cmplx(1:i_GRID,1:j_GRID,1:k_GRID) = dcmplx (array_work(1:i_GRID,1:j_GRID,1:k_GRID)) / (dble(f_GRID**3)) 


    ! FT:
    call dfftw_plan_dft_3d(plan,f_GRID,f_GRID,f_GRID,in_cmplx,out_cmplx,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute(plan)    
    call dfftw_destroy_plan(plan)

!      ! ****
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
!     ! ****

    out_cmplx = out_cmplx * filter


    ! IFFT:
    call dfftw_plan_dft_3d(plan,f_GRID,f_GRID,f_GRID,out_cmplx,in_cmplx,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_execute(plan)
    call dfftw_destroy_plan(plan)

 
    sharpFilter = real(in_cmplx) 
    
  end function sharpFilter

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
