module fourier
  use global
  use, intrinsic :: iso_c_binding
  include '/opt/fftw-3.3.4/include/fftw3.f03' !! FOR KEYWORDS IN FFTW ARGUMENT LIST
contains

subroutine createFilter(filter,scale)
!!$ STATUS : > Passed test for FFTW layout.
!!$          > Passed test for center in spectral space.
!!$          > Check for the axes -- should be in alignment with FFTW output.
!!$ Result : Passed
!!$ Notes  : 1) Create and define variables in a module. 
!!$          2) The filter creation can be parallelized (not needed right now).
!!$          3) Include options for other filters.   
!!$ Layout: The center should lie at an unit offset from the geometric center.
!!$         This is because the DC component occupies the first index. 
!!$ 
!!$  1             129          256 
!!$1  |-------------|------------|
!!$   |             |            |
!!$   |             |            |
!!$   |      A      |     B      |
!!$128|             |            |
!!$129--------------o------------|
!!$   |             |            |
!!$   |      D      |     C      |
!!$   |             |            |
!!$256--------------|------------|

    
implicit none
 integer, intent(in)                                 :: scale
 real(8),dimension(GRID,GRID,GRID),intent(inout)     :: filter

 integer                                   :: center 
 real(8)                                   :: distance
 integer                                   :: i,j,k
 logical                                   :: test=.FALSE.

!!$ Initialize filter :
  filter = 0.d0
  center = 0.5d0*GRID+1.d0

!!$ Create spectrally sharp filter:
  do k = 1,GRID
  do j = 1,GRID
  do i = 1,GRID
          
     distance = sqrt( dble((i-center)**2) &
              +       dble((j-center)**2) &
              +       dble((k-center)**2) )
       if (distance.le.scale) filter(i,j,k) = 1.d0
         
  end do
  end do
  end do

!Sanity checks:
if(test) then
   print*, 'Sanity Checks:'
   write(*,20) filter(center+scale+1,center,center) ! should be 0
   write(*,20) filter(center,center-scale,center)   ! should be 1
   write(*,20) filter(center+scale,center,center)   ! should be 1
   write(*,20) filter(center,center,center+scale)   ! should be 1
   write(*,20) filter(center,center-scale-1,center) ! should be 0
20 format(f8.0)
end if
return
end subroutine createFilter


subroutine fftshift(filter)
!!$ fftshift:
!!$ 
!!$ Layout:                     k 
!!$             E       F       |
!!$           A       B         |
!!$                             /----> j
!!$             H       G      /
!!$           D       C       i
!!$  
!!$
!!$ In clockwise direction from the topleft corner, 
!!$ the front half(i = 129,256) has A,B,C,D
!!$ The rear half (i = 1,128) is labelled as E,F,G,H


!!$ Define arguments:
real(8),dimension(GRID,GRID,GRID),intent(inout):: filter

!!$ Define local variables: 
real(8),allocatable,dimension(:,:,:) :: temp,A,B,C,D,E,F,G,H 
integer                           :: center


center = 0.5d0*GRID+1.d0

allocate(temp(1:(center-1), 1:(center-1) , 1:(center-1) ) )
allocate(A ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
allocate(B ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
allocate(C ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
allocate(D ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
allocate(E ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
allocate(F ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
allocate(G ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )
allocate(H ( 1:(center-1) , 1:(center-1) , 1:(center-1) ) )

!!$Define subarrays:

A = filter( center:GRID  ,  1:(center-1) ,  center:GRID  )
B = filter( center:GRID  ,  center:GRID  ,  center:GRID  )
C = filter( center:GRID  ,  center:GRID  ,  1:(center-1) )
D = filter( center:GRID  ,  1:(center-1) ,  1:(center-1) )
E = filter( 1:(center-1) ,  1:(center-1) ,  center:GRID  )
F = filter( 1:(center-1) ,  center:GRID  ,  center:GRID  )
G = filter( 1:(center-1) ,  center:GRID  ,  1:(center-1) )
H = filter( 1:(center-1) ,  1:(center-1) ,  1:(center-1) )

!!$ For fftshift, the swaps needed are:
!!$ A - G:
temp = A
A = G
G = temp

!!$ C - E:
temp = C
C = E
E = temp

!!$ D - F:
temp = D
D = F
F = temp

!!$ B - H:
temp = B
B = H
H = temp

!!$ Recreate the filter matrix:

filter( center:GRID  , 1:(center-1) , center:GRID  )= A
filter( center:GRID  , center:GRID  , center:GRID  )= B
filter( center:GRID  , center:GRID  , 1:(center-1) )= C
filter( center:GRID  , 1:(center-1) , 1:(center-1) )= D
filter( 1:(center-1) , 1:(center-1) , center:GRID  )= E
filter( 1:(center-1) , center:GRID  , center:GRID  )= F
filter( 1:(center-1) , center:GRID  , 1:(center-1) )= G
filter( 1:(center-1) , 1:(center-1) , 1:(center-1) )= H

deallocate (temp,A,B,C,D,E,F,G,H)
return
end subroutine fftshift


function sharpFilter(array_work,filter)
!!$ STATUS : > Passed test with MATLAB results (10/21/2015).        
!!$ Result : Passed
!!$ Notes  : 1) Runs too slow. 

  implicit none
  
  integer(C_INT)           :: n 
  type(C_PTR)              :: plan

  real(C_DOUBLE),dimension(GRID,GRID,GRID)               :: sharpFilter, array_work,filter
  complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:,:)     :: in_cmplx, out_cmplx
  
  !!$ convert GRID(real*8) to n(C_INT) for FFTW
  !n = GRID
  
  allocate(in_cmplx(GRID,GRID,GRID))
  allocate(out_cmplx(GRID,GRID,GRID))

  in_cmplx = 0.d0; out_cmplx = 0.d0

  in_cmplx = dcmplx(array_work)/(dble(GRID**3)) ! Normalize  ! Create input complex array
  
  
  !!$ Forward Fourier transform
  call dfftw_plan_dft_3d(plan,GRID,GRID,GRID,in_cmplx,out_cmplx,FFTW_FORWARD,FFTW_ESTIMATE)
  call dfftw_execute(plan)
  call dfftw_destroy_plan(plan)
  
  !!$ Apply filter
  out_cmplx = out_cmplx * filter !! double complex * real(kind=8)


  !!$ Inverse Fourier transform
  call dfftw_plan_dft_3d(plan,GRID,GRID,GRID,out_cmplx,in_cmplx,FFTW_BACKWARD,FFTW_ESTIMATE)
  call dfftw_execute(plan)
  call dfftw_destroy_plan(plan)


  sharpFilter = real(in_cmplx)


  deallocate(in_cmplx,out_cmplx)
  return
21 format(a45,5x,2(f30.16))
end function sharpFilter

subroutine computeStress(u,u_f,u_t,tau_ij,T_ij,n_u,n_uu,LES,test,stress)
implicit none

integer,                   intent(in)  :: n_u, n_uu
real(8),dimension(:,:,:,:),intent(in)  :: u, u_f, u_t
real(8),dimension(  :,:,:),intent(in)  :: LES, test
real(8),dimension(:,:,:,:),intent(out) :: tau_ij, T_ij

character(3),intent(in) :: stress 

real(8),allocatable,dimension(:,:,:)   :: dev_t
integer :: i,j,k


! ABSOLUTE STRESS:
if (stress.eq.'abs') then
   k = 1
   do j=1,n_u
   do i=1,n_u
         if (i.ge.j) then
            print *, 'tau(', i, ',', j, ')'
            tau_ij(k,:,:,:) = sharpFilter(u(i,:,:,:)*u(j,:,:,:),LES)       &
                 - u_f(i,:,:,:)*u_f(j,:,:,:)
            print *, 'T(', i, ',', j, ')'
            T_ij(k,:,:,:) = sharpFilter(u_f(i,:,:,:)*u_f(j,:,:,:),test)    &
                 - u_t(i,:,:,:)*u_t(j,:,:,:)
            k = k+1
       end if
   end do
   end do

! DEVIATORIC STRESS:
elseif (stress.eq.'dev')then
   k = 1
   do j=1,n_u
   do i=1,n_u
      if (i.ge.j) then
            print *, 'tau(', i, ',', j, ')'
            tau_ij(k,:,:,:) = sharpFilter(u(i,:,:,:)*u(j,:,:,:),LES)     &
                 - u_f(i,:,:,:)*u_f(j,:,:,:)
            print *, 'T(', i, ',', j, ')'
            T_ij(k,:,:,:) = sharpFilter(u(i,:,:,:)*u(j,:,:,:),test)      &
                 - u_t(i,:,:,:)*u_t(j,:,:,:)
            k = k+1
       end if
   end do
   end do
     allocate(dev_t(GRID,GRID,GRID))

      dev_t = (tau_ij(1,:,:,:) + tau_ij(4,:,:,:) + tau_ij(6,:,:,:)) / 3.d0
      tau_ij(1,:,:,:) = tau_ij(1,:,:,:) - dev_t
      tau_ij(4,:,:,:) = tau_ij(4,:,:,:) - dev_t
      tau_ij(6,:,:,:) = tau_ij(6,:,:,:) - dev_t

      dev_t = (T_ij(1,:,:,:) + T_ij(4,:,:,:) + T_ij(6,:,:,:)) / 3.d0
      T_ij(1,:,:,:) = T_ij(1,:,:,:) - dev_t
      T_ij(4,:,:,:) = T_ij(4,:,:,:) - dev_t
      T_ij(6,:,:,:) = T_ij(6,:,:,:) - dev_t

      ! CHECK DEVIATORIC STRESS: (based on T_ij)                                                                                       
      if (T_ij(1,15,24,10).ne.-5.2544371578038401d-3) then
         print*, 'Precision test for stress: Failed'
         print*, 'Check precision or data for testing is not JHU 256'
         print*, T_ij(1,15,24,10)
         stop
      else
         print*, 'Precision test for deviatoric stress: Passed'
      end if

     deallocate (dev_t)

else
   print*,"Stress must be either abs or dev"
   stop
end if
return
end subroutine computeStress





end module fourier
