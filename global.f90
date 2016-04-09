module global
  real(8), parameter :: pi = 4.d0 * atan(1.d0)
  integer, parameter :: GRID = 256


  integer,parameter :: stride = 1 ! Is the ratio between LES(taken as 1) and test scale
  integer,parameter :: skip = 10
  integer,parameter :: X = 1     ! Number of realizations
  integer,parameter :: n_DAMP = 1  ! Number of lambda's

  ! Stencil parameters:
  integer,parameter :: M = 17576              ! Number of training points 3x3x3x9
  integer,parameter :: N = 3403

  ! Bounding Box parameters:
  real,   parameter :: eps = 1e-3 ! Ensures correct integer values for
  ! the ceiling and floor functions. They can give wrong integer values due to conversion of
  ! an integer value into its machine representation.
  
  integer,parameter :: box       = 252
  integer,parameter :: boxSize   = box**3
  integer,parameter :: bigHalf   = ceiling(0.5*real(box) + eps) ! for 8->5
  integer,parameter :: smallHalf = floor(0.5*real(box) + eps)   ! for 8->4
  integer,parameter :: boxCenter = smallHalf * box*(box + 1) + bigHalf
  integer,parameter :: boxLower  = stride * (bigHalf - 1)
  integer,parameter :: boxUpper  = stride * (box - bigHalf)

  ! Test field parameters: 
  integer,parameter :: testSize = 1
  integer,parameter :: testcutSize = stride * (testSize + box) + 1
  integer,parameter :: testLower = stride * bigHalf + 1
  integer,parameter :: testUpper = stride * (bigHalf - 1 + testSize) + 1

  ! Cutout parameters:
  integer,parameter :: lBound = 0.5*(GRID - testcutSize)
  integer,parameter :: uBound = 0.5*(GRID + testcutSize) - 1

end module global
