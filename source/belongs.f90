function belongs(x, x1, x2)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Test if x belongs to the interval [xlow,xup] or [xup,xlow]
!
! Author    : Stephane Hilaire
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               sgl                ! single precision kind
!
! *** Declaration of local data
!
  logical   :: belongs ! logical function to test if one value is between two others
  real(sgl) :: dnorm   ! difference in distance
  real(sgl) :: epsilon ! tolerance
  real(sgl) :: norm12  ! distance
  real(sgl) :: norm1x  ! distance
  real(sgl) :: pscal   ! scalar product of the vector (x-x1,y-y1) with (x2-x1,y2-y1)
  real(sgl) :: x       ! help variable
  real(sgl) :: x1      ! coordinates of intersection points inside the bin
  real(sgl) :: x1x     ! coordinates of the vector (x-x1,y-y1)
  real(sgl) :: x1x2    ! difference
  real(sgl) :: x2      ! coordinates of the 2nd summit of the triangle
!
! **********************************************************************
!
  belongs = .true.
  x1x = x - x1
  x1x2 = x2 - x1
  norm12 = abs(x1x2)
  norm1x = abs(x1x)
  epsilon = norm12 / 1.0e14
  pscal = x1x * x1x2
  dnorm = norm1x - norm12
!
! test if vector (x-x1,0) and vector (x2-x1,0) have the same direction
!
  if (pscal <  - epsilon) belongs = .false.
!
! test if norm of vector (x-x1,0) lower than norm of vector (x2-x1,0)
!
  if (dnorm > epsilon) belongs = .false.
  return
end function belongs
! Copyright A.J. Koning 2021
