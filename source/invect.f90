function invect(x, y, x1, y1, x2, y2)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Test if (x,y) belongs to the segment defined by the points
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
  logical   :: invect  ! function to test if (x,y) belongs to the segment
  real(sgl) :: dnorm   ! difference in distance
  real(sgl) :: epsilon ! tolerance
  real(sgl) :: norm1   ! norm of vector (x-x1,y-y1)
  real(sgl) :: norm2   ! norm of vector (x2-x1,y2-y1)
  real(sgl) :: pscal   ! scalar product of the vector (x-x1,y-y1) with (x2-x1,y2-y1)
  real(sgl) :: x       ! help variable
  real(sgl) :: x1      ! coordinates of intersection points inside the bin
  real(sgl) :: x12     ! coordinates of the vector (x2-x1,y2-y1)
  real(sgl) :: x1x     ! coordinates of the vector (x-x1,y-y1)
  real(sgl) :: x2      ! coordinates of the 2nd summit of the triangle
  real(sgl) :: y       ! coordinates of the point to test
  real(sgl) :: y1      ! variable for final GOE calculation
  real(sgl) :: y12     ! coordinates of the vector (x2-x1,y2-y1)
  real(sgl) :: y1x     ! coordinates of the vector (x-x1,y-y1)
  real(sgl) :: y2      ! variable for final GOE calculation
!
! *************************** Initialisations **************************
!
  invect = .false.
  x1x = x - x1
  y1x = y - y1
  x12 = x2 - x1
  y12 = y2 - y1
!
! norms and scalar product
!
  norm1 = x1x * x1x + y1x * y1x
  norm2 = x12 * x12 + y12 * y12
  epsilon = norm2 / 1.0e14
  pscal = x12 * x1x + y12 * y1x
!
! test if the two vectors have the same direction (if not return)
!
  if (pscal <  - epsilon) return
!
! test if norm of vector (x-x1,y-y1) lower than norm of segment
!
  dnorm = norm1 - norm2
  if (dnorm > epsilon) return
  invect = .true.
  return
end function invect
! Copyright A.J. Koning 2021
