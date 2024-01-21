subroutine pol1(x1, x2, y1, y2, x, y)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Interpolation of first order
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               sgl            ! single precision kind
!
! *** Declaration of local data
!
  real(sgl) :: fac            ! factor
  real(sgl) :: x              ! help variable
  real(sgl) :: x1             ! coordinates of intersection points inside the bin
  real(sgl) :: x2             ! coordinates of the 2nd summit of the triangle
  real(sgl) :: y              ! coordinates of the point to test
  real(sgl) :: y1             ! variable for final GOE calculation
  real(sgl) :: y2             ! variable for final GOE calculation
!
! ***************************** Interpolation **************************
!
  fac = (x-x1)/(x2-x1)
  y = y1 + fac * (y2 - y1)
  return
end subroutine pol1
! Copyright A.J. Koning 2021
