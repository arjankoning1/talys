subroutine pol2(x1, x2, x3, y1, y2, y3, x, y)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Polynomial interpolation of second order
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
  real(sgl) :: x   ! help variable
  real(sgl) :: x1  ! coordinates of intersection points inside the bin
  real(sgl) :: x2  ! coordinates of the 2nd summit of the triangle
  real(sgl) :: x3  ! coordinates of the 3rd summit of the triangle
  real(sgl) :: y   ! coordinates of the point to test
  real(sgl) :: y1  ! variable for final GOE calculation
  real(sgl) :: y2  ! variable for final GOE calculation
  real(sgl) :: y3  ! variable for final GOE calculation
  real(sgl) :: yy1 ! help variable
  real(sgl) :: yy2 ! help variable
  real(sgl) :: yy3 ! help variable
!
! ******************** Lagrange formula for order 2 ********************
!
!   (x-x2)(x-x3)        (x-x1)(x-x3)        (x-x1)(x-x2)
!    y = --------------.y1 + --------------.y2 + --------------.y3
!   (x1-x2)(x1-x3)      (x2-x1)(x2-x3)      (x3-x1)(x3-x2)
!
  yy1 = (x-x2)*(x-x3)/((x1-x2)*(x1-x3))*y1
  yy2 = (x - x1) * (x - x3) / ((x2 - x1) * (x2 - x3)) * y2
  yy3 = (x - x1) * (x - x2) / ((x3 - x1) * (x3 - x2)) * y3
  y = yy1 + yy2 + yy3
  return
end subroutine pol2
! Copyright A.J. Koning 2021
