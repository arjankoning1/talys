function func1(a, b, c, t1, t2, dab)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculation of the x,x1,x2 terms for GOE
!
! Author    : Stephane Hilaire
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               dbl  , & ! double precision kind
               sgl      ! single precision kind
!
! *** Declaration of local data
!
  implicit none
  real(sgl) :: a     ! fit variables
  real(sgl) :: b     ! parameter for slope
  real(sgl) :: c     ! curvature of neck
  real(sgl) :: dab   ! help variable
  real(sgl) :: func1 ! function for GOE
  real(dbl) :: f1    ! help variable
  real(dbl) :: f2    ! help variable
  real(dbl) :: t1    ! help variable
  real(dbl) :: t2    ! help variable
!
! ****************** Calculation of the x,x1,x2 terms ******************
!
  f1 = dab * (1. - t1) * ((a / (1. + t1 * a) + b / (1. + t1 * b) + 2 * c / (1. - t1 * c)) **2)
  f2 = (1. + dab) * (a * (1. + a) / ((1. + t1 * a) * (1. + t2 * a)) + &
 &  b * (1. + b) / ((1. + t1 * b) * (1. + t2 * b)) + 2 * c * (1. - c) / ((1. - t1 * c) * (1. - t2 * c)))
  func1 = real(f1 + f2)
  return
end function func1
! Copyright A.J. Koning 2021
