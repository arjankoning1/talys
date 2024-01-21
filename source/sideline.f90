function sideline(x, y, xs, ys, xe, ye)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Indicate if (x,y) is on one side (sideline>0.) of the
!
! Author    : Stephane Hilaire
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               sgl      ! single precision kind
!
! *** Declaration of local data
!
  implicit none
  real(sgl) :: a        ! fit variables
  real(sgl) :: b        ! parameter for slope
  real(sgl) :: sideline ! function to indicate if (x,y) is on one side of the segment
  real(sgl) :: x        ! help variable
  real(sgl) :: xdlim    ! tolerance
  real(sgl) :: xe       ! coordinates of the segment second point
  real(sgl) :: xs       ! help variable
  real(sgl) :: y        ! coordinates of the point to test
  real(sgl) :: ydlim    ! tolerance
  real(sgl) :: ye       ! coordinates of the segment second point
  real(sgl) :: ys       ! coordinates of the segment first point
!
! ************************** Initialisation ****************************
!
  xdlim = min(abs(xs), abs(xe))/1.e14
  ydlim = min(abs(ys), abs(ye)) / 1.e14
!
! segment // y-axis2
!
  if (abs(xs - xe) <= xdlim) then
    sideline = x - xs
    if (abs(sideline) <= xdlim) sideline = 0.
    return
  endif
!
! segment // x-axis
!
  if (abs(ys - ye) <= ydlim) then
    sideline = y - ys
    if (abs(sideline) <= ydlim) sideline = 0.
    return
  endif
!
! normal segment
!
  a = (ys - ye) / (xs - xe)
  b = ys - a * xs
  sideline = y - (a * x + b)
  if (abs(sideline) <= 1.e-12) then
    sideline = 0.
    return
  endif
  if (sideline < 0.) then
    sideline = - 1.
    return
  endif
  if (sideline > 0.) then
    sideline = 1.
    return
  endif
end function sideline
! Copyright A.J. Koning 2021
