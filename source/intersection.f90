subroutine intersection(xs, ys, xe, ye, xl, yl, xu, yu, xr, yr, ir)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculation of the intersection points defined by a segment
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
  logical   :: belongs ! logical function to test if one value is between two others
  integer   :: ir      ! loop counter
  real(sgl) :: a       ! fit variables
  real(sgl) :: b       ! parameter for slope
  real(sgl) :: xa      ! intersection of the segment with the vertical axis linking  (xl,yl) with (xl,yu)
  real(sgl) :: xb      ! intersection of the segment with the vertical axis linking  (xu,yl) with (xu,yu)
  real(sgl) :: xc      ! intersection of the segment with the horizontal axis linking  (xl,yl) with (xu,yl)
  real(sgl) :: xd      ! intersection of the segment with the horizontal axis linking  (xl,yu) with (xu,yu)
  real(sgl) :: xdlim   ! tolerance
  real(sgl) :: xe      ! coordinates of the segment second point
  real(sgl) :: xl      ! coordinates defining a bin
  real(sgl) :: xr(4)   ! help array for storing segment intersection points  with the bin sides
  real(sgl) :: xs      ! help variable
  real(sgl) :: xu      ! upper x-bin coordinates
  real(sgl) :: ya      ! intersection of the segment with the vertical axis linking  (xl,yl) with (xl,yu)
  real(sgl) :: yb      ! intersection of the segment with the vertical axis linking  (xu,yl) with (xu,yu)
  real(sgl) :: yc      ! intersection of the segment with the horizontal axis linking  (xl,yl) with (xu,yl)
  real(sgl) :: yd      ! intersection of the segment with the horizontal axis linking  (xl,yu) with (xu,yu)
  real(sgl) :: ydlim   ! tolerance
  real(sgl) :: ye      ! coordinates of the segment second point
  real(sgl) :: yl      ! lower and upper y-bin coordinates
  real(sgl) :: yr(4)   ! help array for storing segment intersection points  with the bin sides
  real(sgl) :: ys      ! coordinates of the segment first point
  real(sgl) :: yu      ! lower and upper y-bin coordinates
!
! **************************** Initialisation **************************
!
  ir = 0
  xdlim = min(abs(xs), abs(xe)) / 1.e14
  ydlim = min(abs(ys), abs(ye)) / 1.e14
!
! segment // y-axis ==> 2 cases
!
  if (abs(xs - xe) <= xdlim) then
    if (belongs(xs, xu, xl)) then
      ir = 2
      xr(1) = xs
      xr(2) = xs
      yr(1) = yu
      yr(2) = yl
    endif
    return
  endif
!
! segment // x-axis
!
  if (abs(ys - ye) <= ydlim) then
    if (belongs(ys, yu, yl)) then
      ir = 2
      xr(1) = xu
      xr(2) = xl
      yr(1) = ys
      yr(2) = ys
    endif
    return
  endif
!
! normal segment crossing the bin
!
  a = (ys - ye) / (xs - xe)
  b = ys - a * xs
!
! intersections with vertical axis
!
  xa = xl
  ya = a * xl + b
  xb = xu
  yb = a * xu + b
!
! intersections with horizontal axis
!
  yc = yl
  xc = (yl - b) / a
  yd = yu
  xd = (yu - b) / a
!
! only keep intersections inside the bin limits
!
  ir = 0
  if (belongs(ya, yu, yl)) then
    ir = ir + 1
    xr(ir) = xa
    yr(ir) = ya
  endif
  if (belongs(yb, yu, yl)) then
    ir = ir + 1
    xr(ir) = xb
    yr(ir) = yb
  endif
  if (belongs(xc, xu, xl)) then
    ir = ir + 1
    xr(ir) = xc
    yr(ir) = yc
  endif
  if (belongs(xd, xu, xl)) then
    ir = ir + 1
    xr(ir) = xd
    yr(ir) = yd
  endif
  return
end subroutine intersection
! Copyright A.J. Koning 2021
