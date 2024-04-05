subroutine splint(xa, ya, y2a, i, x, y)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Spline fit
!
! Author    : Marieke Duijvestijn
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               sgl               ! single precision kind
!
! *** Declaration of local data
!
  integer   :: i      ! level
  integer   :: k      ! designator for particle
  integer   :: khi    ! help variable
  integer   :: klo    ! help variable
  real(sgl) :: a      ! fit variables
  real(sgl) :: b      ! parameter for slope
  real(sgl) :: hsp    ! help variable
  real(sgl) :: x      ! help variable
  real(sgl) :: xa(i)  ! intersection of the segment with the vertical axis linking  (xl,yl) with (xl,yu)
  real(sgl) :: y      ! coordinates of the point to test
  real(sgl) :: y2a(i) ! help variable
  real(sgl) :: ya(i)  ! intersection of the segment with the vertical axis linking  (xl,yl) with (xl,yu)
!
! **********************************************************************
!
  klo = 1
  khi = i
   10 if (khi - klo > 1) then
    k = (khi + klo) / 2
    if (xa(k) > x) then
      khi = k
    else
      klo = k
    endif
    goto 10
  endif
  hsp = xa(khi) - xa(klo)
  if (hsp == 0.) then
    write(*, *) 'bad xa input in splint'
    stop
  endif
  a = (xa(khi) - x) / hsp
  b = (x - xa(klo)) / hsp
  y = a * ya(klo) + b * ya(khi) + ((a **3 - a) * y2a(klo) + (b **3 - b) * y2a(khi)) * (hsp **2) / 6.
  return
end subroutine splint
! Copyright A.J. Koning 2021
