function rtbis(func, x1, x2, xacc)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Search for zero crossings of the function
!
! Author    : Arjan Koning (adapted from Numerical Recipes)
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               sgl      ! single precision kind
!
! *** Declaration of local data
!
  integer, parameter :: jmax=40   ! maximum j-value
  integer            :: j         ! counter
  real(sgl)          :: dx        ! increment
  real(sgl)          :: f         ! E-Ef
  real(sgl)          :: fmid      ! help variable
  real(sgl)          :: func      ! function for which zero crossing is searched
  real(sgl)          :: rtbis     ! function to search for zero crossings of the function
  real(sgl)          :: x1        ! coordinates of intersection points inside the bin
  real(sgl)          :: x2        ! coordinates of the 2nd summit of the triangle
  real(sgl)          :: xacc      ! help variable
  real(sgl)          :: xmid      ! middle x value
!
! ************************** Search ************************************
!
! func: function for which zero crossing is searched
!
  fmid = func(x2)
  f = func(x1)
  if (f < 0.) then
    rtbis = x1
    dx = x2 - x1
  else
    rtbis = x2
    dx = x1 - x2
  endif
  do j = 1, jmax
    dx = dx * 0.5
    xmid = rtbis + dx
    fmid = func(xmid)
    if (fmid <= 0.) rtbis = xmid
    if (abs(dx) < xacc .or. fmid == 0.) return
  enddo
end function rtbis
! Copyright A.J. Koning 2021
