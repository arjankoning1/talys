subroutine zbrak(func, x1, x2, n, xb1, xb2, nb)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Bracket the function
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
  implicit none
  integer   :: i       ! counter
  integer   :: n       ! exciton number
  integer   :: nb      ! help variable
  integer   :: nbb     ! help variable
  real(sgl) :: dx      ! increment
  real(sgl) :: fc      ! function value
  real(sgl) :: fp      ! function value
  real(sgl) :: func    ! function for which zero crossing is searched
  real(sgl) :: x       ! help variable
  real(sgl) :: x1      ! coordinates of intersection points inside the bin
  real(sgl) :: x2      ! coordinates of the 2nd summit of the triangle
  real(sgl) :: xb1(nb) ! help variable
  real(sgl) :: xb2(nb) ! help variable
!
! ************************** Bracket function **************************
!
! The function has a zero between xb1 and xb2.
!
! func: function
! fc  : function value
! fp  : function value
!
  nbb = 0
  x = x1
  dx = (x2 - x1) / n
  fc = 0.
  fp = func(x)
  do i = 1, n
    x = x + dx
    fc = func(x)
    if (fc * fp < 0.) then
      nbb = nbb + 1
      xb1(nbb) = x - dx
      xb2(nbb) = x
      if (nbb == nb) exit
    endif
    fp = fc
  enddo
  nb = nbb
  return
end subroutine zbrak
! Copyright A.J. Koning 2021
