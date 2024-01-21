subroutine unfac(l, rho, rhoc, amun, vl, ps)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Penetrability factor (vl) and phase shift (ps) from NJOY
!
! Author    : Gilles Noguere
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
  integer   :: l    ! multipolarity
  real(sgl) :: amun ! number of degrees of freedom for l,j
  real(sgl) :: ps   ! phase shift
  real(sgl) :: r2   ! value
  real(sgl) :: r4   ! value
  real(sgl) :: rho  ! integrated level density
  real(sgl) :: rhoc ! help variable
  real(sgl) :: vl   ! penetrability factor
!
! **********************************************************************
!
  vl = 1.
  ps = 1.
  r2 = rho * rho
  if (l == 0) then
    vl = amun
    ps = rhoc
  else if (l == 1) then
    vl = amun * r2 / (1. + r2)
    ps = rhoc - atan(rhoc)
  else if (l == 2) then
    r4 = r2 * r2
    vl = amun * r4 / (9. + 3. * r2 + r4)
    ps = rhoc - atan(3. * rhoc / (3. - rhoc * rhoc))
  endif
  return
end subroutine unfac
! Copyright A.J. Koning 2021
