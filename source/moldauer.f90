subroutine moldauer(numtjl, na, nb, st, npmold, xmo, tav, vnu, product, res, numtr, ielas)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Moldauer width fluctuation correction
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
  integer   :: ielas           ! designator for elastic channel
  integer   :: im              ! counter
  integer   :: na              ! help variable
  integer   :: nb              ! help variable
  integer   :: npmold          ! variables for Gauss-Legendre integration
  integer   :: numtjl          ! number of transmission coefficients
  integer   :: numtr           ! number of transmission coefficients
  real(sgl) :: dab             ! help variable
  real(sgl) :: factor          ! multiplication factor
  real(sgl) :: factora         ! help variable
  real(sgl) :: factorb         ! help variable
  real(sgl) :: product(npmold) ! product used in final Moldauer calculation
  real(sgl) :: res             ! width fluctuation factor
  real(sgl) :: vnu(numtr)      ! number of degrees of freedom
  real(sgl) :: x               ! help variable
  real(sgl) :: xmo(npmold)     ! variables for Gauss-Legendre integration
  real(dbl) :: st              ! denominator of compound nucleus formula
  real(dbl) :: tav(numtr)      ! average transmission coefficients
!
! **************** Calculation of Moldauer integral ********************
!
! Initialization
!
  res = 0.
  dab = 0.
  if (ielas == 1) dab = 1.
  factor = 1. + 2. * dab / vnu(na)
  factora = real(2. * tav(na) / (st * vnu(na)))
  factorb = real(2. * tav(nb) / (st * vnu(nb)))
!
! Loop over integration points
!
  do im = 1, npmold
    x = xmo(im)
!
! Final result
!
    if (nb == numtjl + 1) then
!
! Special case for capture
!
      res = res + product(im) / (1. + x * factora)
    else
      res = res + product(im) * factor / ((1. + x * factora) * (1. + x * factorb))
    endif
  enddo
  return
end subroutine moldauer
! Copyright A.J. Koning 2021
