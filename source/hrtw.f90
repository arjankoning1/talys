subroutine hrtw(tjl, na, nb, st, sv, v, w, res, numtr, ielas)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : HRTW width fluctuation correction
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
  integer   :: na              ! help variable
  integer   :: nb              ! help variable
  integer   :: numtr           ! number of transmission coefficients
  real(sgl) :: dab             ! help variable
  real(sgl) :: res             ! width fluctuation factor
  real(sgl) :: sv              ! variable for width fluctuation
  real(sgl) :: v(numtr)        ! real volume potential, radius, diffuseness
  real(sgl) :: w(numtr)        ! imaginary volume potential, radius, diffuseness
  real(dbl) :: res1            ! help variable
  real(dbl) :: st              ! denominator of compound nucleus formula
  real(dbl) :: tjl(0:5, numtr) ! transmission coefficients
  real(dbl) :: tt              ! term with transmission coefficients
!
! ************************** HRTW calculation **************************
!
! We use the model of HRTW (1975), revised in 1980.
!
  dab = 0.
  if (ielas == 1) dab = 1.
!
! Final result
!
  tt = tjl(1, na) * tjl(1, nb) / tjl(0, na)
  if (tt /= 0) then
    res1 = v(na) * v(nb) * tjl(0, nb) / sv * (1. + dab * (w(na) - 1.))
    res = real(res1 * st / tt)
  else
    res = 0.
  endif
  return
end subroutine hrtw
! Copyright A.J. Koning 2021
