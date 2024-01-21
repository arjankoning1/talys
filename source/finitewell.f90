function finitewell(p, h, Eex, Ewell, surfwell)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Correction for finite well depth
!
! Author    : Marieke Duijvestijn and Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl        ! single precision kind
! Constants
!   sgn        ! sign
! Variables for preequilibrium initialization
!   Efermi     ! depth of Fermi well
!   ncomb      ! n! / (n!(n - k)!)
!
! *** Declaration of local data
!
  implicit none
  logical   :: surfwell   ! flag for surface effects in finite well
  integer   :: h          ! help variable
  integer   :: jbin       ! number of bins used in averaging procedure equals 2*jbin+1
  integer   :: jwell      ! help variable in loop over potential well depth bins
  integer   :: k          ! designator for particle
  integer   :: n          ! exciton number
  integer   :: p          ! particle number
  real(sgl) :: E          ! incident energy
  real(sgl) :: Eex        ! excitation energy
  real(sgl) :: Ewell      ! depth of potential well
  real(sgl) :: Ewellj     ! potential well depth used in loop
  real(sgl) :: factor     ! multiplication factor
  real(sgl) :: finitewell ! correction function for finite well depth
  real(sgl) :: fwell      ! help variable
  real(sgl) :: fwtsum     ! sum over all weighted contributions to the finite well  depth correction function
  real(sgl) :: widthdis   ! width of the weighting distribution for averaging the  finite well depth correction function
  real(sgl) :: wt         ! (inverse) weight
  real(sgl) :: wtinv      ! (inverse) weight
  real(sgl) :: wtsum      ! sum over all contributions to the weight
  real(sgl) :: x          ! help variable
!
! ********************** Finite well correction ************************
!
! finitewell: correction function for finite well depth
!
! The general finite well correction comes from Betak and Dobes, Z. Phys. A279 (1976) 319.
! For the first p-h interaction, surface effects can be included according to Kalbach, Phys Rev C32, p. 1157 (1985).
!
  finitewell = 1.
  n = p + h
!
! Kalbach surface effects
!
  if (surfwell .and. Ewell < (Efermi - 0.5)) then
    if (p == 0 .and. h == 1) then
      if (Eex > 1.16 * Efermi) then
        finitewell = 0.
        return
      endif
      if (Eex < Ewell) return
      widthdis = Ewell * (Efermi - Ewell) / (2. * Efermi)
      finitewell = 1. / (1. + exp((Eex - Ewell) / widthdis))
    else
      widthdis = Ewell * (Efermi - Ewell) / (2. * Efermi)
      jbin = 4
      wtsum = 0.
      fwtsum = 0.
      do jwell = - jbin, jbin
        Ewellj = Ewell + jwell * widthdis
        if (Ewellj > Efermi .or. Ewellj < 0.) cycle
        x = (Ewellj - Ewell) / widthdis
        if (x <= 80.) then
          wtinv = (1. + exp(x)) * (1. + exp( - x))
          wt = 1. / wtinv
        else
          wt = 0.
        endif
        fwell = 0.
        do k = 0, h
          E = Eex - k * Ewellj
          if (E > 0.) then
            factor = (E / Eex) **(n - 1)
            fwell = fwell + sgn(k) * ncomb(h, k) * factor
          endif
        enddo
        wtsum = wtsum + wt
        fwtsum = fwtsum + wt * fwell
      enddo
      if (wtsum > 0.) finitewell = fwtsum / wtsum
    endif
  else
!
! Betak and Dobes surface effects
!
    if (Eex <= Ewell) return
    if (h == 1 .and. n == 1) return
    do k = 1, h
      E = Eex - k * Ewell
      if (E > 0.) then
        factor = (E / Eex) **(n - 1)
        finitewell = finitewell + sgn(k) * ncomb(h, k) * factor
      endif
    enddo
  endif
  return
end function finitewell
! Copyright A.J. Koning 2021
