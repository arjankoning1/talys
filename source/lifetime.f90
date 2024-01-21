subroutine lifetime(Zcomp, Ncomp, p, h)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculation of lifetime of exciton state
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl          ! single precision kind
! Variables for exciton model
!   depletion    ! depletion factor at each stage
!   tauexc       ! lifetime of exciton state
!   wemistot     ! total emission rate per exciton number
! Variables for preequilibrium
!   p0           ! initial particle number
!
! *** Declaration of local data
!
  implicit none
  integer   :: h          ! help variable
  integer   :: Ncomp      ! neutron number index for compound nucleus
  integer   :: p          ! particle number
  integer   :: Zcomp      ! proton number index for compound nucleus
  real(sgl) :: lambdaplus ! transition rate for n --> n+2
  real(sgl) :: term1      ! help variable
  real(sgl) :: term2      ! help variable
  real(sgl) :: tplus      ! transition rate for n --> n+2
!
! ******************** Never-come-back solution ************************
!
! The lifetime of the exciton state is calculated, see the manual.
!
  if (p == p0) then
    depletion(p, h) = 1.
  else
    tplus = lambdaplus(Zcomp, Ncomp, p - 1, h - 1)
    term1 = tplus + wemistot(p - 1, h - 1)
    if (term1 /= 0.) then
      depletion(p, h) = depletion(p - 1, h - 1) * tplus / term1
    else
      depletion(p, h) = 0.
    endif
  endif
  term2 = lambdaplus(Zcomp, Ncomp, p, h) + wemistot(p, h)
  if (term2 /= 0.) then
    tauexc(p, h) = depletion(p, h) / term2
  else
    tauexc(p, h) = 0.
  endif
  return
end subroutine lifetime
! Copyright A.J. Koning 2021
