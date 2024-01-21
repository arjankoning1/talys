subroutine lifetime2(ppi, hpi, pnu, hnu)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculation of lifetime of two-component exciton state
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
!   sgl         ! single precision kind
! Variables for preequilibrium initialization
!   maxpar      ! maximal particle number
! Variables for preequilibrium
!   hnu0        ! initial neutron hole number
!   hpi0        ! initial proton hole number
!   p0          ! initial particle number
!   pnu0        ! initial neutron number
!   ppi0        ! initial proton number
!   Spre        ! time - integrated strength of two - component exciton state
! Variables for exciton model
!   Gnupi       ! two - component branching ratio
!   Gnuplus     ! two - component branching ratio
!   Gpinu       ! two - component branching ratio
!   Gpiplus     ! two - component branching ratio
!   Lexc        ! exchange term
!   PP2         ! total strength
!   tauexc2     ! lifetime of two - component exciton state
!
! *** Declaration of local data
!
  implicit none
  integer   :: hnu   ! neutron hole number
  integer   :: hpi   ! proton hole number
  integer   :: p     ! particle number
  integer   :: pnu   ! neutron particle number
  integer   :: ppi   ! proton particle number
  real(sgl) :: term1 ! help variable
  real(sgl) :: term2 ! help variable
  real(sgl) :: term3 ! help variable
  real(sgl) :: term4 ! help variable
  real(sgl) :: term5 ! help variable
  real(sgl) :: term6 ! help variable
  real(sgl) :: term7 ! help variable
  real(sgl) :: term8 ! help variable
!
! **************** Calculation of total strength ***********************
!
! The strength of the exciton state is calculated
!
  p = ppi+pnu
  if (p == p0) then
    PP2(ppi0, hpi0, pnu0, hnu0) = 1.
  else
    term1 = 0.
    term2 = 0.
    term3 = 0.
    term4 = 0.
    term5 = 0.
    term6 = 0.
    term7 = 0.
    term8 = 0.
    if (ppi > ppi0 .and. hpi > hpi0) then
      term1 = PP2(ppi - 1, hpi - 1, pnu, hnu) * Gpiplus(ppi - 1, hpi - 1, pnu, hnu)
      term4 = PP2(ppi - 1, hpi - 1, pnu, hnu) * Gnuplus(ppi - 1, hpi - 1, pnu, hnu)
    endif
    if (pnu > pnu0 .and. hnu > hnu0) then
      term2 = PP2(ppi, hpi, pnu - 1, hnu - 1) * Gnuplus(ppi, hpi, pnu - 1, hnu - 1)
      term6 = PP2(ppi, hpi, pnu - 1, hnu - 1) * Gpiplus(ppi, hpi, pnu - 1, hnu - 1)
    endif
    if (pnu < maxpar .and. hnu < maxpar)  then
      if (ppi > ppi0 .and. hpi > hpi0) term5 = Gnupi(ppi - 1, hpi - 1, pnu + 1, hnu + 1)
      if (ppi - 1 > ppi0 .and. hpi - 1 > hpi0) term3 = PP2(ppi - 2, hpi - 2, pnu + 1, hnu + 1) * &
        Gpiplus(ppi - 2, hpi - 2, pnu + 1, hnu + 1)
    endif
    if (ppi < maxpar .and. hpi < maxpar) then
      if (pnu > pnu0 .and. hnu > hnu0) term8 = Gpinu(ppi + 1, hpi + 1, pnu - 1, hnu - 1)
      if (pnu - 1 > pnu0 .and. hnu - 1 >= hnu0) term7 = PP2(ppi + 1, hpi + 1, pnu - 2, hnu - 2) * &
        Gnuplus(ppi + 1, hpi + 1, pnu - 2, hnu - 2)
    endif
    PP2(ppi, hpi, pnu, hnu) = term1 + term2 + Lexc(ppi, hpi, pnu, hnu) * ((term3 + term4) * term5 + (term6 + term7) * term8)
  endif
  Spre(ppi, hpi, pnu, hnu) = PP2(ppi, hpi, pnu, hnu) * tauexc2(ppi, hpi, pnu, hnu)
  return
end subroutine lifetime2
! Copyright A.J. Koning 2021
