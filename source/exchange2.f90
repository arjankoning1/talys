subroutine exchange2(Zcomp, Ncomp)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculation of two-component exchange terms
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
!   sgl         ! single precision kind
! Variables for preequilibrium initialization
!   maxpar      ! maximal particle number
! Variables for preequilibrium
!   hnu0        ! initial neutron hole number
!   hpi0        ! initial proton hole number
!   pnu0        ! initial neutron number
!   ppi0        ! initial proton number
! Variables for exciton model
!   Gnupi       ! two - component branching ratio
!   Gnuplus     ! two - component branching ratio
!   Gpinu       ! two - component branching ratio
!   Gpiplus     ! two - component branching ratio
!   Lexc        ! exchange term
!   tauexc2     ! lifetime of two - component exciton state
!   wemistot2   ! total two - component emission rate p
!
! *** Declaration of local data
!
  implicit none
  integer   :: h            ! help variable
  integer   :: hnu          ! neutron hole number
  integer   :: hpi          ! proton hole number
  integer   :: Ncomp        ! neutron number index for compound nucleus
  integer   :: p            ! particle number
  integer   :: pnu          ! neutron particle number
  integer   :: ppi          ! proton particle number
  integer   :: Zcomp        ! proton number index for compound nucleus
  real(sgl) :: lambdanupi   ! neutron-proton transition rate for n --> n
  real(sgl) :: lambdanuplus ! neutron transition rate for n --> n+2
  real(sgl) :: lambdapinu   ! proton-neutron transition rate for n --> n
  real(sgl) :: lambdapiplus ! proton transition rate for n --> n+2
  real(sgl) :: tauexc2p     ! lifetime for creation and emission only
  real(sgl) :: tauinv       ! inverse of time
  real(sgl) :: tauinvp      ! help variable
  real(sgl) :: texchange    ! help variable
  real(sgl) :: tnupi        ! help variable
  real(sgl) :: tnuplus      ! help variable
  real(sgl) :: tpinu        ! help variable
  real(sgl) :: tpiplus      ! help variable
  real(sgl) :: tplus        ! transition rate for n --> n+2
  real(sgl) :: wemis        ! total two-component emission rate per exciton number
!
! ************************* Exchange terms *****************************
!
! emissionrate2  : subroutine for two-component emission rate
!
! The strength of the exciton state is calculated
!
  do ppi = ppi0, maxpar
    hpi = hpi0 + ppi - ppi0
    do pnu = pnu0, maxpar
      hnu = hnu0 + pnu - pnu0
      p = ppi + pnu
      h = hpi + hnu
      if (p > maxpar .or. h > maxpar) cycle
      tpiplus = lambdapiplus(Zcomp, Ncomp, ppi, hpi, pnu, hnu)
      tnuplus = lambdanuplus(Zcomp, Ncomp, ppi, hpi, pnu, hnu)
      tpinu = lambdapinu(Zcomp, Ncomp, ppi, hpi, pnu, hnu)
      tnupi = lambdanupi(Zcomp, Ncomp, ppi, hpi, pnu, hnu)
      tplus = tpiplus + tnuplus
      texchange = tpinu + tnupi
      call emissionrate2(Zcomp, Ncomp, ppi, hpi, pnu, hnu)
      wemis = wemistot2(ppi, hpi, pnu, hnu)
      tauinv = tplus + texchange + wemis
      if (tplus == 0.) then
        tauexc2(ppi, hpi, pnu, hnu) = 0.
        Lexc(ppi, hpi, pnu, hnu) = 0.
      else
        tauexc2(ppi, hpi, pnu, hnu) = 1. / tauinv
        tauinvp = tplus + wemis
        if (tauinvp == 0.) then
          tauexc2p = 0.
        else
          tauexc2p = 1. / tauinvp
        endif
        Lexc(ppi, hpi, pnu, hnu) = tauexc2p / tauexc2(ppi, hpi, pnu, hnu)
      endif
      Gpiplus(ppi, hpi, pnu, hnu) = tpiplus * tauexc2(ppi, hpi, pnu, hnu)
      Gnuplus(ppi, hpi, pnu, hnu) = tnuplus * tauexc2(ppi, hpi, pnu, hnu)
      Gpinu(ppi, hpi, pnu, hnu) = tpinu * tauexc2(ppi, hpi, pnu, hnu)
      Gnupi(ppi, hpi, pnu, hnu) = tnupi * tauexc2(ppi, hpi, pnu, hnu)
    enddo
  enddo
  return
end subroutine exchange2
! Copyright A.J. Koning 2021
