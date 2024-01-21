subroutine msdinit
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialization of MSD model parameters
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
! Variables for output
!   flagddx        ! flag for output of double - differential cross sections
! Variables for preequilibrium
!   Emsdmin        ! minimal outgoing energy for MSD calculation
!   msdbins        ! number of energy points for DWBA calculation for MSD
!   flagonestep    ! flag for continuum one - step direct only
! Variables for main input
!   k0             ! index of incident particle
! Variables for energy grid
!   Einc           ! incident energy in MeV
! Constants
!   parN           ! neutron number of particle
!   parZ           ! charge number of particle
! Variables for masses
!   specmass       ! specific mass for residual nucleus
! Variables for preequilibrium initialization
!   maxpar         ! maximal particle number
! Variables for MSD
!   dEmsd          ! energy bin for MSD
!   Emsd           ! minimal outgoing energy for MSD calculation
!   maxJmsd        ! maximal spin for MSD calculation
!   maxmsd         ! number of MSD steps
!   msdbins2       ! number of energy points for MSD calculation
!
! *** Declaration of local data
!
  implicit none
  integer :: nen              ! energy counter
!
! ********* Set parameters and energy grid for MSD calculation *********
!
! interangle : subroutine for intermediate angles by addition theorem for MSD model
!
  maxmsd = maxpar-1
  maxJmsd = 6
  if (Emsdmin == 0 .or. Emsdmin >= Einc) Emsdmin = Einc / 5.
  msdbins2 = msdbins * 2
  dEmsd = (Einc - Emsdmin) / msdbins2
  do nen = 0, msdbins2
    Emsd(nen) = Einc - nen * dEmsd
  enddo
  Emsd(0) = real(Emsd(0) / specmass(parZ(k0), parN(k0), k0))
  if ( .not. flagonestep .and. flagddx) call interangle
  return
end subroutine msdinit
! Copyright A.J. Koning 2021
