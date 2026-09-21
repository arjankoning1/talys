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
  integer :: nenmaxmsd        ! maximum outgoing-energy index for MSD
!
! ********* Set parameters and energy grid for MSD calculation *********
!
! interangle : subroutine for intermediate angles by addition theorem for MSD model
!
  maxmsd = maxpar-1
  maxJmsd = 6
  if (Emsdmin == 0 .or. Emsdmin >= Einc) Emsdmin = Einc / 5.
  msdbins2 = msdbins * 2
  nenmaxmsd = max(eend(1), eend(2))
!
! Main MSD arrays
!
  if (allocated(Emsd)) deallocate(Emsd)
  if (allocated(xsdwin)) deallocate(xsdwin)
  if (allocated(msdstep1)) deallocate(msdstep1)
  if (allocated(msdstep)) deallocate(msdstep)
  if (allocated(msdstepint)) deallocate(msdstepint)
  if (allocated(msdsum)) deallocate(msdsum)
  if (allocated(msdtot)) deallocate(msdtot)

  allocate(Emsd(0:msdbins2))
  allocate(xsdwin(0:msdbins2,0:msdbins2,0:maxJmsd))
  allocate(msdstep1(0:numpar,0:nenmaxmsd))
  allocate(msdstep(0:numpar,1:maxmsd,0:nenmaxmsd))
  allocate(msdstepint(0:numpar,1:maxmsd))
  allocate(msdsum(0:numpar))
  allocate(msdtot(0:numpar,0:nenmaxmsd))

  Emsd = 0.
  xsdwin = 0.
  msdstep1 = 0.
  msdstep = 0.
  msdstepint = 0.
  msdsum = 0.
  msdtot = 0.
!
! Angular MSD arrays
!
  if (allocated(xsdw)) deallocate(xsdw)
  if (allocated(msdstepad1)) deallocate(msdstepad1)
  if (allocated(msdstepad)) deallocate(msdstepad)
  if (allocated(msdstepintad)) deallocate(msdstepintad)
  if (allocated(msdtotad)) deallocate(msdtotad)
  if (allocated(msdtotintad)) deallocate(msdtotintad)

  if (flagddx) then
    allocate(xsdw(0:msdbins2,0:msdbins2,0:maxJmsd,0:nanglecont))
    allocate(msdstepad1(0:numpar,0:nenmaxmsd,0:nanglecont))
    allocate(msdstepad(0:numpar,1:maxmsd,0:nenmaxmsd,0:nanglecont))
    allocate(msdstepintad(0:numpar,1:maxmsd,0:nanglecont))
    allocate(msdtotad(0:numpar,0:nenmaxmsd,0:nanglecont))
    allocate(msdtotintad(0:numpar,0:nanglecont))

    xsdw = 0.
    msdstepad1 = 0.
    msdstepad = 0.
    msdstepintad = 0.
    msdtotad = 0.
    msdtotintad = 0.
  endif
!
! Multi-step MSD arrays
!
  if (allocated(xscont1)) deallocate(xscont1)
  if (allocated(xscont)) deallocate(xscont)
  if (allocated(msdstep0)) deallocate(msdstep0)

  if (.not. flagonestep) then
    allocate(xscont1(0:numpar,0:numpar,0:msdbins2,0:msdbins2))
    allocate(xscont(0:numpar,0:numpar,0:msdbins2,0:msdbins2))
    allocate(msdstep0(0:numpar,1:maxmsd,0:msdbins2))

    xscont1 = 0.
    xscont = 0.
    msdstep0 = 0.
  endif
!
! Multi-step angular arrays
!
  if (allocated(nangleint)) deallocate(nangleint)
  if (allocated(xscontad1)) deallocate(xscontad1)
  if (allocated(xscontad)) deallocate(xscontad)
  if (allocated(msdstepad0)) deallocate(msdstepad0)

  if (.not. flagonestep .and. flagddx) then
    allocate(nangleint(0:nanglecont,0:nanglecont,0:nanglecont))
    allocate(xscontad1(0:numpar,0:numpar,0:msdbins2,0:msdbins2,0:nanglecont))
    allocate(xscontad(0:numpar,0:numpar,0:msdbins2,0:msdbins2,0:nanglecont))
    allocate(msdstepad0(0:numpar,1:maxmsd,0:msdbins2,0:nanglecont))

    nangleint = 0
    xscontad1 = 0.
    xscontad = 0.
    msdstepad0 = 0.
  endif
!
! MSD energy grid
!
  dEmsd = (Einc - Emsdmin) / msdbins2
  do nen = 0, msdbins2
    Emsd(nen) = Einc - nen * dEmsd
  enddo
  Emsd(0) = real(Emsd(0) / specmass(parZ(k0), parN(k0), k0))
  if ( .not. flagonestep .and. flagddx) call interangle
  return
end subroutine msdinit
! Copyright A.J. Koning 2021
