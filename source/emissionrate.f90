subroutine emissionrate(Zcomp, Ncomp, p, h)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Emission rates for exciton model
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
!   sgl            ! single precision kind
! All global variables
!   numen          ! maximum number of outgoing energies
! Variables for preequilibrium
!   flaggshell     ! flag for energy dependence of single particle level density parameter g
!   flagsurface    ! flag for surface effects in exciton model
!   g              ! single - particle level density parameter
!   pairmodel      ! model for preequilibrium pairing energy
!   Rgamma         ! adjustable parameter for pre - equilibrium gamma decay
! Variables for main input
!   Atarget        ! mass number of target nucleus
! Variables for level density
!   alev           ! level density parameter
! Variables for energy grid
!   deltaE         ! energy bin around outgoing energies
!   ebegin         ! first energy point of energy grid
!   egrid          ! outgoing energy grid
! Variables for energies
!   eend           ! last energy point of energy grid
! Variables for inverse channel data
!   xsreac         ! reaction cross section
! Variables for nuclides
!   Nindex         ! neutron number index for residual nucleus
!   parskip        ! logical to skip outgoing particle
!   primary        ! flag to designate primary (binary) reaction
!   Zindex         ! charge number index for residual nucleus
! Constants
!   parA           ! mass number of particle
! Variables for masses
!   S              ! separation energy
! Variables for preequilibrium
!   Ecomp          ! total energy of composite system
!   Esurf          ! well depth for surface interaction
!   wemission      ! emission rate per particle, exciton number
! Variables for exciton model
!   wemispart      ! emission rate per particle and exciton number
!   wemistot       ! total emission rate per exciton number
! Variables for preequilibrium initialization
!   Efermi         ! depth of Fermi well
!   numparx        ! maximum number of particles
! Variables for exciton model initialization
!   Qfactor        ! Q - factor for neutron / proton distinction
!   wfac           ! factor for emission rate
!
! *** Declaration of local data
!
  implicit none
  logical   :: surfgam                                 ! flag for surface effects for photons (always false)
  logical   :: surfwell                                ! flag for surface effects in finite well
  integer   :: aejec                                   ! mass number of leading particle
  integer   :: h                                       ! help variable
  integer   :: n                                       ! exciton number
  integer   :: Ncomp                                   ! neutron number index for compound nucleus
  integer   :: nen                                     ! energy counter
  integer   :: nen1                                    ! energy counter
  integer   :: Nix                                     ! neutron number index for residual nucleus
  integer   :: nres                                    ! exciton number for residual system
  integer   :: p                                       ! particle number
  integer   :: pres                                    ! help variable
  integer   :: type                                    ! particle type
  integer   :: Zcomp                                   ! proton number index for compound nucleus
  integer   :: Zix                                     ! charge number index for residual nucleus
  real(sgl) :: branchplus                              ! branching ratio for n-2 --> n
  real(sgl) :: branchzero                              ! branching ratio for n --> n
  real(sgl) :: dE                                      ! help variable
  real(sgl) :: edepth                                  ! depth of potential well
  real(sgl) :: Emax                                    ! maximal emission energy for particle channel
  real(sgl) :: Eout                                    ! outgoing energy
  real(sgl) :: Eres                                    ! total energy of residual system
  real(sgl) :: Ewell                                   ! depth of potential well
  real(sgl) :: factor                                  ! multiplication factor
  real(sgl) :: g2E                                     ! help variable
  real(sgl) :: gs                                      ! single-particle level density parameter
  real(sgl) :: gsg                                     ! single-particle level density parameter
  real(sgl) :: ignatyuk                                ! function for energy dependent level density parameter a
  real(sgl) :: phcomp                                  ! particle-hole state density for compound system
  real(sgl) :: phcompg                                 ! particle-hole state density for compound system (gamma)
  real(sgl) :: phdens                                  ! function for particle-hole state density
  real(sgl) :: phratio                                 ! ratio between residual and compound particle-hole density
  real(sgl) :: phres                                   ! particle-hole state density for residual system
  real(sgl) :: phres1                                  ! particle-hole state density for residual system
  real(sgl) :: phres2                                  ! particle-hole state density for residual system
  real(sgl) :: preeqpair                               ! pre-equilibrium pairing energy
  real(sgl) :: Qfac                                    ! Q-factor for neutron/proton distinction
  real(sgl) :: U                                       ! excitation energy minus pairing energy
  real(sgl) :: Ures                                    ! excitation energy minus pairing energy
  real(sgl) :: wemissum(0:numparx, 0:numparx, 0:numen) ! emission rate per exciton number and energy
  real(sgl) :: xs                                      ! help variable
!
! *************************** Emission rates ***************************
!
! ignatyuk    : function for energy dependent level density parameter a
! phdens      : function for particle-hole state density
!
! The emission rates are derived from detailed balance, see the manual.
!
  n = p + h
  gs = g(Zcomp, Ncomp)
  gsg = Atarget / 13.
  if (flaggshell) gs = gs * ignatyuk(Zcomp, Ncomp, Ecomp, 0) / alev(Zcomp, Ncomp)
  surfwell = flagsurface .and. h == 1 .and. primary
  if (surfwell) then
    edepth = Esurf
  else
    edepth = Efermi
  endif
  U = Ecomp - preeqpair(Zcomp, Ncomp, n, Ecomp, pairmodel)
  phcomp = phdens(Zcomp, Ncomp, p, h, gs, U, edepth, surfwell)
  phcompg = phdens(Zcomp, Ncomp, p, h, gsg, U, edepth, surfwell)
  wemistot(p, h) = 0.
  do nen = 0, numen
    wemissum(p, h, nen) = 0.
  enddo
  do type = 0, 6
    wemispart(type, p, h) = 0.
    do nen = 0, numen
      wemission(type, p, h, nen) = 0.
    enddo
    if (parskip(type)) cycle
    Zix = Zindex(Zcomp, Ncomp, type)
    Nix = Nindex(Zcomp, Ncomp, type)
    gs = g(Zix, Nix)
    aejec = parA(type)
    if (primary .and. type > 0) then
      Qfac = Qfactor(type, p)
    else
      Qfac = 1.
    endif
    if (type > 2) then
      Ewell = Efermi
    else
      Ewell = edepth
    endif
    if (phcomp <= 1.e-10) cycle
    if (type == 0 .and. phcompg <= 1.e-10) cycle
    do nen = ebegin(type), eend(type)
      Eout = egrid(nen)
      xs = xsreac(type, nen)
      factor = wfac(type) * xs * Eout
      Eres = Ecomp - S(Zcomp, Ncomp, type) - Eout
!
! Avoid anomalies at the drip line
!
      if (type <= 2) Eres = min(Eres, Ecomp)
!
! Check if outgoing energy exceeds maximal possible energy
!
      if (Eres < 0.) then
        nen1 = nen - 1
!
! Correction in integration for last outgoing energy.
!
        Eout = egrid(nen1)
        Emax = Ecomp - S(Zcomp, Ncomp, type)
        dE = Emax - (Eout + 0.5 * deltaE(nen1))
        wemispart(type, p, h) = wemispart(type, p, h) + wemission(type, p, h, nen1) * dE
        exit
      endif
      if (flaggshell) gs = g(Zix, Nix) * ignatyuk(Zix, Nix, Eres, 0) / alev(Zix, Nix)
!
! 1. Gamma emission rates
!
! Gamma emission is only included for primary pre-equilibrium decay and n <= 7.
!
      if (type == 0) then
        if (primary .and. n <= 7) then
          factor = factor * Eout
          U = max(Eres - preeqpair(Zcomp, Ncomp, n, Eres, pairmodel), preeqpair(Zcomp, Ncomp, n, Eres, pairmodel))
          surfgam = .false.
          phres1 = phdens(Zcomp, Ncomp, p - 1, h - 1, gsg, U, Efermi, surfgam)
          phres2 = phdens(Zcomp, Ncomp, p, h, gsg, U, Efermi, surfgam)
          g2E = gsg * gsg * Eout
          if (n >= 2) then
            branchplus = g2E / (gsg * (n - 2) + g2E)
          else
            branchplus = 0.
          endif
          branchzero = gsg * n / (gsg * n + g2E)
          phratio = (branchplus * phres1 + branchzero * phres2) / phcompg
          wemission(type, p, h, nen) = Rgamma * factor * phratio
          wemissum(p, h, nen) = wemission(type, p, h, nen)
        endif
      else
!
! 2. Particle emission rates
!
        pres = p - aejec
        nres = n - aejec
        if (pres < 0 .or. h == 0) cycle
        Ures = max(Eres - preeqpair(Zix, Nix, nres, Eres, pairmodel), preeqpair(Zix, Nix, nres, Eres, pairmodel))
        phres = phdens(Zix, Nix, pres, h, gs, Ures, Ewell, surfwell)
        phratio = phres / phcomp
        wemission(type, p, h, nen) = factor * phratio * Qfac
        wemissum(p, h, nen) = wemissum(p, h, nen) + wemission(type, p, h, nen)
      endif
!
! *** Integration of emission rates over all energies and particles ****
!
      wemispart(type, p, h) = wemispart(type, p, h) + wemission(type, p, h, nen) * deltaE(nen)
    enddo
    wemistot(p, h) = wemistot(p, h) + wemispart(type, p, h)
  enddo
!
! Prevent divergence of pre-equilibrium gamma cross sections in case of absence of particle competition.
!
  if (n > 1) then
    if (wemistot(p, h) == wemispart(0, p, h)) wemispart(0, p, h) = 0.
    do nen = 0, numen
      if (wemissum(p, h, nen) == wemission(0, p, h, nen)) then
        wemissum(p, h, nen) = 0.
        wemission(0, p, h, nen) = 0.
      endif
    enddo
  endif
  return
end subroutine emissionrate
! Copyright A.J. Koning 2021
