subroutine emissionrate2(Zcomp, Ncomp, ppi, hpi, pnu, hnu)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Two-component emission rates for exciton model
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
!   flaggshell     ! flag for energy dependence of single parti
!   flagsurface    ! flag for surface effects in exciton model
!   gn             ! single - particle neutron level density parameter
!   gp             ! single - particle proton level density parameter
!   pairmodel      ! model for preequilibrium pairing energy
!   Rgamma         ! adjustable parameter for pre - equilibrium gamma decay
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
!   parN           ! neutron number of particle
!   parZ           ! charge number of particle
! Variables for masses
!   S              ! separation energy
! Variables for preequilibrium
!   Ecomp          ! total energy of composite system
!   Esurf          ! well depth for surface interaction
!   wemission2     ! two - component emission
! Variables for exciton model
!   wemispart2     ! two - component emission
!   wemistot2      ! total two - component em
! Variables for preequilibrium initialization
!   Efermi         ! depth of Fermi well
!   numparx        ! maximum number of particles
! Variables for exciton model initialization
!   wfac           ! factor for emission rate
!
! *** Declaration of local data
!
  implicit none
  logical   :: surfgam                                                        ! flag for surface effects for photons (alw
  logical   :: surfwell                                                       ! flag for surface effects in finite well
  integer   :: h                                                              ! help variable
  integer   :: hnu                                                            ! neutron hole number
  integer   :: hpi                                                            ! proton hole number
  integer   :: n                                                              ! exciton number
  integer   :: Ncomp                                                          ! neutron number index for compound nucleus
  integer   :: nejec                                                          ! neutron number of leading particle
  integer   :: nen                                                            ! energy counter
  integer   :: nen1                                                           ! energy counter
  integer   :: Nix                                                            ! neutron number index for residual nucleus
  integer   :: nres                                                           ! exciton number for residual system
  integer   :: pnu                                                            ! neutron particle number
  integer   :: pnures                                                         ! help variable
  integer   :: ppi                                                            ! proton particle number
  integer   :: ppires                                                         ! help variable
  integer   :: type                                                           ! particle type
  integer   :: Zcomp                                                          ! proton number index for compound nucleus
  integer   :: zejec                                                          ! charge number of leading particle
  integer   :: Zix                                                            ! charge number index for residual nucleus
  real(sgl) :: branchplus                                                     ! branching ratio for n-2 --> n
  real(sgl) :: branchzero                                                     ! branching ratio for n --> n
  real(sgl) :: damp                                                           ! shell damping factor
  real(sgl) :: dE                                                             ! help variable
  real(sgl) :: edepth                                                         ! depth of potential well
  real(sgl) :: Emax                                                           ! maximal emission energy for particle chan
  real(sgl) :: Eout                                                           ! outgoing energy
  real(sgl) :: Eres                                                           ! total energy of residual system
  real(sgl) :: Ewell                                                          ! depth of potential well
  real(sgl) :: factor                                                         ! multiplication factor
  real(sgl) :: g2E                                                            ! help variable
  real(sgl) :: gs                                                             ! single-particle level density parameter
  real(sgl) :: gsn                                                            ! single-particle neutron level density par
  real(sgl) :: gsp                                                            ! single-particle proton level density para
  real(sgl) :: ignatyuk                                                       ! function for energy dependent level densi
  real(sgl) :: phcomp                                                         ! particle-hole state density for compound
  real(sgl) :: phdens2                                                        ! function for two-component particle-hole
  real(sgl) :: phratio                                                        ! ratio between residual and compound parti
  real(sgl) :: phres                                                          ! particle-hole state density for residual
  real(sgl) :: phres1                                                         ! particle-hole state density for residual
  real(sgl) :: phres2                                                         ! particle-hole state density for residual
  real(sgl) :: preeqpair                                                      ! pre-equilibrium pairing energy
  real(sgl) :: U                                                              ! excitation energy minus pairing energy
  real(sgl) :: Ures                                                           ! excitation energy minus pairing energy
  real(sgl) :: wemissum2(0:numparx, 0:numparx, 0:numparx, 0:numparx, 0:numen) ! two-component emission rate per exciton n
  real(sgl) :: xs                                                             ! help variable
!
! *************************** Emission rates ***************************
!
! ignatyuk   : function for energy dependent level density parameter a
! phdens2    : function for two-component particle-hole state density number
!
  n = ppi + hpi + pnu + hnu
  h = hpi + hnu
  gsp = gp(Zcomp, Ncomp)
  gsn = gn(Zcomp, Ncomp)
  gs = gsp + gsn
  if (flaggshell) then
    damp = ignatyuk(Zcomp, Ncomp, Ecomp, 0) / alev(Zcomp, Ncomp)
    gsp = gsp * damp
    gsn = gsn * damp
    gs = gsp + gsn
  endif
  surfwell = flagsurface .and. h == 1 .and. primary
  if (surfwell) then
    edepth = Esurf
  else
    edepth = Efermi
  endif
  U = Ecomp - preeqpair(Zcomp, Ncomp, n, Ecomp, pairmodel)
  phcomp = phdens2(Zcomp, Ncomp, ppi, hpi, pnu, hnu, gsp, gsn, U, edepth, surfwell)
  wemistot2(ppi, hpi, pnu, hnu) = 0.
  do nen = 0, numen
    wemissum2(ppi, hpi, pnu, hnu, nen) = 0.
  enddo
  do type = 0, 6
    wemispart2(type, ppi, hpi, pnu, hnu) = 0.
    do nen = 0, numen
      wemission2(type, ppi, hpi, pnu, hnu, nen) = 0.
    enddo
    if (parskip(type)) cycle
    Zix = Zindex(Zcomp, Ncomp, type)
    Nix = Nindex(Zcomp, Ncomp, type)
    gsp = gp(Zix, Nix)
    gsn = gn(Zix, Nix)
    zejec = parZ(type)
    nejec = parN(type)
    if (type > 2) then
      Ewell = Efermi
    else
      Ewell = edepth
    endif
    if (phcomp <= 1.e-10) cycle
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
        wemispart2(type, ppi, hpi, pnu, hnu) = wemispart2(type, ppi, hpi, pnu, hnu) + &
 &        wemission2(type, ppi, hpi, pnu, hnu, nen1) * dE
        exit
      endif
      if (flaggshell) then
        damp = ignatyuk(Zix, Nix, Eres, 0) / alev(Zix, Nix)
        gsp = gp(Zix, Nix) * damp
        gsn = gn(Zix, Nix) * damp
        gs = gsp + gsn
      endif
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
          phres1 = 0.5 * (phdens2(Zcomp, Ncomp, ppi - 1, hpi - 1, pnu, hnu, gsp, gsn, U, Efermi, surfgam) + &
            phdens2(Zcomp, Ncomp, ppi, hpi, pnu - 1, hnu - 1, gsp, gsn, U, Efermi, surfgam))
          phres2 = phdens2(Zcomp, Ncomp, ppi, hpi, pnu, hnu, gsp, gsn, U, Efermi, surfgam)
          g2E = gs * gs * Eout
          if (n >= 2) then
            branchplus = g2E / (gs * (n - 2) + g2E)
          else
            branchplus = 0.
          endif
          branchzero = gs * n / (gs * n + g2E)
          phratio = (branchplus * phres1 + branchzero * phres2) / phcomp
          wemission2(type, ppi, hpi, pnu, hnu, nen) = Rgamma * factor * phratio
          wemissum2(ppi, hpi, pnu, hnu, nen) = wemission2(type, ppi, hpi, pnu, hnu, nen)
        endif
      else
!
! 2. Particle emission rates
!
        ppires = ppi - zejec
        pnures = pnu - nejec
        nres = n - zejec - nejec
        if (ppires < 0 .or. pnures < 0 .or. h == 0) cycle
        Ures = max(Eres - preeqpair(Zix, Nix, nres, Eres, pairmodel), preeqpair(Zix, Nix, nres, Eres, pairmodel))
        phres = phdens2(Zix, Nix, ppires, hpi, pnures, hnu, gsp, gsn, Ures, Ewell, surfwell)
        phratio = phres / phcomp
        wemission2(type, ppi, hpi, pnu, hnu, nen) = factor * phratio
        wemissum2(ppi, hpi, pnu, hnu, nen) = wemissum2(ppi, hpi, pnu, hnu, nen) + wemission2(type, ppi, hpi, pnu, hnu, nen)
      endif
!
! *** Integration of emission rates over all energies and particles ****
!
      wemispart2(type, ppi, hpi, pnu, hnu) = wemispart2(type, ppi, hpi, pnu, hnu) + &
 &      wemission2(type, ppi, hpi, pnu, hnu, nen) * deltaE(nen)
    enddo
    wemistot2(ppi, hpi, pnu, hnu) = wemistot2(ppi, hpi, pnu, hnu) + wemispart2(type, ppi, hpi, pnu, hnu)
  enddo
!
! Prevent divergence of pre-equilibrium gamma cross sections in case of absence of particle competition.
!
  if (n > 1) then
    if (wemistot2(ppi, hpi, pnu, hnu) == wemispart2(0, ppi, hpi, pnu, hnu)) then
      wemistot2(ppi, hpi, pnu, hnu) = 0.
      wemispart2(0, ppi, hpi, pnu, hnu) = 0.
    endif
    do nen = 0, numen
      if (wemissum2(ppi, hpi, pnu, hnu, nen) == wemission2(0, ppi, hpi, pnu, hnu, nen)) wemission2(0, ppi, hpi, pnu, hnu, nen) = 0.
    enddo
  endif
  return
end subroutine emissionrate2
! Copyright A.J. Koning 2021
