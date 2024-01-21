function lambdanupi(Zcomp, Ncomp, ppi, hpi, pnu, hnu)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Neutron-proton transition rates for n --> n
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
!   dbl            ! double precision kind
! All global variables
!   numen          ! maximum number of outgoing energies
! Variables for preequilibrium
!   flaggshell     ! flag for energy dependence of single particle level density parameter g
!   flagsurface    ! flag for surface effects in exciton model
!   gn             ! single - particle neutron level density parameter
!   gp             ! single - particle proton level density parameter
!   pairmodel      ! model for preequilibrium pairing energy
!   preeqmode      ! designator for pre - equilibrium model
! Variables for main input
!   Ainit          ! mass number of initial compound nucleus
! Variables for level density
!   alev           ! level density parameter
! Variables for energies
!   nbins          ! number of continuum excitation energy bins
! Variables for nuclides
!   primary        ! flag to designate primary (binary) reaction
! Constants
!   hbar           ! Planck's constant / 2.pi in MeV.s
!   twopihbar      ! 2 * pi / hbar
! Variables for masses
!   S              ! separation energy
! Variables for preequilibrium initialization
!   Apauli2        ! two - component Pauli blocking correction
!   Efermi         ! depth of Fermi well
! Variables for preequilibrium
!   Ecomp          ! total energy of composite system
!   Esurf          ! well depth for surface interaction
! Variables for exciton model
!   M2nupi         ! square of neutron - proton matrix element
!   Wompfac        ! adjustable constant for OMP based transition rates
! Variables for exciton model initialization
!   wvol           ! absorption part of the optical potential averaged over the volume
!
! *** Declaration of local data
!
  implicit none
  logical   :: surfwell     ! flag for surface effects in finite well
  integer   :: A            ! mass number of target nucleus
  integer   :: h            ! help variable
  integer   :: hnu          ! neutron hole number
  integer   :: hpi          ! proton hole number
  integer   :: i            ! level
  integer   :: k            ! designator for particle
  integer   :: n            ! exciton number
  integer   :: Ncomp        ! neutron number index for compound nucleus
  integer   :: nen          ! energy counter
  integer   :: nexcbins     ! number of integration bins
  integer   :: Nix          ! neutron number index for residual nucleus
  integer   :: p            ! particle number
  integer   :: pnu          ! neutron particle number
  integer   :: ppi          ! proton particle number
  integer   :: Zcomp        ! proton number index for compound nucleus
  integer   :: Zix          ! charge number index for residual nucleus
  real(sgl) :: Bfactor      ! help variable
  real(sgl) :: damp         ! shell damping factor
  real(sgl) :: densh        ! help variable
  real(sgl) :: densp        ! help variable
  real(sgl) :: dEx          ! excitation energy bin for population arrays
  real(sgl) :: edepth       ! depth of potential well
  real(sgl) :: eopt         ! incident energy
  real(sgl) :: factor1      ! help variable
  real(sgl) :: factor2      ! help variable
  real(sgl) :: factor23     ! help variable
  real(sgl) :: factor3      ! help variable
  real(sgl) :: factor4      ! help variable
  real(sgl) :: finitewell   ! correction function for finite well depth
  real(sgl) :: gsn          ! single-particle neutron level density parameter
  real(sgl) :: gsp          ! single-particle proton level density parameter
  real(sgl) :: ignatyuk     ! function for energy dependent level density parameter a
  real(sgl) :: L1           ! integration limits
  real(sgl) :: L2           ! integration limits
  real(sgl) :: lambdanupi   ! neutron-proton transition rate for n --> n
  real(sgl) :: lambdanupi1p ! collision probability for neutron-proton particle
  real(sgl) :: phdens2      ! function for two-component particle-hole state density
  real(sgl) :: preeqpair    ! pre-equilibrium pairing energy
  real(sgl) :: ratio        ! if 0<ratio<1 then x is between xl1 and xl2
  real(sgl) :: U            ! excitation energy minus pairing energy
  real(sgl) :: uu           ! residual excitation energy
  real(sgl) :: Weff         ! effective imaginary well depth
  real(dbl) :: lnupi        ! help variable
  real(dbl) :: phtot        ! total particle-hole state density
  real(dbl) :: sumnupi1p    ! help variable
  real(dbl) :: termnupi1p   ! help variable
!
! *************************** Transition rates *************************
!
! matrix       : subroutine for matrix element for exciton model
! ignatyuk     : function for energy dependent level density parameter a
! finitewell   : correction function for finite well depth
!
! A. Analytical solution: The transition rates are taken from Kalbach, PRC33, 818 (1986).
!
  lambdanupi = 0.
  h = hpi + hnu
  p = ppi + pnu
  n = p + h
  if (n == 0) return
  A = Ainit - Zcomp - Ncomp
  call matrix(A, n)
  surfwell = flagsurface .and. h == 1 .and. primary
  if (surfwell) then
    edepth = Esurf
  else
    edepth = Efermi
  endif
  gsp = gp(Zcomp, Ncomp)
  gsn = gn(Zcomp, Ncomp)
  if (flaggshell) then
    damp = ignatyuk(Zcomp, Ncomp, Ecomp, 0) / alev(Zcomp, Ncomp)
    gsp = gsp * damp
    gsn = gsn * damp
  endif
  U = Ecomp - preeqpair(Zcomp, Ncomp, n, Ecomp, pairmodel)
  if ((preeqmode /= 2 .and. preeqmode /= 3) .or. n == 1) then
    factor1 = twopihbar * pnu * hnu * M2nupi / n * gsp * gsp
    Bfactor = Apauli2(ppi, hpi, pnu, hnu)
    if (pnu > 0 .and. hnu > 0) Bfactor = max(Bfactor, Apauli2(ppi + 1, hpi + 1, pnu - 1, hnu - 1))
    factor2 = U - Bfactor
    factor3 = U - Apauli2(ppi, hpi, pnu, hnu)
    if (factor2 <= 0..or.factor3 <= 0.) return
    factor23 = factor2 / factor3
    if (factor23 < 0.01) return
    factor4 = 2. * (U - Bfactor)
    if (pnu > 0 .and. hnu > 0) factor4 = factor4 + n * &
      abs(Apauli2(ppi, hpi, pnu, hnu) - Apauli2(ppi + 1, hpi + 1, pnu - 1, hnu - 1))
    lnupi = factor1 * (factor23 **(n - 1)) * factor4
    lambdanupi = lnupi * finitewell(p, h, U, edepth, surfwell)
  else
!
! B. Numerical solution: Transition rates based on either matrix element (preeqmode=2) or optical model (preeqmode=3).
!
    L1 = Apauli2(ppi, hpi, pnu, hnu) - Apauli2(ppi, hpi, pnu - 1, hnu - 1)
    L2 = U - Apauli2(ppi, hpi, pnu - 1, hnu - 1)
    if (primary) then
      nexcbins = max(nbins / 2, 2)
    else
      nexcbins = max(nbins / 4, 2)
    endif
    dEx = (L2 - L1) / nexcbins
    sumnupi1p = 0.
    do i = 1, nexcbins
      uu = L1 + (i - 0.5) * dEx
      if (flaggshell) then
        damp = ignatyuk(Zcomp, Ncomp, uu, 0) / alev(Zcomp, Ncomp)
        gsp = gp(Zcomp, Ncomp) * damp
        gsn = gn(Zcomp, Ncomp) * damp
      endif
      if (preeqmode == 2) then
        lambdanupi1p = twopihbar * M2nupi * phdens2(Zcomp, Ncomp, 1, 1, 0, 0, gsp, gsn, uu, edepth, surfwell)
      else
        Zix = 0
        Nix = 1
        k = 1
        eopt = max(uu - S(Zix, Nix, k), - 20.)
        nen = min(10 * numen, int(eopt * 10.))
        Weff = Wompfac(0) * wvol(k, nen)
        densh = phdens2(Zcomp, Ncomp, 1, 1, 0, 0, gsp, gsn, uu, edepth, surfwell)
        densp = phdens2(Zcomp, Ncomp, 1, 1, 1, 0, gsp, gsn, uu, edepth, surfwell)
        if (densp > 1.) then
          ratio = densh / densp
        else
          ratio = 1.
        endif
        lambdanupi1p = 2. * Weff / hbar * ratio
      endif
      termnupi1p = lambdanupi1p * phdens2(Zcomp, Ncomp, 0, 0, 1, 1, gsp, gsn, uu, edepth, surfwell) * dEx
      if (flaggshell) then
        damp = ignatyuk(Zcomp, Ncomp, uu, 0) / alev(Zcomp, Ncomp)
        gsp = gp(Zcomp, Ncomp) * damp
        gsn = gn(Zcomp, Ncomp) * damp
      endif
      sumnupi1p = sumnupi1p + termnupi1p * phdens2(Zcomp, Ncomp, ppi, hpi, pnu - 1, hnu - 1, gsp, gsn, U - uu, edepth, surfwell)
    enddo
    if (flaggshell) then
      damp = ignatyuk(Zcomp, Ncomp, U, 0) / alev(Zcomp, Ncomp)
      gsp = gp(Zcomp, Ncomp) * damp
      gsn = gn(Zcomp, Ncomp) * damp
    endif
    phtot = phdens2(Zcomp, Ncomp, ppi, hpi, pnu, hnu, gsp, gsn, U, edepth, surfwell)
    if (phtot > 0.) then
      lnupi = sumnupi1p / phtot
      lambdanupi = real(lnupi)
    endif
  endif
  return
end function lambdanupi
! Copyright A.J. Koning 2021
