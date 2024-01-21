function lambdanuplus(Zcomp, Ncomp, ppi, hpi, pnu, hnu)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Neutron transition rates for n --> n+2
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
!   M2nunu         ! square of neutron - neutron matrix element
!   M2nupi         ! square of neutron - proton matrix element
!   M2pinu         ! square of proton - neutron matrix element
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
  integer   :: j            ! counter
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
  real(sgl) :: damp         ! shell damping factor
  real(sgl) :: densh        ! help variable
  real(sgl) :: densp        ! help variable
  real(sgl) :: dExnuh       ! integration bin width
  real(sgl) :: dExnup       ! integration bin width
  real(sgl) :: dExpih       ! integration bin width
  real(sgl) :: dExpip       ! integration bin width
  real(sgl) :: edepth       ! depth of potential well
  real(sgl) :: eopt         ! incident energy
  real(sgl) :: fac1         ! help variable
  real(sgl) :: factor1      ! help variable
  real(sgl) :: factor2      ! help variable
  real(sgl) :: factor3      ! help variable
  real(sgl) :: finitewell   ! correction function for finite well depth
  real(sgl) :: gsn          ! single-particle neutron level density parameter
  real(sgl) :: gsp          ! single-particle proton level density parameter
  real(sgl) :: ignatyuk     ! function for energy dependent level density parameter a
  real(sgl) :: L1nuh        ! integration limit
  real(sgl) :: L1nup        ! integration limits
  real(sgl) :: L1pih        ! integration limit
  real(sgl) :: L1pip        ! integration limit
  real(sgl) :: L2nuh        ! integration limit
  real(sgl) :: L2nup        ! integration limit
  real(sgl) :: L2pih        ! integration limit
  real(sgl) :: L2pip        ! integration limit
  real(sgl) :: lambdanunu1h ! collision probability for neutron-neutron hole
  real(sgl) :: lambdanunu1p ! collision probability for neutron-neutron particle
  real(sgl) :: lambdanuplus ! neutron transition rate for n --> n+2
  real(sgl) :: lambdapinu1h ! collision probability for proton-neutron hole
  real(sgl) :: lambdapinu1p ! collision probability for proton-neutron particle
  real(sgl) :: phdens2      ! function for two-component particle-hole state density
  real(sgl) :: preeqpair    ! pre-equilibrium pairing energy
  real(sgl) :: ratio        ! if 0<ratio<1 then x is between xl1 and xl2
  real(sgl) :: term1        ! help variable
  real(sgl) :: term12       ! help variable
  real(sgl) :: term2        ! help variable
  real(sgl) :: U            ! excitation energy minus pairing energy
  real(sgl) :: uu           ! residual excitation energy
  real(sgl) :: uunuh        ! residual excitation energy for neutron hole
  real(sgl) :: uunup        ! residual excitation energy
  real(sgl) :: uupih        ! residual excitation energy for proton hole
  real(sgl) :: uupip        ! residual excitation energy for proton particle
  real(sgl) :: Weff         ! effective imaginary well depth
  real(dbl) :: lplus        ! help variable
  real(dbl) :: phtot        ! total particle-hole state density
  real(dbl) :: sumnunu1h    ! help variable
  real(dbl) :: sumnunu1p    ! help variable
  real(dbl) :: sumpinu1h    ! help variable
  real(dbl) :: sumpinu1p    ! help variable
  real(dbl) :: termnunu1h   ! help variable
  real(dbl) :: termnunu1p   ! help variable
  real(dbl) :: termpinu1h   ! help variable
  real(dbl) :: termpinu1p   ! help variable
!
! *************************** Transition rates *************************
!
! matrix      : subroutine for matrix element for exciton model
! ignatyuk    : function for energy dependent level density parameter a
! finitewell  : correction function for finite well depth
!
! A. Analytical solution: The transition rates are taken from Kalbach, PRC33, 818 (1986).
!
  lambdanuplus = 0.
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
    fac1 = 2. * n * (n + 1.)
    factor1 = twopihbar * gsn * gsn / fac1
    term1 = U - Apauli2(ppi, hpi, pnu + 1, hnu + 1)
    term2 = U - Apauli2(ppi, hpi, pnu, hnu)
    if (term1 <= 0 .or. term2 <= 0.) return
    term12 = term1 / term2
    if (term12 < 0.01) return
    factor2 = term1 **2 * (term12 **(n - 1))
    factor3 = (pnu + hnu) * gsn * M2nunu + 2. * (ppi + hpi) * gsp * M2nupi
    lplus = factor1 * factor2 * factor3
    lambdanuplus = lplus * finitewell(p + 1, h + 1, U, edepth, surfwell)
  else
!
! B. Numerical solution: Transition rates based on either matrix element (preeqmode=2) or optical model (preeqmode=3).
!
    L1nup = Apauli2(ppi, hpi, pnu + 1, hnu + 1) - Apauli2(ppi, hpi, pnu - 1, hnu)
    L2nup = U - Apauli2(ppi, hpi, pnu - 1, hnu)
    L1nuh = Apauli2(ppi, hpi, pnu + 1, hnu + 1) - Apauli2(ppi, hpi, pnu, hnu - 1)
    L2nuh = U - Apauli2(ppi, hpi, pnu, hnu - 1)
    L1pip = Apauli2(ppi, hpi, pnu + 1, hnu + 1) - Apauli2(ppi - 1, hpi, pnu, hnu)
    L2pip = U - Apauli2(ppi - 1, hpi, pnu, hnu)
    L1pih = Apauli2(ppi, hpi, pnu + 1, hnu + 1) - Apauli2(ppi, hpi - 1, pnu, hnu)
    L2pih = U - Apauli2(ppi, hpi - 1, pnu, hnu)
    if (primary) then
      nexcbins = max(nbins / 2, 2)
    else
      nexcbins = max(nbins / 4, 2)
    endif
    dExnup = (L2nup - L1nup) / nexcbins
    dExnuh = (L2nuh - L1nuh) / nexcbins
    dExpip = (L2pip - L1pip) / nexcbins
    dExpih = (L2pih - L1pih) / nexcbins
    sumnunu1p = 0.
    sumnunu1h = 0.
    sumpinu1p = 0.
    sumpinu1h = 0.
    do i = 1, nexcbins
      uunup = L1nup + (i - 0.5) * dExnup
      uunuh = L1nuh + (i - 0.5) * dExnuh
      uupip = L1pip + (i - 0.5) * dExpip
      uupih = L1pih + (i - 0.5) * dExpih
      if (flaggshell) then
        damp = ignatyuk(Zcomp, Ncomp, uunup, 0) / alev(Zcomp, Ncomp)
        gsp = gp(Zcomp, Ncomp) * damp
        gsn = gn(Zcomp, Ncomp) * damp
      endif
      if (preeqmode == 2) then
        lambdanunu1p = twopihbar * M2nunu * phdens2(Zcomp, Ncomp, 0, 0, 2, 1, gsp, gsn, uunup, edepth, surfwell)
        lambdanunu1h = twopihbar * M2nunu * phdens2(Zcomp, Ncomp, 0, 0, 1, 2, gsp, gsn, uunuh, edepth, surfwell)
        lambdapinu1p = twopihbar * M2pinu * phdens2(Zcomp, Ncomp, 1, 0, 1, 1, gsp, gsn, uupip, edepth, surfwell)
        lambdapinu1h = twopihbar * M2pinu * phdens2(Zcomp, Ncomp, 0, 1, 1, 1, gsp, gsn, uupih, edepth, surfwell)
      else
        do j = 1, 4
          if (j == 1) uu = uunup
          if (j == 2) uu = uunuh
          if (j == 3) uu = uupip
          if (j == 4) uu = uupih
          if (j > 2) then
            Zix = 1
            Nix = 0
            k = 2
          else
            Zix = 0
            Nix = 1
            k = 1
          endif
          eopt = max(uu - S(Zix, Nix, k), - 20.)
          nen = min(10 * numen, int(eopt * 10.))
          if (j <= 2) then
            Weff = Wompfac(1) * wvol(k, nen)
          else
            Weff = Wompfac(2) * wvol(k, nen)
          endif
          if (j == 1) lambdanunu1p = 2. * Weff / hbar
          if (j == 2) then
            densh = phdens2(Zcomp, Ncomp, 0, 0, 1, 2, gsp, gsn, uu, edepth, surfwell)
            densp = phdens2(Zcomp, Ncomp, 0, 0, 2, 1, gsp, gsn, uu, edepth, surfwell)
            if (densp > 1.) then
              ratio = densh / densp
            else
              ratio = 1.
            endif
            lambdanunu1h = 2. * Weff / hbar * ratio
          endif
          if (j == 3) lambdapinu1p = 2. * Weff / hbar
          if (j == 4) then
            densh = phdens2(Zcomp, Ncomp, 0, 1, 1, 1, gsp, gsn, uu, edepth, surfwell)
            densp = phdens2(Zcomp, Ncomp, 1, 0, 1, 1, gsp, gsn, uu, edepth, surfwell)
            if (densp > 1.) then
              ratio = densh / densp
            else
              ratio = 1.
            endif
            lambdapinu1h = 2. * Weff / hbar * ratio
          endif
        enddo
      endif
      termnunu1p = lambdanunu1p * gsn * dExnup
      termnunu1h = lambdanunu1h * gsn * dExnuh
      termpinu1p = lambdapinu1p * gsp * dExpip
      termpinu1h = lambdapinu1h * gsp * dExpih
      if (flaggshell) then
        damp = ignatyuk(Zcomp, Ncomp, U - uunup, 0) / alev(Zcomp, Ncomp)
        gsp = gp(Zcomp, Ncomp) * damp
        gsn = gn(Zcomp, Ncomp) * damp
      endif
      sumnunu1p = sumnunu1p + termnunu1p * phdens2(Zcomp, Ncomp, ppi, hpi, pnu - 1, hnu, gsp, gsn, U - uunup, edepth, surfwell)
      sumnunu1h = sumnunu1h + termnunu1h * phdens2(Zcomp, Ncomp, ppi, hpi, pnu, hnu - 1, gsp, gsn, U - uunuh, edepth, surfwell)
      sumpinu1p = sumpinu1p + termpinu1p * phdens2(Zcomp, Ncomp, ppi - 1, hpi, pnu, hnu, gsp, gsn, U - uupip, edepth, surfwell)
      sumpinu1h = sumpinu1h + termpinu1h * phdens2(Zcomp, Ncomp, ppi, hpi - 1, pnu, hnu, gsp, gsn, U - uupih, edepth, surfwell)
    enddo
    if (flaggshell) then
      damp = ignatyuk(Zcomp, Ncomp, U, 0) / alev(Zcomp, Ncomp)
      gsp = gp(Zcomp, Ncomp) * damp
      gsn = gn(Zcomp, Ncomp) * damp
    endif
    phtot = phdens2(Zcomp, Ncomp, ppi, hpi, pnu, hnu, gsp, gsn, U, edepth, surfwell)
    if (phtot > 0.) then
      lplus = (sumnunu1p + sumnunu1h + sumpinu1p + sumpinu1h) / phtot
      lambdanuplus = real(lplus)
    endif
  endif
  return
end function lambdanuplus
! Copyright A.J. Koning 2021
