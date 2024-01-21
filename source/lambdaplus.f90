function lambdaplus(Zcomp, Ncomp, p, h)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Transition rates for n --> n+2
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
!   g              ! single - particle level density parameter
!   pairmodel      ! model for preequilibrium pairing energy
!   preeqmode      ! designator for pre - equilibrium model
! Variables for main input
!   Ainit          ! mass number of initial compound nucleus
!   Ninit          ! neutron number of initial compound nucleus
!   Zinit          ! charge number of initial compound nucleus
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
!   Apauli         ! two - component Pauli blocking correction
!   Efermi         ! depth of Fermi well
! Variables for preequilibrium
!   Ecomp          ! total energy of composite system
!   Esurf          ! well depth for surface interaction
! Variables for exciton model
!   M2             ! square of matrix element
!   Wompfac        ! adjustable constant for OMP based transition rates
! Variables for exciton model initialization
!   wvol           ! absorption part of the optical potential averaged over the volume
!
! *** Declaration of local data
!
  implicit none
  logical   :: surfwell   ! flag for surface effects in finite well
  integer   :: A          ! mass number of target nucleus
  integer   :: h          ! help variable
  integer   :: i          ! level
  integer   :: j          ! counter
  integer   :: k          ! designator for particle
  integer   :: n          ! exciton number
  integer   :: Ncomp      ! neutron number index for compound nucleus
  integer   :: nen        ! energy counter
  integer   :: nexcbins   ! number of integration bins
  integer   :: Nix        ! neutron number index for residual nucleus
  integer   :: p          ! particle number
  integer   :: Zcomp      ! proton number index for compound nucleus
  integer   :: Zix        ! charge number index for residual nucleus
  real(sgl) :: densh      ! help variable
  real(sgl) :: densp      ! help variable
  real(sgl) :: dExhh      ! energy bin for holes
  real(sgl) :: dExpp      ! energy bin for particles
  real(sgl) :: edepth     ! depth of potential well
  real(sgl) :: eopt       ! incident energy
  real(sgl) :: factor1    ! help variable
  real(sgl) :: finitewell ! correction function for finite well depth
  real(sgl) :: gs         ! single-particle level density parameter
  real(sgl) :: ignatyuk   ! function for energy dependent level density parameter a
  real(sgl) :: L1h        ! integration limit
  real(sgl) :: L1p        ! integration limit
  real(sgl) :: L2h        ! integration limit
  real(sgl) :: L2p        ! integration limit
  real(sgl) :: lambda1h   ! collision probability for hole
  real(sgl) :: lambda1p   ! collision probability for particle
  real(sgl) :: lambdaplus ! transition rate for n --> n+2
  real(sgl) :: Nratio     ! help variable
  real(sgl) :: phdens     ! function for particle-hole state density
  real(sgl) :: preeqpair  ! pre-equilibrium pairing energy
  real(sgl) :: ratio      ! if 0<ratio<1 then x is between xl1 and xl2
  real(sgl) :: term1      ! help variable
  real(sgl) :: term12     ! help variable
  real(sgl) :: term2      ! help variable
  real(sgl) :: U          ! excitation energy minus pairing energy
  real(sgl) :: uu         ! residual excitation energy
  real(sgl) :: uuh        ! residual excitation energy
  real(sgl) :: uup        ! residual excitation energy
  real(sgl) :: Weff(2)    ! effective imaginary well depth
  real(sgl) :: ww         ! weight
  real(sgl) :: Zratio     ! help variable
  real(dbl) :: lplus      ! help variable
  real(dbl) :: phtot      ! total particle-hole state density
  real(dbl) :: sum1h      ! help variable
  real(dbl) :: sum1p      ! help variable
  real(dbl) :: term1h     ! help variable
  real(dbl) :: term1p     ! help variable
!
! *************************** Transition rates *************************
!
! matrix     : subroutine for matrix element for exciton model
! ignatyuk   : function for energy dependent level density parameter a
!
  lambdaplus = 0.
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
  gs = g(Zcomp, Ncomp)
  if (flaggshell) gs = gs * ignatyuk(Zcomp, Ncomp, Ecomp, 0) / alev(Zcomp, Ncomp)
  U = Ecomp - preeqpair(Zcomp, Ncomp, n, Ecomp, pairmodel)
!
! A. Analytical solution
!    First we calculate the transition density for the infinite well as given by Oblozinsky et al (Nuc Phys A226, 347 (1974))
!    and then multiply it by the finite well correction function.
!
! finitewell   : correction function for finite well depth
!
  if ((preeqmode /= 2 .and. preeqmode /= 3) .or. n == 1) then
    factor1 = twopihbar * M2 * 0.5 / (n + 1.) * gs **3
    term1 = U - Apauli(p + 1, h + 1)
    term2 = U - Apauli(p, h)
    if (term1 <= 0 .or. term2 <= 0.) return
    term12 = term1 / term2
    if (term12 < 0.01) return
    lplus = factor1 * term1 **2 * (term12 **(n - 1))
    lambdaplus = lplus * finitewell(p + 1, h + 1, U, edepth, surfwell)
  else
!
! B. Numerical solution: Transition rates based on either matrix element (preeqmode=2) or optical model (preeqmode=3).
!
    L1p = Apauli(p + 1, h + 1) - Apauli(p - 1, h)
    L2p = U - Apauli(p - 1, h)
    L1h = Apauli(p + 1, h + 1) - Apauli(p, h - 1)
    L2h = U - Apauli(p, h - 1)
    if (primary) then
      nexcbins = max(nbins / 2, 2)
    else
      nexcbins = max(nbins / 4, 2)
    endif
    dExpp = (L2p - L1p) / nexcbins
    dExhh = (L2h - L1h) / nexcbins
    sum1p = 0.
    sum1h = 0.
    do i = 1, nexcbins
      uup = L1p + (i - 0.5) * dExpp
      uuh = L1h + (i - 0.5) * dExhh
      if (flaggshell) gs = g(Zcomp, Ncomp) * ignatyuk(Zcomp, Ncomp, uup, 0) / alev(Zcomp, Ncomp)
      if (preeqmode == 2) then
        lambda1p = twopihbar * M2 * phdens(Zcomp, Ncomp, 2, 1, gs, uup, edepth, surfwell)
        lambda1h = twopihbar * M2 * phdens(Zcomp, Ncomp, 1, 2, gs, uuh, edepth, surfwell)
      else
        Zratio = real(Zinit) / Ainit
        Nratio = real(Ninit) / Ainit
        do j = 1, 2
          if (j == 1) then
            uu = uup
          else
            uu = uuh
          endif
          do k = 1, 2
            if (k == 1) then
              Zix = 0
              Nix = 1
            else
              Zix = 1
              Nix = 0
            endif
            eopt = max(uu - S(Zix, Nix, k), - 20.)
            nen = min(10 * numen, int(eopt * 10.))
            Weff(k) = 0.5 * Wompfac(0) * wvol(k, nen)
          enddo
          ww = Zratio * Weff(2) + Nratio * Weff(1)
          if (j == 1) then
            lambda1p = 2. * ww / hbar
          else
            densh = phdens(Zcomp, Ncomp, 1, 2, gs, uu, edepth, surfwell)
            densp = phdens(Zcomp, Ncomp, 2, 1, gs, uu, edepth, surfwell)
            if (densp > 1.) then
              ratio = densh / densp
            else
              ratio = 1.
            endif
            lambda1h = 2. * ww / hbar * ratio
          endif
        enddo
      endif
      term1p = lambda1p * gs * dExpp
      term1h = lambda1h * gs * dExhh
      if (flaggshell) gs = g(Zcomp, Ncomp) * ignatyuk(Zcomp, Ncomp, U - uup, 0) / alev(Zcomp, Ncomp)
      sum1p = sum1p + term1p * phdens(Zcomp, Ncomp, p - 1, h, gs, U - uup, edepth, surfwell)
      sum1h = sum1h + term1h * phdens(Zcomp, Ncomp, p, h - 1, gs, U - uuh, edepth, surfwell)
    enddo
    if (flaggshell) gs = g(Zcomp, Ncomp) * ignatyuk(Zcomp, Ncomp, U, 0) / alev(Zcomp, Ncomp)
    phtot = phdens(Zcomp, Ncomp, p, h, gs, U, edepth, surfwell)
    if (phtot > 0.) then
      lplus = (sum1p + sum1h) / phtot
      lambdaplus = real(lplus)
    endif
  endif
  return
end function lambdaplus
! Copyright A.J. Koning 2021
