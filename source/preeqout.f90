subroutine preeqout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of pre-equilibrium cross sections
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
!   sgl             ! single precision kind
! Variables for preequilibrium
!   flag2comp       ! flag for two - component pre - equilibrium model
!   flaggshell      ! flag for energy dependence of single particle level den
!   g               ! single - particle level density parameter
!   gn              ! single - particle neutron level density parameter
!   gp              ! single - particle proton level density parameter
!   pairmodel       ! model for preequilibrium pairing energy
! Variables for level density
!   alev            ! level density parameter
! Variables for energy grid
!   ebegin          ! first energy point of energy grid
!   egrid           ! outgoing energy grid
! Variables for energies
!   eend            ! last energy point of energy grid
!   Etotal          ! total energy of compound system (target + projectile)
! Variables for nuclides
!   parskip         ! logical to skip outgoing particle
! Constants
!   parname         ! name of particle
! Variables for incident channel
!   xspreeq         ! preeq. cross section per particle typ and outgoing energye
!   xspreeqsum      ! total preequilibrium cross section summed over particles
!   xspreeqtot      ! preequilibrium cross section per particle type
! Variables for preequilibrium initialization
!   Efermi          ! depth of Fermi well
!   maxexc          ! maximal exciton number
!   maxpar          ! maximal particle number
!   RnJ             ! spin distribution for particle - hole stat
!   RnJsum          ! (2J + 1) * sum over spin distributions
! Variables for preequilibrium
!   Esurf           ! well depth for surface interaction
!   preeqnorm       ! preequilibrium normalizati
!   xspreeqbu       ! preequilibrium cross section per particle type and outgoing energy for brea
!   xspreeqki       ! preequilibrium cross section per particle type and outgoing energy for knoc
!   xspreeqps       ! preequilibrium cross section per particle type and outgoing energy for pick
!   xspreeqtotbu    ! preequilibrium cross section per particle type for breakup
!   xspreeqtotki    ! preequilibrium cross section per particle type for knockout and inelas
!   xspreeqtotps    ! preequilibrium cross section per particle type for pickup and strippin
!   xsstep          ! preeq. cross section per particle type, stage and outgoing E
!   xssteptot       ! preequilibrium cross section per particle type and stage
!
! *** Declaration of local data
!
  implicit none
  logical   :: surfwell  ! flag for surface effects in finite well
  integer   :: h         ! help variable
  integer   :: J         ! spin of level
  integer   :: k         ! designator for particle
  integer   :: n         ! exciton number
  integer   :: nen       ! energy counter
  integer   :: Nix       ! neutron number index for residual nucleus
  integer   :: p         ! particle number
  integer   :: type      ! particle type
  integer   :: Zix       ! charge number index for residual nucleus
  real(sgl) :: damp      ! shell damping factor
  real(sgl) :: Eex       ! excitation energy
  real(sgl) :: gs        ! single-particle level density parameter
  real(sgl) :: gsn       ! single-particle neutron level density parameter
  real(sgl) :: gsp       ! single-particle proton level density parameter
  real(sgl) :: ignatyuk  ! function for energy dependent level density parameter a
  real(sgl) :: nonpski   ! preequilibrium cross section without pickup etc.
  real(sgl) :: phdens    ! function for particle-hole state density
  real(sgl) :: phdens2   ! function for two-component particle-hole state density
  real(sgl) :: preeqpair ! pre-equilibrium pairing energy
!
! ************************ Pre-equilibrium *****************************
!
! 1. Output of particle-hole state densities
!
! ignatyuk  : function for energy dependent level density parameter a
! phdens2   : function for two-component particle-hole state density
!
  Zix = 0
  Nix = 0
  surfwell = .false.
  write(*, '(/" ++++++++++ PARTIAL STATE DENSITIES ++++++++++")')
  if ( .not. flag2comp) then
    write(*, '(/" Particle-hole state densities"/)')
    write(*, '("     Ex    P(n=3)     gs    ", 8(i1, "p", i1, "h", 6x)/)') ((h + k, h, k = 0, 1), h = 1, 4)
    do nen = 1, int(Etotal)
      Eex = real(nen)
      gs = g(0, 0)
      if (flaggshell) gs = gs * ignatyuk(Zix, Nix, Eex, 0) / alev(0, 0)
      write(*, '(1x, 3f8.3, 8es10.3)') Eex, preeqpair(Zix, Nix, 3, Eex, pairmodel), gs, &
 &      ((phdens(Zix, Nix, h + k, h, gs, Eex, Efermi, surfwell), k = 0, 1), h = 1, 4)
    enddo
  else
    write(*, '(/" Particle-hole state densities", /)')
    write(*, '("     Ex    P(n=3)    gp      gn   ", 26x, "Configuration p(p) h(p) p(n) h(n)")')
    write(*, '(28x, 9(2x, 4i2)/)') 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, &
 &    2, 1, 0, 0, 0, 0, 2, 1, 2, 2, 0, 0, 0, 0, 2, 2, 1, 1, 1, 1
    do nen = 1, int(Etotal)
      Eex = real(nen)
      gsp = gp(0, 0)
      gsn = gn(0, 0)
      if (flaggshell) then
        damp = ignatyuk(Zix, Nix, Eex, 0) / alev(0, 0)
        gsp = gsp * damp
        gsn = gsn * damp
      endif
      write(*, '(1x, 4f8.3, 9es10.3)') Eex, preeqpair(Zix, Nix, 3, Eex, pairmodel), gsp, gsn, &
 &      phdens2(Zix, Nix, 1, 1, 0, 0, gsp, gsn, Eex, Efermi, surfwell), &
 &      phdens2(Zix, Nix, 0, 0, 1, 1, gsp, gsn, Eex, Efermi, surfwell), &
        phdens2(Zix, Nix, 1, 1, 1, 0, gsp, gsn, Eex, Efermi, surfwell), &
        phdens2(Zix, Nix, 1, 0, 1, 1, gsp, gsn, Eex, Efermi, surfwell), &
        phdens2(Zix, Nix, 2, 1, 0, 0, gsp, gsn, Eex, Efermi, surfwell), &
        phdens2(Zix, Nix, 0, 0, 2, 1, gsp, gsn, Eex, Efermi, surfwell), &
        phdens2(Zix, Nix, 2, 2, 0, 0, gsp, gsn, Eex, Efermi, surfwell), &
        phdens2(Zix, Nix, 0, 0, 2, 2, gsp, gsn, Eex, Efermi, surfwell), &
        phdens2(Zix, Nix, 1, 1, 1, 1, gsp, gsn, Eex, Efermi, surfwell)
    enddo
  endif
  write(*, '(/" Particle-hole spin distributions"/)')
  write(*, '("   n    ", 9(" J=", i2, "       "), " Sum"/)') (J, J = 0, 8)
  do n = 1, maxexc
    write(*, '(1x, i3, 10es12.4)') n, (RnJ(n, J), J = 0, 8), RnJsum(n)
  enddo
  write(*, '(/" Effective well depth for surface interaction:", f12.5, " MeV")') Esurf
!
! 2. Output of pre-equilibrium cross sections
!
  write(*, '(/" ++++++++++ TOTAL PRE-EQUILIBRIUM CROSS SECTIONS ++++++++++")')
  if (preeqnorm /= 0.) write(*, '(/" Pre-equilibrium normalization factor: ", f8.5/)') preeqnorm
  do type = 0, 6
    if (parskip(type)) cycle
    if (ebegin(type) >= eend(type)) cycle
    write(*, '(/" Pre-equilibrium cross sections for ", a8/)') parname(type)
    write(*, '("     E     Total", 6("       p=", i1), "     Total  Pickup/Strip Knockout Breakup", /)') (p, p = 1, maxpar)
    do nen = ebegin(type), eend(type)
      nonpski = xspreeq(type, nen) - xspreeqps(type, nen) - xspreeqki(type, nen) - xspreeqbu(type, nen)
      write(*, '(1x, f8.3, 11es10.3)') egrid(nen), xspreeq(type, nen), (xsstep(type, p, nen), p = 1, maxpar), nonpski, &
 &      xspreeqps(type, nen), xspreeqki(type, nen), xspreeqbu(type, nen)
    enddo
    nonpski = xspreeqtot(type) - xspreeqtotps(type) - xspreeqtotki(type) - xspreeqtotbu(type)
    write(*, '(/9x, 11es10.3)') xspreeqtot(type), (xssteptot(type, p), p = 1, maxpar), nonpski, xspreeqtotps(type), &
 &    xspreeqtotki(type), xspreeqtotbu(type)
    write(*, '(/" Integrated:", f12.5/)') xspreeqtot(type)
  enddo
  write(*, '(" Total pre-equilibrium cross section:", f12.5)') xspreeqsum
  return
end subroutine preeqout
! Copyright A.J. Koning 2021
