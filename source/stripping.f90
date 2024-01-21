subroutine stripping
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Contribution of stripping and pickup reactions
!
! Author    : Arjan Koning and Vivian Dimitriou
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl          ! single precision kind
! Variables for preequilibrium
!   Cstrip       ! adjustable parameter for stripping / pick - up reactions
!   Kph          ! constant for single - particle level density par. (g = A / Kph)
!   pairmodel    ! model for preequilibrium pairing energy
! Variables for main input
!   Atarget      ! mass number of target nucleus
!   k0           ! index of incident particle
!   Ztarget      ! charge number of target nucleus
! Variables for energy grid
!   ebegin       ! first energy point of energy grid
!   egrid        ! outgoing energy grid
!   Einc         ! incident energy in MeV
! Variables for energies
!   eend         ! last energy point of energy grid
!   eninccm      ! center - of - mass incident energy in MeV
!   Etotal       ! total energy of compound system (target + projectile)
! Variables for inverse channel data
!   xsreac       ! reaction cross section
! Variables for nuclides
!   AA           ! mass number of residual nucleus
!   NN           ! neutron number of residual nucleus
!   parskip      ! logical to skip outgoing particle
!   ZZ           ! charge number of residual nucleus
! Constants
!   parA         ! mass number of particle
!   parmass      ! mass of particle in a.m.u.
!   parN         ! neutron number of particle
!   parspin      ! spin of particle
!   parZ         ! charge number of particle
! Variables for ECIS
!   projmass     ! mass of projectile
! Variables for masses
!   S            ! separation energy
! Variables for preequilibrium
!   xspreeqps    ! preequilibrium cross section per particle type and outgoing energy for pickup and stripping
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: key          ! keyword
  logical           :: surfwell     ! flag for surface effects in finite well
  integer           :: A            ! mass number of target nucleus
  integer           :: hnu          ! neutron hole number
  integer           :: hpi          ! proton hole number
  integer           :: i            ! counter
  integer           :: j            ! counter
  integer           :: k            ! designator for particle
  integer           :: l            ! multipolarity
  integer           :: N            ! neutron number of residual nucleus
  integer           :: ndelta       ! number of transferred particles
  integer           :: ndeltanu     ! number of transferred neutrons
  integer           :: ndeltapi     ! number of transferred protons
  integer           :: nen          ! energy counter
  integer           :: pnu          ! neutron particle number
  integer           :: ppi          ! proton particle number
  integer           :: type         ! particle type
  integer           :: Z            ! charge number of target nucleus
  real(sgl)         :: base         ! help variable
  real(sgl)         :: ejec2sp1     ! 2*spin +1 of ejectile
  real(sgl)         :: ejecmass     ! mass of ejectile
  real(sgl)         :: Eout         ! outgoing energy
  real(sgl)         :: Eres         ! total energy of residual system
  real(sgl)         :: Ewell        ! depth of potential well
  real(sgl)         :: factor       ! multiplication factor
  real(sgl)         :: factor1      ! help variable
  real(sgl)         :: gsn          ! single-particle neutron level density parameter
  real(sgl)         :: gsp          ! single-particle proton level density parameter
  real(sgl)         :: Kap          ! alpha enhancement factor
  real(sgl)         :: omegaNT      ! state density function for pickup and stripping
  real(sgl)         :: omegaph      ! particle-hole state density for compound system
  real(sgl)         :: P            ! pairing energy
  real(sgl)         :: phdens2      ! function for two-component particle-hole state density
  real(sgl)         :: preeqpair    ! pre-equilibrium pairing energy
  real(sgl)         :: proj2sp1     ! 2*spin +1 of projectile
  real(sgl)         :: surface      ! well depth for first hole
  real(sgl)         :: term1        ! help variable
  real(sgl)         :: term2        ! help variable
  real(sgl)         :: term3        ! help variable
  real(sgl)         :: term4        ! help variable
  real(sgl)         :: term5        ! help variable
  real(sgl)         :: termps       ! term for pickup and stripping
  real(sgl)         :: V1well       ! depth of potential well
  real(sgl)         :: Va           ! potential drop
  real(sgl)         :: XNT          ! probability of exciting particle-hole pair
!
! ************************** Kalbach model *****************************
!
! The stripping and pickup model is described in C. Kalbach, "Preequilibrium reactions with complex channels",
! Phys. Rev. C71, 034606 (2005).
!
! Factors for pickup and stripping processes. For reactions involving only neutrons and protons this is thus not considered.
!
  projmass = parmass(k0)
  proj2sp1 = 2. * parspin(k0) + 1.
  do type = 1, 6
    if (parskip(type)) cycle
    if (k0 <= 2 .and. type <= 2) cycle
!
! Calculation of terms independent of emission energy.
!
! Kap=5 for outgoing helions was introduced to better fit (n,h) data, i.e. it does not come from the original Kalbach model.
!
    ejecmass = parmass(type)
    ejec2sp1 = 2. * parspin(type) + 1.
    term1 = ejec2sp1 * ejecmass / (proj2sp1 * projmass)
    Kap = 1.
    if ((k0 == 1 .or. k0 == 2)) then
      if (type == 6) Kap = 12.
      if (type == 3 .and. Einc > 80.) Kap = 80. / Einc
      if (type == 5) Kap = 5.
    endif
    if (k0 == 6 .and. (type == 1 .or. type == 2)) Kap = 12. - 11. * max(eninccm - 20., 0.) / eninccm
!
! Extra adjustment for (a,n) and (a,p) cross sections, TENDL-2021
! Also applied to incident tritons and helions
!
    if (k0 >= 4 .and. (type == 1 .or. type == 2)) Kap = Kap * max(4.-Atarget/80., 1.)
    ndelta = abs(parA(k0) - parA(type))
    ndeltapi = parZ(k0) - parZ(type)
    ndeltanu = parN(k0) - parN(type)
    Va = 12.5 * projmass
    term2 = (Kap / projmass) * (projmass / (Einc + Va)) **(2 * ndelta)
!
! Initial configuration for pickup, stripping or t-h charge exchange
!
    ppi = max(ndeltapi, 0)
    hpi = max( - ndeltapi, 0)
    pnu = max(ndeltanu, 0)
    hnu = max( - ndeltanu, 0)
!
! Further terms
!
! adjust : subroutine for energy-dependent parameter adjustment
!
    A = AA(0, 0, type)
    Z = ZZ(0, 0, type)
    N = NN(0, 0, type)
    if (k0 == 1) then
      term3 = (5500. / real(A)) **ndelta
    else
      term3 = (3800. / real(A)) **ndelta
    endif
    if (parA(k0) < parA(type)) term4 = 1. / (80. * eninccm)
    if (parA(k0) > parA(type)) term4 = 1. / (580. * sqrt(eninccm))
    if (parA(k0) == parA(type)) term4 = 1. / (1160. * sqrt(eninccm))
    base = real(2. * Ztarget) / real(Atarget)
    term5 = base **(2 * (parZ(k0) + 2) * hpi + 2 * pnu)
    key = 'cstrip'
    call adjust(Einc, key, 0, 0, type, 0, factor)
    termps = factor * Cstrip(type) * term1 * term2 * term3 * term4 * term5
!
! XNT function
!
    V1well = 17.
    if (k0 == 1) V1well = surface(1, Einc)
    if (k0 == 2) V1well = surface(2, Einc)
    if (k0 == 5 .or. k0 == 6) V1well = 25.
    XNT = sqrt(min(Einc, 100.) / projmass) * 7. / (V1well * Atarget * Atarget) * (pnu **2 + ppi **2 + hnu **2 + 1.5 * hpi **2)
    gsn = N / Kph
    gsp = Z / Kph
    if (ndeltapi == 0) then
      Ewell = V1well * base
    else
      Ewell = V1well
    endif
!
! Calculation of stripping or pickup spectra.
!
    do nen = ebegin(type), eend(type)
      Eout = egrid(nen)
      factor1 = xsreac(type, nen) * Eout
      P = preeqpair(parZ(type), parN(type), ndelta, Etotal, pairmodel)
      Eres = Etotal - S(0, 0, type) - Eout - P
!
! Check if outgoing energy exceeds maximal possible energy
!
      if (Eres < 0.) cycle
!
! Stripping/pick-up terms that depend on emission energy.
! For reactions in which only one hole is left, e.g. (p,d), the finite well function leads to a discontinuity at the well depth.
! Therefore, for this case the well depth is set equal to the Fermi energy.
!
! omegaNT  : state density function for pickup and stripping
! phdens2  : function for two-component particle-hole state density outgoing energy for pickup and stripping
!
      omegaNT = 0.
      surfwell = .false.
      do i = 0, 3
        do j = 0, 3 - i
          omegaph = phdens2(parZ(type), parN(type), ppi + i, hpi + i, pnu + j, hnu + j, gsp, gsn, Eres, Ewell, surfwell)
          omegaNT = omegaNT + XNT **(i + j) * omegaph
        enddo
      enddo
      do i = 0, ppi
        do j = 0, hpi
          do k = 0, pnu
            do l = 0, hnu
              if (i + j + k + l /= 0) omegaNT = omegaNT + &
                phdens2(parZ(type), parN(type), ppi - i, hpi - j, pnu - k, hnu - l, gsp, gsn, Eres, Ewell, surfwell)
            enddo
          enddo
        enddo
      enddo
      xspreeqps(type, nen) = termps * omegaNT * factor1
    enddo
  enddo
  return
end subroutine stripping
! Copyright A.J. Koning 2021
