subroutine preeqtotal
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Total pre-equilibrium cross sections
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
!   sgl               ! single precision kind
! Variables for preequilibrium
!   breakupmodel      ! model for break - up reaction: 1. Kalbach 2. Avrigeanu
!   flag2comp         ! flag for two - component pre - equilibrium model
!   pespinmodel       ! model for pre - equilibrium or compound spin distribution
! Variables for numerics
!   xseps             ! limit for cross sections
! Variables for main input
!   k0                ! index of incident particle
! Variables for energy grid
!   deltaE            ! energy bin around outgoing energies
!   ebegin            ! first energy point of energy grid
!   Etop              ! top of outgoing energy bin
! Variables for energies
!   eend              ! last energy point of energy grid
!   eoutdis           ! outgoing energy of discrete state reaction
!   nendisc           ! last discrete bin
! Variables for incident channel
!   xsdirdisc         ! direct cross section for discrete state direct cross section
!   xsdirdiscsum      ! total direct cross section
!   xsdirdisctot      ! direct cross section summed over discrete states
!   xsgrsum           ! sum over giant resonance cross sections
!   xspreeq           ! preeq. cross section per particle type and outgoing energy
!   xspreeqsum        ! total preequilibrium cross section summed over particles
!   xspreeqtot        ! preequilibrium cross section per particle type
!   xsreacinc         ! reaction cross section for incident channel
! Variables for nuclides
!   parskip           ! logical to skip outgoing particle
! Constants
!   parN              ! neutron number of particle
!   parZ              ! charge number of particle
! Variables for level density
!   Nlast             ! last discrete level
! Variables for preequilibrium initialization
!   maxJph            ! maximal spin for particle - hole states
!   maxpar            ! maximal particle number
! Variables for preequilibrium
!   p0                ! initial particle number
!   pnu0              ! initial neutron number
!   ppi0              ! initial proton number
!   preeqnorm         ! preequilibrium normalizati
!   xsEB              ! elastic breakup cross section
!   xsflux            ! cross section flux
!   xspreeqbu         ! preequilibrium cross section per particle type and outgoing energy for break-up
!   xspreeqdisc       ! preequilibrium cross section for discrete state
!   xspreeqdiscsum    ! total preequilibrium cross section for discrete states
!   xspreeqdisctot    ! preequilibrium cross section summed over discrete states
!   xspreeqJP         ! preeq. cross section per particle type, outgoing energy, J, P
!   xspreeqki         ! preequilibrium cross section per particle type and outgoing energy for knock-out
!   xspreeqps         ! preequilibrium cross section per particle type and outgoing energy for pick-up
!   xspreeqtotbu      ! preequilibrium cross section per particle type for breakup
!   xspreeqtotki      ! preequilibrium cross section per particle type for knockout and inelastic
!   xspreeqtotps      ! preequilibrium cross section per particle type for pickup and stripping
!   xsstep            ! preeq. cross section per particle type, stage and outgoing E
!   xsstep2           ! two - component preequilibrium cross section
!   xssteptot         ! preequilibrium cross section per particle type and stage
!
! *** Declaration of local data
!
  implicit none
  integer   :: i      ! counter
  integer   :: J      ! spin of level
  integer   :: nen    ! energy counter
  integer   :: NL     ! last discrete level
  integer   :: p      ! particle number
  integer   :: parity ! parity
  integer   :: pnu    ! neutron particle number
  integer   :: ppi    ! proton particle number
  integer   :: type   ! particle type
  real(sgl) :: Elast  ! help variable
  real(sgl) :: frac   ! help variable
  real(sgl) :: norm   ! normalization factor
!
! ************************ Total pre-equilibrium ***********************
!
! The pre-equilibrium spectra and spectra per exciton number are summed to total pre-equilibrium cross sections.
! Special care is taken for the continuum bin with the highest outgoing energy, i.e. the one that overlaps with
! the energy corresponding to the last discrete state.
!
  do type = 0, 6
    if (parskip(type)) cycle
    NL = Nlast(parZ(type), parN(type), 0)
    Elast = eoutdis(type, NL)
    if (Elast > 0.) then
      frac = Etop(nendisc(type)) - Elast
    else
      frac = 0.
    endif
    do p = p0, maxpar
      do nen = ebegin(type), nendisc(type)
        xssteptot(type, p) = xssteptot(type, p) + xsstep(type, p, nen) * deltaE(nen)
      enddo
      xssteptot(type, p) = xssteptot(type, p) - xsstep(type, p, nendisc(type)) * frac
      xspreeqtot(type) = xspreeqtot(type) + xssteptot(type, p)
    enddo
    do nen = ebegin(type), eend(type)
      xspreeqtotps(type) = xspreeqtotps(type) + xspreeqps(type, nen) * deltaE(nen)
      xspreeqtotki(type) = xspreeqtotki(type) + xspreeqki(type, nen) * deltaE(nen)
      xspreeqtotbu(type) = xspreeqtotbu(type) + xspreeqbu(type, nen) * deltaE(nen)
    enddo
    xspreeqtotps(type) = xspreeqtotps(type) - xspreeqps(type, nendisc(type)) * frac
    xspreeqtotki(type) = xspreeqtotki(type) - xspreeqki(type, nendisc(type)) * frac
    xspreeqtotbu(type) = xspreeqtotbu(type) - xspreeqbu(type, nendisc(type)) * frac
    xspreeqtot(type) = xspreeqtot(type) + xspreeqtotps(type) + xspreeqtotki(type) + xspreeqtotbu(type)
    xspreeqsum = xspreeqsum + xspreeqtot(type)
  enddo
  if (breakupmodel == 2 .and. k0 == 3) xspreeqsum = xspreeqsum - xsEB(1)
!
! Prevent divergence of pre-equilibrium gamma cross sections in case of absence of particle competition.
!
  if (xspreeqsum == xspreeqtot(0)) then
    xspreeqtot(0) = 0.
    xspreeqsum = 0.
    do p = p0, maxpar
      xssteptot(0, p) = 0.
      do nen = ebegin(0), nendisc(0)
        xsstep(0, p, nen) = 0.
        xspreeq(0, nen) = 0.
      enddo
    enddo
  endif
!
! ************************* Unitarity condition ************************
!
! In line with unitarity, the summed direct + pre-equilibrium cross section may not exceed the reaction cross section.
! In these cases, we normalize the results.
!
  xsflux = xsreacinc - xsdirdiscsum - xsgrsum
  xsflux = max(xsflux, 0.)
  preeqnorm = 0.
  if (xsflux > xseps .and. xspreeqsum + xspreeqdiscsum > xsflux) then
    norm = xsflux / (xspreeqsum + xspreeqdiscsum)
    preeqnorm = norm
    xspreeqdiscsum = xspreeqdiscsum * norm
    xspreeqsum = xsflux - xspreeqdiscsum
    do type = 0, 6
      if (parskip(type)) cycle
      xspreeqtot(type) = xspreeqtot(type) * norm
      xspreeqtotps(type) = xspreeqtotps(type) * norm
      xspreeqtotki(type) = xspreeqtotki(type) * norm
      xspreeqtotbu(type) = xspreeqtotbu(type) * norm
      xspreeqdisctot(type) = xspreeqdisctot(type) * norm
      do p = p0, maxpar
        xssteptot(type, p) = xssteptot(type, p) * norm
      enddo
      do i = 0, Nlast(parZ(type), parN(type), 0)
        xspreeqdisc(type, i) = xspreeqdisc(type, i) * norm
      enddo
      do nen = ebegin(type), eend(type)
        xspreeq(type, nen) = xspreeq(type, nen) * norm
        xspreeqps(type, nen) = xspreeqps(type, nen) * norm
        xspreeqki(type, nen) = xspreeqki(type, nen) * norm
        xspreeqbu(type, nen) = xspreeqbu(type, nen) * norm
        do p = p0, maxpar
          xsstep(type, p, nen) = xsstep(type, p, nen) * norm
        enddo
        if (flag2comp) then
          do ppi = ppi0, maxpar
            do pnu = pnu0, maxpar
              xsstep2(type, ppi, pnu, nen) = xsstep2(type, ppi, pnu, nen) * norm
            enddo
          enddo
        endif
        if (pespinmodel >= 3) then
          do parity = - 1, 1, 2
            do J = 0, maxJph
              xspreeqJP(type, nen, J, parity) = xspreeqJP(type, nen, J, parity) * norm
            enddo
          enddo
        endif
      enddo
    enddo
  endif
!
! **** Add discrete pre-equilibrium contribution to discrete states ****
!
  xsdirdiscsum = xsdirdiscsum + xspreeqdiscsum
  do type = 0, 6
    if (parskip(type)) cycle
    xsdirdisctot(type) = xsdirdisctot(type) + xspreeqdisctot(type)
    do i = 0, Nlast(parZ(type), parN(type), 0)
      xsdirdisc(type, i) = xsdirdisc(type, i) + xspreeqdisc(type, i)
    enddo
  enddo
  return
end subroutine preeqtotal
! Copyright A.J. Koning 2021
