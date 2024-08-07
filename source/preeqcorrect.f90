subroutine preeqcorrect
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Correct pre-equilibrium cross sections for direct effects
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
! Variables for numerics
!   xseps             ! limit for cross sections
! Variables for preequilibrium
!   flag2comp         ! flag for two - component pre - equilibrium model
! Variables for main
!   k0                ! index of incident particle
!   Ltarget           ! excited level of target
! Variables for energy grid
!   deltaE            ! energy bin around outgoing energies
!   ebegin            ! first energy point of energy grid
!   egrid             ! outgoing energy grid
!   Etop              ! top of outgoing energy bin
! Variables for energies
!   eend              ! last energy point of energy grid
!   eoutdis           ! outgoing energy of discrete state reaction
!   nendisc           ! last discrete bin
! Variables for incident channel
!   dorigin           ! origin of direct cross section (Direct or Preeq)
!   xsdirdisc         ! direct cross section for discrete state direct cross section
!   xspreeq           ! preeq. cross section per particle typ and outgoing energye
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
!   pnu0              ! initial neutron number
!   ppi0              ! initial proton number
!   xspreeqbu         ! preequilibrium cross section per particle type and outgoing energy for brea
!   xspreeqdisc       ! preequilibrium cross section for discrete state
!   xspreeqdiscsum    ! total preequilibrium cross section for discrete states
!   xspreeqdisctot    ! preequilibrium cross section summed over discrete states
!   xspreeqJP         ! preeq. cross section per particle type, outgoing energy, J, P
!   xspreeqki         ! preequilibrium cross section per particle type and outgoing energy for knoc
!   xspreeqps         ! preequilibrium cross section per particle type and outgoing energy for pick
!   xsstep            ! preeq. cross section per particle type, stage and outgoing E
!   xsstep2           ! two - component preequilibrium cross section
!
! *** Declaration of local data
!
  implicit none
  integer   :: i         ! counter
  integer   :: J         ! spin of level
  integer   :: nen       ! energy counter
  integer   :: nen1      ! energy counter
  integer   :: nen2      ! energy counter
  integer   :: NL        ! last discrete level
  integer   :: p         ! particle number
  integer   :: parity    ! parity
  integer   :: pnu       ! neutron particle number
  integer   :: ppi       ! proton particle number
  integer   :: type      ! particle type
  real(sgl) :: Elast     ! help variable
  real(sgl) :: esd       ! outgoing energy of discrete state reaction
  real(sgl) :: esd1      ! help variable
  real(sgl) :: esd2      ! help variable
  real(sgl) :: Rboundary ! factor taking into count first accessible mother bin for discrete state
  real(sgl) :: xs1       ! help variable
!
! ******************** Pre-equilibrium to direct ***********************
!
! 1. Correction of pre-equilibrium cross sections for direct discrete cross sections.
!    If the cross sections for discrete states have NOT been calculated by a direct reaction model, we collapse the
!    continuum pre-equilibrium cross sections in the high energy region on the associated discrete states.
!
  do type = 0, 6
    if (parskip(type)) cycle
    if (ebegin(type) >= eend(type)) cycle
    if (Ltarget /= 0) goto 100
    NL = Nlast(parZ(type), parN(type), 0)
    do i = 0, NL
      if (xsdirdisc(type, i) /= 0.) cycle
      esd = eoutdis(type, i)
      if (esd < 0.) goto 100
      if (i == 0 .or. (type == k0 .and. i == 1)) then
        esd2 = eoutdis(type, 0)
      else
        esd2 = 0.5 * (esd + eoutdis(type, i - 1))
      endif
      if (i == NL) then
        esd1 = eoutdis(type, NL)
      else
        if (eoutdis(type, i + 1) > 0.) then
          esd1 = 0.5 * (esd + eoutdis(type, i + 1))
        else
          esd1 = 0.
        endif
      endif
!
! Find the part of the continuum spectrum that corresponds with the discrete states.
!
! locate        : subroutine to find value in ordered table
!
      if (type == k0) then
        xs1 = xseps
      else
        call locate(egrid, nendisc(type), eend(type), esd1, nen1)
        call locate(egrid, nendisc(type), eend(type), esd2, nen2)
        xs1 = 0.5 * (xspreeq(type, nen1) + xspreeq(type, nen2)) * (esd2 - esd1)
      endif
!
! 2. Set pre-equilibrium cross section in discrete energy region to zero.
!
! The first continuum outgoing energy bin is only partially depleted by the last discrete level.
! This is corrected using Rboundary.
!
      if (i == NL) then
        nen = nendisc(type)
        Elast = eoutdis(type, NL)
        Rboundary = (Etop(nen) - Elast) / deltaE(nen)
        if (abs(Rboundary) > 1.) Rboundary = 0.
        if (type /= k0) xs1 = xspreeq(type, nen) * Rboundary
        xspreeq(type, nen) = xspreeq(type, nen) * (1. - Rboundary)
        xspreeqps(type, nen) = xspreeqps(type, nen) * (1. - Rboundary)
        xspreeqki(type, nen) = xspreeqki(type, nen) * (1. - Rboundary)
        xspreeqbu(type, nen) = xspreeqbu(type, nen) * (1. - Rboundary)
        do p = 1, maxpar
          xsstep(type, p, nen) = xsstep(type, p, nen) * (1. - Rboundary)
        enddo
        if (flag2comp) then
          do ppi = ppi0, maxpar
            do pnu = pnu0, maxpar
              xsstep2(type, ppi, pnu, nen) = xsstep2(type, ppi, pnu, nen) * (1. - Rboundary)
            enddo
          enddo
        endif
        do parity = - 1, 1, 2
          do J = 0, maxJph
            xspreeqJP(type, nen, J, parity) = xspreeqJP(type, nen, J, parity) * (1. - Rboundary)
          enddo
        enddo
      endif
      xspreeqdisc(type, i) = xs1
      xspreeqdisctot(type) = xspreeqdisctot(type) + xs1
      xspreeqdiscsum = xspreeqdiscsum + xs1
      dorigin(type, i) = 'Preeq '
    enddo
!
! The pre-equilibrium spectrum for energies corresponding to discrete transitions is set to zero.
!
100 do nen = nendisc(type) + 1, eend(type)
      xspreeq(type, nen) = 0.
      xspreeqps(type, nen) = 0.
      xspreeqki(type, nen) = 0.
      xspreeqbu(type, nen) = 0.
      do p = 1, maxpar
        xsstep(type, p, nen) = 0.
      enddo
      if (flag2comp) then
        do ppi = ppi0, maxpar
          do pnu = pnu0, maxpar
            xsstep2(type, ppi, pnu, nen) = 0.
          enddo
        enddo
      endif
      do parity = - 1, 1, 2
        do J = 0, maxJph
          xspreeqJP(type, nen, J, parity) = 0.
        enddo
      enddo
    enddo
  enddo
  return
end subroutine preeqcorrect
! Copyright A.J. Koning 2021
