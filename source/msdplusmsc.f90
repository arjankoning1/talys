subroutine msdplusmsc
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Total quantum-mechanical pre-equilibrium cross sections
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
!   flagddx       ! flag for output of double - differential cross sections
! Variables for numerics
!   nanglecont    ! number of angles for continuum
! Variables for energy grid
!   ebegin        ! first energy point of energy grid
! Variables for energies
!   eend          ! last energy point of energy grid
! Variables for nuclides
!   parskip       ! logical to skip outgoing particle
! Variables for preequilibrium initialization
!   maxJph        ! maximal spin for particle - hole states
!   RnJ           ! spin distribution for particle - hole states
!   RnJsum        ! (2J + 1) * sum over spin distributions
! Variables for incident channel
!   xspreeq       ! preeq. cross section per particle typ and outgoing energy
! Variables for preequilibrium
!   xspreeqad     ! preequilibrium angular distribution per particle type and outgoing energy
!   xspreeqJP     ! preeq. cross section per particle type, outgoing energy, J, P
!   xsstep        ! preeq. cross section per particle type, stage and outgoing E
! Variables for MSD
!   maxmsd        ! number of MSD steps
!   msdstep       ! continuum n - step direct cross section
!   msdtot        ! multi - step direct cross section summed over steps
!   msdtotad      ! multi - step direct angular distribution summed over steps
!
! *** Declaration of local data
!
  implicit none
  integer :: ns     ! counter
  integer :: iang   ! running variable for angle
  integer :: J      ! spin of level
  integer :: nen    ! energy counter
  integer :: parity ! parity
  integer :: type   ! particle type
!
! *********** Angle-integrated pre-equilibrium cross sections **********
!
  do type = 1, 2
    if (parskip(type)) cycle
    do nen = ebegin(type), eend(type)
      xspreeq(type, nen) = msdtot(type, nen)
    enddo
    do ns = 1, maxmsd
      do nen = ebegin(type), eend(type)
        xsstep(type, ns, nen) = msdstep(type, ns, nen)
!
! Spin distribution the same as in exciton model. This will be changed in the true MSD model.
!
! Create J-dependent pre-equilibrium cross sections using the spin distribution.
! The result is normalized with the sum of RnJ over J.
!
        do parity = - 1, 1, 2
          do J = 0, maxJph
            xspreeqJP(type, nen, J, parity) = xspreeqJP(type, nen, J, parity) + &
              xsstep(type, ns, nen) * 0.5 * (2. * J + 1.) * RnJ(ns * 2 + 1, J) / RnJsum(ns * 2 + 1)
            enddo
          enddo
      enddo
    enddo
  enddo
!
! **************** Pre-equilibrium angular distributions ***************
!
  if ( .not. flagddx) return
  do type = 1, 2
    if (parskip(type)) cycle
    do iang = 0, nanglecont
      do nen = ebegin(type), eend(type)
        xspreeqad(type, nen, iang) = msdtotad(type, nen, iang)
      enddo
    enddo
  enddo
  return
end subroutine msdplusmsc
! Copyright A.J. Koning 2021
