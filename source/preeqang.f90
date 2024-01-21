subroutine preeqang
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Pre-equilibrium angular distribution
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
! Variables for numerics
!   nangle         ! number of angles
!   nanglecont     ! number of angles for continuum
! Variables for preequilibrium
!   preeqmode      ! designator for pre - equilibrium model
! Variables for energy grid
!   angle          ! angle in degrees
!   anglecont      ! angle in degrees for continuum
!   ebegin         ! first energy point of energy grid
!   egrid          ! outgoing energy grid
!   Einc           ! incident energy in MeV
! Variables for energies
!   eend           ! last energy point of energy grid
!   eoutdis        ! outgoing energy of discrete state reaction
! Variables for incident channel
!   directad       ! direct angular distribution
!   xspreeq        ! preeq. cross section per particle type and outgoing energy
!   xspreeqtot     ! preequilibrium cross section per particle type
! Variables for nuclides
!   parskip        ! logical to skip outgoing particle
! Constants
!   deg2rad        ! conversion factor for degrees to radians
!   parN           ! neutron number of particle
!   parZ           ! charge number of particle
! Variables for level density
!   Nlast          ! last discrete level
! Variables for preequilibrium
!   xspreeqad      ! preequilibrium angular distribution per particle type and outgoing energy
!   xspreeqbu      ! preequilibrium cross section per particle type and outgoing energy for break-up
!   xspreeqdisc    ! preequilibrium cross section for discrete state
!
! *** Declaration of local data
!
  implicit none
  integer   :: i         ! counter
  integer   :: iang      ! running variable for angle
  integer   :: nen       ! energy counter
  integer   :: NL        ! last discrete level
  integer   :: type      ! particle type
  real(sgl) :: ang       ! angle
  real(sgl) :: Eout      ! outgoing energy
  real(sgl) :: kalbach   ! Kalbach function
  real(sgl) :: kalbachBU ! Kalbach function for break-up
  real(sgl) :: xs        ! help variable
  real(sgl) :: xsbu      ! help variable
  real(sgl) :: xspe      ! help variable
!
! ************** Kalbach angular distribution for exciton model ********
!
! kalbach   : Kalbach function
! kalbachBU : Kalbach function for break-up
!
  do type = 0, 6
    if (parskip(type)) cycle
    if (preeqmode == 4 .and. (type == 1 .or. type == 2)) cycle
    if (xspreeqtot(type) == 0.) cycle
    do nen = ebegin(type), eend(type)
      if (xspreeq(type, nen) == 0.) cycle
      Eout = egrid(nen)
      xspe = xspreeq(type, nen) - xspreeqbu(type, nen)
      xsbu = xspreeqbu(type, nen)
      do iang = 0, nanglecont
        ang = anglecont(iang) * deg2rad
        xspreeqad(type, nen, iang) = xspe * kalbach(type, Einc, Eout, ang)
        if (xsbu > 0.) xspreeqad(type, nen, iang) = xspreeqad(type, nen, iang) + xsbu * kalbachBU(type, Einc, ang)
      enddo
    enddo
  enddo
!
! ************ Pre-equilibrium cross sections for direct states ********
!
! Correction of pre-equilibrium cross sections for direct discrete cross sections.
! If the cross sections for discrete states have NOT been calculated by a direct reaction model, we collapse the
! continuum pre-equilibrium cross sections in the high energy region on the associated discrete states.
!
  do type = 0, 6
    if (parskip(type)) cycle
    NL = Nlast(parZ(type), parN(type), 0)
    do i = 0, NL
      xs = xspreeqdisc(type, i)
      if (xs == 0.) cycle
      Eout = eoutdis(type, NL)
      do iang = 0, nangle
        ang = angle(iang) * deg2rad
        directad(type, i, iang) = directad(type, i, iang) + xs * kalbach(type, Einc, Eout, ang)
      enddo
    enddo
  enddo
  return
end subroutine preeqang
! Copyright A.J. Koning 2021
