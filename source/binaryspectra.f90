subroutine binaryspectra
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Creation of binary spectra
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
!   flagddx        ! flag for output of double - differential cross sections
! Variables for basic reaction
!   flagrecoil     ! flag for calculation of recoils
! Variables for numerics
!   nanglecont     ! number of angles for continuum
! Variables for energy grid
!   ebegin         ! first energy point of energy grid
! Variables for energies
!   eend           ! last energy point of energy grid
! Variables for binary emission spectra
!   xsbinemis      ! cross section for emission from first compound nucleus
!   xsbinemisad    ! angular distribution for emission from first compound nucleus
!   xscomp         ! compound elastic cross section
!   xscompad       ! compound emission angular distribution
!   xsemis         ! cross section for emission from compound nucleus
! Variables for nuclides
!   parskip        ! logical to skip outgoing particle
! Variables for giant resonances
!   xsgrad         ! smoothed giant resonance angular distribution
! Constants
!   fourpi         ! 4. * pi
! Variables for incident channel
!   xsgr           ! total smoothed giant resonance cross section
!   xspreeq        ! preeq. cross section per particle typ and outgoing energye
! Variables for preequilibrium
!   xspreeqad      ! preequilibrium angular distribution per particle type and outg
!
! *** Declaration of local data
!
  implicit none
  integer :: iang              ! running variable for angle
  integer :: nen               ! energy counter
  integer :: type              ! particle type
!
! *************** Interpolate decay on emission spectrum ***************
!
! binemission: subroutine for compound emission cross sections for binary reaction
!
  call binemission
  do type = 0, 6
    if (parskip(type)) cycle
    do nen = ebegin(type), eend(type)
      xscomp(type, nen) = xsemis(type, nen)
      xsbinemis(type, nen) = xscomp(type, nen) + xspreeq(type, nen) + xsgr(type, nen)
      if (flagddx .or. flagrecoil) then
        do iang = 0, nanglecont
          xscompad(type, nen, iang) = xscomp(type, nen) / fourpi
          xsbinemisad(type, nen, iang) = xscompad(type, nen, iang) + &
            xspreeqad(type, nen, iang) + xsgrad(type, nen, iang)
        enddo
      endif
    enddo
  enddo
  return
end subroutine binaryspectra
! Copyright A.J. Koning 2021
