subroutine onestepB
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Continuum one-step direct cross sections
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
!   sgl           ! single precision kind
! Variables for output
!   flagddx       ! flag for output of double - differential cross sections
! Variables for numerics
!   nanglecont    ! number of angles for continuum
! Variables for energy grid
!   ebegin        ! first energy point of energy grid
! Variables for energies
!   eend          ! last energy point of energy grid
!   eninccm       ! center - of - mass incident energy in MeV
! Variables for nuclides
!   parskip       ! logical to skip outgoing particle
! Variables for MSD
!   msdstep       ! continuum n - step direct cross section
!   msdstep1      ! continuum one - step direct cross section (unnormalized)
!   msdstepad     ! continuum n - step direct angular distribution
!   msdstepad1    ! continuum one - step direct angular distribution (unnormalized)
!
! *** Declaration of local data
!
  implicit none
  integer   :: iang ! running variable for angle
  integer   :: nen  ! energy counter
  integer   :: type ! particle type
  real(sgl) :: cmsd ! normalization factor for MSD
!
! ************* Calculate continuum one-step cross sections ************
!
  do type = 1, 2
    if (parskip(type)) cycle
    do nen = ebegin(type), eend(type)
      msdstep(type, 1, nen) = cmsd(eninccm) * msdstep1(type, nen)
      if (flagddx) then
        do iang = 0, nanglecont
          msdstepad(type, 1, nen, iang) = cmsd(eninccm) * msdstepad1(type, nen, iang)
        enddo
      endif
    enddo
  enddo
  return
end subroutine onestepB
! Copyright A.J. Koning 2021
