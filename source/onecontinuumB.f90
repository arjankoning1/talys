subroutine onecontinuumB
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : One-step direct cross sections for MSD
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
! Variables for main input
!   k0            ! index of incident particle
! Variables for nuclides
!   Nindex        ! neutron number index for residual nucleus
!   parskip       ! logical to skip outgoing particle
!   Zindex        ! charge number index for residual nucleus
! Variables for masses
!   specmass      ! specific mass for residual nucleus
! Variables for MSD
!   Emsd          ! minimal outgoing energy for MSD calculation
!   Emsdin        ! incident MSD energy
!   msdbins2      ! number of energy points for MSD calculation
!   msdstep0      ! n - step cross section for MSD
!   msdstepad0    ! n - step angular distribution for MSD
!   xscont        ! continuum one - step direct cross section
!   xscont1       ! continuum one - step direct cross section (unnormalized)
!   xscontad      ! continuum one - step direct angular distribution for MSD
!   xscontad1     ! continuum one - step direct angular distribution for MSD (unnormalized)
!
! *** Declaration of local data
!
  implicit none
  integer   :: iang  ! running variable for angle
  integer   :: itype ! help variable
  integer   :: nen1  ! energy counter
  integer   :: nen2  ! energy counter
  integer   :: Nix   ! neutron number index for residual nucleus
  integer   :: type  ! particle type
  integer   :: Zix   ! charge number index for residual nucleus
  real(sgl) :: cmsd  ! normalization factor for MSD
!
! ************* Calculate continuum one-step cross sections ************
!
!   (unnormalized)
!   (unnormalized)
!
  do itype = 1, 2
    if (parskip(itype)) cycle
    do type = 1, 2
      if (parskip(type)) cycle
      Zix = Zindex(0, 0, type)
      Nix = Nindex(0, 0, type)
      do nen1 = 0, msdbins2
        Emsdin = real(specmass(Zix, Nix, itype) * Emsd(nen1))
        do nen2 = nen1, msdbins2
          xscont(itype, type, nen1, nen2) = cmsd(Emsdin) * xscont1(itype, type, nen1, nen2)
          if (flagddx) then
            do iang = 0, nanglecont
              xscontad(itype, type, nen1, nen2, iang) = cmsd(Emsdin) * xscontad1(itype, type, nen1, nen2, iang)
            enddo
          endif
        enddo
      enddo
    enddo
  enddo
!
! ******************* First-step cross sections for MSD ****************
!
  do type = 1, 2
    if (parskip(type)) cycle
    do nen2 = 0, msdbins2
      msdstep0(type, 1, nen2) = xscont(k0, type, 0, nen2)
      if (flagddx) then
        do iang = 0, nanglecont
          msdstepad0(type, 1, nen2, iang) = xscontad(k0, type, 0, nen2, iang)
        enddo
      endif
    enddo
  enddo
  return
end subroutine onecontinuumB
! Copyright A.J. Koning 2021
