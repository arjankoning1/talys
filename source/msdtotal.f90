subroutine msdtotal
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Total multi-step direct cross sections
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
!   flagddx         ! flag for output of double - differential cross sections
! Variables for numerics
!   nanglecont      ! number of angles for continuum
! Variables for energy grid
!   deltaE          ! energy bin around outgoing energies
!   ebegin          ! first energy point of energy grid
! Variables for energies
!   eend            ! last energy point of energy grid
! Variables for MSD
!   maxmsd          ! number of MSD steps
!   msdall          ! total multi - step direct cross section
!   msdstep         ! continuum n - step direct cross section
!   msdstepad       ! continuum n - step direct angular distribution
!   msdstepint      ! n - step direct cross section integrated over energy
!   msdstepintad    ! n - step direct angular distribution integrated over energy
!   msdsum          ! multi - step direct cross section summed over steps and integrated over energy
!   msdtot          ! multi - step direct cross section summed over steps
!   msdtotad        ! multi - step direct angular distribution summed over steps
!   msdtotintad     ! multi - step direct angular distribution summed over steps and integrated over energy
!
! *** Declaration of local data
!
  implicit none
  integer :: ns   ! counter
  integer :: iang ! running variable for angle
  integer :: nen  ! energy counter
  integer :: type ! particle type
!
! ************** Angle-integrated multi-step cross sections ************
!
  msdall = 0.
  do type = 1, 2
    msdsum(type) = 0.
    do nen = ebegin(type), eend(type)
      msdtot(type, nen) = 0.
    enddo
    do ns = 1, maxmsd
      msdstepint(type, ns) = 0.
      do nen = ebegin(type), eend(type)
        msdtot(type, nen) = msdtot(type, nen) + msdstep(type, ns, nen)
        msdstepint(type, ns) = msdstepint(type, ns) + msdstep(type, ns, nen) * deltaE(nen)
      enddo
      msdsum(type) = msdsum(type) + msdstepint(type, ns)
    enddo
    msdall = msdall + msdsum(type)
  enddo
  if ( .not. flagddx) return
!
! ******************* Multi-step angular distributions *****************
!
! Total multi-step angular distributions
!
  do type = 1, 2
    do iang = 0, nanglecont
      msdtotintad(type, iang) = 0.
    enddo
    do nen = ebegin(type), eend(type)
      do iang = 0, nanglecont
        msdtotad(type, nen, iang) = 0.
    enddo
  enddo
    do ns = 1, maxmsd
      do iang = 0, nanglecont
        msdstepintad(type, ns, iang) = 0.
      enddo
      do iang = 0, nanglecont
        do nen = ebegin(type), eend(type)
          msdtotad(type, nen, iang) = msdtotad(type, nen, iang) + msdstepad(type, ns, nen, iang)
          msdstepintad(type, ns, iang) = msdstepintad(type, ns, iang) + msdstepad(type, ns, nen, iang) * deltaE(nen)
        enddo
        msdtotintad(type, iang) = msdtotintad(type, iang) + msdstepintad(type, ns, iang)
      enddo
    enddo
  enddo
  return
end subroutine msdtotal
! Copyright A.J. Koning 2021
