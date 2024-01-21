subroutine dwbaint
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Interpolate DWBA cross sections for MSD
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
! Variables for MSD
!   maxJmsd       ! maximal spin for MSD calculation
!   msdbins2      ! number of energy points for MSD calculation
!   xsdw          ! DWBA angular distribution as a function of incident energy, outgoing energy, ang. mom. and angle
!   xsdwin        ! DWBA cross section as a function of incident energy, outgoing energy and angular momentum
!
! *** Declaration of local data
!
  implicit none
  integer :: iang ! running variable for angle
  integer :: J    ! spin of level
  integer :: nen1 ! energy counter
  integer :: nen2 ! energy counter
!
! ******************* Interpolate DWBA cross sections ******************
!
  do nen1 = 0, msdbins2, 2
    do J = 0, maxJmsd
      xsdwin(nen1, nen1, J, 0) = xsdwin(nen1, nen1 + 2, J, 0)
      xsdwin(nen1 + 1, nen1 + 1, J, 0) = xsdwin(nen1, nen1 + 2, J, 0)
      if (flagddx) then
        do iang = 0, nanglecont
          xsdw(nen1 + 1, nen1 + 1, J, iang, 0) = xsdw(nen1, nen1 + 2, J, iang, 0)
        enddo
      endif
    enddo
    do nen2 = nen1 + 1, msdbins2 - 1, 2
      do J = 0, maxJmsd
        xsdwin(nen1, nen2, J, 0) = 0.5 * (xsdwin(nen1, nen2 - 1, J, 0) + xsdwin(nen1, nen2 + 1, J, 0))
        if (flagddx) then
          do iang = 0, nanglecont
            xsdw(nen1, nen2, J, iang, 0) = 0.5 * (xsdw(nen1, nen2 - 1, J, iang, 0) &
              + xsdw(nen1, nen2 + 1, J, iang, 0))
          enddo
        endif
    enddo
  enddo
  enddo
  do nen1 = 1, msdbins2 - 1, 2
    do nen2 = nen1 + 1, msdbins2, 2
      do J = 0, maxJmsd
        xsdwin(nen1, nen2, J, 0) = 0.5 * (xsdwin(nen1 - 1, nen2, J, 0) + xsdwin(nen1 + 1, nen2, J, 0))
        if (flagddx) then
          do iang = 0, nanglecont
            xsdw(nen1, nen2, J, iang, 0) = 0.5 * (xsdw(nen1 - 1, nen2, J, iang, 0) &
              + xsdw(nen1 + 1, nen2, J, iang, 0))
          enddo
        endif
      enddo
    enddo
  enddo
  do nen1 = 1, msdbins2 - 1, 2
    do nen2 = nen1 + 2, msdbins2 - 1, 2
      do J = 0, maxJmsd
        xsdwin(nen1, nen2, J, 0) = 0.5 * (xsdwin(nen1 - 1, nen2 - 1, J, 0) + xsdwin(nen1 + 1, nen2 + 1, J, 0))
        if (flagddx) then
          do iang = 0, nanglecont
            xsdw(nen1, nen2, J, iang, 0) = 0.5 * (xsdw(nen1 - 1, nen2 - 1, J, iang, 0) + &
              xsdw(nen1 + 1, nen2 + 1, J, iang, 0))
          enddo
        endif
      enddo
    enddo
  enddo
  return
end subroutine dwbaint
! Copyright A.J. Koning 2021
