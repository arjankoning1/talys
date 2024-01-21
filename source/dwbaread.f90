subroutine dwbaread(nen1, nen2)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read ECIS results for DWBA for MSD
!
! Author    : Arjan Koning and Eric Bauge
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   dbl           ! double precision kind
! Variables for output
!   flagddx       ! flag for output of double - differential cross sections
! Variables for numerics
!   nanglecont    ! number of angles for continuum
! Variables for MSD
!   maxJmsd       ! maximal spin for MSD calculation
!   xsdw          ! DWBA angular distribution as a function of incident energy, outgoing energy, ang. mom. and angle
!   xsdwin        ! DWBA cross section as a function of incident energy, outgoing energy and angular momentum
!
! *** Declaration of local data
!
  implicit none
  integer   :: iang  ! running variable for angle
  integer   :: istat ! logical for file access
  integer   :: itype ! help variable
  integer   :: J     ! spin of level
  integer   :: k     ! designator for particle
  integer   :: nen1  ! energy counter
  integer   :: nen2  ! energy counter
  integer   :: nS    ! number of states
  real(dbl) :: xs    ! help variable
!
! ********************** Read DWBA cross sections **********************
!
  read(10, '()')
  do J = 0, maxJmsd
    read(10, * ) xs
    xsdwin(nen1, nen2, J, 0) = real(xs)
  enddo
  if (flagddx) then
    read(8, '()')
    read(8, '(12x, i3)') nS
    do iang = 0, nanglecont
      do k = 1, nS
        read(8, '()')
    enddo
  enddo
    do J = 0, maxJmsd
      read(8, '(12x, i3)', iostat = istat) nS
      if (istat /= 0) cycle
      do iang = 0, nanglecont
        do k = 1, nS
          read(8, '(i3, 12x, e12.5)', iostat = istat) itype, xs
          if (istat /= 0) cycle
          if (itype == 0) xsdw(nen1, nen2, J, iang, 0) = real(xs)
        enddo
      enddo
    enddo
  endif
  return
end subroutine dwbaread
! Copyright A.J. Koning 2021
