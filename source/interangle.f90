subroutine interangle
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Intermediate angles by addition theorem for MSD model
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
! Variables for numerics
!   nanglecont    ! number of angles for continuum
! Variables for energy grid
!   anglecont     ! angle in degrees for continuum
! Constants
!   deg2rad       ! conversion factor for degrees to radians
!   rad2deg       ! conversion factor for radians to degrees
!   twopi         ! 2 * pi
! Variables for MSD
!   nangleint     ! number of possibilities to link intermediate angle to
!
! *** Declaration of local data
!
  implicit none
  integer   :: iang    ! running variable for angle
  integer   :: iangint ! intermediate angle index
  integer   :: iangout ! outgoing angle index
  integer   :: iphi    ! counter for angle
  real(sgl) :: angint  ! intermediate angle
  real(sgl) :: angout  ! outgoing angle
  real(sgl) :: angstep ! angle step in degrees
  real(sgl) :: cosgam  ! cosine of gamma angle
  real(sgl) :: cosine  ! help variable
  real(sgl) :: gam     ! Brosa parameter
  real(sgl) :: phi     ! help variable
  real(sgl) :: sine    ! help variable
!
! ************************ Addition theorem ****************************
!
! This is necessary to transform the ingoing angle of the second step and link it with the outgoing angle of the first step.
!
  nangleint = 0
  angstep = 180./nanglecont
  do iangout = 0, nanglecont
    angout = anglecont(iangout) * deg2rad
    do iangint = 0, nanglecont
      angint = anglecont(iangint) * deg2rad
      cosine = cos(angout) * cos(angint)
      sine = sin(angout) * sin(angint)
      do iphi = 0, 2 * nanglecont - 1
        phi = real(iphi) * 0.5 / nanglecont * twopi
        cosgam = cosine + cos(phi) * sine
        if (abs(cosgam) > 1.0) cosgam = 1.0
        gam = acos(cosgam) * rad2deg
        iang = int((gam + angstep / 2.) / angstep)
        nangleint(iangout, iangint, iang) = nangleint(iangout, iangint, iang) + 1
      enddo
    enddo
  enddo
  return
end subroutine interangle
! Copyright A.J. Koning 2021
