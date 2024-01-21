subroutine initial_adjust
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialize parameters for adjustment
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
! Variables for adjustment
!   adjustfile        ! file for local adjustment
!   adjustix          ! local adjustment index
!   adjustkey         ! keyword for local adjustment
!   adjustpar         ! local adjustment parameters
!   Dadjust           ! tabulated depth of local adjustment
!   Eadjust           ! tabulated energy of local adjustment
!   Nadjust           ! number of adjustable parameters
!   nenadjust         ! number of tabulated energies of local adjustment
!
! *** Declaration of local data
!
  implicit none
!
! ************** Defaults for adjustment variables *************
!
  Nadjust = 0
  adjustkey = ' '
  adjustfile = ' '
  nenadjust = 0
  adjustix = 0
  adjustpar = 0.
  Eadjust = 0.
  Dadjust = 0.
  return
end subroutine initial_adjust
! Copyright A.J. Koning 2021
