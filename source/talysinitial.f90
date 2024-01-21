subroutine talysinitial
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialization of nuclear structure
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
!   flagmain        ! flag for main output
! Variables for basic reaction
!   flagreaction    ! flag for calculation of nuclear reactions
!
! *** Declaration of local data
!
  implicit none
!
! ********** Initialization of constants and nuclear structure *********
!
! particles   : subroutine to determine included light particles
! nuclides    : subroutine for properties of nuclides
! grid        : subroutine for energy and angle grid
! mainout     : subroutine for main output
!
  call particles
  call nuclides
  if (flagreaction) call grid
  if (flagmain) call mainout
  return
end subroutine talysinitial
! Copyright A.J. Koning 2021
