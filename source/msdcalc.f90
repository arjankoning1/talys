subroutine msdcalc
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : MSD calculation
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
! Variables for preequilibrium
!   flagonestep    ! flag for continuum one - step direct only
!
! *** Declaration of local data
!
  implicit none
!
! ****************************** MSD model *****************************
!
! onestepB     : subroutine for continuum one-step direct cross sections
! onecontinuumB: subroutine for one-step direct cross sections for MSD
! multistepA   : subroutine for multi-step direct cross sections
! multistepB   : subroutine for multi-step direct cross sections on
!   outgoing energy grid
!
  call onestepB
  if ( .not. flagonestep) then
    call onecontinuumB
    call multistepA
    call multistepB
  endif
  return
end subroutine msdcalc
! Copyright A.J. Koning 2021
