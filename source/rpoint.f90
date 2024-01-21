function rpoint(z)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : rpoint=0 gives for a given nucleon number hidden in the
!
! Author    : Marieke Duijvestijn
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl     ! single precision kind
! Variables for Brosa model
!   di      ! nucleon number density
!   rest    ! parameter for neck rupture
!   z1      ! parameter for neck rupture
!   z2      ! parameter for neck rupture
!
! *** Declaration of local data
!
  implicit none
  real(sgl) :: rpoint  ! function for rupture point
  real(sgl) :: vr2     ! function for volume of the neck
  real(sgl) :: z       ! charge number
!
! **********************************************************************
!
! rpoint: function for rupture point
!
  rpoint = di * vr2(z1, z2, z) + rest
  return
end function rpoint
! Copyright A.J. Koning 2021
