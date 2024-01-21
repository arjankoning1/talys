function vr1(z)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Volume of the projectile-like section.
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
!   sgl    ! single precision kind
! Constants
!   pi     ! pi
! Variables for Brosa model
!   r1     ! parameter for neck rupture
!
! *** Declaration of local data
!
  implicit none
  real(sgl) :: vr1  ! function for volume of the projectile-like section
  real(sgl) :: z    ! charge number
!
! **********************************************************************
!
! vr1: function for volume of the projectile-like section
!
  vr1 = pi*((2.*r1**3-z**3)/3.+r1**2*z)
  return
end function vr1
! Copyright A.J. Koning 2021
