function vr2(zz1, zz2, zz3)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Volume of the neck
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
!   aaa    ! parameter for neck rupture
!   cur    ! parameter for neck rupture
!   r2     ! parameter for neck rupture
!
! *** Declaration of local data
!
  implicit none
  real(sgl) :: b    ! parameter for slope
  real(sgl) :: vr2  ! function for volume of the neck
  real(sgl) :: z21  ! asymmetry value
  real(sgl) :: z32  ! asymmetry value
  real(sgl) :: zz1  ! Z value
  real(sgl) :: zz2  ! Z value
  real(sgl) :: zz3  ! Z value
!
! **********************************************************************
!
! vr2: function for volume of the neck
!
  z21 = (zz2  -zz1) / aaa
  z32 = (zz3 - zz2) / aaa
  if (z21 > 30.) z21 = 30.
  if (z32 > 30.) z32 = 30.
  b = aaa * aaa * cur
  vr2 = pi * (((r2 - b) **2 + .5 * b **2) * (zz3 - zz1) + aaa * b * (2. * (r2 - b) * (sinh(z32) + sinh(z21)) + &
    0.25 * b * (sinh(2. * z32) + sinh(2. * z21))))
  return
end function vr2
! Copyright A.J. Koning 2021
