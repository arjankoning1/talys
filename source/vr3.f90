function vr3(z)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Volume of the target-like section.
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
! Constants
!   pi      ! pi
! Variables for Brosa model
!   r1      ! parameter for neck rupture
!   r3      ! parameter for neck rupture
!   totl    ! parameter for neck rupture
!
! *** Declaration of local data
!
  implicit none
  real(sgl) :: d              ! parameter for energy smoothing
  real(sgl) :: vr3            ! help value
  real(sgl) :: z              ! charge number
!
! **********************************************************************
!
  d = totl - r1 - r3
  vr3 = pi * ((2. * r3 **3 - (d - z) **3) / 3. + r3 **2 * (d - z))
  return
end function vr3
! Copyright A.J. Koning 2021
