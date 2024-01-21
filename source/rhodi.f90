function rhodi(z)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Returns the shape of the dinuclear system.
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
!   aaa     ! parameter for neck rupture
!   cur     ! parameter for neck rupture
!   r1      ! parameter for neck rupture
!   r2      ! parameter for neck rupture
!   r3      ! parameter for neck rupture
!   totl    ! parameter for neck rupture
!   z1      ! parameter for neck rupture
!   z2      ! parameter for neck rupture
!   z3      ! parameter for neck rupture
!
! *** Declaration of local data
!
  implicit none
  real(sgl) :: d     ! parameter for energy smoothing
  real(sgl) :: rhodi ! function that returns the shape of the dinuclear system
  real(sgl) :: z     ! charge number
!
! **********************************************************************
!
! rhodi: function that returns the shape of the dinuclear system
!
  d = totl-r1-r3
  if (z <  - r1 .or. z > d+r3) goto 3
  if (z > z1) goto 1
  rhodi = sqrt(r1 **2 - z **2)
  return
 1    if (z > z3) goto 2
  rhodi = r2 + aaa * aaa * cur * (cosh((z - z2) / aaa) - 1.)
  return
 2    rhodi = sqrt(r3 **2 - (z - d) **2)
  return
 3    rhodi = 0.
  return
end function rhodi
! Copyright A.J. Koning 2021
