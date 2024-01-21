function radius(a)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Radius function for Tripathi formula
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               sgl      ! single precision kind
!
! *** Declaration of local data
!
  implicit none
  integer   :: i       ! counter
  integer   :: ia      ! mass number from abundance table
  integer   :: na(23)  ! help variable
  real(sgl) :: a       ! fit variables
  real(sgl) :: fact    ! reaction rate factor for particles
  real(sgl) :: radius  ! radius function
  real(sgl) :: rms(23) ! root-mean-square radius of the potential
!
! *************************** Radius function **************************
!
  na = (/1, 2, 3, 4, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26 /)
  rms = (/ 0.85, 2.095, 1.976, 1.671, 2.57, 2.41, 2.519, 2.45, 2.42, &
    2.471, 2.440, 2.58, 2.611, 2.730, 2.662, 2.727, 2.900, 3.040, 2.969, 2.94, 3.075, 3.11, 3.06 /)
  fact = sqrt(5. / 3.)
  ia = int(a + 0.4)
  radius = fact * (0.84 * a **(1. / 3.) + 0.55)
  do i = 1, 23
    if (ia == na(i)) radius = fact * rms(i)
  enddo
  return
end function radius
! Copyright A.J. Koning 2021
