function fidi(x, i)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : fidi=0 determines the exact shape of the dinuclear complex.
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
!   c0      ! curvature of neck
!   cur     ! parameter for neck rupture
!   r1      ! parameter for neck rupture
!   r2      ! parameter for neck rupture
!   r3      ! parameter for neck rupture
!   rp      ! parameter for neck rupture
!   rpt     ! parameter for neck rupture
!   rt      ! parameter for neck rupture
!   totl    ! parameter for neck rupture
!   vtot    ! total volume of the complex
!   z1      ! parameter for neck rupture
!   z2      ! parameter for neck rupture
!   z3      ! parameter for neck rupture
!
! *** Declaration of local data
!
  implicit none
  integer   :: i     ! counter
  real(sgl) :: al    ! help variable
  real(sgl) :: be    ! help variable
  real(sgl) :: d     ! parameter for energy smoothing
  real(sgl) :: de    ! single-nucleon exchange term J00
  real(sgl) :: ds    ! help variable
  real(sgl) :: epsc  ! help variable
  real(sgl) :: eq    ! help variable
  real(sgl) :: fidi  ! function for dinuclear shape
  real(sgl) :: ga    ! help variable
  real(sgl) :: sq1   ! help variable
  real(sgl) :: sq3   ! help variable
  real(sgl) :: vr1   ! function for volume of the projectile-like section
  real(sgl) :: vr2   ! function for volume of the neck
  real(sgl) :: vr3   ! help value
  real(sgl) :: wu1   ! help variable
  real(sgl) :: wu3   ! help variable
  real(sgl) :: x(i)  ! help variable
  real(sgl) :: z21   ! asymmetry value
  real(sgl) :: z32   ! asymmetry value
!
! **********************************************************************
!
! fidi: function for dinuclear shape
!
  d = totl - r1 - r3
  aaa = .8 + x(1) * 2.
  z2 = d * (.5 - x(2))
  z1 = r1 * x(3)
  z3 = d - r3 * x(4)
  r1 = rp * x(5)
  r3 = rt * x(6)
  cur = c0 * abs(x(7))
  wu1 = r1 **2 - z1 **2
  if (wu1 < 1.e-8) wu1 = 1.e-8
  wu3 = r3 **2 - (d - z3) **2
  if (wu3 < 1.e-8) wu3 = 1.e-8
  sq1 = sqrt(wu1)
  sq3 = sqrt(wu3)
  z21 = (z2 - z1) / aaa
  z32 = (z3 - z2) / aaa
  if (z21 > 30.) z21 = 30.
  if (z32 > 30.) z32 = 30.
  al = (vtot - vr1(z1) - vr2(z1, z2, z3) - vr3(z3)) **2
  be = ((rp **3 - r1 **3) - rpt * (rt **3 - r3 **3)) **2 * 1.e1
  ga = (sq1 - r2 - aaa * aaa * cur * (cosh(z21) - 1.)) **2 * 1.e5
  de = (sq3 - r2 - aaa * aaa * cur * (cosh(z32) - 1.)) **2 * 1.e5
  eq = (z1 / sq1 - aaa * cur * sinh(z21)) **2 * 1.e5
  epsc = ((d - z3) / sq3 - aaa * cur * sinh(z32)) **2 * 1.e5
  ds = x(7) **2 * 1.e-4
  fidi = al + be + ga + de + eq + epsc + ds
  return
end function fidi
! Copyright A.J. Koning 2021
