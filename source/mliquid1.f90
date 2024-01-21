function mliquid1(Z, A)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Myers-Swiatecki liquid drop mass
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
! Definition of single and double precision variables
!   dbl         ! double precision kind
! Constants
!   amu         ! atomic mass unit in MeV
!   onethird    ! 1 / 3
!   twothird    ! 2 / 3
!
! *** Declaration of local data
!
  implicit none
  integer   :: A        ! mass number of target nucleus
  integer   :: N        ! neutron number of residual nucleus
  integer   :: oddN     ! help variable
  integer   :: oddZ     ! help variable
  integer   :: Z        ! charge number of target nucleus
  real(dbl) :: a1       ! Myers-Swiatecki parameter
  real(dbl) :: a2       ! Myers-Swiatecki parameter
  real(dbl) :: c3       ! character
  real(dbl) :: c4       ! character
  real(dbl) :: deltaP   ! pairing energy
  real(dbl) :: Ecoul    ! Coulomb energy
  real(dbl) :: Eldm     ! liquid drop energy
  real(dbl) :: Esur     ! surface energy
  real(dbl) :: Ev       ! volume energy
  real(dbl) :: factor   ! multiplication factor
  real(dbl) :: kappa    ! Myers-Swiatecki parameter
  real(dbl) :: Mh       ! mass excess of neutron and hydrogen
  real(dbl) :: mliquid1 ! function for liquid drop mass
  real(dbl) :: Mn       ! mass excess of neutron and hydrogen
  real(dbl) :: rA       ! mass number of residual nucleus
  real(dbl) :: rN       ! neutron number of residual nucleus
  real(dbl) :: rZ       ! charge number of residual nucleus
!
! ************** Spherical Myers-Swiatecki parameters ******************
!
! mliquid1 : function for liquid drop mass
!
! Myers-Swiatecki model: Nucl. Phys. 81 (1966) 1.
! We use the original M-S parameters, see Mengoni and Nakajima, J. Nucl. Sci. Technol.  31[2], p 151 (1994).
!
  a1 = 15.677
  a2 = 18.56
  kappa = 1.79
  c3 = 0.717
  c4 = 1.21129
  Mn = 8.07144
  Mh = 7.28899
  N = A - Z
  rA = real(A)
  rZ = real(Z)
  rN = real(N)
  factor = 1. - kappa * ((rN - rZ) / rA) **2
  Ev = - a1 * factor * rA
  Esur = a2 * factor * rA **twothird
  Ecoul = c3 * rZ **2 / (rA **onethird) - c4 * rZ **2 / rA
  oddZ = mod(Z, 2)
  oddN = mod(N, 2)
  if (oddZ /= oddN) deltaP = 0.
  if (oddZ == 0 .and. oddN == 0) deltaP = - 11. / sqrt(rA)
  if (oddZ == 1 .and. oddN == 1) deltaP = 11. / sqrt(rA)
  Eldm = Z * Mh + N * Mn + Ev + Esur + Ecoul + deltaP
  mliquid1 = A + Eldm / amu
  return
end function mliquid1
! Copyright A.J. Koning 2021
