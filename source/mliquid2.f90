function mliquid2(Z, A)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Goriely liquid drop mass
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
  integer   :: Z        ! charge number of target nucleus
  real(dbl) :: ac       ! channel radius in ENDF-6 convention
  real(dbl) :: as       ! diffuseness of the functional form of the imaginary volume integral
  real(dbl) :: ass      ! Goriely parameter
  real(dbl) :: asym     ! asymmetry parameter
  real(dbl) :: avol     ! Goriely parameter
  real(dbl) :: Ecoul    ! Coulomb energy
  real(dbl) :: Eldm     ! liquid drop energy
  real(dbl) :: Esur     ! surface energy
  real(dbl) :: Esym     ! symmetry energy
  real(dbl) :: Ev       ! volume energy
  real(dbl) :: factor   ! multiplication factor
  real(dbl) :: Mh       ! mass excess of neutron and hydrogen
  real(dbl) :: mliquid2 ! function for liquid drop mass (Goriely)
  real(dbl) :: Mn       ! mass excess of neutron and hydrogen
  real(dbl) :: rA       ! mass number of residual nucleus
  real(dbl) :: rN       ! neutron number of residual nucleus
  real(dbl) :: rZ       ! charge number of residual nucleus
!
! ****************** Spherical Goriely parameters **********************
!
! mliquid2: function for liquid drop mass
!
! Goriely model: S. Goriely, ND2001, Tsukuba, Japan, J. Nuc. Sci Techn. August 2002, suppl 2. p.536 (2002)
!
  avol = - 15.6428
  as = 17.5418
  asym = 27.9418
  ass = - 25.3440
  ac = 0.70
  Mn = 8.07132281
  Mh = 7.2889694
  N = A - Z
  rA = real(A)
  rZ = real(Z)
  rN = real(N)
  factor = (rN - rZ) / rA
  Ev = avol * rA
  Esur = as * (rA **twothird)
  Esym = (asym + ass * (rA **( - onethird))) * rA * factor **2
  Ecoul = ac * rZ **2 / (rA **onethird)
  Eldm = Z * Mh + N * Mn + Ev + Esur + Esym + Ecoul - 1.4333e-5 * rZ **2.39
  mliquid2 = A + Eldm / amu
  return
end function mliquid2
! Copyright A.J. Koning 2021
