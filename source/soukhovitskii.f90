subroutine soukhovitskii(k, Z, A, eopt)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Global optical model parameters for actinides by
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
!   sgl         ! single precision kind
! Constants
!   onethird    ! 1 / 3
! Variables for OMP
!   av          ! real volume diffuseness
!   avd         ! real surface diffuseness
!   avso        ! real spin - orbit diffuseness
!   aw          ! imaginary volume diffuseness
!   awd         ! imaginary surface diffuseness
!   awso        ! imaginary spin - orbit diffuseness
!   Fav         ! adjustable factor for OMP (default 1.)
!   Favd        ! adjustable factor for OMP (default 1.)
!   Favso       ! adjustable factor for OMP (default 1.)
!   Faw         ! adjustable factor for OMP (default 1.)
!   Fawd        ! adjustable factor for OMP (default 1.)
!   Fawso       ! adjustable factor for OMP (default 1.)
!   Fd1         ! adjustable factor for OMP (default 1.)
!   Fd2         ! adjustable factor for OMP (default 1.)
!   Fd3         ! adjustable factor for OMP (default 1.)
!   Frc         ! adjustable factor for OMP (default 1.)
!   Frv         ! adjustable factor for OMP (default 1.)
!   Frvd        ! adjustable factor for OMP (default 1.)
!   Frvso       ! adjustable factor for OMP (default 1.)
!   Frw         ! adjustable factor for OMP (default 1.)
!   Frwd        ! adjustable factor for OMP (default 1.)
!   Frwso       ! adjustable factor for OMP (default 1.)
!   Fv1         ! adjustable factor for OMP (default 1.)
!   Fvso1       ! adjustable factor for OMP (default 1.)
!   Fvso2       ! adjustable factor for OMP (default 1.)
!   Fw1         ! adjustable factor for OMP (default 1.)
!   Fw2         ! adjustable factor for OMP (default 1.)
!   Fwso1       ! adjustable factor for OMP (default 1.)
!   Fwso2       ! adjustable factor for OMP (default 1.)
!   rc          ! Coulomb radius
!   rv          ! real volume radius
!   rvd         ! real surface radius
!   rvso        ! real spin - orbit radius
!   rw          ! imaginary volume radius
!   rwd         ! imaginary surface radius
!   rwso        ! imaginary spin - orbit radius
!   v           ! real volume depth
!   vd          ! real surface depth
!   vso         ! real spin - orbit depth
!   w           ! imaginary volume depth
!   wd          ! imaginary surface depth
!   wso         ! imaginary spin - orbit depth
! Variables for masses
!   S           ! separation energy
!
! *** Declaration of local data
!
  implicit none
  integer   :: A       ! mass number of target nucleus
  integer   :: k       ! designator for particle
  integer   :: Z       ! charge number of target nucleus
  real(sgl) :: asym    ! asymmetry parameter
  real(sgl) :: Ccoul   ! optical potential parameter
  real(sgl) :: Crr     ! optical potential parameter
  real(sgl) :: Cviso   ! optical potential parameter
  real(sgl) :: Cwiso   ! optical potential parameter
  real(sgl) :: d1loc   ! help variable
  real(sgl) :: d2loc   ! help variable
  real(sgl) :: d3loc   ! help variable
  real(sgl) :: eferm   ! Fermi energy
  real(sgl) :: eopt    ! incident energy
  real(sgl) :: f       ! E-Ef
  real(sgl) :: lambdaR ! Soukhovitskii OMP parameter
  real(sgl) :: phicoul ! Soukhovitskii OMP parameter
  real(sgl) :: rr      ! running variable in integration over the radius
  real(sgl) :: V0r     ! Soukhovitskii OMP parameter
  real(sgl) :: v1r     ! Soukhovitskii OMP parameter
  real(sgl) :: v2r     ! Soukhovitskii OMP parameter
  real(sgl) :: Var     ! Soukhovitskii OMP parameter
  real(sgl) :: viso    ! Soukhovitskii OMP parameter
  real(sgl) :: vrdisp  ! Soukhovitskii OMP parameter
  real(sgl) :: vso1loc ! help variable
  real(sgl) :: vso2loc ! help variable
  real(sgl) :: w1loc   ! help variable
  real(sgl) :: w2loc   ! help variable
  real(sgl) :: Wad     ! Soukhovitskii OMP parameter
  real(sgl) :: wddisp  ! Soukhovitskii OMP parameter
  real(sgl) :: widr    ! Soukhovitskii OMP parameter
  real(sgl) :: wso1loc ! help variable
  real(sgl) :: wso2loc ! help variable
!
! *** Parameters of Soukhovitskii et al, J. Phys. G30, p. 905 (2004) ***
!
  asym = (A - 2. * Z) / real(A)
  if (k == 1) then
    eferm = - 0.5 * (S(0, 1, 1) + S(0, 0, 1))
  else
    eferm = - 0.5 * (S(1, 0, 2) + S(0, 0, 2))
  endif
  f = max(eopt - eferm, -20.)
  Cviso = 10.5
  V0r = - 41.45
  Var = - 0.06667
  vrdisp = 92.44
  v1r = 0.03
  v2r = 2.05e-4
  lambdaR = 3.9075e-3
  viso = 1. + (( - 1) **k) * Cviso * asym / (V0r + Var * (A - 232.) + vrdisp)
  v = (V0r + Var * (A - 232.) + v1r * f + v2r * (f **2) + vrdisp * exp( - lambdaR * f)) * viso
  if (k == 2) then
    Ccoul = 0.9
    phicoul = (lambdaR * vrdisp * exp( - lambdaR * f) - v1r - 2. * v2r * f) * viso
    v = v + Ccoul * Z / (A **onethird) * phicoul
  endif
  v = Fv1 * v
  rr = 1.245
  Crr = 0.05
  widr = 100.
  rv = Frv * rr * (1. - Crr * f **2 / (f **2 + widr **2))
  av = Fav * (0.660 + 2.53e-4 * eopt)
  w1loc = Fw1 * 14.74
  w2loc = Fw2 * 81.63
  w = w1loc * f **2 / (f **2 + w2loc **2)
  rw = Frw * 1.2476
  aw = Faw * 0.594
  vd = 0.
  rvd = Frvd * 1.2080
  avd = Favd * 0.614
  wddisp = 17.38
  Cwiso = 24.
  Wad = 0.03833
  d1loc = Fd1 * (wddisp + Wad * (A - 232.) + (( - 1) **k) * Cwiso * asym)
  d2loc = Fd2 * 0.01759
  d3loc = Fd3 * 11.79
  wd = d1loc * f **2 * exp( - d2loc * f) / (f **2 + d3loc **2)
  rwd = Frwd * 1.2080
  awd = Fawd * 0.614
  vso1loc = Fvso1 * 5.86
  vso2loc = Fvso2 * 0.0050
  vso = vso1loc * exp( - vso2loc * f)
  rvso = Frvso * 1.1213
  avso = Favso * 0.59
  wso1loc = - 3.1 * Fwso1
  wso2loc = Fwso2 * 160.
  wso = wso1loc * f **2 / (f **2 + wso2loc **2)
  rwso = Frwso * 1.1213
  awso = Fawso * 0.59
  if (k == 1) then
    rc = 0.
  else
    rc = Frc * 1.2643
  endif
  return
end subroutine soukhovitskii
! Copyright A.J. Koning 2021
