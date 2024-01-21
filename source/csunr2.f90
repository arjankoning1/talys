subroutine csunr2(ay, l, j)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Unresolved resonance cross section from NJOY
!
! Author    : Gilles Noguere
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
! Variables for main input
!   Atarget     ! mass number of target nucleus
! Variables for energy grid
!   Einc        ! incident energy in MeV
! Variables for energies
!   wavenum     ! wave number
! Variables for compound nucleus from target
!   nulj        ! (l, j) number of degrees of freedom for URR calculation
! Variables for nuclides
!   tarmass     ! mass of target nucleus
! Constants
!   onethird    ! 1 / 3
!   pi          ! pi
!   twopi       ! 2 * pi
! Variables for resonance parameters
!   Dlj         ! mean resonance spacing per J, l value
! Variables for levels
!   jdis        ! spin of level
! Variables for URR
!   sigurrc     ! (l, j) capture cross section for URR
!   sigurrf     ! (l, j) fission cross section for URR
!   sigurrs     ! (l, j) scattering cross section for URR
!   spot        ! potential scattering contribution
!   urrwidth    ! channel width in URR
!
! *** Declaration of local data
!
  implicit none
  integer   :: j      ! counter
  integer   :: l      ! multipolarity
  integer   :: lamda  ! number of degress of freedom per l,j
  integer   :: mu     ! help variable
  integer   :: nu     ! number of radial grid point
  integer   :: odd    ! odd (1) or even (0) nucleus
  real(sgl) :: ac     ! channel radius in ENDF-6 convention
  real(sgl) :: add    ! help variable
  real(sgl) :: aj     ! spin
  real(sgl) :: amun   ! number of degrees of freedom for l,j
  real(sgl) :: ay     ! help variable
  real(sgl) :: const  ! constant
  real(sgl) :: den    ! help variable
  real(sgl) :: diff   ! difference
  real(sgl) :: dx     ! increment
  real(sgl) :: e      ! energy
  real(sgl) :: galpha ! URR width
  real(sgl) :: gam    ! Brosa parameter
  real(sgl) :: gbeta  ! URR width
  real(sgl) :: gc     ! help variable
  real(sgl) :: gff    ! help variable
  real(sgl) :: gfx    ! URR width
  real(sgl) :: ggx    ! URR width
  real(sgl) :: gj     ! spin factor
  real(sgl) :: gnox   ! URR width
  real(sgl) :: gnx    ! URR width
  real(sgl) :: gs     ! single-particle level density parameter
  real(sgl) :: gxx    ! URR width
  real(sgl) :: k      ! designator for particle
  real(sgl) :: ll     ! angular momentum
  real(sgl) :: ps     ! phase shift
  real(sgl) :: rho    ! integrated level density
  real(sgl) :: rhoc   ! help variable
  real(sgl) :: spi    ! spin
  real(sgl) :: sqrte  ! square root of energy
  real(sgl) :: temp   ! nuclear temperature
  real(sgl) :: terf   ! help variable
  real(sgl) :: terg   ! help variable
  real(sgl) :: ters   ! help variable
  real(sgl) :: vl     ! penetrability factor
!
!   link with NJOY parameters
!
  spot(l) = 0.
  gnox = urrwidth(3, l, j)
  gxx = urrwidth(1, l, j)
  ggx = urrwidth(0, l, j)
  gfx = urrwidth( - 1, l, j)
  dx = Dlj(l, j)
  if (dx == 0.) return
  amun = real(nulj(0, l, j))
  mu = min(max(nulj(0, l, j), 1), 4)
  nu = min(max(nulj( - 1, l, j), 1), 4)
  lamda = min(max(nulj(1, l, j), 1), 4)
  e = Einc * 1.e6
  sqrte = sqrt(e)
  k = wavenum * 10.
  const = 2. * e * pi **2 / k **2
  odd = mod(Atarget + 1, 2)
  aj = j + 0.5 * odd
  ll = real(l)
  spi = jdis(0, 1, 0)
  gj = (2. * aj + 1.) / (4. * spi + 2.)
  ac = 0.123 * tarmass **onethird + 0.08
  rho = k * ac
  rhoc = k * ay
!
!   calculate penetrability (vl) and phase shift(ps)
!
  call unfac(l, rho, rhoc, amun, vl, ps)
  vl = vl * sqrte
!
!   calculate potential scattering
!
  spot(l) = 2. * twopi * (2. * ll + 1.) * (sin(ps) / k) **2
!
!   compute cross section contributions
!
  gnx = gnox * vl
  gam = ggx
  galpha = gnx
  gbeta = gfx
  diff = gxx
  den = e * dx
  temp = const * gj * gnx / den
  terg = temp * gam
  ters = temp * gnx
  terf = temp * gbeta
!
!   calculate fluctuation integrals
!
  call gnrl(galpha, gbeta, gam, mu, nu, lamda, gs , diff, 1)
  call gnrl(galpha, gbeta, gam, mu, nu, lamda, gc , diff, 2)
  call gnrl(galpha, gbeta, gam, mu, nu, lamda, gff, diff, 3)
  gc = gc * terg
  gff = gff * terf
  gs = gs * ters
!
! add interference correction term
!
  add = const * gj * 2. * gnx * (sin(ps)) **2
  add = add / (e * dx)
  gs = gs - add
!
! cross sections
!
  sigurrs(l, j) = gs
  sigurrf(l, j) = gff
  sigurrc(l, j) = gc
  return
end subroutine csunr2
! Copyright A.J. Koning 2021
