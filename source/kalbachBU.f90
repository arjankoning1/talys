function kalbachBU(type, Ein, ang)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Kalbach systematics for break-up
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
!   sgl       ! single precision kind
! Variables for main input
!   k0        ! index of incident particle
! Constants
!   fourpi    ! 4. * pi
!   parA      ! mass number of particle
!   parN      ! neutron number of particle
!   parZ      ! charge number of particle
!   twopi     ! 2 * pi
! Variables for preequilibrium
!   Ca        ! effective Coulomb barrier
!   Ecent     ! centroid energy for emission spectrum
!   Sab       ! separation energy for projectile
!
! *** Declaration of local data
!
  implicit none
  integer   :: Ab        ! help variable
  integer   :: Ap        ! two-component Pauli blocking correction factor
  integer   :: Npr       ! neutron number of projectile
  integer   :: type      ! particle type
  integer   :: Za        ! charge of projectile
  integer   :: Zb        ! charge of ejectile
  real(sgl) :: abu       ! Kalbach parameter
  real(sgl) :: ang       ! angle
  real(sgl) :: ang0      ! angular barrier
  real(sgl) :: Ein       ! incident energy
  real(sgl) :: K1        ! constant of Kalbach systematics
  real(sgl) :: K2        ! constant of Kalbach systematics
  real(sgl) :: K3        ! constant of Kalbach systematics
  real(sgl) :: kalbachBU ! Kalbach function for break-up
  real(sgl) :: Tc        ! function for Coulomb dip
  real(sgl) :: term      ! help variable
  real(sgl) :: wang      ! angular width
!
! ************************ Kalbach systematics *************************
!
! kalbachBU : Kalbach function for break-up
! Tc        : function for Coulomb dip
!
!  Systematics of Kalbach: Phys. Rev. C95, 014606, (2017)
!
  kalbachBU = 1./fourpi
  Ap = parA(k0)
  Npr = parN(k0)
  Za = parZ(k0)
  Ab = parA(type)
  Zb = parZ(type)
  abu = 4.7 + Ab
  if (Ab == Ap - 1) then
    if (Sab > 0.) then
      term = 1. + exp((12. * Sab - Ecent) / (0.84 * Sab))
      abu = 4. * Ab + Zb - 2. + 0.029 * Ecent + 7.6 / Ap / term
    endif
  endif
  K1 = 1.8
  K2 = 0.26
  K3 = 0.035
  if (Npr > 0.) then
    term = 1 + exp((K2 - Ca / Ein) / (K3 * Za / Npr))
    ang0 = K1 / (Ap **1.5) / term
    wang = min(0.09, ang0 / 3.)
    term = 1. + exp((ang0 - ang) / wang)
    Tc = 1. / term
    kalbachBU = (abu * abu + 1.) / twopi * exp( - abu * cos(ang)) * Tc
  endif
  return
end function kalbachBU
! Copyright A.J. Koning 2021
