function ignatyuk(Zix, Nix, Eex, ibar)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Energy dependent level density parameter a
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
!   sgl             ! single precision kind
! Variables for level density
!   alimit          ! asymptotic level density parameter
!   cfermi          ! width of Fermi distribution for damping
!   deltaW          ! shell correction in nuclear mass
!   flagcolldamp    ! flag for damping of coll. effects in eff. level density (without explicit coll. enh.)
!   gammald         ! gamma - constant for asymptotic level density parameter
!   Ufermi          ! energy of Fermi distribution for damping
! Variables for nuclides
!   AA              ! mass number of residual nucleus
! Variables for level density
!   delta           ! energy shift
!
! *** Declaration of local data
!
  implicit none
  integer   :: A        ! mass number of target nucleus
  integer   :: ibar     ! fission barrier
  integer   :: Nix      ! neutron number index for residual nucleus
  integer   :: Zix      ! charge number index for residual nucleus
  real(sgl) :: aldlim   ! asymptotic level density parameter
  real(sgl) :: aldlow   ! lower limit of a
  real(sgl) :: damp     ! shell damping factor
  real(sgl) :: Eex      ! excitation energy
  real(sgl) :: expo     ! help variable
  real(sgl) :: fU       ! help variable
  real(sgl) :: ignatyuk ! function for energy dependent level density parameter a
  real(sgl) :: qfermi   ! Fermi distribution
  real(sgl) :: U        ! excitation energy minus pairing energy
!
! *********************** Level density formula ************************
!
! ignatyuk: function for energy dependent level density parameter a
!
! Formalism from Ignatyuk et al. Sov. Jour. Nuc. Phys. 21 (1975), 255.
!
  U = Eex-delta(Zix, Nix, ibar)
  aldlim = alimit(Zix, Nix)
!
! 1. For very low Eex, i.e. U < 0, we use the first order Taylor
!    expansion
!
  if (U <= 0.) then
    damp = (1. + deltaW(Zix, Nix, ibar) * gammald(Zix, Nix))
  else
!
! 2. Higher energies
!
    expo = gammald(Zix, Nix) * U
    fU = 1.
    if (abs(expo) <= 80.) fU = 1. - exp( - expo)
    damp = 1. + fU * deltaW(Zix, Nix, ibar) / U
  endif
!
! Only for fission models with damping of collective effects in
! effective level density model. (Re-installed from TALYS-0.64).
! Fermi distribution for asymptotic level density parameter. This takes
! the damping of collective effects into account in a phenomenological
! way. The level density parameters are then equal to A/13 in the
! high-energy limit rather than A/8.
!
!   level density (without explicit collective enhancement)
!   Only used for Bruyeres-le-Chatel (Pascal Romain) fission
!   model
!
  if (flagcolldamp) then
    A = AA(Zix, Nix, 0)
    aldlow = A / 13.
    expo = (U - Ufermi(Zix,Nix,0)) / cfermi(Zix,Nix,0)
    qfermi = 0.
    if (expo >  - 80.) qfermi = 1. / (1. + exp( - expo))
    aldlim = aldlow * qfermi + aldlim * (1. - qfermi)
  endif
  ignatyuk = max(aldlim * damp, 1.)
  return
end function ignatyuk
! Copyright A.J. Koning 2021

