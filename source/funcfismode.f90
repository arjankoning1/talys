function funcfismode(Zix, Nix, Epar)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Transmission coefficient per fission mode
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
!   sgl        ! single precision kind
!   dbl        ! double precision kind
! Constants
!   pi         ! pi
! Variables for mass distribution
!   excfis     ! excitation energy at fission
! Variables for Brosa model
!   bf         ! barrier heigt
!   bfsplin    ! barrier height splin fit parameters
!   hw         ! barrier width
!   hwsplin    ! barrier width splin fit parameters
!   numtemp    ! running variable denoting the numtempth temperature
!
! *** Declaration of local data
!
  implicit none
  integer   :: Nix          ! neutron number index for residual nucleus
  integer   :: Zix          ! charge number index for residual nucleus
  real(sgl) :: ald          ! level density parameter
  real(sgl) :: bft          ! barrier
  real(sgl) :: Epar         ! energy
  real(sgl) :: hwt          ! width
  real(sgl) :: ignatyuk     ! function for energy dependent level density parameter a
  real(sgl) :: temps(9)     ! temperature
  real(sgl) :: tmp          ! temperature
  real(dbl) :: density      ! level density
  real(dbl) :: expo         ! help variable
  real(dbl) :: funcfismode  ! function for transmission coefficient per fission mode
  temps =   (/ 0.0, 0.3, 0.6, 0.9, 1.2, 1.6, 2.0, 2.5, 3.0 /)
!
! **********************************************************************
!
  ald = ignatyuk(Zix, Nix, Epar, 0)
  tmp = sqrt(Epar / ald)
  call splint(temps, bf, bfsplin, numtemp, tmp, bft)
  call splint(temps, hw, hwsplin, numtemp, tmp, hwt)
  expo = 2. * pi * (bft + Epar - excfis) / hwt
  if (expo <= 80.) then
    funcfismode = 1. / (1. + exp(expo))
  else
    funcfismode = exp( - 80.) / (exp( - 80.) + exp(expo - 80.))
  endif
  if (Epar > 0.) funcfismode = funcfismode * density(Zix, Nix, Epar, 0., 1, 0, 1)
  return
end function funcfismode
! Copyright A.J. Koning 2021
