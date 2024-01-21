function gilcam(Zix, Nix, ald, Eex, P, ibar)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Gilbert-Cameron level density formula
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
!   sgl        ! single precision kind
!   dbl        ! double precision kind
! Variables for level density
!   E0         ! particle constant of temperature formula
!   Exmatch    ! matching point for Ex
!   T          ! temperature
!
! *** Declaration of local data
!
  implicit none
  integer   :: ibar   ! fission barrier
  integer   :: Nix    ! neutron number index for residual nucleus
  integer   :: Zix    ! charge number index for residual nucleus
  real(sgl) :: ald    ! level density parameter
  real(sgl) :: Eex    ! excitation energy
  real(sgl) :: P      ! pairing energy
  real(dbl) :: expo   ! help variable
  real(dbl) :: fermi  ! function for Fermi gas level density formula
  real(dbl) :: gilcam ! Gilbert-Cameron level density formula
!
! *********************** Level density formula ************************
!
! fermi  : function for Fermi gas level density formula
!
  if (Eex > Exmatch(Zix, Nix, ibar)) then
!
! 1. Fermi Gas
!
    gilcam = fermi(Zix, Nix, ald, Eex, P, ibar)
  else
!
! 2. Constant temperature
!
    expo = min((Eex - E0(Zix, Nix, ibar)) / T(Zix, Nix, ibar), 300.)
    gilcam = exp(expo) / T(Zix, Nix, ibar)
  endif
  return
end function gilcam
! Copyright A.J. Koning 2021
