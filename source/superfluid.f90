function superfluid(Zix, Nix, ald, Eex, P, ibar)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Superfluid model level density formula
!
! Author    : Arjan Koning and Stephane Hilaire
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl          ! single precision kind
!   dbl          ! double precision kind
! Variables for level density
!   pair         ! pairing energy
!   Pshift       ! adjustable pairing shift
! Constants
!   sqrttwopi    ! sqrt(2. * pi)
! Variables for level density
!   Dcrit        ! critical determinant
!   Scrit        ! critical entropy
!   Tcrit        ! critical temperature
!   Ucrit        ! critical U
!
! *** Declaration of local data
!
  implicit none
  integer   :: ibar       ! fission barrier
  integer   :: Nix        ! neutron number index for residual nucleus
  integer   :: Zix        ! charge number index for residual nucleus
  real(sgl) :: ald        ! level density parameter
  real(sgl) :: Df         ! determinant
  real(sgl) :: Eex        ! excitation energy
  real(sgl) :: P          ! pairing energy
  real(sgl) :: phi1       ! phi function of superfluid model
  real(sgl) :: phi2       ! help variable
  real(sgl) :: Sf         ! entropy
  real(sgl) :: sigma      ! help variable
  real(sgl) :: spincut    ! spin cutoff factor
  real(sgl) :: Tf         ! temperature
  real(sgl) :: U          ! excitation energy minus pairing energy
  real(dbl) :: fermi      ! function for Fermi gas level density formula
  real(dbl) :: superfluid ! Superfluid model level density formula
!
! *********************** Level density formula ************************
!
! fermi     : function for Fermi gas level density formula
! phi1      : phi function of superfluid model
!
! Superfluid model
!
  U = Eex + pair(Zix, Nix) + Pshift(Zix, Nix, ibar)
  if (U > 0.) then
    if (U > Ucrit(Zix, Nix, ibar)) then
      superfluid = fermi(Zix, Nix, ald, Eex, P, ibar)
    else
      phi2 = 1. - U / Ucrit(Zix, Nix, ibar)
      Df = Dcrit(Zix, Nix, ibar) * (1. - phi2) * (1. + phi2) * (1. + phi2)
      phi1 = sqrt(phi2)
      Tf = 2. * Tcrit(Zix, Nix) * phi1 / log((phi1 + 1.) / (1. - phi1))
      Sf = Scrit(Zix, Nix, ibar) * Tcrit(Zix, Nix) / Tf * (1. - phi2)
      sigma = sqrt(spincut(Zix, Nix, ald, Eex, ibar, 0))
      superfluid = exp(dble(Sf)) / sqrt(Df) / sqrttwopi / sigma
    endif
  else
    superfluid = 1.
  endif
  return
end function superfluid
! Copyright A.J. Koning 2021
