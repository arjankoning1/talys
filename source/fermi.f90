function fermi(Zix, Nix, ald, Eex, P, ibar)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Fermi gas level density formula
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               dbl  , &           ! double precision kind
               sgl                ! single precision kind
!
! *** Declaration of local data
!
  integer   :: ibar    ! fission barrier
  integer   :: Nix     ! neutron number index for residual nucleus
  integer   :: Zix     ! charge number index for residual nucleus
  real(sgl) :: ald     ! level density parameter
  real(sgl) :: denom   ! help variable
  real(sgl) :: Eex     ! excitation energy
  real(sgl) :: P       ! pairing energy
  real(sgl) :: sigma   ! help variable
  real(sgl) :: spincut ! spin cutoff factor
  real(sgl) :: U       ! excitation energy minus pairing energy
  real(dbl) :: factor  ! multiplication factor
  real(dbl) :: fermi   ! function for Fermi gas level density formula
!
! *********************** Total level density **********************
!
  U = Eex - P
  if (U > 0.) then
    factor = min(2. * sqrt(ald * U), 700.)
    sigma = sqrt(spincut(Zix, Nix, ald, Eex, ibar, 0))
    denom = 12. * sqrt(2.) * sigma * (ald **0.25) * (U **1.25)
    fermi = exp(factor) / denom
  else
    fermi = 1.
  endif
  return
end function fermi
! Copyright A.J. Koning 2021
