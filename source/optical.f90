subroutine optical(Zix, Nix, kopt, eopt)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Determination of optical potential
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               sgl      ! single precision kind
!
! *** Declaration of local data
!
  implicit none
  integer   :: kopt ! optical model index for fast particle
  integer   :: Nix  ! neutron number index for residual nucleus
  integer   :: Zix  ! charge number index for residual nucleus
  real(sgl) :: eopt ! incident energy
!
! ********************** Call optical model module *********************
!
  eopt = max(eopt, 0.)
  if (kopt <= 2) then
    call opticalnp(Zix, Nix, kopt, eopt)
  else
    call opticalcomp(Zix, Nix, kopt, eopt)
  endif
  return
end subroutine optical
! Copyright A.J. Koning 2021
