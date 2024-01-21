subroutine direct
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculation of direct inelastic cross sections
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
! Variables for direct reactions
!   flagdirect     ! flag for output of direct reaction results
! Variables for OMP
!   flagomponly    ! flag to execute ONLY an optical model calculation
! Variables for main input
!   k0             ! index of incident particle
!   Ltarget        ! excited level of target
! Variables for energies
!   flaggiant      ! flag for collective contribution from giant resonances
!
! ****************** Direct cross section calculation ******************
!
! directecis : subroutine for ECIS calculation of direct cross section
! directread : subroutine to read ECIS results for direct cross section
! giant      : subroutine for giant resonance contribution
! directout  : subroutine for output of direct reaction cross sections
!
  if (k0 == 0) return
  if (Ltarget /= 0) return
  if (flagomponly) return
  call directecis
  call directread
  if (flaggiant) call giant
  if (flagdirect) call directout
  return
end subroutine direct
! Copyright A.J. Koning 2021
