subroutine endf
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Cross sections and information for ENDF-6 file
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
! Variables for input energies
!   Ninc         ! number of incident energies
! Variables for main input
!   k0             ! index of incident particle
! Variables for ECIS
!   flagecisinp    ! flag for existence of ecis input file
!
! ************************** ECIS calculation **************************
!
! endfinfo    : subroutine for info for ENDF-6 file
! endfenergies: subroutine for energy grid for ENDF-6 file
! endfecis    : subroutine for ECIS calculation for incident particle on ENDF-6 energy grid
! endfread    : subroutine to read ECIS results for incident particle on ENDF-6 energy grid
!
  if (Ninc == 1) return
  call endfinfo
  call endfenergies
  if (k0 /= 0) then
    call endfecis
    if (flagecisinp) call endfread
  endif
  return
end subroutine endf
! Copyright A.J. Koning 2021
