subroutine gamma(Zcomp, Ncomp)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Gamma cross section and transmission coefficients
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
! Variables for gamma rays
!   flaggamma    ! flag for output of gamma - ray information
!
! *** Declaration of local data
!
  implicit none
  integer :: Ncomp ! neutron number index for compound nucleus
  integer :: Zcomp ! proton number index for compound nucleus
!
! ********** Gamma cross section and transmission coefficients *********
!
! tgamma   : subroutine for photon transmission coefficients
! gammaout : subroutine for output of gamma-ray strength functions, transmission coefficients and cross sections
!
  call tgamma(Zcomp, Ncomp)
  if (flaggamma .and. nin == 1) call gammaout(Zcomp, Ncomp)
  return
end subroutine gamma
! Copyright A.J. Koning 2022
