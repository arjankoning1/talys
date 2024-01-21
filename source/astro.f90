subroutine astro
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Astrophysical reaction rates
!
! Author    : Stephane Goriely and Stephane Hilaire
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
! stellarrate: subroutine for calculation of reaction rate for a Maxwell-Boltzmann distribution
! astroout   : subroutine for output of astrophysical results
!
  call stellarrate
  call astroout
  return
end subroutine astro
! Copyright A.J. Koning 2021
