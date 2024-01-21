subroutine basicxs(Zcomp, Ncomp)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Basic cross sections and transmission coefficients
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
! Variables for compound reactions
!   flagcomp       ! flag for compound angular distribution calculation
! Variables for OMP
!   flagomponly    ! flag to execute ONLY an optical model calculation
!
! *** Declaration of local data
!
  implicit none
  integer :: Ncomp              ! neutron number index for compound nucleus
  integer :: Zcomp              ! proton number index for compound nucleus
!
! ******************* ECIS calculations and output *********************
!
! inverse     : subroutine for ECIS calculation of total, reaction and elastic cross sections and transmission coefficients
!               for outgoing energy grid.
! basicinitial: subroutine for initialization of arrays for basic cross sections
!
! The transmission coefficients and inverse reaction cross sections for the outgoing energy grid need to be calculated only once,
! for the maximal incident energy.
!
  if (flagomponly .and. .not. flagcomp) return
  call basicinitial
  call inverse(Zcomp, Ncomp)
  return
end subroutine basicxs
! Copyright A.J. Koning 2021
