subroutine basicinitial
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialization of arrays for various cross sections
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
! All global variables
!   numen         ! maximum number of outgoing energies
! Variables for inverse channel data
!   Tjl           ! transmission coefficient per particle, energy, spin and l-value
!   Tl            ! transmission coefficients per particle, energy and l-value
!   xselas        ! total elastic cross section (shape + compound)
!   xsopt         ! optical model reaction cross section
!   xsreac        ! reaction cross section
!   xstot         ! total cross section (neutrons only)
! Variables for gamma rays
!   gammax        ! number of l - values for gamma multipolarity
!   lmax          ! maximal l - value for transmission coefficients
!
! *** Declaration of local data
!
  implicit none
  integer           :: nen                  ! energy counter
!
! Initialization
!
  Tjl = 0.
  Tl = 0.
  xselas = 0.
  xsopt = 0.
  xsreac = 0.
  xstot = 0.
  lmax = 0
  do nen = 0, numen
     lmax(0, nen) = gammax
  enddo
  return
end subroutine basicinitial
! Copyright A.J. Koning 2021
