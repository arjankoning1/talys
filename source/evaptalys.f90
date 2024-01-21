subroutine evaptalys
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Execute TALYS subroutines for secondary evaporation
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
! Variables for output
!   flagddx         ! flag for output of double - differential cross sections
!   flagspec        ! flag for output of spectra
! Variables for basic reaction
!   flagchannels    ! flag for exclusive channels calculation
! Variables for compound reactions
!   flagcomp        ! flag for compound angular distribution calculation
! Variables for input energies
!   enincmax        ! maximum incident energy
! Variables for preequilibrium
!   epreeq          ! on - set incident energy for preequilibrium calculation
! Variables for energies
!   flagpreeq       ! flag for pre - equilibrium calculation
! Variables for nuclides
!   parinclude      ! logical to include outgoing particle
!
! *** Declaration of local data
!
  implicit none
!
! These subroutines are executed for each nuclide as a fission fragment or as a remaining resdiual product of a
! high-energy calculation.
!
! talysinput  : subroutine for user input and defaults
! talysinitial: subroutine for initialization of nuclear structure
!
  call talysinput
  call talysinitial
!
! basicxs     : subroutine for basic cross sections and transmission coefficients
! gamma       : subroutine for gamma cross section and transmission coefficients
! preeqinit   : subroutine for initialization of general pre-equilibrium parameters
! excitoninit : subroutine for initialization of exciton model parameters
! compoundinit: subroutine for initialization of compound model parameters
!
  call basicxs(0, 0)
  if (parinclude(0)) call gamma(0, 0)
  if (enincmax >= epreeq) call preeqinit
  if (flagcomp) call compoundinit
!
! energies   : subroutine for energies
! reacinitial: subroutine for initialization of arrays for various cross sections
!
  call energies
  call reacinitial
!
! Optical model
!
! incident: subroutine for main settings and basic cross sections for incident energy
! exgrid  : subroutine to set excitation energy grid
!
  call incident
  call exgrid(0, 0)
!
! Pre-equilibrium reactions
!
! preeq     : subroutine for preequilibrium reactions
! population: subroutine for processing of pre-equilibrium spectra into population bins
!
  if (flagpreeq) then
    call preeq
    call population
  endif
!
! Multiple emission
!
! multiple: subroutine for multiple emission
!
  call multiple
!
! Exclusive channels
!
! channels    : subroutine for exclusive reaction channels
!
  if (flagchannels) call channels
!
! Collecting total cross sections, spectra, angular distributions, etc.
!
! totalxs : subroutine for total cross sections
! spectra : subroutine for creation of particle spectra
! residual: subroutine for residual production cross sections
! output  : subroutine for output
!
  call totalxs
  if (flagspec .or. flagddx) call spectra
  call residual
  call output
  return
end subroutine evaptalys
! Copyright A.J. Koning 2021
