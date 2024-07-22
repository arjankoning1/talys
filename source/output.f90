subroutine output
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output
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
!   flaggamdis      ! flag for output of discrete gamma - ray intensities
!   flagmain        ! flag for main output
!   flagspec        ! flag for output of spectra
! Variables for basic reaction
!   flagffruns      ! flag to designate subsequent evaporation of fission products
!   flagrpruns      ! flag to designate that run is for residual product
! Variables for basic reaction
!   flagang         ! flag for output of angular distributions
!   flagchannels    ! flag for exclusive channels calculation
!   flagrecoil      ! flag for calculation of recoils
! Variables for direct reactions
!   flagdisc        ! flag for output of discrete state cross sections
! Variables for fission
!   flagfission     ! flag for fission
! Variables for gamma rays
!   flagracap       ! flag for radiative capture model
! Variables for input energies
!   flaginitpop     ! flag for initial population distribution
! Variables for main input
!   Atarget         ! mass number of target nucleus
!
! ******************************* Output *******************************
!
! totalout     : subroutine for output of total cross sections
! binaryout    : subroutine for output of binary cross sections
! productionout: subroutine for output of particle production cross sections
! residualout  : subroutine for output of residual production cross sections
! fissionout   : subroutine for output of fission cross sections
! discreteout  : subroutine for output of cross sections for discrete states
! channelsout  : subroutine for output of exclusive reaction channels
! spectraout   : subroutine for output of particle spectra
! recoilout    : subroutine for output of recoils
! angleout     : subroutine for output of discrete angular distributions
! ddxout       : subroutine for output of double-differential cross sections
! gamdisout    : subroutine for output of discrete gamma-ray intensities
! racapout     : subroutine for output of radiative capture model
!
  if (flagmain) then
    call totalout
    if ( .not. flaginitpop) call binaryout
    call productionout
    call residualout
    if (flagfission) call fissionout
  endif
  if (flagdisc) call discreteout
  if (flagchannels) call channelsout
  if (flagspec) call spectraout
  if (flagrecoil) call recoilout
  if (flagang .or. fileelastic) call angleout
  if (flagddx) call ddxout
  if (flaggamdis) call gamdisout
  if (flagracap) call racapout
  if (flagffruns .or. flagrpruns) write(*, '(/" End of calculation for ", a/)') trim(targetnuclide)
  return
end subroutine output
! Copyright A.J. Koning 2021
