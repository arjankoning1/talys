subroutine incident
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Main settings and basic cross sections for incident energy
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
!   flaginverse    ! flag for output of transmission coefficients and inverse cross sections
!   flagmain       ! flag for main output
! Variables for gamma rays
!   strength       ! E1 strength function model
! Variables for compound reactions
!   flagcomp    ! flag for compound angular distribution calculation
!   flagurr     ! flag for output of unresolved resonance parameters
! Variables for input energies
!   flaginitpop    ! flag for initial population distribution
! Variables for main input
!   k0             ! index of incident particle
! Variables for energy grid
!   coullimit      ! energy limit for charged particle OMP calculation
!   Einc           ! incident energy in MeV
! Variables for energies
!   flagmulpre     ! flag for multiple pre - equilibrium calculation
!   flagpreeq      ! flag for pre - equilibrium calculation
!   flagwidth      ! flag for width fluctuation calculation
!   nbins          ! number of continuum excitation energy bins
! Variables for nuclides
!   parinclude     ! logical to include outgoing particle
!
! *** Declaration of local data
!
  implicit none
  character(len=1) :: yesno     ! function to assign y or n to logical value
!
! ************ Setting of some parameters for incident channel *********
!
! yesno     : function to assign y or n to logical value
!
! Multiple pre-equilibrium emission is always set off if there is no primary pre-equilibrium emission.
!

  if (.not. flagpreeq) flagmulpre = .false.
!
! Write the energy dependent flags to the output file.
!
  if (flagmain) then
    if (Einc >= 0.001) then
      write(*, '(/, " ########## RESULTS FOR E=", f10.5, " ##########"/)') Einc
    else
      write(*, '(/, " ########## RESULTS FOR E=", es12.5, " ##########"/)') Einc
    endif
    write(*, '(" Energy dependent input flags"/)')
    write(*, '(" Width fluctuations (flagwidth)            : ", a1)') yesno(flagwidth)
    write(*, '(" Unresolved resonance parameters (flagurr) : ", a1)') yesno(flagurr)
    write(*, '(" Preequilibrium (flagpreeq)                : ", a1)') yesno(flagpreeq)
    write(*, '(" Multiple preequilibrium (flagmulpre)      : ", a1)') yesno(flagmulpre)
    write(*, '(" Number of continuum excitation energy bins:", i3)') nbins
  endif
!
! *** Calculation of total, reaction, elastic cross section, transmission coefficients and elastic angular distribution for
!   incident energy ***
!
! incidentecis : subroutine for ECIS calculation for incident energy
! incidentread : subroutine to read ECIS results for incident energy
! incidentnorm : subroutine for normalization of reaction cross sections and transmission coefficients for incident channel
! incidentgamma: subroutine for incident photons
! strength     : model for E1 gamma-ray strength function
! radwidtheory : subroutine for theoretical calculation of total radiative width
! tgamma       : subroutine for photon transmission coefficients
! spr          : subroutine for S, P and R' resonance parameters inverse reaction cross sections
! incidentout  : subroutine for reaction output for incident channel
!
  if (k0 > 1 .and. Einc < coullimit(k0)) return
  if (k0 > 0) then
    call incidentecis
    call incidentread
    call incidentnorm
  else
    if ( .not. flaginitpop) call incidentgamma
  endif
  if ((parinclude(0) .or. flagcomp) .and. Einc <= 100.) then
    call radwidtheory(0, 0, Einc)
    if (strength == 1) call tgamma(0, 0)
  endif
  if (k0 == 1 .and. (parinclude(0) .or. flagcomp)) call spr
  if (flaginverse) call incidentout
  return
end subroutine incident
! Copyright A.J. Koning 2021
