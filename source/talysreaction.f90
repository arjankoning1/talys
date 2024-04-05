subroutine talysreaction
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Reaction models
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
! Variables for numerics
!   xseps           ! limit for cross sections
! Variables for main input
!   k0            ! index of incident particle
! Variables for output
!   flagddx         ! flag for output of double - differential cross sections
!   flagintegral    ! flag for calculation of effective cross section using integral spectrum
!   flagmain        ! flag for main output
!   flagsacs        ! flag for statistical analysis of cross sections
!   flagspec        ! flag for output of spectra
! Variables for basic reaction
!   flagang         ! flag for output of angular distributions
!   flagastro       ! flag for calculation of astrophysics reaction rate
!   flagchannels    ! flag for exclusive channels calculation
!   flagendf        ! flag for information for ENDF - 6 file
!   flagmassdis     ! flag for calculation of fission fragment mass yields
!   flagreaction    ! flag for calculation of nuclear reactions
!   flagrecoil      ! flag for calculation of recoils
!   flagrpevap      ! flag for evaporation of residual products at high incident energies
! Variables for best files
!   flagrescue      ! flag for final rescue: normalization to data
! Variables for basic parameters
!   eninclow        ! minimal incident energy for nuclear model calculations
! Variables for input energies
!   eninc           ! incident energy in MeV
!   enincmax        ! maximum incident energy
!   flaginitpop     ! flag for initial population distribution
!   nin             ! counter for incident energy
!   Ninc          ! number of incident energies
! Variables for compound reactions
!   flagcomp        ! flag for compound angular distribution calculation
!   flagres         ! flag for output of low energy resonance cross sections
!   flagurr         ! flag for output of unresolved resonance parameters
! Variables for preequilibrium
!   breakupmodel    ! model for break-up reaction: 1. Kalbach 2. Avrigeanu
!   epreeq          ! on - set incident energy for preequilibrium calculation
! Variables for medical isotope production
!   flagprod        ! flag for isotope production
! Variables for fission
!   flagfission     ! flag for fission
!   fymodel         ! fission yield model, 1: Brosa 2: GEF
! Variables for gamma rays
!   flagracap       ! flag for radiative capture model
! Variables for OMP
!   flagompall      ! flag for new optical model calculation for all residual
! Variables for energy grid
!   Einc            ! incident energy in MeV
! Variables for energies
!   flagpreeq       ! flag for pre - equilibrium calculation
!   Ninclow       ! number of incident energies below Elow
! Variables for incident channel
!   lmaxinc         ! maximal l - value for transm. coeff. for incident channel
!   xsreacinc       ! reaction cross section for incident channel
! Variables for nuclides
!   parinclude      ! logical to include outgoing particle
!
! *** Declaration of local data
!
  implicit none
!
! ************************* Reaction calculation ***********************
!
! Basic cross sections (optical model) and initialisation
!
! basicxs     : subroutine for basic cross sections and transmission coefficients
! gamma       : subroutine for gamma cross section and transmission coefficients
! preeqinit   : subroutine for initialization of general pre-equilibrium parameters
! excitoninit : subroutine for initialization of exciton model parameters
! racapinit   : subroutine for initialization of radiative capture model
! compoundinit: subroutine for initialization of compound model parameters
! astroinit   : subroutine for initialization of astrophysics quantities
!
  if (flagreaction) then
    if (.not. flagompall) call basicxs(0, 0)
    if (parinclude(0)) call gamma(0, 0)
    if (enincmax >= epreeq .or. flagracap) call preeqinit
    if (flagracap) call racapinit
    if (flagcomp) call compoundinit
    if (flagastro) call astroinit
!
! Loop over incident energies
!
! Initialisation
!
! energies   : subroutine for energies
! reacinitial: subroutine for initialization of arrays for various cross sections
!
    do nin = 1, Ninc
      if (flagompall) call basicxs(0, 0)
      Einc = eninc(nin)
      call energies
      call reacinitial
      if (Einc < eninclow) cycle
!
! Optical model
!
! incident   : subroutine for main settings and basic cross sections for incident energy
! exgrid     : subroutine to set excitation energy grid
! recoilinit : subroutine for calculation of initial recoil velocity and direction
!
      call incident
      call exgrid(0, 0)
      if (flagrecoil) call recoilinit
!
! In certain cases, there will be no nuclear reaction calculation
!
      if ((lmaxinc /=  - 1 .and. xsreacinc > xseps) .or. flaginitpop) then
!
! Direct reactions
!
! direct   : subroutine for calculation of direct inelastic cross sections
! racap    : subroutine for radiative capture model
!
        call direct
        if (flagracap) call racap
!
! Pre-equilibrium reactions
!
! preeqspindis: subroutine for preequilibrium spin distribution
! preeq       : subroutine for preequilibrium reactions
! population  : subroutine for processing of pre-equilibrium spectra into population bins
!
        if (flagpreeq) then
          call preeqspindis
          call preeq
          call population
        endif
!
! Binary compound reactions
!
! compnorm  : subroutine for normalization of compound nucleus cross section
! comptarget: subroutine for compound reaction for initial compound nucleus
!
        if (flagcomp) then
          call compnorm
          call comptarget
        endif
!
! Collecting binary reaction results
!
! binary : subroutine for binary reaction results
! angdis : subroutine for calculation of angular distributions for discrete states
!
        call binary
        if (flagang .or. flagddx .or. flagrecoil) call angdis
!
! Multiple emission
!
! multiple:preeesubroutine for multiple emission
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
! totalxs      : subroutine for total cross sections
! spectra      : subroutine for creation of particle spectra
! massdis      : subroutine for fission fragment yields
! massdisout   : subroutine for output of fission fragment yields
! nubarout     : subroutine for output of number of fission neutrons and spectra
! residual     : subroutine for residual production cross sections
! totalrecoil  : subroutine for recoil results
! normalization: subroutine to normalize cross sections to experimental or evaluated data
! thermal      : subroutine for estimate of thermal cross sections
! urr          : subroutine for unresolved resonance range parameters
! output       : subroutine for output
!
        call totalxs
        if (flagspec .or. flagddx) call spectra
        if (flagfission .and. flagmassdis) then
          call massdis
          if (fymodel <= 2) call massdisout
          if (fymodel == 2) then
            call nubarout
            call nudisout
          endif
        endif
        if (breakupmodel == 2 .and. k0 == 3) call residualBU
        call residual
        if (flagrecoil) call totalrecoil
        if (flagrescue) call normalization
        if (nin == Ninclow + 1 .and. Ninclow > 0) call thermal
      endif
      if (flagurr) call urr
      if ( .not. flagastro) call output
      if (flagfission .and. flagmassdis .and. fymodel >= 3) call ffevap
      if (flagrpevap) call rpevap
    enddo
!
! Final output
!
! astro       : subroutine for astrophysical reaction rates
! finalout    : subroutine for output of final results integral spectrum
! integral    : subroutine to calculate effective cross section for integral spectrum
! endf        : subroutine for cross sections and information for ENDF-6 file
! isoprod     : subroutine for isotope production
! timer       : subroutine for output of execution time
!
    if (flagastro) then
      call astro
    else
      call finalout
    endif
    if (flagres) call resonance
    if (flagintegral) call integral
    if (flagsacs) call sacs
    if (flagendf) call endf
  endif
  if (flagprod) call isoprod
  if (flagmain) call timer
  return
end subroutine talysreaction
! Copyright A.J. Koning 2021
