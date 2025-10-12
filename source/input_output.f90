subroutine input_output
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for output variables
!
! Author    : Arjan Koning
!
! 2025-10-08: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
  use A1_error_handling_mod
!
! Definition of single and double precision variables
!   sgl               ! single precision kind
! Variables for output
!   flagang         ! flag for output of angular distributions
!   flagbinspec     ! flag for output of emission spectrum per excitation bin
!   flagblock       ! flag to block spectra, angle and gamma files
!   flagblockddx    ! flag to block DDX files
!   flagblockspectra! flag to block spectra files
!   flagblockangle  ! flag to block angle files
!   flagblocklevels ! flag to block discrete level files
!   flagblockomp    ! flag to block optical model files
!   flagblockpreeq  ! flag to block preequilibrium files
!   flagblockbin    ! flag to block binary files
!   flagblockdirect ! flag to block direct reaction files
!   flagblockastro  ! flag to block astro reaction rate files
!   flagblockyield  ! flag to block isotopic yield files
!   flagcheck       ! flag for output of numerical checks
!   flagcompo       ! flag for output of cross section components
!   flagddx         ! flag for output of double-differential cross sections
!   flagdecay       ! flag for output of decay of each population bin
!   flagexc         ! flag for output of excitation functions
!   flaggamdis      ! flag for output of discrete gamma-ray intensities
!   flagintegral    ! flag for calc. of effective cross section using integral data
!   flaginverse     ! flag for output of transmission coeff. and inverse cross sections
!   flagmain        ! flag for main output
!   flagpop         ! flag for output of population
!   flagsacs        ! flag for statistical analysis of cross sections
!   flagspec        ! flag for output of spectra
!   fluxname        ! name of integral spectrum
!   xsfluxfile      ! TALYS cross section file for integral data
!   ddxmode         ! mode for DDX: 0: None, 1: Ang. distr., 2: Spect. per angle, 3 both
!   Nflux           ! number of reactions with integral data
!   integralexp     ! experimental effective cross section
! All global variables
!   numflux    ! number of integral experiments
! Variables for reading input lines
!   inline                ! input line
!   nlines                ! number of input lines
! Variables for basic reaction
!   flagastro       ! flag for calculation of astrophysics reaction rate
!   flagbasic       ! flag for output of basic information and results
!   flagchannels    ! flag for exclusive channel calculation
!   flagendf        ! flag for information for ENDF - 6 file
!   flagendfdet     ! flag for detailed ENDF - 6 information per channel
!   flagrecoil      ! flag for calculation of recoils
! Variables for fission
!   flagfission     ! flag for fission
! Variables for direct reactions
!   flagdisc        ! flag for output of discrete state cross sections
! Variables for input energies
!   Ninc       ! number of incident energies
! Variables for main input
!   flagnatural    ! flag for calculation of natural element
!   k0             ! index of incident particle
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=1)   :: ch         ! character
  character(len=132) :: key        ! keyword
  character(len=132) :: value      ! value or string
  character(len=132) :: word(40)   ! words on input line
  character(len=132) :: line       ! input line
  integer            :: i          ! counter
  integer            :: istat      ! logical for file access
!
! ************** Defaults *************
!
  flagbinspec = .false.
  flagblock = .false.
  flagblockddx = .false.
  flagblockspectra = .false.
  flagblockangle = .false.
  flagblocklevels = .false.
  flagblockomp = .true.
  flagblockpreeq = .false.
  flagblockbin = .false.
  flagblockdirect = .false.
  flagblockastro = .false.
  flagblockyield = .false.
  flagblockZA = .true.
  flagcheck = flagbasic
  flagang = .false.
  flagcompo = .false.
  flagddx = .false.
  flagdecay = .false.
  if (Ninc == 1) then
    flagexc = .false.
  else
    flagexc = .true.
  endif
  if (flagnatural) flagexc = .true.
  flaggamdis = .false.
  flagintegral = .false.
  flaginverse = flagbasic
  flagmain = .true.
  flagpop = flagbasic
  flagsacs = .false.
  flagspec = .false.
  if (flagrecoil) flagspec = .true.
  fluxname = ' '
  xsfluxfile = ' '
  ddxmode = 0
  Nflux = 0
  integralexp = 0.
!
! If the results of TALYS are used to create ENDF-6 data files, several output flags are automatically set.
!
  if (flagendf) then
    flagang = .true.
    flagcheck = .true.
    flagblock = .true.
    flagblockddx = .true.
    flagblockspectra = .true.
    flagblockangle = .true.
    if (k0 == 3) then
      ddxmode = 2
      flagddx = .true.
    endif
    flagspec = .true.
    flagexc = .true.
    if (flagendfdet) then
      flaggamdis = .true.
      flagchannels = .true.
    endif
  endif
!
! **************** Read input variables *******************
!
! getkeywords: subroutine to retrieve keywords and values from input line
!
! The keyword is identified and the corresponding values are read.
! Erroneous input is immediately checked.
! The keywords and number of values on each line are retrieved from the input.
!
  do i = 1, nlines
    line = inline(i)
    call getkeywords(line, word)
    key = word(1)
    value = word(2)
    ch = word(2)(1:1)
!
! Test for keywords
!
! Here, the various model parameters can be set to overrule the default values.
! Most default values will be computed later on, since they require more computation (e.g. level density parameters).
!
! Each keyword is characterized by a certain order of parameter and value input.
! They are distinguished by different classes.
!
! getvalues : subroutine to assign values to keywords
!
    if (key == 'outmain') then
      if (ch == 'n') flagmain = .false.
      if (ch == 'y') flagmain = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'block') then
      if (ch == 'n') flagblock = .false.
      if (ch == 'y') flagblock = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'blockddx') then
      if (ch == 'n') flagblockddx = .false.
      if (ch == 'y') flagblockddx = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'blockspectra') then
      if (ch == 'n') flagblockspectra = .false.
      if (ch == 'y') flagblockspectra = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'blockangle') then
      if (ch == 'n') flagblockangle = .false.
      if (ch == 'y') flagblockangle = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'blockdirect') then
      if (ch == 'n') flagblockdirect = .false.
      if (ch == 'y') flagblockdirect = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'blocklevels') then
      if (ch == 'n') flagblocklevels = .false.
      if (ch == 'y') flagblocklevels = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'blockomp') then
      if (ch == 'n') flagblockomp = .false.
      if (ch == 'y') flagblockomp = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'blockpreeq') then
      if (ch == 'n') flagblockpreeq = .false.
      if (ch == 'y') flagblockpreeq = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'blockbin') then
      if (ch == 'n') flagblockbin = .false.
      if (ch == 'y') flagblockbin = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'blockastro') then
      if (ch == 'n') flagblockastro = .false.
      if (ch == 'y') flagblockastro = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'blockyield') then
      if (ch == 'n') flagblockyield = .false.
      if (ch == 'y') flagblockyield = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'blockza') then
      if (ch == 'n') flagblockZA = .false.
      if (ch == 'y') flagblockZA = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outpopulation') then
      if (ch == 'n') flagpop = .false.
      if (ch == 'y') flagpop = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outcheck') then
      if (ch == 'n') flagcheck = .false.
      if (ch == 'y') flagcheck = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outinverse') then
      if (ch == 'n') flaginverse = .false.
      if (ch == 'y') flaginverse = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outdecay') then
      if (ch == 'n') flagdecay = .false.
      if (ch == 'y') flagdecay = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outspectra') then
      if (ch == 'n') flagspec = .false.
      if (ch == 'y') flagspec = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'ddxmode') then
      read(value, * , iostat = istat) ddxmode
      if (istat /= 0) call read_error(line, istat)
      if (ddxmode == 0) flagddx = .false.
      if (ddxmode > 0) then
        flagddx = .true.
        flagspec = .true.
      endif
      cycle
    endif
    if (key == 'outgamdis') then
      if (ch == 'n') flaggamdis = .false.
      if (ch == 'y') flaggamdis = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outexcitation') then
      if (ch == 'n') flagexc = .false.
      if (ch == 'y') flagexc = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outbinspectra') then
      if (ch == 'n') flagbinspec = .false.
      if (ch == 'y') flagbinspec = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outangle') then
      if (ch == 'n') flagang = .false.
      if (ch == 'y') flagang = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'components') then
      if (ch == 'n') flagcompo = .false.
      if (ch == 'y') flagcompo = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'integral') then
      if (ch == 'n') flagintegral = .false.
      if (ch == 'y') flagintegral = .true.
      if (ch /= 'y' .and. ch /= 'n') then
        call range_integer_error(key, k0, 0, 1)
        Nflux = Nflux + 1
        call range_integer_error(key, Nflux, 1, numflux)
        xsfluxfile(Nflux) = value
        fluxname(Nflux) = word(3)
        flagintegral = .true.
        read(word(4), * , iostat = istat) integralexp(Nflux)
        if (istat /= 0) cycle
      endif
      cycle
    endif
    if (key == 'sacs') then
      if (ch == 'n') flagsacs = .false.
      if (ch == 'y') flagsacs = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
  enddo
  return
end subroutine input_output
! Copyright A.J. Koning 2021
