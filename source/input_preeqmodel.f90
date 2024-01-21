subroutine input_preeqmodel
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for preequilibrium models
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
!
  use A0_talys_mod
  use A1_error_handling_mod
!
! Variables for preequilibrium
!   flag2comp    ! flag for two-component pre-equilibrium model
!   flagecisdwba ! flag for new ECIS calculation for DWBA for MSD
!   flagoutdwba  ! flag for output of DWBA cross sections for MSD
!   flaggshell   ! flag for energy dependence of p.h. level density parameter
!   flagonestep  ! flag for continuum one-step direct only
!   flagpecomp   ! flag for Kalbach complex particle emission model
!   flagpeout    ! flag for output of pre-equilibrium results
!   flagsurface  ! flag for surface effects in exciton model
!   preeqadjust  ! logical for energy-dependent pre-eq adjustment
!   breakupmodel ! model for break-up reaction: 1. Kalbach 2. Avrigeanu
!   mpreeqmode   ! designator for multiple pre-equilibrium model
!   pairmodel    ! model for preequilibrium pairing energy
!   pespinmodel  ! model for pre-equilibrium or compound spin
!   phmodel      ! particle-hole state density model
!   preeqmode    ! designator for pre-equilibrium model
!   emulpre      ! on-set incident energy for multiple preequilibrium
!   epreeq       ! on-set incident energy for preequilibrium
! Constants
!   Emaxtalys           ! maximum acceptable energy for TALYS
! Variables for reading input lines
!   inline         ! input line
!   nlines            ! number of input lines
! Error handling
!   read_error ! Message for file reading error
! Variables for main input
!   k0             ! index of incident particle
!   ptype0            ! type of incident particle
!
! *** Declaration of local data
!
  implicit none
  character(len=1)   :: ch                   ! character
  character(len=132) :: key                  ! keyword
  character(len=132) :: value                ! value or string
  character(len=132) :: word(40)             ! words on input line
  character(len=132) :: line                 ! input line
  integer            :: i                    ! counter
  integer            :: istat                ! logical for file access
!
! ************* Defaults ************************
!
  flag2comp = .true.
  flagecisdwba = .true.
  flagoutdwba = .false.
  flaggshell = .false.
  flagonestep = .false.
  flagpeout = .false.
  if (k0 >= 1) then
    flagpecomp = .true.
    flagsurface = .true.
  else
    flagpecomp = .false.
    flagsurface = .false.
  endif
  preeqadjust = .false.
  breakupmodel = 1
  mpreeqmode = 2
  pairmodel = 2
  if (k0 <= 1) then
    pespinmodel = 1
  else
    pespinmodel = 2
  endif
  phmodel = 1
  preeqmode = 2
  emulpre = 20.
  if (ptype0 == '0') then
    epreeq = Emaxtalys
  else
    epreeq = - 1.
  endif
  if (flagomponly) then
    epreeq = Emaxtalys
    emulpre = Emaxtalys
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
    if (key == 'preequilibrium') then
      if (ch == 'y') then
        epreeq = 0.
        cycle
      endif
      if (ch == 'n') then
        epreeq = Emaxtalys
        cycle
      endif
      read(value, * , iostat = istat) epreeq
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'multipreeq') then
      if (ch == 'y') then
        emulpre = 0.
        cycle
      endif
      if (ch == 'n') then
        emulpre = Emaxtalys
        cycle
      endif
      read(value, * , iostat = istat) emulpre
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'preeqspin') then
      if (ch == 'n') then
        pespinmodel = 1
        cycle
      endif
      if (ch == 'y') then
        pespinmodel = 3
        cycle
      endif
      read(value, * , iostat = istat) pespinmodel
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'preeqsurface') then
      if (ch == 'n') flagsurface = .false.
      if (ch == 'y') flagsurface = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'preeqcomplex') then
      if (ch == 'n') flagpecomp = .false.
      if (ch == 'y') flagpecomp = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'twocomponent') then
      if (ch == 'n') flag2comp = .false.
      if (ch == 'y') flag2comp = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'onestep') then
      if (ch == 'n') flagonestep = .false.
      if (ch == 'y') flagonestep = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'gshell') then
      if (ch == 'n') flaggshell = .false.
      if (ch == 'y') flaggshell = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'preeqmode') then
      read(value, * , iostat = istat) preeqmode
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'mpreeqmode') then
      read(value, * , iostat = istat) mpreeqmode
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'phmodel') then
      read(value, * , iostat = istat) phmodel
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'ecisdwba') then
      if (ch == 'n') flagecisdwba = .false.
      if (ch == 'y') flagecisdwba = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outdwba') then
      if (ch == 'n') flagoutdwba = .false.
      if (ch == 'y') flagoutdwba = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outpreequilibrium') then
      if (ch == 'n') flagpeout = .false.
      if (ch == 'y') flagpeout = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'pairmodel') then
      read(value, * , iostat = istat) pairmodel
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'breakupmodel') then
      read(value, * , iostat = istat) breakupmodel
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
  enddo
  return
end subroutine input_preeqmodel
! Copyright A.J. Koning 2021
