subroutine input_directpar
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for direct reaction parameters
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
  use A1_error_handling_mod
!
! Variables for direct reactions
!   eadd        ! on-set incident energy for addition of discrete states
!   eaddel      ! on-set incident energy for addition of elastic peak
!   elwidth     ! width of elastic peak in MeV
!   maxband     ! highest vibrational band added to rotational model
!   maxrot      ! number of included excited rotational levels
!   core        ! even-even core for weakcoupling (-1 or 1)
!   soswitch    ! switch for deformed spin-orbit calculation
! Constants
!   Emaxtalys       ! maximum acceptable energy for TALYS
!   parsym          ! symbol of particle
! All global variables
!   numpar      ! number of particles
! Variables for main input
!   k0           ! index of incident particle
! Variables for basic reaction
!   flagendf     ! flag for information for ENDF - 6 file
! Variables for reading input lines
!   inline            ! input line
!   nlines            ! number of input lines
! Error handling
!   read_error ! Message for file reading error
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
! ************** Defaults *************
!
  maxband = 0
  maxrot = 4
  core = -1
  eadd = 0.
  eaddel = 0.
  if (flagendf) then
    if (k0 == 1) then
      eadd = 30.
      eaddel = Emaxtalys
    endif
  endif
  elwidth = 0.5
  soswitch = 3.
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
    if (key == 'maxband') then
      read(value, * , iostat = istat) maxband
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'maxrot')  then
      read(value, * , iostat = istat) maxrot
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'core') then
      read(value, * , iostat = istat) core
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'adddiscrete') then
      if (ch == 'y') then
        eadd = 0.
        cycle
      endif
      if (ch == 'n') then
        eadd = Emaxtalys
        cycle
      endif
      read(value, * , iostat = istat) eadd
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'addelastic') then
      if (ch == 'y') then
        eaddel = 0.
        cycle
      endif
      if (ch == 'n') then
        eaddel = Emaxtalys
        cycle
      endif
      read(value, * , iostat = istat) eaddel
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'elwidth')  then
      read(value, * , iostat = istat) elwidth
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'soswitch') then
      read(value, * , iostat = istat) soswitch
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
  enddo
  return
end subroutine input_directpar
! Copyright A.J. Koning 2021
