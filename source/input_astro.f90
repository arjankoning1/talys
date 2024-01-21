subroutine input_astro
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for astro variables
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
! Variables for input of astrophysics
!   astroE          ! energy, in MeV, for Maxwellian average
!   astroT9         ! temperature, in 10^9 K, for Maxwellian average
!   flagastroex     ! flag for calculation of astrophysics reaction rate to f
!   flagastrogs     ! flag for calculation of astrophysics reaction rate with
!   nonthermlev     ! non - thermalized level in the calculation of astrophysic
!   nTmax           ! effective number of temperatures for Maxwellian average
! Variables for basic reactions
!   flagastro        ! flag for calculation of astrophysics reaction rate
! All global variables
!   numT              ! number of temperatures
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
  astroE = 0.
  astroT9 = 0.
  flagastroex = .false.
  flagastrogs = .false.
  nonthermlev = -1
  nTmax = numT
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
    if (key == 'astroe')  then
      read(value, * , iostat = istat) astroE
      if (istat /= 0) call read_error(line, istat)
      nTmax = 1
      flagastro = .true.
      flagastrogs = .true.
      cycle
    endif
    if (key == 'astrot')  then
      read(value, * , iostat = istat) astroT9
      if (istat /= 0) call read_error(line, istat)
      nTmax = 1
      flagastro = .true.
      cycle
    endif
    if (key == 'astrogs') then
      if (ch == 'n') flagastrogs = .false.
      if (ch == 'y') flagastrogs = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'astroex') then
      if (ch == 'n') flagastroex = .false.
      if (ch == 'y') flagastroex = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'nonthermlev') then
      read(value, * , iostat = istat) nonthermlev
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
  enddo
  return
end subroutine input_astro
! Copyright A.J. Koning 2021
