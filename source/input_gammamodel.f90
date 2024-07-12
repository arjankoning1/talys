subroutine input_gammamodel
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for gamma models
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
! Definition of single and double precision variables
!   sgl               ! single precision kind
! Variables for gamma rays
!   Exlfile           ! tabulated gamma strength function
!   fiso              ! correction factor for isospin forbidden transitions
!   flagracap         ! flag for radiative capture model
!   flagupbend        ! flag for low-energy upbend of photon strength function
!   flagpsfglobal     ! flag for global photon strength functions only
!   flaggnorm         ! flag to normalize PSF to average radiative width
!   ldmodelracap      ! level density model for direct radiative capture
!   strength          ! E1 strength function model
!   strengthM1        ! model for M1 gamma - ray strength function
! Variables for reading input lines
!   inline                ! input line
!   nlines                ! number of input lines
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical            :: flagassign           ! flag to assign value or not
  character(len=1)   :: ch                   ! character
  character(len=132) :: cval                 ! character value
  character(len=132) :: key                  ! keyword
  character(len=132) :: value                ! value or string
  character(len=132) :: word(40)             ! words on input line
  character(len=132) :: line                 ! input line
  integer            :: class                ! input class
  integer            :: i                    ! counter
  integer            :: ibar                 ! fission barrier
  integer            :: igr                  ! giant resonance
  integer            :: irad                 ! variable to indicate M(=0) or E(=1) radiation
  integer            :: istat                ! logical for file access
  integer            :: ival                 ! integer value
  integer            :: lval                 ! multipolarity
  integer            :: Nix                  ! neutron number index for residual nucleus
  integer            :: type                 ! particle type
  integer            :: Zix                  ! charge number index for residual nucleus
  real(sgl)          :: val                  ! real value
!
! ************** Defaults *************
!
  Exlfile = ' '
  flagracap = .false.
  filepsf = .false.
  flaggamma = flagbasic
  if (k0 >= 1) then
    flagupbend = .true.
  else
    flagupbend = .false.
  endif
  flagpsfglobal = .false.
  flaggnorm = .false.
  strengthM1 = 3
  if (strength == 8) strengthM1 = 8
  if (strength == 10) strengthM1 = 10
  if (strength <= 2) strengthM1 = 2
  ldmodelracap = 3
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
    Zix = 0
    Nix = 0
    type = 0
    lval = 0
    ibar = 0
    igr = 1
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
    if (key == 'e1file') then
      class = 11
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Exlfile(Zix, Nix, 1, 1) = cval
      cycle
    endif
    if (key == 'outgamma') then
      if (ch == 'n') flaggamma = .false.
      if (ch == 'y') flaggamma = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'filepsf') then
      if (ch == 'n') filepsf = .false.
      if (ch == 'y') filepsf = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'ldmodelracap') then
      read(value, * , iostat = istat) ldmodelracap
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'racap') then
      if (ch == 'n') flagracap = .false.
      if (ch == 'y') flagracap = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'strengthm1') then
      read(value, * , iostat = istat) strengthM1
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'upbend') then
      if (ch == 'n') flagupbend = .false.
      if (ch == 'y') flagupbend = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'psfglobal') then
      if (ch == 'n') flagpsfglobal = .false.
      if (ch == 'y') flagpsfglobal = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'gnorm') then
      if (ch == 'n') flaggnorm = .false.
      if (ch == 'y') flaggnorm = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'm1file') then
      class = 11
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Exlfile(Zix, Nix, 0, 1) = cval
      cycle
    endif
  enddo
  return
end subroutine input_gammamodel
! Copyright A.J. Koning 2021
