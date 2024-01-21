subroutine input_mass
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for mass variables
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
! Variables for input of masses
!   beta2          ! deformation parameter
!   flagexpmass    ! flag for using experimental nuclear mass if available
!   massdir        ! directory with mass tables
!   massexcess     ! mass excess in MeV as read from user input file
!   massmodel      ! model for theoretical nuclear mass
!   massnucleus    ! mass of nucleus in amu as read from user input file
! All global variables
!   numbar            ! number of fission barriers
!   numN              ! maximum number of neutrons from initial compound nucleus
!   numZ              ! maximum number of protons from initial compound nucleus
! Variables for reading input lines
!   inline            ! input line
!   nlines            ! number of input lines
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
  beta2 = 0.
  do Zix = 0, numZ + 4
    do Nix = 0, numN + 4
      beta2(Zix, Nix, 1) = 0.6
      beta2(Zix, Nix, 2) = 0.8
      do ibar = 3, numbar
        beta2(Zix, Nix, ibar) = 1.
      enddo
    enddo
  enddo
  flagexpmass = .true.
  massdir = ' '
  massexcess = 0.
  massmodel = 2
  massnucleus = 0.
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
    if (key == 'beta2') then
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) beta2(Zix, Nix, ibar) = val
      cycle
    endif
    if (key == 'expmass') then
      if (ch == 'n') flagexpmass = .false.
      if (ch == 'y') flagexpmass = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'massdir') then
      read(value, * , iostat = istat) massdir
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'massexcess') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) massexcess(Zix, Nix) = val
      cycle
    endif
    if (key == 'massmodel') then
      read(value, * , iostat = istat) massmodel
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'massnucleus') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) massnucleus(Zix, Nix) = val
      cycle
    endif
  enddo
  return
end subroutine input_mass
! Copyright A.J. Koning 2021
