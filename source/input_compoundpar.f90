subroutine input_compoundpar
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for compound reaction parameters
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
! Variables for compound reactions
!   adjustTJ        ! logical for energy-dependent TJ adjustment
!   reslib          ! library with resonance parameters
!   lenreslib       ! length of library name with resonance parameters
!   flagurr         ! flag for output of unresolved resonance parameters
!   eurr            ! off-set incident energy for URR calculation
!   lurr            ! maximal orbital angular momentum for URR
!   TJadjust        ! adjustable factor for TJ (default 1.)
!   Tres            ! temperature for broadening low energy cross sections
!   xsalphatherm    ! thermal (n,a) cross section
!   xscaptherm      ! thermal capture cross section
!   xsptherm        ! thermal (n,p) cross section
! Variables for reading input lines
!   inline            ! input line
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
  reslib = 'tendl.2025'
  lenreslib = 7
  if (k0 /= 1 .or. .not. flagcomp) then
    eurr = 0.
  else
    eurr = -1.
  endif
  flagurr = .false.
  if (flagendf .and. k0 == 1 .and. Atarget > 20) flagurr = .true.
  lurr = 2
  Tres = 293.16
  xsalphatherm = 0.
  xscaptherm = 0.
  xsptherm = 0.
  adjustTJ = .false.
  TJadjust = 1.
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
    if (key == 'reslib') then
      reslib = value
      cycle
    endif
    if (key == 'urr') then
      flagurr = .true.
      if (ch == 'y') cycle
      if (ch == 'n') then
        eurr = 0
        flagurr = .false.
        cycle
      endif
      read(value, * , iostat = istat) eurr
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'lurr') then
      read(value, * , iostat = istat) lurr
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'tres') then
      read(value, * , iostat = istat) Tres
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'xscaptherm')  then
      read(value, * , iostat = istat) xscaptherm(-1)
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'xsptherm')  then
      read(value, * , iostat = istat) xsptherm(-1)
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'xsalphatherm')  then
      read(value, * , iostat = istat) xsalphatherm(-1)
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'tjadjust') then
      class = 12
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        TJadjust(Zix, Nix, type) = val
        adjustTJ(Zix, Nix, type) = .true.
      endif
      cycle
    endif
  enddo
  return
end subroutine input_compoundpar
! Copyright A.J. Koning 2021
