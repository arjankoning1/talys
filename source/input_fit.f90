subroutine input_fit
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for fitted model parameters 
!
! Author    : Arjan Koning
!
! 2024-06-25: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
!
  use A0_talys_mod
  use A1_error_handling_mod
!
! Variables for basic reactions
!   flagastro        ! flag for calculation of astrophysics reaction rate
!   flagfit          ! flag for using fitted nuclear model parameters
!   flagngfit        ! flag for using fitted (n,g) nuclear model parameters
!   flagnffit        ! flag for using fitted (n,f) nuclear model parameters
!   flagnnfit        ! flag for using fitted (n,n'), (n,2n) and (n,p) nuclear model parameters
!   flagnafit        ! flag for using fitted (n,a) nuclear model parameters
!   flagndfit        ! flag for using fitted (n,d) nuclear model parameters
!   flagpnfit        ! flag for using fitted (p,n) nuclear model parameters
!   flagdnfit        ! flag for using fitted (g,n) nuclear model parameters
!   flaggnfit        ! flag for using fitted (d,n) nuclear model parameters
!   flaganfit        ! flag for using fitted (a,n) nuclear model parameters
!   flaggamgamfit    ! flag for using fitted Gamma_gamma nuclear model parameters
!   flagmacsfit      ! flag for using fitted MACS nuclear model parameters
! Variables for reading input lines
!   inline            ! input line
!   nlines                ! number of input lines
! Variables for main input
!   k0             ! index of incident particle
! Error handling
!   range_index_error    ! Test if index is out of range
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
  flagngfit = (k0 == 1 .and. flagfit .and. .not.flagastro)
  flagnnfit = (k0 == 1 .and. flagfit)
  flagnffit = (k0 == 1 .and. flagfit)
  flagnafit = (k0 == 1 .and. flagfit)
  flagndfit = (k0 == 1 .and. flagfit)
  flagpnfit = (k0 == 2 .and. flagfit)
  flaggnfit = (k0 == 0 .and. flagfit)
  flagdnfit = (k0 == 3 .and. flagfit)
  flaganfit = (k0 == 6 .and. flagfit)
  flagmacsfit = (k0 == 1 .and. flagfit .and. flagastro)
  flaggamgamfit = .false.
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
    if (key == 'ngfit') then
      if (ch == 'n') flagngfit = .false.
      if (ch == 'y') then
        flagngfit = .true.
        flagfit = .true.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'nffit') then
      if (ch == 'n') flagnffit = .false.
      if (ch == 'y') then
        flagnffit = .true.
        flagfit = .true.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'nnfit') then
      if (ch == 'n') flagnnfit = .false.
      if (ch == 'y') then
        flagnnfit = .true.
        flagfit = .true.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'nafit') then
      if (ch == 'n') flagnafit = .false.
      if (ch == 'y') then
        flagnafit = .true.
        flagfit = .true.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'ndfit') then
      if (ch == 'n') flagndfit = .false.
      if (ch == 'y') then
        flagndfit = .true.
        flagfit = .true.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'pnfit') then
      if (ch == 'n') flagpnfit = .false.
      if (ch == 'y') then
        flagpnfit = .true.
        flagfit = .true.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'dnfit') then
      if (ch == 'n') flagdnfit = .false.
      if (ch == 'y') then
        flagdnfit = .true.
        flagfit = .true.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'gnfit') then
      if (ch == 'n') flaggnfit = .false.
      if (ch == 'y') then
        flaggnfit = .true.
        flagfit = .true.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'anfit') then
      if (ch == 'n') flaganfit = .false.
      if (ch == 'y') then
        flaganfit = .true.
        flagfit = .true.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'gamgamfit') then
      if (ch == 'n') flaggamgamfit = .false.
      if (ch == 'y') then
        flaggamgamfit = .true.
        flagfit = .true.
        flagngfit = .false.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'macsfit') then
      if (ch == 'n') flagmacsfit = .false.
      if (ch == 'y') then
        flagmacsfit = .true.
        flagfit = .true.
        flagngfit = .false.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
  enddo
  return
end subroutine input_fit
! Copyright A.J. Koning 2021
