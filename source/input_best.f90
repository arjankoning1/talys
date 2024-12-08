subroutine input_best
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for best files
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
! All global variables
!   numisom     ! number of isomers
!   numlines    ! number of input lines
!   nummt       ! number of MT numbers
! Variables for best files
!   flagbest         ! flag to use best set of adjusted parameters
!   flagbestend      ! flag to put best set of parameters at end of input file
!   flagrescue       ! flag for final rescue: normalization to data
!   rescuefile       ! file with incident energy dependent adjustment factors
!   grescue          ! global multiplication factor for incident energy dependence
! Variables for files
!   path        ! directory containing files to be read
! Constants
!   iso            ! counter for isotope
!   nuc            ! symbol of nucleus
! Variables for reading input lines
!   inline         ! input line
!   nlines         ! number of input lines
! Error handling
!   range_index_error    ! Test if index is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=1)   :: ch                   ! character
  character(len=132) :: key                  ! keyword
  character(len=132) :: word(40)             ! words on input line
  character(len=132) :: line                 ! input line
  integer            :: i                    ! counter
  integer            :: istat                ! logical for file access
  integer            :: is                   ! counter for isomer
  integer            :: k                    ! counter
  integer            :: mt                   ! MT number
  real(sgl)          :: val                  ! real value
!
! ************** Defaults *************
!
  flagbest = .false.
  flagbestend = .false.
  flagrescue = .false.
  rescuefile = ' '
  grescue = 1.
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
    if (key == 'best') then
      if (ch == 'n') flagbest = .false.
      if (ch == 'y') flagbest = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'bestend') then
      if (ch == 'n') flagbestend = .false.
      if (ch == 'y') flagbestend = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'rescuefile') then
      read(word(2), * , iostat = istat) mt
      if (istat /= 0) call read_error(line, istat)
      call range_index_error(key, 'MT', mt, 1, nummt)
      is = - 1
      do k = 1, 71
        if (word(3)(k:k+1) == '_g') then
          is = 0
          exit
        endif
        if (word(3)(k:k+1) == '_m') then
          is = 1
          exit
        endif
        if (word(3)(k:k+1) == '_n') then
          is = 2
          exit
        endif
        if (word(3)(k:k+1) == '_o') then
          is = 3
          exit
        endif
        if (word(3)(k:k+1) == '_p') then
          is = 4
          exit
        endif
        if (word(3)(k:k+1) == '_q') then
          is = 5
          exit
        endif
        if (word(3)(k:k+1) == '_r') then
          is = 6
          exit
        endif
        if (word(3)(k:k+1) == '_s') then
          is = 7
          exit
        endif
        if (word(3)(k:k+1) == '_t') then
          is = 8
          exit
        endif
        if (word(3)(k:k+1) == '_u') then
          is = 9
          exit
        endif
        if (word(3)(k:k+1) == '_v') then
          is = 10
          exit
        endif
      enddo
      rescuefile(mt, is) = word(3)
      val = 1.
      read(word(4), * , iostat = istat) val
      if (istat > 0) call read_error(line, istat)
      grescue(mt, is) = val
      flagrescue = .true.
      cycle
    endif
  enddo
  return
end subroutine input_best
! Copyright A.J. Koning 2021
