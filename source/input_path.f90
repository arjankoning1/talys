subroutine input_path
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for paths
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
! Variables for files
!   bestpath    ! alternative directory for best values
!   nulldev     ! null device
!   path        ! directory containing files to be read
! Variables for reading input lines
!   inline    ! input line
!   nlines    ! number of input lines
!
! *** Declaration of local data
!
  implicit none
  character(len=132) :: key      ! keyword
  character(len=132) :: value    ! value or string
  character(len=132) :: word(40) ! words on input line
  character(len=132) :: line     ! input line
  integer            :: i        ! counter
!
! ************** Defaults *************
!
  bestpath = ' '
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
!
! Test for keywords
!
! Here, the various model parameters can be set to overrule the default values.
! Most default values will be computed later on, since they require more computation (e.g. level density parameters).
!
! Each keyword is characterized by a certain order of parameter and value input.
! They are distinguished by different classes.
!
    if (key == 'strucpath') then
      path = value
      cycle
    endif
    if (key == 'nulldev') then
      nulldev = value
      cycle
    endif
    if (key == 'bestpath') then
      bestpath = value
      cycle
    endif
  enddo
  return
end subroutine input_path
! Copyright A.J. Koning 2021
