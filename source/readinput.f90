subroutine readinput
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input
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
! All global variables
!   numlines   ! number of input lines
! Variables for reading input lines
!   inline     ! input line
!   nlines     ! number of input lines
!   nlines0    ! number of input lines
! Constants
!   iso        ! counter for isotope
! Error handling
!   range_integer_error ! Test if integer variable is out of range
!   read_error          ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=132) :: profile
  logical           :: lexist
  integer           :: i     ! counter
  integer           :: istat ! logical for file access
!
! ************************** User Input ********************************
!
! We read the complete input file first as a set of character strings.
! The actual keywords will be read from these later on.
! For natural elements, the input file only needs to be read once.
!
  if (iso /= 1) return
  if (nlines > 0) return
  i = 1
  do
    read(*, '(a132)', iostat = istat) inline(i)
    if (istat ==  -1) exit
    if (istat /= 0) call read_error(inline(i), istat)
    i = i + 1
    call range_integer_error('inline', i, 1, numlines)
  enddo
!
! ***** Add profile to input file ******
!
! You may edit structure/profile to insert keywords which are always included in your input files
!
  profile = trim(path)//'profile'
  inquire (file = profile, exist = lexist)
  if (lexist) then
    open (unit = 1, file = profile, status = 'unknown')
    do
      read(1, '(a132)', iostat = istat) inline(i)
      if (istat ==  -1) exit
      if (istat /= 0) call read_error(inline(i), istat)
      i = i + 1
      call range_integer_error('inline', i, 1, numlines)
    enddo
    close (unit = 1)
  endif
  nlines = i - 1
!
! ************** Convert uppercase to lowercase characters *************
!
! For easy handling of all the input parameters, the whole input is converted to lowercase characters,
! with the exception of filenames or other character strings.
!
! convert: subroutine to convert input line from upper case to lowercase
!
  do i = 1, nlines
    call convert(i)
  enddo
  nlines0 = nlines
  return
end subroutine readinput
! Copyright A.J. Koning 2021
