subroutine input_basicpar
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for basic parameter variables
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
! Variables for basic parameters
!   eninclow     ! minimal incident energy for nuclear model calculation
!   flagequi     ! flag to use equidistant excitation instead of logarithmic bins
!   flagequispec ! flag to use equidistant bins for emission spectra
!   isomer       ! definition of isomer in seconds
!   Lisoinp      ! user assignment of target isomer number
!   outtype      ! type of outgoing particles
!   flagfit          ! flag to use automatically fitted parameters
! Constants
!   parsym          ! symbol of particle
! Variables for reading input lines
!   inline            ! input line
!   nlines                ! number of input lines
! Variables for basic reaction
!   flagendf    ! flag for information for ENDF - 6 file
! Variables for main input
!   k0            ! index of incident particle
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
  integer            :: ix                   ! index
  integer            :: i                    ! counter
  integer            :: i2                   ! counter
  integer            :: ip                   ! counter
  integer            :: istat                ! logical for file access
  integer            :: type                 ! particle type
!
! ************** Defaults *************
!
  if (k0 == 1 .and. flagendf) then
    eninclow = 0.
  else
    eninclow = 1.e-6
  endif
  flagequi = .true.
  flagequispec = .false.
  isomer = 1.
  Lisoinp = -1
  outtype = ' '
  source = 'TALYS-2.1'
  oformat = 'YANDF-0.2'
  flagfit = .false.
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
    if (key == 'elow') then
      read(value, * , iostat = istat) eninclow
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'ejectiles') then
      ip = - 1
Loop1: do i2 = 2, 40
        ch = word(i2)(1:1)
        do type = 0, 6
          if (ch == parsym(type)) then
            ip = ip + 1
            if (ip <= 6) outtype(ip) = ch
            cycle Loop1
          endif
        enddo
        if (ip ==  -1)  call read_error('ejectiles', istat)
      enddo Loop1
      cycle
    endif
    if (key == 'liso') then
      read(value, * , iostat = istat) Lisoinp
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'isomer') then
      read(value, * , iostat = istat) isomer
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'equidistant') then
      if (ch == 'n') flagequi = .false.
      if (ch == 'y') flagequi = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'equispec') then
      if (ch == 'n') flagequispec = .false.
      if (ch == 'y') flagequispec = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'source') then
      ix=index(line,'source')+7
      source=trim(adjustl(line(ix:132)))
      cycle
    endif
    if (key == 'user') then
      ix=index(line,'user')+5
      user=trim(adjustl(line(ix:132)))
      cycle
    endif
    if (key == 'format') then
      ix=index(line,'format')+7
      oformat=trim(adjustl(line(ix:132)))
      cycle
    endif
    if (key == 'fit') then
      if (ch == 'n') flagfit = .false.
      if (ch == 'y') flagfit = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif 
  enddo
  return
end subroutine input_basicpar
! Copyright A.J. Koning 2021
