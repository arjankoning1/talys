subroutine input_medical
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for medical variables
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
! Variables for medical isotope production
!   Area          ! target area in cm^2
!   Eback         ! lower end of energy range in MeV for isotope
!   Ebeam         ! incident energy in MeV for isotope production
!   flagprod      ! flag for isotope production
!   Ibeam         ! beam current in mA for isotope production
!   radiounit     ! unit for radioactivity: Bq, kBq, MBq, Gbq, mCi,
!   rhotarget     ! target material density
!   Tcool         ! cooling time per unit cooling time unit (y,d,h,m,s)
!   Tirrad        ! irradiation time per unit irradiation time unit
!   unitTcool     ! cooling time unit (y,d,h,m,s)
!   unitTirrad    ! irradiation time unit (y,d,h,m,s)
!   yieldunit     ! unit for isotope yield: num (number), mug, mg, g, or kg
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
  integer            :: k                    ! counter
!
! ************** Defaults *************
!
  Area = 1.
  Eback = -1.
  Ebeam = -1.
  flagprod = .false.
  Ibeam = 1.
  radiounit = 'gbq'
  rhotarget = -1.
  Tcool = 0
  unitTcool = ' '
  Tcool(1) = 1
  unitTcool(1) = 'd'
  Tirrad = 0
  unitTirrad = ' '
  Tirrad(1) = 1
  unitTirrad(1) = 'd'
  yieldunit = 'num'
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
    if (key == 'area') then
      read(value, * , iostat = istat) Area
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'eback') then
      read(value, * , iostat = istat) Eback
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'ebeam') then
      read(value, * , iostat = istat) Ebeam
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'ibeam') then
      read(value, * , iostat = istat) Ibeam
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'production') then
      if (ch == 'n') flagprod = .false.
      if (ch == 'y') flagprod = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'radiounit') then
      read(value, * , iostat = istat) radiounit
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'rho') then
      read(value, * , iostat = istat) rhotarget
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'tcool') then
      Tcool(1) = 0
      unitTcool(1) = ' '
      do k = 1, 5
        read(word(2 * k), '(i9)' , iostat = istat) Tcool(k)
        if (istat > 0) cycle
        if (istat /= 0) call read_error(line, istat)
        read(word(2 * k + 1), '(a1)', iostat = istat) unitTcool(k)
        if (istat /= 0) call read_error(line, istat)
      enddo
      cycle
    endif
    if (key == 'tirrad') then
      Tirrad(1) = 0
      unitTirrad(1) = ' '
      do k = 1, 5
        read(word(2*k), '(i9)', iostat = istat) Tirrad(k)
        if (istat > 0) cycle
        if (istat /= 0) call read_error(line, istat)
        read(word(2 * k + 1), '(a1)', iostat = istat) unitTirrad(k)
        if (istat /= 0) call read_error(line, istat)
      enddo
      cycle
    endif
    if (key == 'yieldunit') then
      read(value, * , iostat = istat) yieldunit
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
  enddo
  return
end subroutine input_medical
! Copyright A.J. Koning 2021
