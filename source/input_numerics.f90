subroutine input_numerics
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for numerics variables
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
! Variables for numerics
!   maxenrec      ! number of recoil energies
!   maxN          ! maximal number of neutrons away from initial compound nucleus
!   maxNrp        ! maximal number of neutrons away from the initial compound nucleus
!   maxZ          ! maximal number of protons away from initial compound nucleus
!   maxZrp        ! maximal number of protons away from the initial compound nucleus
!   nangle        ! number of angles
!   nanglecont    ! number of angles for continuum
!   nbins0        ! number of continuum excitation energy bins
!   segment       ! help array for storing segment intersection points
!   maxchannel    ! maximal number of outgoing particles in individual channel description
!   nanglerec     ! number of recoil angles
!   transpower    ! power for transmission coefficient limit
!   popeps        ! limit for population cross sections
!   transeps      ! absolute limit for transmission coefficient
!   xseps         ! limit for cross sections
! All global variables
!   numang       ! maximum number of angles
!   numangrec    ! maximum number of recoil angles
!   numenrec     ! maximum number of recoil energies
!   numN         ! maximum number of neutrons from initial compound nucleus
!   numZ         ! maximum number of protons from initial compound nucleus
! Variables for basic reaction
!   flagffruns        ! flag to designate subsequent evaporation of fission products
!   flagrpruns        ! flag to designate that run is for residual product
! Variables for basic reaction
!   flagastro      ! flag for calculation of astrophysics reaction rate
!   flaglabddx     ! flag for calculation of DDX in LAB system
!   flagmassdis    ! flag for calculation of fission fragment mass yields
! Variables for reading input lines
!   inline            ! input line
!   nlines            ! number of input lines
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=132) :: key                  ! keyword
  character(len=132) :: value                ! value or string
  character(len=132) :: word(40)             ! words on input line
  character(len=132) :: line                 ! input line
  integer            :: i                    ! counter
  integer            :: istat                ! logical for file access
!
! ************** Defaults *************
!
  maxenrec = numenrec
  maxN = numN - 2
  maxNrp = numN - 2
  maxZ = numZ - 2
  maxZrp = numZ - 2
  nbins0 = 40
  segment = 1
  nangle = numang
  nanglecont = 18
  maxchannel = 4
  if (flaglabddx) then
    nanglerec = numangrec
  else
    nanglerec = 1
  endif
  transpower = 5
  transeps = 1.e-8
  xseps = 1.e-7
  popeps = 1.e-3
  if (flagmassdis) then
    transpower = 10
    transeps = 1.e-12
    xseps = 1.e-12
    popeps = 1.e-10
  endif
  if (flagastro) then
    transpower = 15
    transeps = 1.e-18
    xseps = 1.e-17
    popeps = 1.e-13
  endif
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
! getvalues : subroutine to assign values to keywords
!
    if (key == 'anglesrec')  then
      read(value, * , iostat = istat) nanglerec
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'angles') then
      read(value, * , iostat = istat) nangle
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'anglescont') then
      read(value, * , iostat = istat) nanglecont
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'maxchannel') then
      read(value, * , iostat = istat) maxchannel
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'maxenrec') then
      read(value, * , iostat = istat) maxenrec
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'maxzrp') then
      read(value, * , iostat = istat) maxZrp
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'maxnrp') then
      read(value, * , iostat = istat) maxNrp
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'maxz') then
      read(value, * , iostat = istat) maxZ
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'maxn') then
      read(value, * , iostat = istat) maxN
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'bins') then
      read(value, * , iostat = istat) nbins0
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'segment') then
      read(value, * , iostat = istat) segment
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
   if (key == 'transpower') then
      read(value, * , iostat = istat) transpower
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'transeps') then
      read(value, * , iostat = istat) transeps
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'xseps') then
      read(value, * , iostat = istat) xseps
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'popeps') then
      read(value, * , iostat = istat) popeps
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
  enddo
  return
end subroutine input_numerics
! Copyright A.J. Koning 2021
