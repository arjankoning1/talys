subroutine input_preeqpar
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for preequilibrium variables
!
! Author    : Arjan Koning
!
! 2022-12-07: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
!
  use A0_talys_mod
  use A1_error_handling_mod
!
! All global variables
!   numN         ! maximum number of neutrons from initial compound nucleus
!   numpar       ! number of particles
!   numZ         ! maximum number of protons from initial compound nucleus
! Variables for preequilibrium
!   preeqadjust     ! logical for energy-dependent pre-eq adjustment
!   msdbins         ! number of energy points for DWBA calculation for MSD
!   Cbreak          ! adjustable parameter for break-up reactions
!   Cknock          ! adjustable parameter for knockout reactions
!   Cstrip          ! adjustable parameter for stripping/pick-up reactions
!   Emsdmin         ! minimal outgoing energy for MSD calculation
!   Esurf0          ! well depth for surface interaction
!   g               ! single-particle level density parameter
!   gadjust         ! adjustable factor for single - particle particle-hole states
!   gn              ! single-particle neutron level density parameter
!   gnadjust        ! adjustable factor for single-particle proton parameter
!   gp              ! single-particle proton level density parameter
!   gpadjust        ! adjustable factor for single-particle neutron parameter
!   Kph             ! constant for single-particle level density parameter
!   M2constant      ! constant for matrix element in exciton model
!   M2limit         ! constant for asymptotic value for matrix element
!   M2shift         ! constant for energy shift for matrix element
!   Rgamma          ! adjustable parameter for pre-equilibrium gamma decay
!   Rnunu           ! ratio for two-component matrix element
!   Rnupi           ! ratio for two-component matrix element
!   Rpinu           ! ratio for two-component matrix element
!   Rpipi           ! ratio for two-component matrix element
!   Rspincutpreeq   ! adjustable constant (global) for preequilibrium spin cutoff factor
!   GMRadjustE      ! adjustable factor for GMR energy
!   GQRadjustE      ! adjustable factor for GQR energy
!   LEORadjustE     ! adjustable factor for LEOR energy
!   HEORadjustE     ! adjustable factor for HEOR energy
!   GMRadjustG      ! adjustable factor for GMR width
!   GQRadjustG      ! adjustable factor for GQR width
!   LEORadjustG     ! adjustable factor for LEOR width
!   HEORadjustG     ! adjustable factor for HEOR width
!   GMRadjustD      ! adjustable factor for GMR deformation parameter
!   GQRadjustD      ! adjustable factor for GQR deformation parameter
!   LEORadjustD     ! adjustable factor for LEOR deformation parameter
!   HEORadjustD     ! adjustable factor for HEOR deformation parameter
! Variables for reading input lines
!   inline            ! input line
!   nlines            ! number of input lines
! Error handling
!   read_error ! Message for file reading error
! Variables for main input
!   k0                ! index of incident particle
!
! *** Declaration of local data
!
  implicit none
  logical            :: flagassign           ! flag to assign value or not
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
! ************* Defaults ************************
!
  msdbins = 6
  Cbreak = 1.
  if (k0 == 6) Cbreak = 0.
  Cknock = 1.
  Cstrip = 1.
  Emsdmin = 0.
  Esurf0 = -1.
  g = 0.
  gadjust = 1.
  gn = 0.
  gnadjust = 1.
  gp = 0.
  gpadjust = 1.
  Kph = 15.
  M2constant = 1.
  M2limit = 1.
  M2shift = 1.
  Rgamma = 2.
  Rnunu = 1.5
  Rnupi = 1.
  Rpinu = 1.
  Rpipi = 1.
  Rspincutpreeq = 1.
  GMRadjustE = 1.
  GQRadjustE = 1.
  LEORadjustE = 1.
  HEORadjustE = 1.
  GMRadjustG = 1.
  GQRadjustG = 1.
  LEORadjustG = 1.
  HEORadjustG = 1.
  GMRadjustD = 1.
  GQRadjustD = 1.
  LEORadjustD = 1.
  HEORadjustD = 1.
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
    if (key == 'g') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) g(Zix, Nix) = val
      cycle
    endif
    if (key == 'gp') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) gp(Zix, Nix) = val
      cycle
    endif
    if (key == 'gadjust') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) gadjust(Zix, Nix) = val
      cycle
    endif
    if (key == 'gn') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) gn(Zix, Nix) = val
      cycle
    endif
    if (key == 'gnadjust') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) gnadjust(Zix, Nix) = val
      cycle
    endif
    if (key == 'gpadjust') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) gpadjust(Zix, Nix) = val
      cycle
    endif
    if (key == 'cstrip') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Cstrip(type) = val
      cycle
    endif
    if (key == 'cknock') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Cknock(type) = val
      cycle
    endif
    if (key == 'cbreak') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Cbreak(type) = val
      cycle
    endif
    if (key == 'kph') then
      read(value, * , iostat = istat) Kph
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'm2constant') then
      class = 9
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        M2constant = val
        preeqadjust = .true.
      endif
      cycle
    endif
    if (key == 'm2limit') then
      read(value, * , iostat = istat) M2limit
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'm2shift') then
      read(value, * , iostat = istat) M2shift
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'rpipi') then
      read(value, * , iostat = istat) Rpipi
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'rnunu') then
      read(value, * , iostat = istat) Rnunu
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'rpinu') then
      read(value, * , iostat = istat) Rpinu
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'rnupi') then
      read(value, * , iostat = istat) Rnupi
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'rspincutpreeq') then
      read(value, * , iostat = istat) Rspincutpreeq
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'esurf') then
      read(value, * , iostat = istat) Esurf0
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'rgamma') then
      read(value, * , iostat = istat) Rgamma
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'msdbins')  then
      read(value, * , iostat = istat) msdbins
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'emsdmin')  then
      read(value, * , iostat = istat) Emsdmin
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'gmradjuste')  then
      read(value, * , iostat = istat) GMRadjustE
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'gqradjuste')  then
      read(value, * , iostat = istat) GQRadjustE
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'leoradjuste')  then
      read(value, * , iostat = istat) LEORadjustE
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'heoradjuste')  then
      read(value, * , iostat = istat) HEORadjustE
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'gmradjustg')  then
      read(value, * , iostat = istat) GMRadjustG
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'gqradjustg')  then
      read(value, * , iostat = istat) GQRadjustG
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'leoradjustg')  then
      read(value, * , iostat = istat) LEORadjustG
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'heoradjustg')  then
      read(value, * , iostat = istat) HEORadjustG
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'gmradjustd')  then
      read(value, * , iostat = istat) GMRadjustD
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'gqradjustd')  then
      read(value, * , iostat = istat) GQRadjustD
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'leoradjustd')  then
      read(value, * , iostat = istat) LEORadjustD
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'heoradjustd')  then
      read(value, * , iostat = istat) HEORadjustD
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
  enddo
  return
end subroutine input_preeqpar
! Copyright A.J. Koning 2021
