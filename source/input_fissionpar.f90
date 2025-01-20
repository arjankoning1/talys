subroutine input_fissionpar
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for fission variables
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
! Variables for fission
!   axtype            ! type of axiality of barrier
!   bdamp             ! fission partial damping parameter
!   bdampadjust       ! correction for fission partial damping parameter
!   betafiscor        ! adjustable factor for fission path width
!   betafiscoradjust  ! adjustable factor for fission path width
!   Cnubar1           ! adjustable parameter for nubar constant value
!   Cnubar2           ! adjustable parameter for nubar energy slope
!   Cbarrier          ! global multiplier for fission barrier for Sierk model
!   Tmadjust          ! adjustable parameter for PFNS temperature
!   Fsadjust          ! adjustable parameter for PFNS scission fraction
!   fbaradjust        ! adjustable factor for fission parameters
!   fbarrier          ! height of fission barrier
!   fisadjust         ! logical for energy-dependent fission adjustment
!   fismodel          ! fission model alternative fission model for default barriers
!   fismodelx         ! fission model
!   flagfisout        ! flag for output of fission information
!   flagoutfy         ! flag for output detailed fission yield calculation
!   fismodel          ! fission model alternative fission model for default barriers
!   fwidth            ! width of fission barrier
!   fwidthadjust      ! adjustable factor for fission parameters
!   fymodel           ! fission yield model, 1: Brosa 2: GEF
!   gefran            ! number of random events for GEF calculation
!   Rfiseps           ! ratio for limit for fission cross section per nucleus
!   vfiscor           ! adjustable factor for fission path height
!   vfiscoradjust     ! adjustable factor for fission path height
!   widthc2           ! width of class2 states
! All global variables
!   numbar            ! number of fission barriers
!   numN              ! maximum number of neutrons from initial compound nucleus
!   numZ              ! maximum number of protons from initial compound nucleus
! Constants
!   fislim          ! mass above which nuclide fissions
! Variables for basic reaction
!   flagastro    ! flag for calculation of astrophysics reaction rate
!   flagbasic    ! flag for output of basic information and results
!   flagendf     ! flag for information for ENDF - 6 file
!   flagmassdis     ! flag for calculation of fission fragment mass yields
!   flagmicro    ! flag for completely microscopic TALYS calculation
! Variables for main input
!   Atarget      ! mass number of target nucleus
!   k0           ! index of incident particle
!   Ninit        ! neutron number of initial compound nucleus
!   Zinit        ! charge number of initial compound nucleus
! Variables for reading input lines
!   inline                ! input line
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
  integer            :: A                    ! mass number of residual nucleus
  integer            :: Aact                 ! help variable
  integer            :: N                    ! neutron number of residual nucleus
  integer            :: oddN                 ! help variable
  integer            :: oddZ                 ! help variable
  integer            :: Z                    ! charge number of target nucleus
  real(sgl)          :: rfiscor              ! helpvariable
  real(sgl)          :: Vf0                  ! helpvariable
  real(sgl)          :: val                  ! real value
!
! ************** Defaults *************
!
  axtype = 1
  bdamp = 0.01
  bdampadjust = 1.
  betafiscor = 1.
  betafiscoradjust = 1.
  fbaradjust = 1.
  fbarrier = 0.
  fisadjust = .false.
  flagfisout = flagbasic
  if (flagfission) then
    if (flagendf) flagfisout = .true.
  else
    flagfisout = .false.
  endif
  fismodelx = fismodel
  fwidth = 0.
  fwidthadjust = 1.
  gefran = 50000
  Rfiseps = 1.e-3
  if (flagmassdis) Rfiseps = 1.e-9
  if (flagastro) Rfiseps = 1.e-6
  vfiscor = 1.
  vfiscoradjust = 1.
  Cnubar1 = 1.
  Cnubar2 = 1.
  Tmadjust = 1.
  Fsadjust = 1.
  if (Atarget <= fislim) then
    if (ldmodel(0,0) <= 3) then
      Cbarrier = 0.85
    else
      Cbarrier = 2.00
    endif
  else
    if (ldmodel(0,0) <= 3) then
      Cbarrier = 1.20
    else
      Cbarrier = 2.00
    endif
  endif
  do Zix = 0, numZ
    do Nix = 0, numN
      Z = Zinit - Zix
      N = Ninit - Nix
      oddZ = mod(Z, 2)
      oddN = mod(N, 2)
      Vf0 = 0.86
      if (oddZ == 0 .and. oddN == 0) vf0 = 0.83
      if (oddZ == 1 .and. oddN == 0) vf0 = 0.91
      if (oddZ == 0 .and. oddN == 1) vf0 = 0.86
      if (oddZ == 1 .and. oddN == 1) vf0 = 0.85
      A = Z + N
      Aact = max(min(A, 255), 225)
      rfiscor = 0.005
      vfiscor(Zix, Nix) = Vf0 - rfiscor * (Aact - 240)
      if (Ninit - Nix > 144 .or. fismodel == 5) axtype(Zix, Nix, 1) = 3
      if (fismodel < 5) axtype(Zix, Nix, 2) = 2
    enddo
  enddo
  widthc2 = 0.2
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
    if (key == 'fisbar') then
      ibar = 1
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        fbarrier(Zix, Nix, ibar) = val
        fisadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'fishw') then
      ibar = 1
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        fwidth(Zix, Nix, ibar) = val
        fisadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'fisbaradjust') then
      ibar = 1
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        fbaradjust(Zix, Nix, ibar) = val
        fisadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'fishwadjust') then
      ibar = 1
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        fwidthadjust(Zix, Nix, ibar) = val
        fisadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'bdamp') then
      ibar = 1
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        bdamp(Zix, Nix, ibar) = val
        fisadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'bdampadjust') then
      ibar = 1
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        bdampadjust(Zix, Nix, ibar) = val
        fisadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'class2width') then
      ibar = 1
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) widthc2(Zix, Nix, ibar) = val
      cycle
    endif
    if (key == 'axtype') then
      ibar = 1
      class = 4
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) axtype(Zix, Nix, ibar) = ival
      cycle
    endif
    if (key == 'vfiscor') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) vfiscor(Zix, Nix) = val
      cycle
    endif
    if (key == 'vfiscoradjust') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) vfiscoradjust(Zix, Nix) = val
      cycle
    endif
    if (key == 'betafiscor') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) betafiscor(Zix, Nix) = val
      cycle
    endif
    if (key == 'cbarrier') then
      class = 9
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Cbarrier = val
      cycle
    endif
    if (key == 'betafiscoradjust') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) betafiscoradjust(Zix, Nix) = val
      cycle
    endif
    if (key == 'gefran') then
      read(value, * , iostat = istat) gefran
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'outfission') then
      if (ch == 'n') flagfisout = .false.
      if (ch == 'y') flagfisout = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'rfiseps') then
      read(value, * , iostat = istat) Rfiseps
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'cnubar1') then
      read(value, * , iostat = istat) Cnubar1
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'cnubar2') then
      read(value, * , iostat = istat) Cnubar2
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'tmadjust') then
      class = 9
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Tmadjust = val
      cycle
    endif
    if (key == 'fsadjust') then
      class = 9
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Fsadjust = val
      cycle
    endif
  enddo
  return
end subroutine input_fissionpar
! Copyright A.J. Koning 2021
