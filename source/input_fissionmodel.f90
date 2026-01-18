subroutine input_fissionmodel
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for fission models
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
!   clas2file         ! file with class 2 transition states
!   fismodel          ! fission model alternative fission model for default barriers
!   fismodelalt       ! alternative fission model for default barriers
!   flagfispartdamp   ! flag for fission partial damping
!   flagsffactor      ! flag to constrain vfiscor and rmiufiscor by sp. fis. half-life
!   flagclass2        ! flag for class2 states in fission
!   flagffevap        ! flag for calculation of particle evaporation per FF
!   flagffspin        ! flag to use spin distribution in initial FF population
!   flagfisfeed       ! flag for output of fission per excitation energy bin
!   flagfission       ! flag for fission
!   flaghbstate       ! flag for head band states in fission
!   flagfisout        ! flag for output of fission information
!   flagoutfy         ! flag for output detailed fission yield calculation
!   fymodel           ! fission yield model, 1: Brosa 2: GEF
!   ffmodel           ! fission fragment model, 1: GEF 2: HF3D (Okumura) 3: SPY 4: Langevin-4D
!   pfnsmodel         ! PFNS  model, 1: Iwamoto 2: from FF decay
!   hbtransfile       ! file with head band transition states
!   yieldfile         ! fission yield file
! Constants
!   fislim          ! mass above which nuclide fissions
! Variables for basic reaction
!   flagffruns  ! flag to designate subsequent evaporation of fission products
! Variables for basic reaction
!   flagmassdis     ! flag for calculation of fission fragment mass yields
!   flagmicro    ! flag for completely microscopic TALYS calculation
! Variables for main input
!   Atarget      ! mass number of target nucleus
!   k0           ! index of incident particle
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
  real(sgl)          :: val                  ! real value
!
! ************** Defaults *************
!
  flagfission = .false.
  if (Atarget > fislim) flagfission = .true.
  if (flagffruns) flagfission = .false.
  clas2file = ' '
  yieldfile = ' '
  hbtransfile = ' '
  fisadjust = .false.
  fismodel = 6
  fismodelalt = 3
  flaghbstate = .false.
  flagclass2 = .false.
  flagfispartdamp = .false.
  flagsffactor = .false.
  flagffevap = .true.
  flagfisfeed = .false.
  flagffspin = .false.
  flagoutfy = .false.
  fymodel = 2
  ffmodel = 1
  if (flagmassdis) then
    pfnsmodel = 2
  else
    pfnsmodel = 1
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
    if (key == 'fission') then
      if (ch == 'n') flagfission = .false.
      if (ch == 'y') flagfission = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'hbstate') then
      if (ch == 'n') flaghbstate = .false.
      if (ch == 'y') flaghbstate = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'class2') then
      if (ch == 'n') flagclass2 = .false.
      if (ch == 'y') flagclass2 = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'ffevaporation') then
      if (ch == 'n') flagffevap = .false.
      if (ch == 'y') flagffevap = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'fisfeed') then
      if (ch == 'n') flagfisfeed = .false.
      if (ch == 'y') flagfisfeed = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'fispartdamp') then
      if (ch == 'n') flagfispartdamp = .false.
      if (ch == 'y') flagfispartdamp = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'sffactor') then
      if (ch == 'n') flagsffactor = .false.
      if (ch == 'y') flagsffactor = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'ffspin') then
      if (ch == 'n') flagffspin = .false.
      if (ch == 'y') flagffspin = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'fymodel') then
      read(value, * , iostat = istat) fymodel
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'ffmodel') then
      read(value, * , iostat = istat) ffmodel
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'fismodel') then
      read(value, * , iostat = istat) fismodel
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'fismodelalt') then
      read(value, * , iostat = istat) fismodelalt
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'pfnsmodel') then
      read(value, * , iostat = istat) pfnsmodel
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'gefran') then
      read(value, * , iostat = istat) gefran
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'outfy') then
      if (ch == 'n') flagoutfy = .false.
      if (ch == 'y') flagoutfy = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'hbtransfile') then
      class = 11
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) hbtransfile(Zix, Nix) = cval
      cycle
    endif
    if (key == 'class2file') then
      class = 11
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) clas2file(Zix, Nix) = cval
      cycle
    endif
    if (key == 'yieldfile') then
      yieldfile = value
      cycle
    endif
  enddo
  return
end subroutine input_fissionmodel
! Copyright A.J. Koning 2021
