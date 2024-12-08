subroutine input_densitymodel
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for level density variables
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
! Definition of single and double precision variables
!   sgl               ! single precision kind
! Variables for input of level density parameters
!   cglobal              ! global constant to adjust tabulated level densities
!   col                  ! flag for collective enhancement of level density
!   colall               ! flag for collective enhancement of level density
!   filedensity          ! flag for level densities on separate files
!   flagasys             ! flag for all level density parameters a from systematic
!   flagcol              ! flag for collective enhancement of level density
!   flagcolall           ! flag for collective enhancement of level density
!   flagcolldamp         ! flag for damping of coll. effects in eff. level density (without explicit coll. enh.)
!   flagctmglob          ! flag for global CTM model (no discrete level info)
!   flagdensity          ! flag for output of level densities
!   kvibmodel            ! model for vibrational enhancement
!   ldmodel              ! level density model
!   ldmodelCN            ! level density model for compound nucleus
!   ldmodelall           ! level density model for all nuclides
!   pglobal              ! global constant to adjust tabulated level densities
!   Rspincutff           ! adjustable parameter (global) for FF spin cutoff factor
!   shellmodel           ! model for shell correction energies
!   spincutmodel         ! model for spin cutoff factor for ground state
!   strength             ! model for E1 gamma-ray strength function
! All global variables
!   numbar    ! number of fission barriers
!   numN      ! maximum number of neutrons from initial compound nucleus
!   numZ      ! maximum number of protons from initial compound nucleus
! Variables for reading input lines
!   inline            ! input line
!   nlines            ! number of input lines
! Constants
!   fislim
! Variables for main input
!   Atarget                  ! mass number of target nucleus
! Variables for basic reaction
!   flagmicro    ! flag for completely microscopic TALYS calculation
!   flagbasic    ! flag for output of basic information and results
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical            :: flagassign            ! flag to assign value or not
  logical            :: fcol                  ! flag for collective enhancement
  character(len=1)   :: ch                    ! character
  character(len=132) :: cval                  ! character value
  character(len=132) :: key                   ! keyword
  character(len=132) :: value                 ! value or string
  character(len=132) :: word(40)              ! words on input line
  character(len=132) :: line                  ! input line
  integer            :: class                 ! input class
  integer            :: col(0:numZ,0:numN)    ! help variable for collective enhancement
  integer            :: i                     ! counter
  integer            :: ibar                  ! fission barrier
  integer            :: igr                   ! giant resonance
  integer            :: irad                  ! variable to indicate M(=0) or E(=1) radiation
  integer            :: istat                 ! logical for file access
  integer            :: ival                  ! integer value
  integer            :: lval                  ! multipolarity
  integer            :: Nix                   ! neutron number index for residual nucleus
  integer            :: type                  ! particle type
  integer            :: Zix                   ! charge number index for residual nucleus
  real(sgl)          :: val                   ! real value
!
! ************** Defaults *************
!
  if (flagmicro) then
    ldmodelall = 5
  else
    ldmodelall = 1
  endif
  strength= 8
  if (k0 <= 1 .and. Atarget > fislim) ldmodelall = 5
  ldmodelCN = 0
  shellmodel = 1
  spincutmodel = 1
  kvibmodel = 2
  if (Atarget > fislim) then
    flagcolall = .true.
  else
    flagcolall = .false.
  endif
  flagcol = flagcolall
  flagcolldamp = .false.
  ldmodel = 0
  col = 0.
  flagasys = .false.
  flagctmglob = .false.
  flagdensity = flagbasic
  filedensity = .false.
  cglobal = 1.e-20
  pglobal = 1.e-20
  Rspincutff = 4.
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
    if (key == 'spincutmodel') then
      read(value, * , iostat = istat) spincutmodel
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'shellmodel') then
      read(value, * , iostat = istat) shellmodel
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'kvibmodel') then
      read(value, * , iostat = istat) kvibmodel
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'colldamp') then
      if (ch == 'n') flagcolldamp = .false.
      if (ch == 'y') flagcolldamp = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'ctmglobal') then
      if (ch == 'n') flagctmglob = .false.
      if (ch == 'y') flagctmglob = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'asys') then
      if (ch == 'n') flagasys = .false.
      if (ch == 'y') flagasys = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'ldmodelcn') then
      read(value, * , iostat = istat) ldmodelCN
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'ldmodel') then
      if (word(3) == ' ') then
        read(value, * , iostat = istat) ldmodelall
        if (istat /= 0) call read_error(line, istat)
      else
        class = 14
        call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
        if (flagassign) ldmodel(Zix, Nix) = ival
      endif
      cycle
    endif
    if (key == 'colenhance') then
      if (ch == 'n') fcol = .false.
      if (ch == 'y') fcol = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      if (word(3) == ' ') then
        flagcolall = fcol
      else
        class = 13
        call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
        if (flagassign) then
          flagcol(Zix, Nix) = fcol
          col(Zix, Nix) = 1
        endif
      endif
      cycle
    endif
    if (key == 'cglobal') then
      read(value, * , iostat = istat) cglobal
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'pglobal') then
      read(value, * , iostat = istat) pglobal
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'rspincutff') then
      class = 9
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Rspincutff = val
      cycle
    endif
    if (key == 'filedensity') then
      if (ch == 'n') filedensity = .false.
      if (ch == 'y') filedensity = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outdensity') then
      if (ch == 'n') flagdensity = .false.
      if (ch == 'y') flagdensity = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'strength') then
      read(value, * , iostat = istat) strength
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
  enddo
!
! Set level density models and spectroscopic factors per nucleus
!
  if (ldmodelCN > 0) then
    ldmodel(0,0) = ldmodelCN
  else
    ldmodelCN = ldmodelall
  endif
  do Nix = 0, numN
    do Zix = 0, numZ
      if (ldmodel(Zix, Nix) == 0) ldmodel(Zix, Nix) = ldmodelall
      if (col(Zix, Nix) == 0) flagcol(Zix, Nix) = flagcolall
    enddo
  enddo
  return
end subroutine input_densitymodel
! Copyright A.J. Koning 2021
