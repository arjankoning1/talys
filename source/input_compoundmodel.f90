subroutine input_compoundmodel
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for compound reaction models
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
! Variables for compound reactions
!   flagcomp        ! flag for compound nucleus calculation
!   flageciscomp    ! flag for compound nucleus calculation by ECIS
!   flagfullhf      ! flag for full spin dependence of transmission coefficients
!   flaggroup       ! flag for output of low energy groupwise cross sections
!   flagres         ! flag for output of low energy resonance cross sections
!   flagurrnjoy     ! normalization of URR parameters with NJOY method
!   lenreslib       ! length of library name with resonance parameters
!   wmode           ! designator for width fluctuation model
!   skipCN          ! flag to skip compound nucleus in evaporation chain
!   ewfc            ! off-set incident energy for width fluctuation correction
! Variables for basic reactions
!   flagendf      ! flag for information for ENDF-6 file
! Variables for main input
!   Atarget    ! mass number of target nucleus
!   k0         ! index of incident particle
! Variables for OMP
!   flagomponly     ! flag to execute ONLY an optical model calculation
! Variables for reading input lines
!   inline            ! input line
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
  if (flagomponly) then
    flagcomp = .false.
  else
    flagcomp = .true.
  endif
  flageciscomp = .false.
  flagfullhf = .false.
  flaggroup = .false.
  flagres = .false.
  flagurrnjoy = .false.
  skipCN = 0
  if (k0 == 1) then
    wmode = 1
  else
    wmode = 2
  endif
  wfcfactor = 1
  ewfc = -1.
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
    if (key == 'compound') then
      if (ch == 'n') flagcomp = .false.
      if (ch == 'y') flagcomp = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'widthfluc') then
      if (ch == 'y') then
        if (k0 > 1) ewfc = 10.
        cycle
      endif
      if (ch == 'n') then
        ewfc = 0.
        cycle
      endif
      read(value, * , iostat = istat) ewfc
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'eciscompound') then
      if (ch == 'n') flageciscomp = .false.
      if (ch == 'y') flageciscomp = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'fullhf') then
      if (ch == 'n') flagfullhf = .false.
      if (ch == 'y') flagfullhf = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'widthmode') then
      read(value, * , iostat = istat) wmode
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'wfcfactor') then
      read(value, * , iostat = istat) WFCfactor
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'skipcn') then
      class = 0
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) skipCN(Zix, Nix) = 1
      cycle
    endif
    if (key == 'resonance') then
      if (ch == 'n') flagres = .false.
      if (ch == 'y') flagres = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'group') then
      if (ch == 'n') flaggroup = .false.
      if (ch == 'y') flaggroup = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'urrnjoy') then
      if (ch == 'n') flagurrnjoy = .false.
      if (ch == 'y') flagurrnjoy = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
  enddo
  return
end subroutine input_compoundmodel
! Copyright A.J. Koning 2021
