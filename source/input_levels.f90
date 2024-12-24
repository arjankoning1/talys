subroutine input_levels
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for discrete level variables
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
! Variables for discrete levels
!   branchlevel     ! level to which branching takes place
!   branchratio     ! gamma-ray branching ratio to level
!   deformfile      ! deformation parameter file
!   levelfile       ! discrete level file
!   disctable       ! table with discrete levels
!   flagbestbr      ! flag to use only best set of branching ratios
!   flagelectron    ! flag for application of electron conversion coefficient
!   flaglevels      ! flag for output of discrete level information
!   nbranch         ! number of branching levels
!   nlev            ! number of levels for nucleus
!   nlevbin         ! number of excited levels for binary nucleus
!   nlevmax         ! maximum number of included discrete levels for target nucleus
!   nlevmaxres      ! maximum number of included discrete levels for residual nucleus
! All global variables
!   numlev    ! maximum number of discrete levels
!   numN      ! maximum number of neutrons from initial compound nucleus
!   numpar    ! number of particles
!   numZ      ! maximum number of protons from initial compound nucleus
! Variables for reading input lines
!   inline    ! input line
!   nlines    ! number of input lines
! Variables for main input
!   Ninit     ! neutron number of initial compound nucleus
!   Zinit     ! charge number of initial compound nucleus
! Variables for main input
!   k0         ! index of incident particle
!   Ltarget    ! excited level of target
! Variables for basic reaction
!   flagbasic    ! flag for output of basic information and results
!   flagendf     ! flag for information for ENDF-6 file
! Error handling
!   range_index_error    ! Test if index is out of range
!   range_real_error    ! Test if real variable is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical            :: flagassign           ! flag to assign value or not
  logical            :: flagoutside          ! flag for value outside range
  character(len=1)   :: ch                   ! character
  character(len=132) :: cval                 ! character value
  character(len=132) :: key                  ! keyword
  character(len=132) :: value                ! value or string
  character(len=132) :: word(40)             ! words on input line
  character(len=132) :: line                 ! input line
  integer            :: class                ! input class
  integer            :: i                    ! counter
  integer            :: ia                   ! mass number from abundance table
  integer            :: ibar                 ! fission barrier
  integer            :: igr                  ! giant resonance
  integer            :: ilev0                ! counter for level
  integer            :: ilev1                ! counter for level
  integer            :: in                   ! neutron number of residual nucleus
  integer            :: irad                 ! variable to indicate M(=0) or E(=1) radiation
  integer            :: istat                ! logical for file access
  integer            :: ival                 ! integer value
  integer            :: iword                ! word counter
  integer            :: iz                   ! charge number of residual nucleus
  integer            :: k                    ! counter
  integer            :: lval                 ! multipolarity
  integer            :: nbr                  ! number of branching levels
  integer            :: Nix                  ! neutron number index for residual nucleus
  integer            :: type                 ! particle type
  integer            :: Zix                  ! charge number index for residual nucleus
  real(sgl)          :: br                   ! branching ratio multiplied by initial flux
  real(sgl)          :: sum                  ! help variable
  real(sgl)          :: val                  ! real value
!
! ************** Defaults *************
!
  branchlevel = 0
  branchratio = 0.
  Risomer = 1.
  flagpseudores = .false.
  flagelectron = .true.
  if (flagendf) flagelectron = .true.
  deformfile = ' '
  levelfile = ' '
  disctable = 1
  flagbestbr = .false.
  flaglevels = flagbasic
  nbranch = 0
  nlev = 0
  nlevmax = max(30, Ltarget)
  nlevbin = 30
  nlevmaxres = 30
! do type = 0, 6
!   if (type <= 2 .or. type == 6) then
!     nlevbin(type) = 10
!   else
!     nlevbin(type) = 5
!   endif
! enddo
  nlevbin(k0) = nlevmax
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
    if (key == 'branch') then
      read(word(2), * , iostat = istat) iz
      if (istat /= 0) call read_error(line, istat)
      read(word(3), * , iostat = istat) ia
      if (istat /= 0) call read_error(line, istat)
      read(word(4), * , iostat = istat) ilev0
      if (istat /= 0) call read_error(line, istat)
      read(word(5), * , iostat = istat) nbr
      if (istat /= 0) call read_error(line, istat)
      Zix = Zinit - iz
      in = ia - iz
      Nix = Ninit - in
      call range_index_error(key, 'Z', iz, Zinit - numZ, Zinit, error = 'continue', flagoutside = flagoutside)
      if (.not. flagoutside) call range_index_error(key, 'N', in, Ninit - numN, Ninit, error = 'continue', &
 &      flagoutside = flagoutside)
      if (flagoutside) write(*,'("Z= ",i3," A= ",i3," out of range")')  iz, ia
      if (flagoutside) cycle
      call range_index_error(key, 'Level', ilev0, 0, numlev)
      call range_index_error(key, '# branchings', nbr, 0, numlev)
      iword = 5
      sum = 0.
      do k = 1, nbr
        iword = iword + 1
        read(word(iword), * , iostat = istat) ilev1
        if (istat /= 0) call read_error(line, istat)
        call range_index_error(key, 'Level', ilev1, 0, numlev)
        iword = iword + 1
        read(word(iword), * , iostat = istat) br
        if (istat /= 0) call read_error(line, istat)
        call range_real_error(key, br, 0., 1.)
        branchlevel(Zix, Nix, ilev0, k) = ilev1
        branchratio(Zix, Nix, ilev0, k) = br
        sum = sum + br
      enddo
      if (sum > 0.) then
        do k = 1, nbr
          branchratio(Zix, Nix, ilev0, k) = branchratio(Zix, Nix, ilev0, k) / sum
        enddo
      endif
      nbranch(Zix, Nix, ilev0) = nbr
      cycle
    endif
    if (key == 'disctable') then
      read(value, * , iostat = istat) disctable
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'bestbranch') then
      if (ch == 'n') flagbestbr = .false.
      if (ch == 'y') flagbestbr = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'electronconv') then
      if (ch == 'n') flagelectron = .false.
      if (ch == 'y') flagelectron = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'pseudoresonances') then
      if (ch == 'n') flagpseudores = .false.
      if (ch == 'y') flagpseudores = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'levelfile') then
      class = 10
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) levelfile(Zix) = cval
      cycle
    endif
    if (key == 'deformfile') then
      class = 10
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) deformfile(Zix) = cval
      cycle
    endif
    if (key == 'outlevels') then
      if (ch == 'n') flaglevels = .false.
      if (ch == 'y') flaglevels = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'nlevels') then
      class = 2
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) nlev(Zix, Nix) = ival
      cycle
    endif
    if (key == 'maxlevelsbin') then
      class = 7
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) nlevbin(type) = ival
      cycle
    endif
    if (key == 'maxlevelstar') then
      read(value, * , iostat = istat) nlevmax
      if (istat /= 0) call read_error(line, istat)
      nlevmax = max(nlevmax, Ltarget)
      nlevbin(k0) = nlevmax
      cycle
    endif
    if (key == 'maxlevelsres') then
      read(value, * , iostat = istat) nlevmaxres
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'risomer') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Risomer(Zix, Nix) = val
      cycle
    endif
  enddo
  return
end subroutine input_levels
! Copyright A.J. Koning 2021
