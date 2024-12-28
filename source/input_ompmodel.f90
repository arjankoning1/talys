subroutine input_ompmodel
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for OMP variables
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
! Variables for input 5
!   alphaomp         ! alpha optical model
!   altomp           ! flag for alternative optical model
!   deuteronomp      ! deuteron optical model
!   flagdisp         ! flag for dispersive optical model
!   flagincadj       ! flag for OMP adjustment on incident channel also
!   flagsoukho       ! flag for Soukhovitskii OMP for actinides
!   flaglocalomp     ! flag for local (y) or global (n) optical model
!   flaglocalomp     ! flag for local (y) or global (n) optical model
!   flagompall       ! flag for new optical model calculation for all residual
!   flagomponly      ! flag to execute ONLY an optical model calculation
!   flagoutkd        ! flag for output of KD03 OMP parameters
!   flagoutomp       ! flag for output of optical model parameters
!   flagjlm          ! flag for using semi - microscopic JLM OMP
!   flagriplrisk     ! flag for going outside RIPL mass validity range
!   jlmmode          ! option for JLM imaginary potential normalization
!   optmodfileN      ! optical model parameter file for neutrons
!   optmodfileP      ! optical model parameter file for protons
!   optmod           ! file with optical model parameters
!   radialfile       ! radial matter density file
!   radialmodel      ! model for radial matter densities (JLM OMP only)
!   riplomp          ! RIPL OMP
! Constants
!   parsym            ! symbol of particle
! Variables for reading input lines
!   inline            ! input line
!   nlines            ! number of input lines
! Variables for main input
!   Ninit                 ! neutron number of initial compound nucleus
!   Zinit                 ! charge number of initial compound nucleus
! All global variables
!   numNph       ! maximum number of protons away from the initial compound nucleus
!   numompadj    ! number of adjustable ranges for OMP
!   numpar       ! number of particles
!   numrange     ! number of energy ranges
!   numZ         ! maximum number of protons from initial compound nucleus
!   numZph       ! maximum number of protons away from the initial compound nucleus
! Constants
!   fislim         ! mass above which nuclide fissions
! Variables for basic reaction
!   flagbasic    ! flag for output of basic information and results
!   flagmicro    ! flag for completely microscopic TALYS calculation
! Variables for main input
!   Atarget          ! mass number of target nucleus
!   Ztarget          ! charge number of target nucleus
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
  integer            :: ia                   ! mass number from abundance table
  integer            :: ibar                 ! fission barrier
  integer            :: igr                  ! giant resonance
  integer            :: irad                 ! variable to indicate M(=0) or E(=1) radiation
  integer            :: istat                ! logical for file access
  integer            :: ival                 ! integer value
  integer            :: iz                   ! charge number of residual nucleus
  integer            :: lval                 ! multipolarity
  integer            :: Nix                  ! neutron number index for residual nucleus
  integer            :: nr                   ! number of radial grid point
  integer            :: omptype              ! type of optical model (spherical or coupled)
  integer            :: type                 ! particle type
  integer            :: type2                ! particle type
  integer            :: Zix                  ! charge number index for residual nucleus
  real(sgl)          :: val                  ! real value
!
! ********************************** Defaults **************************
!
  flagdisp = .false.
  flagincadj = .true.
  if (flagmicro) then
    flagjlm = .true.
  else
    flagjlm = .false.
  endif
  pruitt = 'n'
  flaglocalomp = .true.
  flagompall = .false.
  flagomponly = .false.
  flagoutomp = flagbasic
  flagoutkd = .false.
  flagriplrisk = .false.
  flagsoukho = .true.
  flagsoukhoinp = .false.
  optmod = ' '
  optmodfileN = ' '
  optmodfileP = ' '
  radialfile = ' '
  radialmodel = 2
  jlmmode = 0
  alphaomp = 6
  deuteronomp = 4
  altomp = .false.
  altomp(6) = .true.
!
! *********************** Read input variables *************************
!
! getkeywords: subroutine to retrieve keywords and values from input line
!
! The keyword is identified and the corresponding values are read.
! Erroneous input is immediately checked.
! The keywords and number of values on each line are retrieved from the input.
!
loop1:  do i = 1, nlines
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
    if (key == 'optmodfilen') then
      class = 10
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) optmodfilen(Zix) = cval
      cycle
    endif
    if (key == 'optmodfilep') then
      class = 10
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) optmodfilep(Zix) = cval
      cycle
    endif
    if (key == 'optmod') then
      read(word(2), * , iostat = istat) iz
      if (istat /= 0) call read_error(line, istat)
      read(word(3), * , iostat = istat) ia
      if (istat /= 0) call read_error(line, istat)
      Zix = Zinit - iz
      Nix = Ninit - ia + iz
      if (Zix < 0 .or. Zix > numZph .or. Nix < 0 .or. Nix > numNph) then
        write(*, '(" TALYS-warning: Z, N index out of range, keyword ignored: ", a)') trim(line)
        cycle
      else
        ch = word(5)(1:1)
        if (ch == ' ') ch = 'n'
        do type = 1, 6
          if (ch == parsym(type)) then
            optmod(Zix, Nix, type) = word(4)
            cycle
          endif
        enddo
      endif
      cycle
    endif
    if (key == 'radialfile') then
      class = 10
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) radialfile(Zix) = cval
      cycle
    endif
    if (key == 'localomp') then
      if (ch == 'n') flaglocalomp = .false.
      if (ch == 'y') flaglocalomp = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'dispersion') then
      if (ch == 'n') flagdisp = .false.
      if (ch == 'y') flagdisp = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'optmodall') then
      if (ch == 'n') flagompall = .false.
      if (ch == 'y') flagompall = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'jlmomp') then
      if (ch == 'n') flagjlm = .false.
      if (ch == 'y') flagjlm = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'soukho') then
      if (ch == 'n') flagsoukho = .false.
      if (ch == 'y') flagsoukho = .true.
      flagsoukhoinp = flagsoukho
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'riplrisk') then
      if (ch == 'n') flagriplrisk = .false.
      if (ch == 'y') flagriplrisk = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'incadjust') then
      if (ch == 'n') flagincadj = .false.
      if (ch == 'y') flagincadj = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'pruitt') then
      pruitt = ch
      if (ch /= 'y' .and. ch /= 'n' .and. ch /= 'd' .and. ch /= 'f') call read_error(line, istat)
      cycle
    endif
    if (key == 'pruitt') then
      read(value, * , iostat = istat) pruitt
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'rvadjustf' .or. key == 'avadjustf' .or. key == 'rwadjustf' .or. key == 'awadjustf' .or. &
      key == 'rvdadjustf' .or. key == 'avdadjustf' .or. key == 'rwdadjustf' .or. key == 'awdadjustf' .or. &
      key == 'rvsoadjustf' .or. key == 'avsoadjustf' .or. key == 'rwsoadjustf' .or. key == 'awsoadjustf') then
      if (key == 'rvadjustf') omptype = 1
      if (key == 'avadjustf') omptype = 2
      if (key == 'rwadjustf') omptype = 3
      if (key == 'awadjustf') omptype = 4
      if (key == 'rvdadjustf') omptype = 5
      if (key == 'avdadjustf') omptype = 6
      if (key == 'rwdadjustf') omptype = 7
      if (key == 'awdadjustf') omptype = 8
      if (key == 'rvsoadjustf') omptype = 9
      if (key == 'avsoadjustf') omptype = 10
      if (key == 'rwsoadjustf') omptype = 11
      if (key == 'awsoadjustf') omptype = 12
      if (key == 'rcadjustf') omptype = 13
      do type = 1, 6
        if (ch == parsym(type)) then
          type2 = type
          ompadjustF(type2) = .true.
          ompadjustN(type2, omptype) = ompadjustN(type2, omptype) + 1
          nr = ompadjustN(type2, omptype)
          read(word(3), * , iostat = istat) ompadjustE1(type2, omptype, nr)
          if (istat /= 0) call read_error(line, istat)
          read(word(4), * , iostat = istat) ompadjustE2(type2, omptype, nr)
          if (istat /= 0) call read_error(line, istat)
          read(word(5), * , iostat = istat) ompadjustD(type2, omptype, nr)
          if (istat /= 0) call read_error(line, istat)
          read(word(6), * , iostat = istat) ompadjusts(type2, omptype, nr)
          if (istat > 0) cycle
          if (istat /= 0) call read_error(line, istat)
          cycle Loop1
        endif
      enddo
      call read_error(line, istat)
    endif
    if (key == 'jlmmode') then
      read(value, * , iostat = istat) jlmmode
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'radialmodel') then
      read(value, * , iostat = istat) radialmodel
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'omponly') then
      if (ch == 'n') flagomponly = .false.
      if (ch == 'y') flagomponly = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outomp') then
      if (ch == 'n') flagoutomp = .false.
      if (ch == 'y') flagoutomp = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outkd') then
      if (ch == 'n') flagoutkd = .false.
      if (ch == 'y') flagoutkd = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'alphaomp') then
      read(value, * , iostat = istat) alphaomp
      if (istat /= 0) call read_error(line, istat)
      if (alphaomp == 1) then
        altomp(6) = .false.
      else
        altomp(6) = .true.
      endif
      cycle
    endif
    if (key == 'deuteronomp') then
      read(value, * , iostat = istat) deuteronomp
      if (istat /= 0) call read_error(line, istat)
      altomp(3) = .true.
      cycle
    endif
  enddo Loop1
  return
end subroutine input_ompmodel
! Copyright A.J. Koning 2021
