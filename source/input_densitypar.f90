subroutine input_densitypar
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
  use A0_talys_mod
  use A1_error_handling_mod
!
! Definition of single and double precision variables
!   sgl               ! single precision kind
! Variables for input of level density parameters
!   aadjust              ! adjustable factor for level density parameter
!   alev                 ! level density parameter
!   alimit               ! asymptotic level density parameter
!   alphald              ! alpha - constant for asymptotic level density parameter
!   alphaldall           ! variable for level density
!   betald               ! beta - constant for asymptotic level density parameter
!   betaldall            ! variable for level density
!   cfermi               ! width of Fermi distribution for damping
!   cglobal              ! global constant to adjust tabulated level densities
!   col                  ! flag for collective enhancement of level density
!   ctable               ! constant to adjust tabulated level densities
!   ctableadjust         ! adjustable correction to adjust tabulated level densities
!   D0                   ! s - wave resonance spacing in eV
!   deltaW               ! shell correction in nuclear mass
!   E0                   ! particle constant of temperature formula
!   E0adjust             ! adjustable factor for E0
!   Exmatch              ! matching point for Ex
!   Exmatchadjust        ! adjustable factor for matching energy
!   filedensity          ! flag for level densities on separate files
!   flagasys             ! flag for all level density parameters a from systematic
!   flagcol              ! flag for collective enhancement of level density
!   flagcolall           ! flag for collective enhancement of level density
!   flagcolldamp         ! flag for damping of coll. effects in eff. level density (without explicit coll. enh.)
!   flagctmglob          ! flag for global CTM model (no discrete level info)
!   flagdensity          ! flag for output of level densities
!   flagparity           ! flag for non - equal parity distribution
!   gammald              ! gamma - constant for asymptotic level density parameter
!   gammashell1          ! gamma - constant for asymptotic level density parameter
!   gammashell1all       ! variable for level density
!   gammashell2          ! gamma - constant for asymptotic level density parameter
!   Krotconstant         ! normalization constant for rotational enhancement
!   ldadjust             ! logical for energy - dependent level density adjustment
!   ldmodel              ! level density model
!   ldmodelall           ! level density model for all nuclides
!   Nlow                 ! lowest discrete level for temperature matching
!   Ntop                 ! highest discrete level for temperature matching
!   pair                 ! pairing energy
!   pairconstant         ! constant for pairing energy systematics
!   pglobal              ! global constant to adjust tabulated level densities
!   Pshift               ! adjustable pairing shift
!   Pshiftadjust         ! adjustable correction to pairing shift
!   Pshiftconstant       ! global constant for pairing shift
!   Pshiftconstantall    ! variable for level density
!   ptable               ! constant to adjust tabulated level densities
!   ptableadjust         ! adjustable correction to adjust tabulated level densities
!   Rclass2mom           ! norm. constant for moment of inertia for class 2 states
!   Rspincut             ! adjustable constant (global) for spin cutoff factor
!   Rtransmom            ! norm. constant for moment of inertia for transition states
!   s2adjust             ! adjustable constant (Z, A, barrier - dependent) for spin
!   T                    ! temperature
!   Tadjust              ! adjustable factor for temperature
!   Ufermi               ! energy of Fermi distribution for damping
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
  ldadjust = .false.
  densfile = ' '
  do Zix = 0, numZ
    do Nix = 0, numN
      if (ldmodel(Zix, Nix) == 1 .or. ldmodel(Zix, Nix) >= 4) then
        if (flagcol(Zix, Nix)) then
          alphald(Zix, Nix) = 0.0207305
          betald(Zix, Nix) = 0.229537
          gammashell1(Zix, Nix) = 0.473625
          Pshiftconstant(Zix, Nix) = 0.
        else
          alphald(Zix, Nix) = 0.0692559
          betald(Zix, Nix) = 0.282769
          gammashell1(Zix, Nix) = 0.433090
          Pshiftconstant(Zix, Nix) = 0.
        endif
        if (flagcolldamp) then
          alphald(Zix, Nix) = 0.0666
          betald(Zix, Nix) = 0.258
          gammashell1(Zix, Nix) = 0.459
          Pshiftconstant(Zix, Nix) = 0.
        endif
      endif
      if (ldmodel(Zix, Nix) == 2) then
        if (flagcol(Zix, Nix)) then
          alphald(Zix, Nix) = 0.0381563
          betald(Zix, Nix) = 0.105378
          gammashell1(Zix, Nix) = 0.546474
          Pshiftconstant(Zix, Nix) = 0.743229
        else
          alphald(Zix, Nix) = 0.0722396
          betald(Zix, Nix) = 0.195267
          gammashell1(Zix, Nix) = 0.410289
          Pshiftconstant(Zix, Nix) = 0.173015
        endif
      endif
      if (ldmodel(Zix, Nix) == 3) then
        if (flagcol(Zix, Nix)) then
          alphald(Zix, Nix) = 0.0357750
          betald(Zix, Nix) = 0.135307
          gammashell1(Zix, Nix) = 0.699663
          Pshiftconstant(Zix, Nix) = - 0.149106
        else
          alphald(Zix, Nix) = 0.110575
          betald(Zix, Nix) = 0.0313662
          gammashell1(Zix, Nix) = 0.648723
          Pshiftconstant(Zix, Nix) = 1.13208
        endif
      endif
    enddo
  enddo
  aadjust = 1.
  alev = 0.
  alimit = 0.
  ctableadjust = 0.
  ptableadjust = 0.
  deltaW = 0.
  D0 = 0.
  E0 = 1.e-20
  E0adjust = 1.
  Exmatch = 0.
  Exmatchadjust = 1.
  if (ldmodelall >= 5) then
    flagparity = .true.
  else
    flagparity = .false.
  endif
  alphaldall=-99.
  betaldall=-99.
  gammashell1all=-99.
  Pshiftconstantall=-99.
  gammald = -1.
  gammashell2 = 0.
  Krotconstant = 1.
  Nlow = - 1
  Ntop = - 1
  pair = 1.e-20
  pairconstant = 12.
  Pshift = 1.e-20
  Pshiftadjust = 0.
  Rclass2mom = 1.
  Rspincut = 1.
  Rtransmom = 1.
  do Zix = 0, numZ
    do Nix = 0, numN
      Rtransmom(Zix, Nix, 1) = 0.6
      Ufermi(Zix, Nix, 0) = 30.
      cfermi(Zix, Nix, 0) = 5.
      do ibar = 1, numbar
        Ufermi(Zix, Nix, ibar) = 45.
        cfermi(Zix, Nix, ibar) = 5.
      enddo
    enddo
  enddo
  s2adjust = 1.
  T = 0.
  Tadjust = 1.
  ctable = cglobal
  ptable = pglobal
  filedensity = flagdensity
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
    if (key == 'a') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) alev(Zix, Nix) = val
      cycle
    endif
    if (key == 'aadjust') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) aadjust(Zix, Nix) = val
      cycle
    endif
    if (key == 'alimit') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) alimit(Zix, Nix) = val
      cycle
    endif
    if (key == 'gammald') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) gammald(Zix, Nix) = val
      cycle
    endif
    if (key == 'pair') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) pair(Zix, Nix) = val
      cycle
    endif
    if (key == 'pshift') then
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) pshift(Zix, Nix, ibar) = val
      cycle
    endif
    if (key == 'pshiftadjust') then
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) pshiftadjust(Zix, Nix, ibar) = val
      cycle
    endif
    if (key == 'deltaw') then
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) deltaW(Zix, Nix, ibar) = val
      cycle
    endif
    if (key == 'exmatch') then
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Exmatch(Zix, Nix, ibar) = val
      cycle
    endif
    if (key == 't') then
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) T(Zix, Nix, ibar) = val
      cycle
    endif
    if (key == 'e0') then
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) E0(Zix, Nix, ibar) = val
      cycle
    endif
    if (key == 'd0') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) D0(Zix, Nix) = val * 1000.
      cycle
    endif
    if (key == 'nlow') then
      class = 4
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Nlow(Zix, Nix, ibar) = ival
      cycle
    endif
    if (key == 'ntop') then
      class = 4
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Ntop(Zix, Nix, ibar) = ival
      cycle
    endif
    if (key == 's2adjust') then
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) s2adjust(Zix, Nix, ibar) = val
      cycle
    endif
    if (key == 'krotconstant') then
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        Krotconstant(Zix, Nix, ibar) = val
        ldadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'ctable') then
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        ctable(Zix, Nix, ibar) = val
        ldadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'ptable') then
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) ptable(Zix, Nix, ibar) = val
      cycle
    endif
    if (key == 'ctableadjust') then
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        ctableadjust(Zix, Nix, ibar) = val
        ldadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'ptableadjust') then
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) ptableadjust(Zix, Nix, ibar) = val
      cycle
    endif
    if (key == 'tadjust') then
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Tadjust(Zix, Nix, ibar) = val
      cycle
    endif
    if (key == 'e0adjust') then
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) E0adjust(Zix, Nix, ibar) = val
      cycle
    endif
    if (key == 'exmatchadjust') then
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Exmatchadjust(Zix, Nix, ibar) = val
      cycle
    endif
    if (key == 'rtransmom') then
      ibar = 1
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Rtransmom(Zix, Nix, ibar) = val
      cycle
    endif
    if (key == 'rclass2mom') then
      ibar = 1
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Rclass2mom(Zix, Nix, ibar) = val
      cycle
    endif
    if (key == 'densfile') then
      class = 11
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) densfile(Zix, Nix) = cval
      cycle
    endif
    if (key == 'rspincut') then
      class = 9
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Rspincut = val
      cycle
    endif
    if (key == 'alphald') then
      read(value, * , iostat = istat) alphaldall
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'betald') then
      read(value, * , iostat = istat) betaldall
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'gammashell1') then
      read(value, * , iostat = istat) gammashell1all
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'gammashell2') then
      read(value, * , iostat = istat) gammashell2
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'pairconstant') then
      read(value, * , iostat = istat) pairconstant
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'pshiftconstant') then
      read(value, * , iostat = istat) Pshiftconstantall
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'ufermi') then
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        Ufermi(Zix, Nix, ibar) = val
        ldadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'cfermi') then
      class = 3
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        cfermi(Zix, Nix, ibar) = val
        ldadjust(Zix, Nix) = .true.
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
        if (flagassign) flagcol(Zix, Nix) = fcol
      endif
      cycle
    endif
    if (key == 'parity') then
      if (ch == 'n') flagparity = .false.
      if (ch == 'y') flagparity = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'filedensity') then
      if (ch == 'n') filedensity = .false.
      if (ch == 'y') filedensity = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
  enddo
!
! Set level density parameters per nucleus
!
  if (alphaldall /= -99.) alphald = alphaldall
  if (betaldall /= -99.) betald = betaldall
  if (gammashell1all /= -99.) gammashell1 = gammashell1all
  if (Pshiftconstantall /= -99.) Pshiftconstant = Pshiftconstantall
  return
end subroutine input_densitypar
! Copyright A.J. Koning 2021
