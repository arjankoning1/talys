subroutine getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Assign values to keywords
!
! Author    : Arjan Koning
!
! 2024-12-08: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
  use A1_error_handling_mod
!
! Definition of single and double precision variables
!   sgl           ! single precision kind
! All global variables
!   numadj        ! maximum number of adjustable parameters
!   numbar        ! number of fission barriers
!   numenadj      ! maximum number of energies for adjustable parameters
!   numgam        ! maximum number of l - values for gamma multipolarity
!   numl          ! number of l values
!   numN          ! maximum number of neutrons from initial compound nucleus
!   numZ          ! maximum number of protons from initial compound nucleus
! Variables for adjustment
!   adjustfile    ! file for local adjustment
!   adjustix      ! local adjustment index
!   adjustkey     ! keyword for local adjustment
!   adjustpar     ! local adjustment parameters
!   Dadjust       ! tabulated depth of local adjustment
!   Eadjust       ! tabulated energy of local adjustment
!   Nadjust       ! number of adjustable parameters
!   nenadjust     ! number of tabulated energies of local adjustment
! Variables for main input
!   Ninit         ! neutron number of initial compound nucleus
!   Zinit         ! charge number of initial compound nucleus
! Constants
!   parsym        ! symbol of particle
! Error handling
!   range_index_error    ! Test if index is out of range
!   range_real_error    ! Test if real variable is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: numEkey=87         ! number of keywords
  logical            :: flagassign         ! flag to assign value or not
  logical            :: flagoutside        ! flag for value outside range
  logical            :: lexist             ! logical to determine existence
  character(len=1)   :: ch                 ! character
  character(len=132) :: adfile             ! file with tabulated adjustments
  character(len=132) :: cval               ! character value
  character(len=132) :: key                ! keyword
  character(len=132) :: keyword(numEkey)   ! keyword
  character(len=132) :: word(40)           ! words on input line
  integer            :: class              ! input class
  integer            :: i                  ! level
  integer            :: ia                 ! mass number from abundance table
  integer            :: ibar               ! fission barrier
  integer            :: igr                ! giant resonance
  integer            :: in                 ! counter for neutrons
  integer            :: irad               ! variable to indicate M(=0) or E(=1) radiation
  integer            :: istat              ! logical for file access
  integer            :: ival               ! integer value
  integer            :: iz                 ! charge number of residual nucleus
  integer            :: j                  ! counter
  integer            :: k                  ! designator for particle
  integer            :: lval               ! multipolarity
  integer            :: Nix                ! neutron number index for residual nucleus
  integer            :: type               ! particle type
  integer            :: type2              ! particle type
  integer            :: Zix                ! charge number index for residual nucleus
  real(sgl)          :: D                  ! depth of local adjustment
  real(sgl)          :: Dadj(numenadj)     ! depth of local adjustment
  real(sgl)          :: Ea                 ! start energy of local adjustment
  real(sgl)          :: Eadj(0:numenadj)   ! energy of local adjustment
  real(sgl)          :: Eb                 ! end energy of local adjustment
  real(sgl)          :: Em                 ! intermediate energy of local adjustment
  real(sgl)          :: val                ! real value
!
! Several keywords may be altered over local energy ranges
!
  data (keyword(i), i = 1, numEkey) / 'adepthcor', 'aradialcor', &
 &  'avadjust', 'avdadjust', 'avsoadjust', 'awadjust', 'awdadjust', 'awsoadjust', 'bdamp', 'bdampadjust', &
 &  'betafiscor', 'betafiscoradjust', 'cbreak', 'cfermi', 'cknock', &
 &  'cstrip', 'ctable', 'ctableadjust', 'd1adjust', 'd2adjust', 'd3adjust', 'egr', 'egradjust', 'epr', 'epradjust', 'etable', &
 &  'fisbar', 'fisbaradjust', 'fishw', 'fishwadjust', 'fsadjust', 'ftable', 'ftableadjust', &
 &  'ggr', 'ggradjust', 'gpr', 'gpradjust', 'krotconstant', 'lv1adjust', 'lvadjust', 'lvsoadjust', &
 &  'lw1adjust', 'lwadjust', 'lwsoadjust', 'm2constant', 'ptable', 'ptableadjust', 'rcadjust', 'rspincut', 'rspincutff', &
 &  'rspincutpreeq', 'rvadjust', 'rvdadjust', 'rvsoadjust', 'rwadjust', 'rwdadjust', 'rwsoadjust', 's2adjust', 'sgr', 'sgradjust', &
 &  'spr', 'spradjust', 'tjadjust', 'tmadjust', 'ufermi', 'upbendc', 'upbendcadjust', 'upbende', 'upbendeadjust', 'upbendf', &
 &  'upbendfadjust','v1adjust', 'v2adjust', 'v3adjust', 'v4adjust', 'vfiscor', &
 &  'vfiscoradjust', 'vso1adjust', 'vso2adjust', 'w1adjust', 'w2adjust', 'w3adjust', 'w4adjust', 'wso1adjust', 'wso2adjust', &
 &  'wtable', 'wtableadjust'/
!
! ************************ Read values for keywords ********************
!
! Each keyword is characterized by a certain order of parameter and value input.
! They are distinguished by different classes.
!
! Classes:
!
!  0: keyword Z A
!  1: keyword Z A real-value [optional: local adjustment]
!  2: keyword Z A integer-value
!  3: keyword Z A real-value barrier [optional: local adjustment]
!  4: keyword Z A integer-value barrier [optional: local adjustment]
!  5: keyword Z A real-value rad-type l-val [optional: local adjustment]
!  6: keyword particle-type real-value [optional: local adjustment]
!  7: keyword particle-type integer-value
!  8: keyword particle-type real value integer-value
!   [optional: local adjustment]
!  9: keyword value [optional: local adjustment]
! 10: keyword Z filename
! 11: keyword Z A filename
! 12: keyword Z A particle-type real-value filename
! 13: keyword character-value Z A
! 14: keyword integer-value Z A
! 15: keyword real-value Z A [optional: integer value]
! 16: keyword particle-type Z A real-value
!
  flagoutside = .false.
  flagassign = .false.
  key = word(1)
  ch = word(2)(1:1)
  if (ch == 'e') ch = 'g'
  val = 1.
  ival = -1
  cval = '                                                           '
!
! Z,A dependent keywords
!
  if (class <= 5) then
    read(word(2), * , iostat = istat) iz
    if (istat > 0) call read_error(key, istat, ival = iz)
    read(word(3), * , iostat = istat) ia
    if (istat > 0) call read_error(key, istat, ival = ia)
    if (class >= 1) then
      if (class == 2 .or. class == 4) then
        read(word(4), * , iostat = istat) ival
        if (istat > 0) call read_error(key, istat, ival = ival)
      else
        read(word(4), * , iostat = istat) val
        if (istat > 0) call read_error(key, istat, xval = val)
      endif
    endif
    if (iz <= 2 .and. ia <= 4) then
      Zix = iz
      Nix = ia - iz
      iz = Zinit - Zix
      in = Ninit - Nix
      ia = iz + in
    else
      in = ia - iz
      Zix = Zinit - iz
      Nix = Ninit - in
      ia = iz + in
    endif
    call range_index_error(key, 'Z', iz, max(Zinit - numZ,0), Zinit, error = 'continue', flagoutside = flagoutside)
    if (.not. flagoutside) &
 &    call range_index_error(key, 'N', in, max(Ninit - numN,0), Ninit, error = 'continue', flagoutside = flagoutside)
    if (flagoutside) write(*,'(a,": Z= ",i3," A= ",i3," out of range")') trim(key), iz, ia
    i = 4
!
! Z,A dependent keywords with possible fission barriers
!
    if (class == 0) i = 3
    if (class == 3 .or. class == 4) then
      read(word(5), * , end = 10, err = 40) ibar
      call range_index_error(key, 'barrier', ibar, 0, numbar)
      i = 5
    endif
!
! Z,A dependent keywords with gamma parameters
!
    if (class == 5) then
      ch = word(5)(1:1)
      irad = 1
      lval = 1
      igr = 1
      if (ch == 'm' .or. ch == 'e') then
        if (ch == 'm') irad = 0
        if (ch == 'e') irad = 1
        read(word(5)(2:2), * , iostat = istat) lval
        if (istat /= 0) call read_error(key, istat, ival = lval)
        call range_index_error(key, 'l-value', lval, 1, numgam)
        read(word(6), * , end = 10, err = 100) igr
        call range_index_error(key, 'igr', igr, 1, 2)
        i = 6
      else
        read(word(5), * , iostat = istat) lval
        if (istat /= 0) call read_error(key, istat, ival = lval)
        call range_index_error(key, 'l-value', lval, 0, numl)
        i = 5
      endif
    endif
  endif
!
! Particle type dependent keywords
!
  if (class >= 6 .and. class <= 8) then
    do type2 = 0, 6
      if (ch == parsym(type2)) then
        type = type2
        goto 30
      endif
    enddo
    goto 180
   30   if (class == 7) then
      read(word(3), * , iostat = istat) ival
      if (istat > 0) call read_error(key, istat, ival = ival)
    else
      read(word(3), * , iostat = istat) val
      if (istat > 0) call read_error(key, istat, xval = val)
    endif
    i = 3
    if (class == 8) then
      read(word(4), * , end = 10, err = 40) lval
      call range_index_error(key, 'l-value', lval, 0, numl)
      i = 4
    endif
  endif
!
! Simple dependent keywords
!
  if (class == 9) then
    read(word(2), * , iostat = istat) val
    if (istat /= 0) call read_error(key, istat, xval = val)
    i = 2
  endif
!
! Keywords with input files
!
  if (class >= 10 .and. class <= 12) then
    read(word(2), * , iostat = istat) iz
    if (istat /= 0) call read_error(key, istat, ival = iz)
    Zix = Zinit - iz
    if (Zix < 0 .or. Zix > numZ) then
      goto 200
    else
      if (class >= 11) then
        read(word(3), * , iostat = istat) ia
        if (istat /= 0) call read_error(key, istat, ival = ia)
        in = ia - iz
        Nix = Ninit - in
        if (Nix < 0 .or. Nix > numN) then
          goto 200
        else
          if (class == 12) then
            read(word(4), * , iostat = istat) ch
            if (istat /= 0) call read_error(key, istat, cval = ch)
            do type2 = -1, 6
              if (ch == parsym(type2)) then
                type = type2
                goto 34
              endif
            enddo
            goto 180
 34         read(word(5), * , iostat = istat) val
            if (istat /= 0) call read_error(key, istat, xval = val)
            i = 5
          else
            cval = word(4)
            i = 4
          endif
        endif
      else
        cval = word(3)
        i = 3
      endif
    endif
  endif
  if (class >= 13) then
    if (class == 13) then
      read(word(2), * , iostat = istat) cval
      if (istat /= 0) call read_error(key, istat, cval = cval)
    endif
    if (class == 14) then
      read(word(2), * , iostat = istat) ival
      if (istat /= 0) call read_error(key, istat, ival = ival)
    endif
    if (class == 15) then
      read(word(2), * , iostat = istat) val
      if (istat /= 0) call read_error(key, istat, xval = val)
      read(word(5), * , iostat = istat) ival
      if (istat > 0) call read_error(key, istat, ival = ival)
    endif
    if (class == 16) then
      do type2 = 0, 6
        if (ch == parsym(type2)) then
          type = type2
          goto 330
        endif
      enddo
      goto 180
330   read(word(5), * , iostat = istat) val
      if (istat /= 0) call read_error(key, istat, xval = val)
    endif
    read(word(3), * , iostat = istat) iz
    if (istat > 0) call read_error(key, istat, ival = iz)
    read(word(4), * , iostat = istat) ia
    if (istat > 0) call read_error(key, istat, ival = ia)
    Zix = Zinit - iz
    in = ia - iz
    Nix = Ninit - in
    call range_index_error(key, 'Z', iz, Zinit - numZ, Zinit, error = 'continue', flagoutside = flagoutside)
    if (.not. flagoutside) &
 &    call range_index_error(key, 'N', in, Ninit - numN, Ninit, error = 'continue', flagoutside = flagoutside)
    if (flagoutside) write(*,'("Z= ",i3," A= ",i3," out of range")')  iz, ia
  endif
!
! Local energy-dependent adjustment of parameters
!
   40 Ea = 0.
  Eb = 0.
  Em = 0.
  D = 0.
  adfile = '                                                         '
  do k = 1, numEkey
    if (key == keyword(k)) goto 60
  enddo
  goto 10
   60 Eadj(0) = 0.
  k = 1
  if ((word(i+1)(1:1) >= 'a' .and. word(i+1)(1:1) <= 'z') .or. (word(i+1)(1:1) >= 'A' .and. word(i+1)(1:1) <= 'Z')) then
    read(word(i + 1), * , end = 10, err = 100) adfile
    inquire (file = adfile, exist = lexist)
    if ( .not. lexist) goto 140
    open (unit = 1, file = adfile, status = 'old')
   70   read(1, * , end = 80, err = 150) Eadj(k), Dadj(k)
    if (Eadj(k) == Eadj(k - 1)) goto 70
    if (Eadj(k) < Eadj(k - 1)) goto 160
    k = k + 1
    if (k > numenadj) goto 170
    goto 70
   80   close (1)
  else
    read(word(i + 1), * , end = 10, err = 10) Ea
    read(word(i + 2), * , end = 10, err = 10) Eb
    read(word(i + 3), * , end = 10, err = 10) Em
    read(word(i + 4), * , end = 10, err = 10) D
  endif
  if (Nadjust < numadj) Nadjust = Nadjust + 1
  nenadjust(Nadjust) = k - 1
  do j = 1, nenadjust(Nadjust)
    Eadjust(Nadjust, j) = Eadj(j)
    Dadjust(Nadjust, j) = val * Dadj(j)
  enddo
  adjustkey(Nadjust) = key
  adjustfile(Nadjust) = adfile
  adjustix(Nadjust, 1) = Zix
  adjustix(Nadjust, 2) = Nix
  adjustix(Nadjust, 3) = type
  if (class == 5) then
    adjustix(Nadjust, 4) = lval
  else
    adjustix(Nadjust, 4) = ibar
  endif
  adjustpar(Nadjust, 1) = Ea
  adjustpar(Nadjust, 2) = Eb
  adjustpar(Nadjust, 3) = Em
  adjustpar(Nadjust, 4) = D
   10 if (.not. flagoutside) flagassign = .true.
  return
!
! Error and warning messages
!
  100 write(*, '(" TALYS-error: Wrong input for: ", a)') trim(key)
  stop
  140 write(*, '(" TALYS-error: parameter file ", a, " does not exist for keyword ", a)') trim(adfile), trim(key)
  stop
  150 write(*, '(" TALYS-error: parameter file ", a, " has wrong format for keyword ", a)') trim(adfile), trim(key)
  stop
  160 write(*, '(" TALYS-error: parameter file ", a, " must have energies in increasing order for keyword ", a)') &
 &  trim(adfile), trim(key)
  stop
  170 write(*, '(" TALYS-error: parameter file ", a, " has more than ", i6, " energies for keyword ", a)') &
 &  trim(adfile), numenadj, trim(key)
  stop
  180 write(*, '(" TALYS-error: wrong particle symbol for: ", a)') trim(key)
  stop
  200 write(*, '(" TALYS-warning: Z, N index out of range, keyword ignored: ", a)') trim(key)
  return
end subroutine getvalues
! Copyright A.J. Koning 2021
