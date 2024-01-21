subroutine input_omppar
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
!   adepthcor        ! adjustable parameter for depth of DF alpha potential
!   aradialcor       ! adjustable parameter for shape of DF alpha potential
!   avadjust         ! adjustable factor for OMP (default 1.)
!   avdadjust        ! adjustable factor for OMP (default 1.)
!   avsoadjust       ! adjustable factor for OMP (default 1.)
!   awadjust         ! adjustable factor for OMP (default 1.)
!   awdadjust        ! adjustable factor for OMP (default 1.)
!   awsoadjust       ! adjustable factor for OMP (default 1.)
!   d1adjust         ! adjustable factor for OMP (default 1.)
!   d2adjust         ! adjustable factor for OMP (default 1.)
!   d3adjust         ! adjustable factor for OMP (default 1.)
!   Ejoin            ! joining energy for high energy OMP
!   flagjlm          ! flag for using semi - microscopic JLM OMP
!   flagsoukho       ! flag for Soukhovitskii OMP for actinides
!   flagspher        ! flag to force spherical optical model
!   lv1adjust        ! adjustable parameter for JLM OMP
!   lvadjust         ! adjustable parameter for JLM OMP
!   lvsoadjust       ! adjustable parameter for JLM OMP
!   lw1adjust        ! adjustable parameter for JLM OMP
!   lwadjust         ! adjustable parameter for JLM OMP
!   lwsoadjust       ! adjustable parameter for JLM OMP
!   ompadjustD       ! depth of local OMP adjustment
!   ompadjustE1      ! start energy of local OMP adjustment
!   ompadjustE2      ! end energy of local OMP adjustment
!   ompadjustF       ! logical for local OMP adjustment
!   ompadjustN       ! number of energy ranges for local OMP adjustment
!   ompadjustp       ! flag for local optical model parameter adjustment
!   ompadjusts       ! variance of local OMP adjustment
!   RprimeU          ! potential scattering radius
!   rcadjust         ! adjustable factor for OMP (default 1.)
!   riplomp          ! RIPL OMP
!   rvadjust         ! adjustable factor for OMP (default 1.)
!   rvdadjust        ! adjustable factor for OMP (default 1.)
!   rvsoadjust       ! adjustable factor for OMP (default 1.)
!   rwadjust         ! adjustable factor for OMP (default 1.)
!   rwdadjust        ! adjustable factor for OMP (default 1.)
!   rwsoadjust       ! adjustable factor for OMP (default 1.)
!   v1adjust         ! adjustable factor for OMP (default 1.)
!   v2adjust         ! adjustable factor for OMP (default 1.)
!   v3adjust         ! adjustable factor for OMP (default 1.)
!   v4adjust         ! adjustable factor for OMP (default 1.)
!   Vinfadjust       ! adj. factor for high energy limit of real central potential
!   vso1adjust       ! adjustable factor for OMP (default 1.)
!   vso2adjust       ! adjustable factor for OMP (default 1.)
!   w1adjust         ! adjustable factor for OMP (default 1.)
!   w2adjust         ! adjustable factor for OMP (default 1.)
!   w3adjust         ! adjustable factor for OMP (default 1.)
!   w4adjust         ! adjustable factor for OMP (default 1.)
!   wso1adjust       ! adjustable factor for OMP (default 1.)
!   wso2adjust       ! adjustable factor for OMP (default 1.)
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
  integer            :: ibar                 ! fission barrier
  integer            :: igr                  ! giant resonance
  integer            :: irad                 ! variable to indicate M(=0) or E(=1) radiation
  integer            :: istat                ! logical for file access
  integer            :: ival                 ! integer value
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
  RprimeU = 0.
  v1adjust = 1.
  v2adjust = 1.
  v3adjust = 1.
  v4adjust = 1.
  rvadjust = 1.
  avadjust = 1.
  w1adjust = 1.
  w2adjust = 1.
  w3adjust = 1.
  w4adjust = 1.
  rwadjust = -1.
  awadjust = -1.
  rvdadjust = -1.
  avdadjust = -1.
  d1adjust = 1.
  d2adjust = 1.
  d3adjust = 1.
  rwdadjust = 1.
  awdadjust = 1.
  vso1adjust = 1.
  vso2adjust = 1.
  rvsoadjust = 1.
  avsoadjust = 1.
  wso1adjust = 1.
  wso2adjust = 1.
  rwsoadjust = -1.
  awsoadjust = -1.
  rcadjust = 1.
  Ejoin = 200.
  Vinfadjust = 1.
  ompadjustN = 0
  ompadjustE1 = 0.
  ompadjustE2 = 0.
  ompadjustD = 1.
  ompadjusts = 1.
  ompadjustF = .false.
  ompadjustp = .false.
  lvadjust = 1.
  lwadjust = 1.
  lv1adjust = 1.
  lw1adjust = 1.
  lvsoadjust = 1.
  lwsoadjust = 1.
  aradialcor = 1.
  adepthcor = 1.
  if (flagjlm) then
    flagspher = .true.
  else
    flagspher = .false.
  endif
  riplomp = 0
  flagriplomp= .false.
  if (.not.flagsoukhoinp .and. Atarget > fislim) then
    if ((Ztarget >= 90 .and. Ztarget <= 97 .and. Atarget >= 228 .and. Atarget <= 249) .or. flagriplrisk) then
      flagriplrisk = .true.
      flagriplomp = .true.
      flagsoukho = .false.
      riplomp(1) = 2408
    endif
  endif
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
    if (key == 'spherical') then
      if (ch == 'n') flagspher = .false.
      if (ch == 'y') flagspher = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'soukho') then
      if (ch == 'n') flagsoukho = .false.
      if (ch == 'y') flagsoukho = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'rprime')  then
      read(value, * , iostat = istat) RprimeU
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'v1adjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        v1adjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'v2adjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        v2adjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'v3adjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        v3adjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'v4adjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        v4adjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'rvadjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        rvadjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'avadjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        avadjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'w1adjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        w1adjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'w2adjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        w2adjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'w3adjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        w3adjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'w4adjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        w4adjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'rwadjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        rwadjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'awadjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        awadjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'rvdadjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        rvdadjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'avdadjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) avdadjust(type) = val
      ompadjustp(type) = .true.
      cycle
    endif
    if (key == 'rwdadjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        rwdadjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'awdadjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        awdadjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'd1adjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        d1adjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'd2adjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        d2adjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'd3adjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        d3adjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'rvsoadjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        rvsoadjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'avsoadjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        avsoadjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'vso1adjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        vso1adjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'vso2adjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        vso2adjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'rwsoadjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        rwsoadjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'awsoadjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        awsoadjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'wso1adjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        wso1adjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'wso2adjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        wso2adjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'rcadjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        rcadjust(type) = val
        ompadjustp(type) = .true.
      endif
      cycle
    endif
    if (key == 'riplomp') then
      if (ch == 'n' .and. trim(word(3)) == '') then
        riplomp = 0
      else
        class = 7
        call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
        if (flagassign) riplomp(type) = ival
        flagriplomp = .true.
      endif
      cycle
    endif
    if (key == 'ejoin') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Ejoin(type) = val
      cycle
    endif
    if (key == 'vinfadjust') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) Vinfadjust(type) = val
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
    if (key == 'lvadjust') then
      class = 9
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        lvadjust = val
        ompadjustp(1) = .true.
      endif
      cycle
    endif
    if (key == 'lwadjust') then
      class = 9
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        lwadjust = val
        ompadjustp(1) = .true.
      endif
      cycle
    endif
    if (key == 'lv1adjust') then
      class = 9
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        lv1adjust = val
        ompadjustp(1) = .true.
      endif
      cycle
    endif
    if (key == 'lw1adjust') then
      class = 9
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        lw1adjust = val
        ompadjustp(1) = .true.
      endif
      cycle
    endif
    if (key == 'lvsoadjust') then
      class = 9
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        lvsoadjust = val
        ompadjustp(1) = .true.
      endif
      cycle
    endif
    if (key == 'lwsoadjust') then
      class = 9
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        lwsoadjust = val
        ompadjustp(1) = .true.
      endif
      cycle
    endif
    if (key == 'aradialcor') then
      class = 9
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        aradialcor = val
        ompadjustp(6) = .true.
      endif
      cycle
    endif
    if (key == 'adepthcor') then
      class = 9
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        adepthcor = val
        ompadjustp(6) = .true.
      endif
      cycle
    endif
  enddo Loop1
!
! Apply consistent OMP adjustment factors for Koning-Delaroche and other potentials.
!
  do type=1,6
    if ((type == 3 .and. deuteronomp >= 2) .or. (type == 6 .and. alphaomp >= 2)) then
      if (rwadjust(type) == -1.) rwadjust(type) = 1.
      if (awadjust(type) == -1.) awadjust(type) = 1.
      if (rvdadjust(type) == -1.) rvdadjust(type) = 1.
      if (avdadjust(type) == -1.) avdadjust(type) = 1.
      if (rwsoadjust(type) == -1.) rwsoadjust(type) = 1.
      if (awsoadjust(type) == -1.) awsoadjust(type) = 1.
    else
      if (rwadjust(type) == -1.) rwadjust(type) = rvadjust(type)
      if (awadjust(type) == -1.) awadjust(type) = avadjust(type)
      if (rvdadjust(type) == -1.) rvdadjust(type) = rwdadjust(type)
      if (avdadjust(type) == -1.) avdadjust(type) = awdadjust(type)
      if (rwsoadjust(type) == -1.) rwsoadjust(type) = rvsoadjust(type)
      if (awsoadjust(type) == -1.) awsoadjust(type) = avsoadjust(type)
    endif
  enddo
  return
end subroutine input_omppar
! Copyright A.J. Koning 2021
