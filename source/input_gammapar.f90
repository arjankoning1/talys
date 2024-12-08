subroutine input_gammapar
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for gamma parameters
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
! Variables for gamma rays
!   egr               ! energy of GR
!   egradjust         ! adjustable factor for energy of GR
!   epr               ! energy of PR
!   epradjust         ! adjustable factor for energy of PR
!   etable            ! constant to adjust tabulated strength functions
!   etableadjust      ! adjustable correction to adjust tabulated strength functions
!   filepsf           ! flag for photon strength functions on separate files
!   fiso              ! correction factor for isospin forbidden transitions
!   fisom             ! correction factor for isospin forbidden transitions for multiple emission
!   ftable            ! constant to adjust tabulated strength functions
!   ftableadjust      ! adjustable correction to adjust tabulated strength functions
!   gamadjust         ! logical for energy-dependent gamma adjustment
!   gamgam            ! total radiative width in eV
!   gamgamadjust      ! adjustable factor for radiative parameter
!   gammax            ! number of l - values for gamma multipolarity
!   ggr               ! width of GR
!   ggradjust         ! adjustable factor for width of GR
!   gpr               ! width of PR
!   gpradjust         ! adjustable factor for width of PR
!   sfexpall          ! variable for spectrocopic factor
!   sfthall           ! variable for spectrocopic factor
!   sgr               ! strength of GR
!   sgradjust         ! adjustable factor for strength of GR
!   spectfacexp       ! experimental spectroscopic factor
!   spectfacth        ! theoretical spectroscopic factor
!   tpr               ! strength of PR
!   tpradjust         ! adjustable factor for strength of PR
!   wtable            ! constant to adjust tabulated strength functions
!   wtableadjust      ! adjustable correction to adjust tabulated strength functions
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
  logical            :: wtadjust             ! flag for wtable adjustment 
  character(len=132) :: cval                 ! character value
  character(len=132) :: key                  ! keyword
  character(len=132) :: value                ! value or string
  character(len=132) :: word(40)             ! words on input line
  character(len=132) :: line                 ! input line
  integer            :: class                ! input class
  integer            :: i                    ! counter
  integer            :: l                    ! counter
  integer            :: ibar                 ! fission barrier
  integer            :: igr                  ! giant resonance
  integer            :: irad                 ! variable to indicate M(=0) or E(=1) radiation
  integer            :: istat                ! logical for file access
  integer            :: ival                 ! integer value
  integer            :: lval                 ! multipolarity
  integer            :: nex                  ! counter
  integer            :: Nix                  ! neutron number index for residual nucleus
  integer            :: type                 ! particle type
  integer            :: Zix                  ! charge number index for residual nucleus
  real(sgl)          :: sfexp                ! variable for spectrocopic factor
  real(sgl)          :: sfth                 ! variable for spectrocopic factor
  real(sgl)          :: val                  ! real value
!
! ************** Defaults *************
!
  egr = 0.
  ggr = 0.
  sgr = 0.
  epr = 0.
  gpr = 0.
  tpr = 0.
  egradjust = 1.
  ggradjust = 1.
  sgradjust = 1.
  epradjust = 1.
  gpradjust = 1.
  tpradjust = 1.
  etable = 0.
  ftable = 1.
  wtable = 1.
  wtadjust = .false.
  etableadjust = 0.
  ftableadjust = 1.
  wtableadjust = 1.
  upbend = 0.
  fiso = -1.
  fisom = -1.
  fisominit = 1.
  gamadjust = .false.
  gamgam = 0.
  gamgamadjust = 1.
  gammax = 2
  spectfacexp = 0.
  spectfacth = 0.
  upbendadjust = 1.
  do Zix = 0, numZ
    do Nix = 0, numN
      if (k0 <= 1) then
        if (strength == 8) then
          if (ldmodel(Zix,Nix) == 1 .and. flagcol(Zix, Nix)) wtable(Zix,Nix,1,1) = 1.052
          if (ldmodel(Zix,Nix) == 1 .and. .not.flagcol(Zix, Nix)) wtable(Zix,Nix,1,1) = 1.076
          if (ldmodel(Zix,Nix) == 2 .and. flagcol(Zix, Nix)) wtable(Zix,Nix,1,1) = 0.942
          if (ldmodel(Zix,Nix) == 2 .and. .not.flagcol(Zix, Nix)) wtable(Zix,Nix,1,1) = 0.943
          if (ldmodel(Zix,Nix) == 3 .and. flagcol(Zix, Nix)) wtable(Zix,Nix,1,1) = 0.905
          if (ldmodel(Zix,Nix) == 3 .and. .not.flagcol(Zix, Nix)) wtable(Zix,Nix,1,1) = 0.911
          if (ldmodel(Zix,Nix) == 4) wtable(Zix,Nix,1,1) = 0.918
          if (ldmodel(Zix,Nix) == 5) wtable(Zix,Nix,1,1) = 1.017
          if (ldmodel(Zix,Nix) == 6) wtable(Zix,Nix,1,1) = 0.942
        endif
        if (strength == 9) then
          if (ldmodel(Zix,Nix) == 1 .and. flagcol(Zix, Nix)) wtable(Zix,Nix,1,1) = 1.048
          if (ldmodel(Zix,Nix) == 1 .and. .not.flagcol(Zix, Nix)) wtable(Zix,Nix,1,1) = 1.081
          if (ldmodel(Zix,Nix) == 2 .and. flagcol(Zix, Nix)) wtable(Zix,Nix,1,1) = 0.911
          if (ldmodel(Zix,Nix) == 2 .and. .not.flagcol(Zix, Nix)) wtable(Zix,Nix,1,1) = 0.934
          if (ldmodel(Zix,Nix) == 3 .and. flagcol(Zix, Nix)) wtable(Zix,Nix,1,1) = 0.884
          if (ldmodel(Zix,Nix) == 3 .and. .not.flagcol(Zix, Nix)) wtable(Zix,Nix,1,1) = 0.919
          if (ldmodel(Zix,Nix) == 4) wtable(Zix,Nix,1,1) = 0.921
          if (ldmodel(Zix,Nix) == 5) wtable(Zix,Nix,1,1) = 1.021
          if (ldmodel(Zix,Nix) == 6) wtable(Zix,Nix,1,1) = 0.936
        endif
      endif
      if (strengthM1 == 8 .or. strengthM1 == 10) then
        if (Ainit >= 105) then
          upbend(Zix, Nix, 0, 1, 1) = 1.e-8
          upbend(Zix, Nix, 0, 1, 3) = 0.
        else
          upbend(Zix, Nix, 0, 1, 1) = 3.e-8
          upbend(Zix, Nix, 0, 1, 3) = 4.
        endif
      endif
      if (strengthM1 == 3) then
        upbend(Zix, Nix, 0, 1, 1) = 3.5e-8
        upbend(Zix, Nix, 0, 1, 3) = 6.
      endif
      upbend(Zix, Nix, 0, 1, 2) = 0.8
      if (strength == 8) then
        upbend(Zix, Nix, 1, 1, 1) = 1.e-10
        upbend(Zix, Nix, 1, 1, 2) = 3.0
      else
        upbend(Zix, Nix, 1, 1, 1) = 0.
        upbend(Zix, Nix, 1, 1, 2) = 0.
      endif
    enddo
  enddo
  sfthall = 1.
  sfexpall = 0.347
  if (mod(Atarget,2) /= 0) sfexpall=1.
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
    if (key == 'egr') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        egr(Zix, Nix, irad, lval, igr) = val
        gamadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'ggr') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        ggr(Zix, Nix, irad, lval, igr) = val
        gamadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'sgr') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        sgr(Zix, Nix, irad, lval, igr) = val
        gamadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'epr') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        epr(Zix, Nix, irad, lval, igr) = val
        gamadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'gpr') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        gpr(Zix, Nix, irad, lval, igr) = val
        gamadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'spr') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        tpr(Zix, Nix, irad, lval, igr) = val
        gamadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'egradjust') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        egradjust(Zix, Nix, irad, lval, igr) = val
        gamadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'ggradjust') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        ggradjust(Zix, Nix, irad, lval, igr) = val
        gamadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'sgradjust') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        sgradjust(Zix, Nix, irad, lval, igr) = val
        gamadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'epradjust') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        epradjust(Zix, Nix, irad, lval, igr) = val
        gamadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'gpradjust') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        gpradjust(Zix, Nix, irad, lval, igr) = val
        gamadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'spradjust') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        tpradjust(Zix, Nix, irad, lval, igr) = val
        gamadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'etable') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        etable(Zix, Nix, irad, lval) = val
        gamadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'ftable') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        ftable(Zix, Nix, irad, lval) = val
        gamadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'wtable') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        wtable(Zix, Nix, irad, lval) = val
        gamadjust(Zix, Nix) = .true.
        wtadjust = .true.
      endif
      cycle
    endif
    if (key == 'etableadjust') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        etableadjust(Zix, Nix, irad, lval) = val
        gamadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'ftableadjust') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        ftableadjust(Zix, Nix, irad, lval) = val
        gamadjust(Zix, Nix) = .true.
      endif
      cycle
    endif
    if (key == 'wtableadjust') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        wtableadjust(Zix, Nix, irad, lval) = val
        gamadjust(Zix, Nix) = .true.
        wtadjust = .true.
      endif
      cycle
    endif
    if (key == 'fiso') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) fiso(type) = val
      cycle
    endif
    if (key == 'fisom') then
      class = 6
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) then
        fisom(type) = val
        fisominit(type) = val
      endif
      cycle
    endif
    if (key == 'gamgam') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) gamgam(Zix, Nix) = val
      cycle
    endif
    if (key == 'gamgamadjust') then
      class = 1
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) gamgamadjust(Zix, Nix) = val
      cycle
    endif
    if (key == 'gammax') then
      read(value, * , iostat = istat) gammax
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'upbendc') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) upbend(Zix, Nix, irad, lval, 1) = val
      cycle
    endif
    if (key == 'upbende') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) upbend(Zix, Nix, irad, lval, 2) = val
      cycle
    endif
    if (key == 'upbendf') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) upbend(Zix, Nix, irad, lval, 3) = val
      cycle
    endif
    if (key == 'upbendcadjust') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) upbendadjust(Zix, Nix, irad, lval, 1) = val
      cycle
    endif
    if (key == 'upbendeadjust') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) upbendadjust(Zix, Nix, irad, lval, 2) = val
      cycle
    endif
    if (key == 'upbendfadjust') then
      class = 5
      call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
      if (flagassign) upbendadjust(Zix, Nix, irad, lval, 3) = val
      cycle
    endif
    if (key == 'sfexp') then
      class = 15
      read(value, * , iostat = istat) sfexp
      if (istat /= 0) call read_error(line, istat)
      if (word(3) == ' ') then
        sfexpall = sfexp
      else
        call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
        if (flagassign) then
          if (ival >= 0) then
            spectfacexp(Zix, Nix, ival) = val
          else
            sfexpall = sfexp
          endif
        endif
      endif
      cycle
    endif
    if (key == 'sfth') then
      class = 15
      read(value, * , iostat = istat) sfth
      if (istat /= 0) call read_error(line, istat)
      if (word(3) == ' ') then
        sfthall = sfth
      else
        call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
        if (flagassign) spectfacth(Zix, Nix) = val
      endif
      cycle
    endif
  enddo
!
! In case of external PSF files, no default normalization of width
!
  if (.not.wtadjust) then
    do Zix=0,numZ
      do Nix=0,numN
        do irad=0,1
          do l=1,numgam
            if (Exlfile(Zix,Nix,irad,l)(1:1).ne.' ') wtable(Zix,Nix,irad,l)=1.
          enddo
        enddo
      enddo
    enddo
  endif
!
! Set spectroscopic factors per nucleus
!
  do Nix = 0, numN
    do Zix = 0, numZ
      do nex = 0, numlev
        if (spectfacexp(Zix, Nix, nex) == 0.) spectfacexp(Zix, Nix, nex) = sfexpall
      enddo
      if (spectfacth(Zix, Nix) == 0.) spectfacth(Zix, Nix) = sfthall
    enddo
  enddo
  return
end subroutine input_gammapar
! Copyright A.J. Koning 2021
