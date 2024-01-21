subroutine input_directmodel
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for direct reaction models
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
! Variables for direct reactions
!   flagautorot     ! flag for automatic rotational coupled channels calculation
!   flagcoulomb     ! flag for Coulomb excitation calculation with ECIS
!   flagcpang       ! flag for compound angular distribution calculation
!   flagdirect      ! flag for output of direct reaction results
!   flagdisc        ! flag for output of discrete state cross sections
!   flageciscalc    ! flag for new ECIS calculation for outgoing channels
!   flagecissave    ! flag for saving ECIS input and output files
!   flaggiant0      ! flag for collective contribution from giant resonances
!   flaginccalc     ! flag for new ECIS calculation for incident channel
!   flaglegendre    ! flag for output of Legendre coefficients
!   flagoutecis     ! flag for output of ECIS results
!   flagrot         ! flag for use of rotational optical model per outgoing particle
!   flagstate       ! flag for optical model potential for each excited state
!   flagsys         ! flag for reaction cross section from systematics
!   flagtransen     ! flag for output of transmission coefficients per energy
! Constants
!   parsym          ! symbol of particle
! All global variables
!   numpar      ! number of particles
! Variables for main input
!   k0           ! index of incident particle
! Variables for basic reaction
!   flagffruns        ! flag to designate subsequent evaporation of fission products
!   flagrpruns        ! flag to designate that run is for residual product
! Variables for basic reaction
!   flagbasic    ! flag for output of basic information and results
!   flagendf     ! flag for information for ENDF - 6 file
!   flagmicro    ! flag for completely microscopic TALYS calculation
! Variables for OMP
!   flagompall    ! flag for new optical model calculation for all residual
!   flagomponly   ! flag to execute ONLY an optical model calculation
! Variables for compound reactions
!   flageciscomp ! flag for compound nucleus calculation by ECIS
! Variables for reading input lines
!   inline            ! input line
!   nlines            ! number of input lines
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=1)   :: ch                   ! character
  character(len=132) :: key                  ! keyword
  character(len=132) :: word(40)             ! words on input line
  character(len=132) :: line                 ! input line
  integer            :: i                    ! counter
  integer            :: i2                   ! counter
  integer            :: istat                ! logical for file access
  integer            :: type                 ! particle type
!
! ************** Defaults *************
!
  if (flagmicro) then
    flagautorot = .true.
  else
    flagautorot = .false.
  endif
  flagrot = .false.
  flagcoulomb = .true.
  flagcpang = .false.
  flagdirect = flagbasic
  flagdisc = flagbasic
  if (flagendf) flagdisc = .true.
  if (flagffruns .or. flagrpruns) flagdisc = .false.
  flageciscalc = .true.
  if (flagffruns) flageciscalc = .true.
  if (flagrpruns) flageciscalc = .true.
  flaginccalc = .true.
  if (flagffruns) flaginccalc = .true.
  if (flagrpruns) flaginccalc = .true.
  if (flagompall) then
    flagecissave = .true.
  else
    flagecissave = .false.
  endif
  if (k0 == 1 .or. k0 == 2) then
    flaggiant0 = .true.
  else
    flaggiant0 = .false.
  endif
  if (flagomponly) flaggiant0 = .false.
  if (flagendf) then
    flaglegendre = .true.
  else
    flaglegendre = .false.
  endif
  flagoutecis = flageciscomp
  flagstate = .false.
  flagsys = .false.
  flagtransen = .true.
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
    ch = word(2)(1:1)
    type = 0
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
    if (key == 'autorot') then
      if (ch == 'n') flagautorot = .false.
      if (ch == 'y') flagautorot = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'eciscalc') then
      if (ch == 'n') flageciscalc = .false.
      if (ch == 'y') flageciscalc = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'inccalc') then
      if (ch == 'n') flaginccalc = .false.
      if (ch == 'y') flaginccalc = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'cpang') then
      if (ch == 'n') flagcpang = .false.
      if (ch == 'y') flagcpang = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'statepot') then
      if (ch == 'n') flagstate = .false.
      if (ch == 'y') flagstate = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'rotational') then
Loop2: do i2 = 2, 40
        ch = word(i2)(1:1)
        do type = 1, 6
          if (ch == parsym(type)) then
            flagrot(type) = .true.
            cycle Loop2
          endif
        enddo
      enddo Loop2
      cycle
    endif
    if (key == 'sysreaction') then
      do type = 0, 6
        flagsys(type) = .false.
      enddo
Loop1: do i2 = 2, 40
        ch = word(i2)(1:1)
        do type = 0, 6
          if (ch == parsym(type)) then
            flagsys(type) = .true.
            cycle Loop1
          endif
        enddo
      enddo Loop1
      cycle
    endif
    if (key == 'outecis') then
      if (ch == 'n') flagoutecis = .false.
      if (ch == 'y') flagoutecis = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'ecissave') then
      if (ch == 'n') flagecissave = .false.
      if (ch == 'y') flagecissave = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outdirect') then
      if (ch == 'n') flagdirect = .false.
      if (ch == 'y') flagdirect = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outlegendre') then
      if (ch == 'n') flaglegendre = .false.
      if (ch == 'y') flaglegendre = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'coulomb') then
      if (ch == 'n') flagcoulomb = .false.
      if (ch == 'y') flagcoulomb = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'giantresonance') then
      if (ch == 'n') flaggiant0 = .false.
      if (ch == 'y') flaggiant0 = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outtransenergy') then
      if (ch == 'n') flagtransen = .false.
      if (ch == 'y') flagtransen = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outdiscrete') then
      if (ch == 'n') flagdisc = .false.
      if (ch == 'y') flagdisc = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
  enddo
  return
end subroutine input_directmodel
! Copyright A.J. Koning 2021
