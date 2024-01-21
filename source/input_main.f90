subroutine input_main
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for main variables
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
! All global variables
!   numelem        ! number of elements
! Variables for main
!   ptype0     ! type of incident particle
!   Atarget    ! mass number of target nucleus
!   Ltarget    ! excited level of target
!   Starget    ! symbol of target nucleus
!   Ztarget    ! charge number of target nucleus
!   deninc     ! incident energy increment
!   energyfile ! file with energies for OMP calculation
!   eninc      ! incident energy in MeV
!   enincF     ! final incident energy
!   Estop      ! incident energy above which TALYS stops
!   flagastro  ! flag for calculation of astrophysics reaction rate
! Variables for reading input lines
!   inline    ! input line
!   nlines    ! number of input lines
! Variables for abundance
!   isotope       ! isotope of natural element
! Constants
!   iso           ! counter for isotope
!   nuc           ! symbol of nucleus
!   parsym        ! symbol of particle
!   parN          ! neutron number of particle
!   parZ          ! charge number of particle
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical            :: elemexist    ! logical for existence of element
  logical            :: massexist    ! logical for existence of mass
  logical            :: projexist    ! logical for existence of projectile
  logical            :: enerexist    ! logical for existence of energy
  character(len=1)   :: ch           ! character
  character(len=132) :: key          ! keyword
  character(len=132) :: value        ! value or string
  character(len=132) :: word(40)     ! words on input line
  character(len=132) :: line         ! input line
  integer            :: i            ! counter
  integer            :: istat        ! logical for file access
  integer            :: iz           ! counter
  integer            :: type         ! particle type
!
! ************** Defaults *************
!
  projexist = .false.
  massexist = .false.
  elemexist = .false.
  enerexist = .false.
  ptype0 = ' '
  Starget = '  '
  Atarget = 0
  Ztarget = 0
  Ltarget = 0
  flagnatural = .false.
  Estop = Emaxtalys
  energyfile = ' '
  eninc = 0.
  enincF = 0.
  deninc = 0.
  flagastro = .false.
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
    if (key == 'projectile') then
      projexist = .true.
      ptype0 = ch
      cycle
    endif
    if (key == 'mass') then
      massexist = .true.
      read(value, * , iostat = istat) Atarget
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'element') then
      elemexist = .true.
      if (ch >= '0' .and. ch <= '9') then
        read(value, * , iostat = istat) Ztarget
        if (istat /= 0) call read_error(line, istat)
        call range_integer_error(key, Ztarget, 1, numelem)
        cycle
      else
        read(value, * , iostat = istat) Starget
        if (istat /= 0) call read_error(line, istat)
        Starget(1:1) = achar(iachar(Starget(1:1)) - 32)
        cycle
      endif
    endif
    if (key == 'ltarget') then
      read(value, * , iostat = istat) Ltarget
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'energy') then
      enerexist = .true.
      if ((ch >= '0' .and. ch <= '9') .or. ch == '.') then
        read(value, * , iostat = istat) eninc(1)
        if (istat /= 0) call read_error(line, istat)
        read(word(3), * , iostat = istat) enincF
        if (istat /= 0) cycle
        read(word(4), * , iostat = istat) deninc
        if (istat /= 0) cycle
        cycle
      else
        eninc(1) = 0.
        energyfile = value
        cycle
      endif
    endif
    if (key == 'estop') then
      read(value, * , iostat = istat) Estop
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'astro') then
      if (ch == 'n') flagastro = .false.
      if (ch == 'y') flagastro = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
  enddo
!
! These main keywords MUST be present in the input file.
!
  if ( .not. projexist) then
    write(*, '(" TALYS-error: projectile must be given")')
    stop
  endif
  if ( .not. massexist) then
    write(*, '(" TALYS-error: mass must be given")')
    stop
  endif
  if ( .not. elemexist) then
    write(*, '(" TALYS-error: element must be given")')
    stop
  endif
  if ( .not. enerexist) then
    write(*, '(" TALYS-error: energy must be given")')
    stop
  endif
!
! Identification of target and nucleus
!
  if (Ztarget == 0) then
    do iz = 1, numelem
      if (nuc(iz) == Starget) then
        Ztarget = iz
        exit
      endif
    enddo
  else
    Starget = nuc(Ztarget)
  endif
  call initial_loop
!
! A calculation for a natural element is specified by target mass 0
!
! abundance: subroutine for natural abundances
!
  if (Atarget == 0) then
    flagnatural = .true.
    if (iso == 1) call abundance
    Atarget = isotope(iso)
  endif
  Ntarget = Atarget - Ztarget
!
! Assignment of index k0 to incident particle
!
! Throughout TALYS, the initial particle can always be identified by the index k0.
!
  k0 = 0
  do type = 0, 6
    if (ptype0 == parsym(type)) then
      k0 = type
      exit
    endif
  enddo
!
! Identification of initial compound nucleus
!
  Zinit = Ztarget + parZ(k0)
  Ninit = Ntarget + parN(k0)
  Ainit = Zinit + Ninit
  return
end subroutine input_main
! Copyright A.J. Koning 2021
