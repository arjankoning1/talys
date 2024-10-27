subroutine input_outfiles
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for output files
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
! All global variables
!   numfile         ! maximum number of separate output files
!   numflux         ! number of integral experiments
!   numlev          ! maximum number of discrete levels
! Variables for output
!   fileangle       ! designator for angular distributions on separate file
!   filechannels    ! flag for exclusive channel cross sections on separate file
!   filediscrete    ! flag for discrete level cross sections on separate file
!   fileelastic     ! flag for elastic angular distribution on separate file
!   filefission     ! flag for fission cross sections on separate file
!   filegamdis      ! flag for gamma-ray intensities on separate file
!   filerecoil      ! flag for recoil spectra on separate file
!   fileresidual    ! flag for residual production cross sections on separate file
!   filespectrum    ! designator for spectrum on separate file
!   filetotal       ! flag for total cross sections on separate file
!   flagexc         ! flag for output of excitation functions
!   ddxacount       ! counter for double-differential cross section files
!   ddxecount       ! counter for double-differential cross section files
!   fileddxa        ! designator for double-differential cross sections on separate file
!   fileddxe        ! designator for double-differential cross sections on separate file
! Constants
!   parsym          ! symbol of particle
! All global variables
!   numfile    ! maximum number of separate output files
!   numflux    ! number of integral experiments
!   numlev     ! maximum number of discrete levels
!   numpar     ! number of particles
! Variables for reading input lines
!   inline                ! input line
!   nlines                ! number of input lines
! Variables for basic reaction
!   flagastro       ! flag for calculation of astrophysics reaction rate
!   flagchannels    ! flag for exclusive channel calculation
!   flagendf        ! flag for information for ENDF - 6 file
!   flagendfdet     ! flag for detailed ENDF - 6 information per channel
!   flagrecoil      ! flag for calculation of recoils
! Variables for fission
!   flagfission     ! flag for fission
! Variables for direct reactions
!   flagdisc        ! flag for output of discrete state cross sections
! Variables for input energies
!   Ninc       ! number of incident energies
! Variables for main input
!   flagnatural    ! flag for calculation of natural element
!   k0             ! index of incident particle
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=1)   :: ch         ! character
  character(len=132) :: key        ! keyword
  character(len=132) :: value      ! value or string
  character(len=132) :: word(40)   ! words on input line
  character(len=132) :: line       ! input line
  integer            :: i          ! counter
  integer            :: i2         ! counter
  integer            :: istat      ! logical for file access
  integer            :: ivalue     ! counter
  integer            :: type       ! particle type
  real(sgl)          :: val        ! real value
!
! ************** Defaults *************
!
  fileangle = .false.
  filechannels = .false.
  filediscrete = .false.
  fileelastic = .false.
  filegamdis = .false.
  filerecoil = .false.
  fileresidual = .false.
  if (flagspec) then
    filespectrum = .true.
  else
    filespectrum = .false.
  endif
  filetotal = .false.
  ddxacount = 0
  ddxecount = 0
  fileddxa = 0.
  fileddxe = 0.
!
! If the results of TALYS are used to create ENDF-6 data files, several output flags are automatically set.
!
  if (flagendf) then
    fileelastic = .true.
    filetotal = .true.
    fileresidual = .true.
    if (flagrecoil) filerecoil = .true.
    filespectrum = .true.
    if (flagendfdet) then
      filechannels = .true.
      filegamdis = .true.
      fileangle = .true.
      filediscrete = .true.
    endif
  endif
!
! If the results of TALYS are written as excitation functions in the output file, several output flags are automatically set.
!
  if (flagexc) then
    filetotal = .true.
    fileresidual = .true.
    if (flagchannels) filechannels = .true.
    if (flaggamdis) filegamdis = .true.
    if (flagdisc) filediscrete = .true.
  endif
  if (flagastro) fileresidual = .true.
  filefission = .false.
  if (flagfission) then
    if (flagendf) filefission = .true.
    if (flagexc) filefission = .true.
  endif
!
! Explicit double-differential cross sections for deuteron ENDF files
!
  if (flagendf .and. k0 == 3) then
    do type = 1, 2
      filespectrum(type) = .true.
      ddxacount(type) = 18
      do i = 1, 7
        fileddxa(type, i) = 5. * (i - 1)
      enddo
      do i = 8, 14
        fileddxa(type, i) = 30. + 10. * (i - 7)
      enddo
      do i = 15, 18
        fileddxa(type, i) = 100. + 20. * (i - 14)
      enddo
    enddo
  endif
!
! **************** Read input variables *******************
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
    if (key == 'filespectrum') then
      filespectrum = .false.
Loop2: do i2 = 2, 40
        ch = word(i2)(1:1)
        do type = 0, 6
          if (ch == parsym(type)) then
            filespectrum(type) = .true.
            cycle Loop2
          endif
        enddo
      enddo Loop2
      cycle
    endif
    if (key == 'fileelastic') then
      if (ch == 'n') fileelastic = .false.
      if (ch == 'y') fileelastic = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'filetotal') then
      if (ch == 'n') filetotal = .false.
      if (ch == 'y') filetotal = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'fileresidual') then
      if (ch == 'n') fileresidual = .false.
      if (ch == 'y') fileresidual = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'filefission') then
      if (ch == 'n') filefission = .false.
      if (ch == 'y') filefission = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'filechannels') then
      if (ch == 'n') filechannels = .false.
      if (ch == 'y') filechannels = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'filerecoil') then
      if (ch == 'n') filerecoil = .false.
      if (ch == 'y') filerecoil = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'fileddxe') then
      do type = 0, 6
        if (ch == parsym(type)) then
          ddxecount(type) = ddxecount(type) + 1
          call range_integer_error(key, ddxecount(type), 1, numfile)
          read(word(3), * , iostat = istat) val
          fileddxe(type, ddxecount(type)) = val
          cycle Loop1
        endif
      enddo
      call read_error(line, istat)
    endif
    if (key == 'fileddxa') then
      do type = 0, 6
        if (ch == parsym(type)) then
          ddxacount(type) = ddxacount(type) + 1
          call range_integer_error(key, ddxacount(type), 1, numfile)
          read(word(3), * , iostat = istat) val
          if (istat /= 0) call read_error(line, istat)
          fileddxa(type, ddxacount(type)) = val
          cycle Loop1
        endif
      enddo
      call read_error(line, istat)
    endif
    if (key == 'fileangle') then
      read(value, * , iostat = istat) ivalue
      if (istat /= 0) call read_error(line, istat)
      call range_integer_error(key, ivalue, 0, numlev)
      fileangle(ivalue) = .true.
      cycle
    endif
    if (key == 'filediscrete') then
      read(value, * , iostat = istat) ivalue
      if (istat /= 0) call read_error(line, istat)
      call range_integer_error(key, ivalue, 0, numlev)
      filediscrete(ivalue) = .true.
      cycle
    endif
    if (key == 'filegamdis') then
      if (ch == 'n') filegamdis = .false.
      if (ch == 'y') filegamdis = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
  enddo Loop1
  return
end subroutine input_outfiles
! Copyright A.J. Koning 2021
