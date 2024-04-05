subroutine input_basicreac
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for basic reaction variables
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
! 2023-03-07: Added nffit keyword and flag
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
!
  use A0_talys_mod
  use A1_error_handling_mod
!
! Variables for basic reactions
!   flagastro        ! flag for calculation of astrophysics reaction rate
!   flagbasic        ! flag for output of basic information and results
!   flagoutall       ! flag for output of all data in main output file
!   flagchannels     ! flag for exclusive channel calculation
!   flagEchannel     ! flag for channel energy for emission spectrum
!   flagendf         ! flag for information for ENDF-6 file
!   flagendfdet      ! flag for detailed ENDF-6 information per channel
!   flagendfecis     ! flag for new ECIS calculation for ENDF-6 files
!   flagffruns       ! flag to designate subsequent evap. of fission products
!   flaglabddx       ! flag for calculation of DDX in LAB system
!   flagmassdis      ! flag for calculation of fission fragment mass yields
!   flagmicro        ! flag for completely microscopic TALYS calculation
!   flagngfit        ! flag for using fitted (n,g) nuclear model parameters
!   flagnffit        ! flag for using fitted (n,f) nuclear model parameters
!   flagnnfit        ! flag for using fitted (n,n'), (n,2n) and (n,p) nuclear model parameters
!   flagnafit        ! flag for using fitted (n,a) nuclear model parameters
!   flagpnfit        ! flag for using fitted (p,n) nuclear model parameters
!   flagdnfit        ! flag for using fitted (g,n) nuclear model parameters
!   flaggnfit        ! flag for using fitted (d,n) nuclear model parameters
!   flaganfit        ! flag for using fitted (a,n) nuclear model parameters
!   flaggamgamfit    ! flag for using fitted Gamma_gamma nuclear model parameters
!   flagmacsfit      ! flag for using fitted MACS nuclear model parameters
!   flagpartable     ! flag for output of model parameters on separate file
!   flagpopMeV       ! flag to use initial population per MeV instead of histogram
!   flagreaction     ! flag for calculation of nuclear reactions
!   flagrecoil       ! flag for calculation of recoils
!   flagrecoilav     ! flag for average velocity in recoil calculation
!   flagrel          ! flag for relativistic kinematics
!   flagrpevap       ! flag for evaporation of residual products at high incident energies
!   flagrpruns       ! flag to designate that run is for residual product
!   ompenergyfile    ! file with energies for OMP calculation (ENDF files only)
! Variables for reading input lines
!   inline            ! input line
!   nlines                ! number of input lines
! Variables for best files
!   flagbest      ! flag to use best set of adjusted parameters
! Variables for main input
!   k0             ! index of incident particle
! Error handling
!   range_index_error    ! Test if index is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=1)   :: ch                   ! character
  character(len=132) :: key                  ! keyword
  character(len=132) :: value                ! value or string
  character(len=132) :: word(40)             ! words on input line
  character(len=132) :: line                 ! input line
  integer            :: i                    ! counter
  integer            :: istat                ! logical for file access
!
! ************** Defaults *************
!
  flagrpruns = .false.
  flagbasic = .false.
  flagoutall = .false.
  flagEchannel = .false.
  flagendf = .false.
  if (k0 <= 1) then
    flagendfdet = .true.
  else
    flagendfdet = .false.
  endif
  flagendfecis = .true.
  flagchannels = .false.
  flaglabddx = .false.
  flagpopMeV = .false.
  flagmassdis = .false.
  flagmicro = .false.
  flagpartable = .false.
  flagreaction = .true.
  flagngfit = (k0 == 1 .and. flagfit .and. .not.flagastro)
  flagnnfit = (k0 == 1 .and. flagfit)
  flagnffit = (k0 == 1 .and. flagfit)
  flagnafit = (k0 == 1 .and. flagfit)
  flagpnfit = (k0 == 2 .and. flagfit)
  flaggnfit = (k0 == 0 .and. flagfit)
  flagdnfit = (k0 == 3 .and. flagfit)
  flaganfit = (k0 == 6 .and. flagfit)
  flagmacsfit = (k0 == 1 .and. flagfit .and. flagastro)
  flaggamgamfit = .false.
  flagrecoil = .false.
  flagrecoilav = .false.
  flagrel = .true.
  flagrpevap = .false.
  ompenergyfile = ' '
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
    if (key == 'endfecis') then
      if (ch == 'n') flagendfecis = .false.
      if (ch == 'y') flagendfecis = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'relativistic') then
      if (ch == 'n') flagrel = .false.
      if (ch == 'y') flagrel = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'channels') then
      if (ch == 'n') flagchannels = .false.
      if (ch == 'y') flagchannels = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'endf') then
      if (ch == 'n') flagendf = .false.
      if (ch == 'y') flagendf = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'endfdetail') then
      if (ch == 'n') flagendfdet = .false.
      if (ch == 'y') flagendfdet = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'micro') then
      if (ch == 'n') flagmicro = .false.
      if (ch == 'y') flagmicro = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outbasic') then
      if (ch == 'n') flagbasic = .false.
      if (ch == 'y') flagbasic = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'outall') then
      if (ch == 'n') flagoutall = .false.
      if (ch == 'y') flagoutall = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'recoil') then
      if (ch == 'n') flagrecoil = .false.
      if (ch == 'y') flagrecoil = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'labddx') then
      if (ch == 'n') flaglabddx = .false.
      if (ch == 'y') flaglabddx = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'popmev') then
      if (ch == 'n') flagpopMeV = .false.
      if (ch == 'y') flagpopMeV = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'recoilaverage') then
      if (ch == 'n') flagrecoilav = .false.
      if (ch == 'y') flagrecoilav = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'channelenergy') then
      if (ch == 'n') flagEchannel = .false.
      if (ch == 'y') flagEchannel = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'reaction') then
      if (ch == 'n') flagreaction = .false.
      if (ch == 'y') flagreaction = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'ngfit') then
      if (ch == 'n') flagngfit = .false.
      if (ch == 'y') then
        flagngfit = .true.
        flagfit = .true.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'nffit') then
      if (ch == 'n') flagnffit = .false.
      if (ch == 'y') then
        flagnffit = .true.
        flagfit = .true.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'nnfit') then
      if (ch == 'n') flagnnfit = .false.
      if (ch == 'y') then
        flagnnfit = .true.
        flagfit = .true.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'nafit') then
      if (ch == 'n') flagnafit = .false.
      if (ch == 'y') then
        flagnafit = .true.
        flagfit = .true.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'pnfit') then
      if (ch == 'n') flagpnfit = .false.
      if (ch == 'y') then
        flagpnfit = .true.
        flagfit = .true.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'dnfit') then
      if (ch == 'n') flagdnfit = .false.
      if (ch == 'y') then
        flagdnfit = .true.
        flagfit = .true.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'gnfit') then
      if (ch == 'n') flaggnfit = .false.
      if (ch == 'y') then
        flaggnfit = .true.
        flagfit = .true.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'anfit') then
      if (ch == 'n') flaganfit = .false.
      if (ch == 'y') then
        flaganfit = .true.
        flagfit = .true.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'gamgamfit') then
      if (ch == 'n') flaggamgamfit = .false.
      if (ch == 'y') then
        flaggamgamfit = .true.
        flagfit = .true.
        flagngfit = .false.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'macsfit') then
      if (ch == 'n') flagmacsfit = .false.
      if (ch == 'y') then
        flagmacsfit = .true.
        flagfit = .true.
        flagngfit = .false.
      endif
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'ompenergyfile') then
      ompenergyfile = value
      cycle
    endif
    if (key == 'massdis') then
      if (ch == 'n') flagmassdis = .false.
      if (ch == 'y') flagmassdis = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'partable') then
      if (ch == 'n') flagpartable = .false.
      if (ch == 'y') flagpartable = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == 'rpevap') then
      if (ch == 'n') flagrpevap = .false.
      if (ch == 'y') flagrpevap = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
  enddo Loop1
  return
end subroutine input_basicreac
! Copyright A.J. Koning 2021
