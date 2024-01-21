subroutine particles
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Determine included light particles
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Variables for compound reactions
!   eurr           ! off - set incident energy for URR calculation
!   flagurr        ! flag for output of unresolved resonance parameters
! Variables for fission
!   flagfission    ! flag for fission
! Variables for basic parameters
!   outtype        ! type of outgoing particles
! Variables for main input
!   k0             ! index of incident particle
! Variables for OMP
!   flagomponly    ! flag to execute ONLY an optical model calculation
! Variables for nuclides
!   parinclude     ! logical to include outgoing particle
!   parskip        ! logical to skip outgoing particle
! Constants
!   parsym         ! symbol of particle
!
! *** Declaration of local data
!
  implicit none
  logical :: default              ! logical to determine default outgoing particles
  integer :: type                 ! particle type
  integer :: type2                ! particle type
!
! ************ Determine outgoing particles to be included *************
!
! The default is to include all particles from photons to alpha's as competing channels.
! If specific outgoing particles in the input are given, only those will be included as competing channels.
! The logicals parinclude and parskip will be referred to throughout TALYS.
! They are always each others' opposite.
!
! 1. Default: All competing channels included
!
  parinclude = .true.
  parskip = .false.
  if ( .not. flagfission) then
    parinclude( - 1) = .false.
    parskip( - 1) = .true.
  endif
  if (flagomponly) then
    parinclude = .false.
    parskip = .true.
  endif
!
! Check if default is to be used
!
  default = .true.
  do type = 0, 6
    if (outtype(type) /= ' ') then
      default = .false.
      exit
    endif
  enddo
  if (default) then
    do type = 0, 6
      outtype(type) = parsym(type)
    enddo
!
! 2. No default, but specific competing outgoing particles are requested.
!    Now, we first turn all competing channels off, and then determine which are to be included and which are to be skipped.
!
  else
    do type = 0, 6
      parinclude(type) = .false.
    enddo
    do type = 0, 6
      do type2 = 0, 6
        if (outtype(type) == parsym(type2)) parinclude(type2) = .true.
      enddo
    enddo
    do type = 0, 6
      parskip(type) = .not. parinclude(type)
    enddo
  endif
!
! The incident particle is always included as outgoing particle
!
  parinclude(k0) = .true.
  parskip(k0) = .false.
!
! Setting for URR
!
  if (parskip(0)) then
    eurr = 0.
    flagurr = .false.
  endif
  return
end subroutine particles
! Copyright A.J. Koning 2021
