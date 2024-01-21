subroutine isotrans(Z, N)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Correction factors for isospin forbidden transitions
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
! Definition of single and double precision variables
!   sgl        ! single precision kind
! Variables for gamma rays
!   fiso       ! correction factor for isospin forbidden transitions
!   fisom      ! correction factor for isospin forbidden transitions for multiple emission
! Variables for main input
!   k0         ! index of incident particle
! All global variables
!   numpar     ! number of particles
! Variables for nuclides
!   primary    ! flag to designate primary (binary) reaction
!
!
! *** Declaration of local data
!
  implicit none
  logical   :: flagiso       !
  integer   :: N             ! neutron number of nucleus
  integer   :: type          ! particle type
  integer   :: Z             ! charge number of nucleus
  real(sgl) :: ff(-1:numpar) !
!
! ******** Correction factors for isospin forbidden transitions ********
!
  ff = 1.
  flagiso = .false.
  if (Z == N) then
    flagiso = .true.
    if (k0 == 0) then
      ff(1) = 2.
      ff(2) = 2.
      ff(6) = 5.
    endif
    if (k0 == 1) ff(0) = 2.
    if (k0 == 2) ff(0) = 2.
    if (k0 == 6) ff(0) = 5.
  endif
  if (Z == N - 1 .or. Z == N + 1) then
    flagiso = .true.
    if (k0 == 0) then
      ff(1) = 1.5
      ff(2) = 1.5
      ff(6) = 1.5
    endif
    if (k0 == 1) ff(0) = 1.5
    if (k0 == 2) ff(0) = 1.5
    if (k0 == 6) ff(0) = 1.5
  endif
  if (flagiso) then
    if (primary) then
      do type = - 1, 6
        if (fiso(type) == -1.) fiso(type) = ff(type)
      enddo
    else
      do type = - 1, 6
        if (fisom(type) == -1.) fisom(type) = ff(type)
      enddo
    endif
  else
    if (primary) then
      do type = - 1, 6
        if (fiso(type) == -1.) fiso(type) = 1.
      enddo
    else
      do type = - 1, 6
        if (fisom(type) == -1.) fisom(type) = 1.
      enddo
    endif
  endif
end subroutine isotrans
! Copyright A.J. Koning 2021
