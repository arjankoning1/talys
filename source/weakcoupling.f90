subroutine weakcoupling(Zix, Nix, type)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Weak coupling model
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
!   sgl           ! single precision kind
! All global variables
!   numlev2       ! maximum number of levels
! Variables for direct reactions
!   core          ! even - even core for weakcoupling ( - 1 or 1)
! Variables for main input
!   k0            ! index of incident particle
!   Ninit         ! neutron number of initial compound nucleus
!   Zinit         ! charge number of initial compound nucleus
! Variables for nuclides
!   AA            ! mass number of residual nucleus
!   NN            ! neutron number of residual nucleus
!   strucexist    ! flag to state whether structure info for this nucleus exists
!   ZZ            ! charge number of residual nucleus
! Constants
!   onethird      ! 1 / 3
!   parN          ! neutron number of particle
!   parZ          ! charge number of particle
! Variables for levels
!   edis          ! energy of level
!   jdis          ! spin of level
!   nlevmax2      ! maximum number of levels
!   parlev        ! parity of level
! Variables for deformation parameters
!   deform        ! deformation parameter
!   deftype       ! deformation length (D) or parameter (B)
!   jcore         ! spin of level of core nucleus
!   leveltype     ! type of level (rotational (R) or vibrational (V))
!   pcore         ! parity of level of core nucleus
!
! *** Declaration of local data
!
  implicit none
  integer   :: A        ! mass number of target nucleus
  integer   :: i        ! counter
  integer   :: i2       ! value
  integer   :: i3       ! value
  integer   :: i4       ! value
  integer   :: i5       ! help variable
  integer   :: j        ! counter
  integer   :: middle   ! help variable
  integer   :: N        ! neutron number of residual nucleus
  integer   :: Ncore    ! neutron number of core nucleus
  integer   :: Nix      ! neutron number index for residual nucleus
  integer   :: Nixcore  ! neutron number index for core nucleus
  integer   :: type     ! particle type
  integer   :: Z        ! charge number of target nucleus
  integer   :: Zcore    ! charge number of core nucleus
  integer   :: Zix      ! charge number index for residual nucleus
  integer   :: Zixcore  ! charge number index for core nucleus
  real(sgl) :: factor   ! multiplication factor
  real(sgl) :: spinbeg  ! begin of possible spin values
  real(sgl) :: spinend  ! end of possible spin values
!
! ************************* Weak coupling model ************************
!
! The even-even core is determined by subtracting or adding the extra nucleon.
!
  Z = ZZ(Zix, Nix, 0)
  N = NN(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  if (mod(Z, 2) == 1) then
    Zcore = Z + core
    Ncore = N
  else
    Zcore = Z
    Ncore = N + core
  endif
  Zix = abs(Zinit - Zcore)
  Nix = abs(Ninit - Ncore)
  Z = Zcore
!
! The levels and deformation parameters for the core nucleus are retrieved.
!
! levels    : subroutine for discrete levels
! deformpar : subroutine for deformation parameters
!
  if ( .not. strucexist(Zix, Nix)) call levels(Zix, Nix)
  call deformpar(Zix, Nix)
  Zixcore = Zix
  Nixcore = Nix
  Zix = parZ(type)
  Nix = parN(type)
!
! The levels of the core nucleus are distributed over the odd nucleus by angular momentum splitting.
!
  do i = 0, numlev2
    if (i == 0 .and. type == k0) cycle
    if (leveltype(Zix, Nix, i) == 'R' .or. leveltype(Zix, Nix, i) == 'V') cycle
    if (deform(Zixcore, Nixcore, i) /= 0.) then
!
! Determine corresponding energy region in odd nucleus
!
      middle = 0
      do i2 = 1, numlev2
        if (edis(Zix, Nix, i2) > edis(Zixcore, Nixcore, i)) then
          middle = i2 + 1
          exit
        endif
      enddo
!
! Loop over possible split levels
!
      spinbeg = abs(jdis(Zix, Nix, 0) - jdis(Zixcore, Nixcore, i))
      spinend = jdis(Zix, Nix, 0) + jdis(Zixcore, Nixcore, i)
Loop1:  do j = int(spinbeg), int(spinend)
!
! Search for odd-spin level that corresponds with j.
! Start from the level 'middle found above, and work outwards until level is found.
! Next, assign parameters to this level.
!
Loop2:  do i3 = 0, numlev2
          do i4 = 1, - 1, - 2
            i5 = middle + i4 * i3
            if (i3 == 0 .and. i4 ==  - 1) cycle Loop2
            if (i5 < 1 .or. i5 > nlevmax2(Zix, Nix)) cycle
            if (deform(Zix, Nix, i5) /= 0.) cycle
            if (leveltype(Zix, Nix, i5) /= 'D') cycle
            if (int(jdis(Zix, Nix, i5)) == j) then
              jcore(Zix, Nix, i5) = jdis(Zixcore, Nixcore, i)
              pcore(Zix, Nix, i5) = parlev(Zixcore, Nixcore, i)
              factor = sqrt((2. * (j + 0.5) + 1.) / ((2. * jdis(Zixcore, Nixcore, i) + 1.) * (2. * jdis(Zix, Nix, 0) + 1.)))
              deform(Zix, Nix, i5) = factor * deform(Zixcore, Nixcore, i)
              if (deftype(Zix, Nix) == 'D' .and. deftype(Zixcore, Nixcore) == 'B') deform(Zix, Nix, i5) = &
                deform(Zix, Nix, i5) * 1.24 * (A **onethird)
              if (deftype(Zix, Nix) == 'B' .and. deftype(Zixcore, Nixcore) == 'D') deform(Zix, Nix, i5) = &
                deform(Zix, Nix, i5) / (1.24 * (A **onethird))
              cycle Loop1
            endif
          enddo
        enddo Loop2
      enddo Loop1
    endif
  enddo
  return
end subroutine weakcoupling
! Copyright A.J. Koning 2021
