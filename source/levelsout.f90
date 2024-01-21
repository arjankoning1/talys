subroutine levelsout(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of discrete levels
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
! Variables for discrete levels
!   nlev           ! number of levels for nucleus
! Variables for nuclides
!   AA             ! mass number of residual nucleus
!   NN             ! neutron number of residual nucleus
!   parskip        ! logical to skip outgoing particle
!   ZZ             ! charge number of residual nucleus
! Constants
!   cparity        ! parity (character)
!   nuc            ! symbol of nucleus
!   parname        ! name of particle
! Variables for levels
!   bassign        ! flag for assignment of branching ratio
!   branchlevel    ! level to which branching takes place
!   branchratio    ! gamma-ray branching ratio to level
!   edis           ! energy of level
!   ENSDF          ! string from original ENSDF discrete level file
!   jassign        ! flag for assignment of spin
!   jdis           ! spin of level
!   nbranch        ! number of branching levels
!   parlev         ! parity of level
!   passign        ! flag for assignment of parity
!   tau            ! lifetime of state in seconds
! Variables for masses
!   nucmass        ! mass of nucleus
!   S              ! separation energy
!
! *** Declaration of local data
!
  implicit none
  integer :: A     ! mass number of target nucleus
  integer :: i     ! level
  integer :: k     ! designator for particle
  integer :: N     ! neutron number of residual nucleus
  integer :: Nix   ! neutron number index for residual nucleus
  integer :: type  ! particle type
  integer :: Z     ! charge number of target nucleus
  integer :: Zix   ! charge number index for residual nucleus
!
! *************************** Discrete levels **************************
!
!
  Z = ZZ(Zix, Nix, 0)
  N = NN(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  write(*, '(/" NUCLEAR STRUCTURE INFORMATION FOR Z=", i3, " N=", i3, " (", i3, a2, ") ")') Z, N, A, nuc(Z)
  write(*, '(/" Mass in a.m.u.     : ", f10.6)') nucmass(Zix, Nix)
  write(*, '(/" Separation energies:")')
  write(*, '(/" Particle        S         "/)')
  do type = 1, 6
    if (parskip(type)) cycle
    write(*, '(1x, a8, f12.5)') parname(type), S(Zix, Nix, type)
  enddo
  write(*, '(/" Discrete levels of Z=", i3, " N=", i3, " (", i3, a2, ") ")') Z, N, A, nuc(Z)
  write(*, '(/" Number  Energy Spin Parity  Branching Ratio (%) Lifetime(sec) Assignment        ENSDF"/)')
  do i = 0, nlev(Zix, Nix)
    if (tau(Zix, Nix, i) /= 0.) then
      write(*, '(1x, i3, 4x, f7.4, 1x, f4.1, 3x, a1, 24(" "), 2x, es10.3, 7x, 2a1, a18)') levnum(Zix, Nix, i), &
 &      edis(Zix, Nix, i), jdis(Zix, Nix, i), cparity(parlev(Zix, Nix, i)), tau(Zix, Nix, i), jassign(Zix, Nix, i), &
 &      passign(Zix, Nix, i), ENSDF(Zix, Nix, i)
    else
      write(*, '(1x, i3, 4x, f7.4, 1x, f4.1, 3x, a1, 36(" "), 7x, 2a1, a18)') i, edis(Zix, Nix, i), jdis(Zix, Nix, i), &
 &      cparity(parlev(Zix, Nix, i)), jassign(Zix, Nix, i), passign(Zix, Nix, i), ENSDF(Zix, Nix, i)
    endif
    do k = 1, nbranch(Zix, Nix, i)
      write(*, '(31x, "--->", i3, 2x, f8.4, 18x, a1)') branchlevel(Zix, Nix, i, k), branchratio(Zix, Nix, i, k) * 100., &
 &      bassign(Zix, Nix, i, k)
    enddo
  enddo
  return
end subroutine levelsout
! Copyright A.J. Koning 2021
