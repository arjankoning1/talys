subroutine mainout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Main output
!
! Author    : Arjan Koning
!
! 2023-12-30: Original code
! 2025-02-21: Current revision
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Variables for compound reactions
!   flagcomp       ! flag for compound angular distribution calculation
! Variables for input energies
!   EdistE         ! excitation energy of population distribution
!   eninc          ! incident energy in MeV
!   flaginitpop    ! flag for initial population distribution
!   npopE          ! number of energies for population distribution
!   npopJ          ! number of spins for population distribution
!   Ninc           ! number of incident energies
!   PdistE         ! population distribution, spin - independent
!   PdistJP        ! population distribution per spin and parity
! Variables for main input
!   Atarget        ! mass number of target nucleus
!   k0             ! index of incident particle
!   Ltarget        ! excited level of target
! Variables for fission
!   flagfisout     ! flag for output of fission information
! Variables for discrete levels
!   flaglevels     ! flag for output of discrete level information
! Variables for level density
!   filedensity    ! flag for level densities on separate files
!   flagdensity    ! flag for output of level densities
! Variables for OMP
!   flagomponly    ! flag to execute ONLY an optical model calculation
! Variables for nuclides
!   Nindex         ! neutron number index for residual nucleus
!   parskip        ! logical to skip outgoing particle
!   Q              ! Q - value
!   strucwrite     ! flag for output of nuclear structure info
!   tarmass        ! mass of target nucleus
!   Zindex         ! charge number index for residual nucleus
! Constants
!   cparity        ! parity (character)
!   parmass        ! mass of particle in a.m.u.
!   parname        ! name of particle
!   parsym         ! symbol of particle
! Variables for levels
!   edis           ! energy of level
!   jdis           ! spin of level
!   parlev         ! parity of level
!   tau            ! lifetime of state in seconds
!
! *** Declaration of local data
!
  implicit none
  integer :: i                   ! counter
  integer :: J                   ! spin of level
  integer :: Ncomp               ! neutron number index for compound nucleus
  integer :: Nix                 ! neutron number index for residual nucleus
  integer :: parity              ! parity
  integer :: type                ! particle type
  integer :: Zcomp               ! proton number index for compound nucleus
  integer :: Zix                 ! charge number index for residual nucleus
!
! *************************** Code and version *************************
!
  write(*, '(/"    TALYS-2.1 (Version: February 20, 2025)"/)')
  write(*, '(" Copyright (C) 2025  A.J. Koning, S. Hilaire and S. Goriely"/)')
  write(*, '(" Dimensions - Cross sections: mb, Energies: MeV, Angles: degrees")')
  write(*, '(/" User: ",a)') trim(user)
  write(*, '(" Date: ",a)') trim(date)
!
! ***************** Write input file and default parameters ************
!
! inputout: subroutine to write the input parameters
!
  call inputout
!
! ********************** Main nuclear parameters ***********************
!
  Zcomp = 0
  Ncomp = 0
  Zix = Zindex(Zcomp, Ncomp, k0)
  Nix = Nindex(Zcomp, Ncomp, k0)
  write(*, '(/" ########## BASIC REACTION PARAMETERS ##########"/)')
  write(*, '(" Projectile           : ", a8, t37, "Mass in a.m.u.      : ", f10.6)') parname(k0), parmass(k0)
  write(*, '(" Target               : ", a, t37, "Mass in a.m.u.      : ", f10.6)') trim(targetnuclide), tarmass
  if (Ltarget /= 0) then
    write(*, '(/" Excited target level : Number  Energy  Spin Parity Lifetime(sec)")')
    write(*, '(24x, i3, 4x, f7.4, 2x, f4.1, 3x, a1, 4x, es10.3)') Ltarget, edis(Zix, Nix, Ltarget), jdis(Zix, Nix, Ltarget), &
 &    cparity(parlev(Zix, Nix, Ltarget)), tau(Zix, Nix, Ltarget)
  endif
  write(*, '(/" Included channels:")')
  do type = -1, 6
    if (parskip(type)) cycle
    write(*, '(21x, a8)') parname(type)
  enddo
  if (flagomponly .and. .not. flagcomp) return
!
! Projectile
!
  if ( .not. flaginitpop) then
    if (Ninc == 1) then
      write(*, '(/, "     1 incident energy (LAB):"/)')
    else
      write(*, '(/, i6, " incident energies (LAB):"/)') Ninc
    endif
    do i = 1, Ninc
      if (eninc(i) < 0.001) then
        write(*, '(1x, es10.3)') eninc(i)
      else
        write(*, '(1x, f10.3)') eninc(i)
      endif
    enddo
  else
!
! Initial population distribution
!
    write(*, '(/, " Initial population distribution - Bins: ", i4, " Spins: ", i3, " Maximum excitation energy:", f12.5/)') &
 &    npopE, npopJ, eninc(1)
    if (npopJ == 0) then
      write(*, '("    Ex     Population "/)')
      do i = 1, npopE
        write(*, '(2es10.3)') EdistE(i), PdistE(i)
      enddo
    else
      write(*, '(" Parity   Ex ", 11("      J=", i2)/)') (J, J = 0, 10)
      do parity = -1, 1, 2
        do i = 1, npopE
          write(*, '(i6,12es10.3)') parity, EdistE(i), (PdistJP(i, J, parity), J = 0, 10)
        enddo
      enddo
    endif
  endif
  write(*, '(/" Q-values for binary reactions:"/)')
  do type = 0, 6
    if (parskip(type)) cycle
    write(*, '(" Q(", a1, ",", a1, "):", f9.5)') parsym(k0), parsym(type), Q(type)
  enddo
  if (k0 > 1) write(*, '(/," Coulomb barrier:",f9.5)') coulbar(k0)
  if (flagcheck) call arraysize
!
! * Write nuclear structure parameters for target and compound nucleus *
!
! levelsout    : subroutine for output of discrete levels
! densityout   : subroutine for output of level density parameters
! fissionparout: subroutine for output for fission parameters
!
  if (flaglevels) call levelsout(Zix, Nix)
  if (flagdensity) call densityout(Zix, Nix)
  if (flagfisout) call fissionparout(Zix, Nix)
  strucwrite(Zix, Nix) = .true.
  if (parskip(0)) return
  if (k0 /= 0) then
    if (flaglevels) call levelsout(Zcomp, Ncomp)
    if (flagdensity .or. filedensity) call densityout(Zcomp, Ncomp)
    if (flagfisout) call fissionparout(Zcomp, Ncomp)
    strucwrite(Zcomp, Ncomp) = .true.
  endif
  return
end subroutine mainout
! Copyright A.J. Koning 2025
