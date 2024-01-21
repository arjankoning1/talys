subroutine endfinfo
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Info for ENDF-6 file
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
! Variables for URR
!   flagurrendf    ! flag for URR info to ENDF
! Variables for basic reaction
!   flagendfdet    ! flag for detailed ENDF - 6 information per channel
!   flagrecoil     ! flag for calculation of recoils
! Variables for output
!   flagblock      ! flag to block spectra, angle and gamma files
! Variables for input energies
!   eninc          ! incident energy in MeV
!   Ninc           ! number of incident energies
! Variables for main input
!   Atarget        ! mass number of target nucleus
!   k0             ! index of incident particle
!   Ltarget        ! excited level of target
!   Ztarget        ! charge number of target nucleus
! Variables for discrete levels
!   nlevmax        ! maximum number of included discrete levels for target
! Variables for energy grid
!   Nrescue        ! number of energies for adjustment factors
! Variables for nuclides
!   targetE        ! excitation energy of target
!   targetspin     ! spin of target
!   tarmass        ! mass of target nucleus
! Variables for levels
!   Liso           ! isomeric number of target
!
! *** Declaration of local data
!
  implicit none
  logical            :: flagel            ! flag for elastic scattering
  character(len=1)   :: yesno             ! y or n function
  character(len=132) :: endfenergyfile    ! file with energies for ENDF file
  integer            :: nen               ! energy counter
!
! ****************** Reaction information for ENDF-6 file **************
!
! yesno      : y or n function
!
  open (unit = 1, file = 'tefal.inf', status = 'replace')
  write(1, '(i3, "       : projectile type")') k0
  write(1, '(i3, "       : Z of target")') Ztarget
  write(1, '(i3, "       : A of target")') Atarget
  write(1, '(f10.6, ": mass of target in a.m.u.")') tarmass
  write(1, '(f4.1, "      : spin of target")') targetspin
  write(1, '(f10.6, ": energy of target")') targetE
  write(1, '(2i4, "  : level and isomeric number of target")') Ltarget, Liso
  endfenergyfile = 'energies.endf'
  open (unit = 2, file = endfenergyfile, status = 'replace')
  do nen = 1, Ninc
    write(2, '(es12.5)') eninc(nen)
  enddo
  close (unit = 2)
  write(1, '(i4, "      : number of incident energies")') Ninc
  write(1, '(i4, "      : number of discrete levels")') nlevmax
  write(1, '(a, " : file with incident energies")') trim(endfenergyfile)
  write(1, '(a1, "         : detailed ENDF-6 information", " per channel")') yesno(flagendfdet)
  if (Nrescue(1, - 1) /= 0) then
    flagel = .false.
  else
    flagel = .true.
  endif
  write(1, '(a1, "         : keep elastic cross section", " in normalization")') yesno(flagel)
  write(1, '(a1, "         : recoils")') yesno(flagrecoil)
  write(1, '(a1, "         : urr")') yesno(flagurrendf)
  write(1, '(a1, "         : block")') yesno(flagblock)
  close (unit = 1)
  return
end subroutine endfinfo
! Copyright A.J. Koning 2021
