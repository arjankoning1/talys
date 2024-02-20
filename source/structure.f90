subroutine structure(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Nuclear structure parameters
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
! 2022-02-15: Added flagnffit keyword
! 2022-04-18: Added flagngfit keyword
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! All global variables
!   numNph          ! maximum number of neutrons away from the initial compound nucleus
!   numZph          ! maximum number of protons away from the initial compound nucleus
! Variables for basic reaction
!   flagendf        ! flag for information for ENDF - 6 file
!   flagfit         ! flag for using fitted nuclear model parameters
!   flagpartable    ! flag for output of model parameters on separate file
! Variables for compound reactions
!   flagcomp        ! flag for compound angular distribution calculation
! Variables for preequilibrium
!   phmodel         ! particle - hole state density model
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
!   Ztarget         ! charge number of target nucleus
! Variables for fission
!   flagfission     ! flag for fission
! Variables for level density
!   ldmodel         ! level density model
! Variables for OMP
!   alphaomp        ! alpha optical model
!   flagomponly     ! flag to execute ONLY an optical model calculation
!   flagjlm         ! flag for using semi - microscopic JLM OMP
!   flagompall      ! flag for new optical model calculation for all residual
! Variables for nuclides
!   parinclude      ! logical to include outgoing particle
!   primary         ! flag to designate primary (binary) reaction
! Constants
!   parN            ! neutron number of particle
!   parZ            ! charge number of particle
!
! *** Declaration of local data
!
  implicit none
  integer :: Nix              ! neutron number index for residual nucleus
  integer :: Zix              ! charge number index for residual nucleus
!
! ********** Calculate and read various nuclear parameters *************
!
! levels        : subroutine for discrete levels
! gammadecay    : subroutine for scheme for discrete gamma decay
! deformpar     : subroutine for deformation parameters
! resonancepar  : subroutine for s-wave resonance parameters
! gammapar      : subroutine for gamma ray parameters residual nuclides
! omppar        : subroutine for optical model parameters
! radialtable   : subroutine for tabulated radial matter densities
! fissionpar    : subroutine for fission parameters
! densitypar    : subroutine for level density parameters
! densitytable  : subroutine for tabulated level densities
! densitymatch  : subroutine for level density matching solution
! phdensitytable: subroutine for tabulated particle-hole state densities
! thermalxs     : subroutine for cross sections at thermal energies
! partable      : subroutine to write model parameters per nucleus to separate file
!
! All the nuclear structure info is read and/or calculated for the nucleus under consideration.
!
  call levels(Zix, Nix)
  if (flagendf .or. flagbasic) then
    call levelsout(Zix, Nix)
    call gammadecay(Zix, Nix)
  endif
  call deformpar(Zix, Nix)
  if (flagfit .and. Zix == 0 .and. Nix == 0) call xsfit(Ztarget, Atarget)
  if (parinclude(0) .or. flagcomp) then
    call resonancepar(Zix, Nix)
    call gammapar(Zix, Nix)
  endif
  if ((Zix <= 2 .and. Nix <= 2) .or. flagompall) call omppar(Zix, Nix)
  if (flagjlm .or. alphaomp >= 3 .and. alphaomp <= 5) call radialtable(Zix, Nix)
  if (flagomponly .and. .not. flagcomp) return
  if (flagfission) call fissionpar(Zix, Nix)
  call densitypar(Zix, Nix)
  if (ldmodel(Zix, Nix) >= 4) call densitytable(Zix, Nix)
  call densitymatch(Zix, Nix)
  if (phmodel == 2 .and. Zix <= numZph .and. Nix <= numNph) call phdensitytable(Zix, Nix)
  if (k0 == 1 .and. Zix == parZ(k0) .and. Nix == parN(k0)) call thermalxs
  if (flagpartable) call partable(Zix, Nix)
  return
end subroutine structure
! Copyright A.J. Koning 2021
