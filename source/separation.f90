subroutine separation
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Separation energies
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
! Variables for numerics
!   maxN       ! maximal number of neutrons away from initial compound nucleus
!   maxZ       ! maximal number of protons away from initial compound nucleus
! Constants
!   amu        ! atomic mass unit in MeV
!   excmass    ! mass excess of particle in a.m.u.
!   parN       ! neutron number of particle
!   parZ       ! charge number of particle
! Variables for masses
!   dumexc     ! theoretical mass excess from Duflo - Zuker formula
!   expmexc    ! experimental mass excess
!   S          ! separation energy
!   thmexc     ! theoretical mass excess
!
! *** Declaration of local data
!
  implicit none
  integer :: Nix               ! neutron number index for residual nucleus
  integer :: Nr                ! neutron number index for residual nucleus
  integer :: type              ! particle type
  integer :: Zix               ! charge number index for residual nucleus
  integer :: Zr                ! charge number index for residual nucleus
!
! For consistency, separation energies are always calculated using two nuclear masses of the same type,
! i.e. both experimental or both theoretical.
! Hence if (Zix,Nix) is in the AME2020 table but (Zr,Nr) is not, for both nuclides the theoretical masses are used.
! For the calculation of separation energies, Zix and Nix act as compound nucleus indices.
!
  do Zix = 0, maxZ+2
    do Nix = 0, maxN + 2
      do type = 1, 6
        Zr = Zix + parZ(type)
        Nr = Nix + parN(type)
        if (expmexc(Zix, Nix) /= 0. .and. expmexc(Zr, Nr) /= 0.) then
          S(Zix, Nix, type) = expmexc(Zr, Nr) - expmexc(Zix, Nix) + excmass(type) * amu
          cycle
        endif
        if (thmexc(Zix, Nix) /= 0. .and. thmexc(Zr, Nr) /= 0.) then
          S(Zix, Nix, type) = thmexc(Zr, Nr) - thmexc(Zix, Nix) + excmass(type) * amu
        else
          S(Zix, Nix, type) = dumexc(Zr, Nr) - dumexc(Zix, Nix) + excmass(type) * amu
        endif
      enddo
    enddo
  enddo
  return
end subroutine separation
! Copyright A.J. Koning 2021
