function twkbphaseint(efis, ibar, Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Interpolation of direct WKB penetrability
!
! Author    : Guillaume Scamps
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl         ! single precision kind
! Variables for WKB
!   nbinswkb     ! integration step for WKB calculation
!   Twkbphase    ! transmission coefficient of WKB potential
!   Uwkb         ! energy of WKB potential
!
! *** Declaration of local data
!
  implicit none
  integer   :: ibar             ! fission barrier
  integer   :: nen              ! energy counter
  integer   :: Nix              ! neutron number index for residual nucleus
  integer   :: Zix              ! charge number index for residual nucleus
  real(sgl) :: Ea               ! begin energy
  real(sgl) :: Eb               ! end energy
  real(sgl) :: efis             ! energy of fission
  real(sgl) :: Ewkb(0:nbinswkb) ! energy
  real(sgl) :: Ta               ! transmission coefficient
  real(sgl) :: Tb               ! transmission coefficient
  real(sgl) :: Tf               ! transmission coefficient
  real(sgl) :: twkbphaseint     ! function for interpolation of direct WKB penetrability
!
! ************************* Interpolation ******************************
!
  do nen = 0, nbinswkb
    Ewkb(nen) = Uwkb(Zix, Nix, nen)
  enddo
  if (efis > Ewkb(nbinswkb)) then
    twkbphaseint = 1.
  else
    call locate(Ewkb, 0, nbinswkb, efis, nen)
    Ea = Ewkb(nen)
    Eb = Ewkb(nen + 1)
    Ta = Twkbphase(Zix, Nix, nen, ibar)
    Tb = Twkbphase(Zix, Nix, nen + 1, ibar)
    call pol1(Ea, Eb, Ta, Tb, efis, Tf)
    twkbphaseint = Tf
  endif
  return
end function twkbphaseint
! Copyright A.J. Koning 2021
