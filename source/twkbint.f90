function twkbint(efis, ibar, Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Interpolation of WKB penetrability
!
! Author    : Stephane Hilaire and Guillaume Scamps
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
!   nbinswkb    ! integration step for WKB calculation
!   Twkb        ! transmission coefficient of WKB potential
!   Uwkb        ! energy of WKB potential
!
! *** Declaration of local data
!
  implicit none
  integer   :: ibar             ! fission barrier
  integer   :: nen              ! energy counter
  integer   :: Nix              ! neutron number index for residual nucleus
  integer   :: Zix              ! charge number index for residual nucleus
  real(sgl) :: Ea               ! start energy of local adjustment
  real(sgl) :: Eb               ! end energy of local adjustment
  real(sgl) :: efis             ! parameter for fission
  real(sgl) :: Ewkb(0:nbinswkb) ! energy
  real(sgl) :: Ta               ! transmission coefficient
  real(sgl) :: Tb               ! temperature for QRPA
  real(sgl) :: Tf               ! temperature
  real(sgl) :: twkbint          ! WKB penetrability
!
! ************************* Interpolation ******************************
!
  do nen = 0, nbinswkb
    Ewkb(nen) = Uwkb(Zix, Nix, nen)
  enddo
  if (efis > Ewkb(nbinswkb)) then
    twkbint = 1.
  else
    call locate(Ewkb, 0, nbinswkb, efis, nen)
    Ea = Ewkb(nen)
    Eb = Ewkb(nen + 1)
    Ta = Twkb(Zix, Nix, nen, ibar)
    Tb = Twkb(Zix, Nix, nen + 1, ibar)
    if (Ea == Eb) then
      twkbint = Ta
      return
    endif
    if (Ta > 0 .and. Tb > 0) then
      Ta = log(Ta)
      Tb = log(Tb)
      call pol1(Ea, Eb, Ta, Tb, efis, Tf)
      twkbint = exp(Tf)
    else
      call pol1(Ea, Eb, Ta, Tb, efis, Tf)
      twkbint = Tf
    endif
  endif
  return
end function twkbint
! Copyright A.J. Koning 2021
