subroutine cm2lab(type)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Recoil for binary reaction
!
! Author    : Stephane Hilaire
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl            ! single precision kind
! Variables for recoil
!   cosejcm1       ! CM recoil energy cosine corresponding to (E1, ang1)
!   cosejcm2       ! CM recoil energy cosine corresponding to (E1, ang2)
!   cosejlab11     ! LAB ejectile angle cosine corresponding to (E1, ang1)
!   cosejlab12     ! LAB ejectile angle cosine corresponding to (E1, ang2)
!   cosejlab21     ! LAB ejectile angle cosine corresponding to (E2, ang1)
!   cosejlab22     ! LAB ejectile angle cosine corresponding to (E2, ang2)
!   cosreclab11    ! LAB recoil angle cosine corresponding to (E1, ang1)
!   cosreclab12    ! LAB recoil angle cosine corresponding to (E1, ang2)
!   cosreclab21    ! LAB recoil angle cosine corresponding to (E2, ang1)
!   cosreclab22    ! LAB recoil angle cosine corresponding to (E2, ang2)
!   Eejcm1         ! Lower limit of CM ejectile energy bin
!   Eejcm2         ! Upper limit of CM ejectile energy bin
!   Eejlab11       ! LAB ejectile energy corresponding to (E1, ang1)
!   Eejlab12       ! LAB ejectile energy corresponding to (E1, ang2)
!   Eejlab21       ! LAB ejectile energy corresponding to (E2, ang1)
!   Eejlab22       ! LAB ejectile energy corresponding to (E2, ang2)
!   ejectmass      ! Ejectile mass
!   Ereclab11      ! LAB recoil energy corresponding to (E1, ang1)
!   Ereclab12      ! LAB recoil energy corresponding to (E1, ang2)
!   Ereclab21      ! LAB recoil energy corresponding to (E2, ang1)
!   Ereclab22      ! LAB recoil energy corresponding to (E2, ang2)
!   recoilmass     ! Recoil mass
!   sinejcm1       ! CM recoil energy cosine corresponding to (E2, ang1)
!   sinejcm2       ! CM recoil energy cosine corresponding to (E2, ang2)
!   sinejlab11     ! LAB ejectile angle sine corresponding to (E1, ang1)
!   sinejlab12     ! LAB ejectile angle sine corresponding to (E1, ang2)
!   sinejlab21     ! LAB ejectile angle sine corresponding to (E2, ang1)
!   sinejlab22     ! LAB ejectile angle sine corresponding to (E2, ang2)
!   sinreclab11    ! LAB recoil angle sine corresponding to (E1, ang1)
!   sinreclab12    ! LAB recoil angle sine corresponding to (E1, ang2)
!   sinreclab21    ! LAB recoil angle sine corresponding to (E2, ang1)
!   sinreclab22    ! LAB recoil angle sine corresponding to (E2, ang2)
!   vcm            ! Compound nucleus velocity
!   vejcm1         ! velocity corresponding to Eejcm1
!   vejcm2         ! velocity corresponding to Eejcm2
!   vreccm1        ! Recoil velocity corresponding to Eejcm1
!   vreccm2        ! Recoil velocity corresponding to Eejcm2
!
! *** Declaration of local data
!
  implicit none
  integer   :: type             ! particle type
  real(sgl) :: coeff            ! help variable
  real(sgl) :: vlab             ! LAB recoil velocity
  real(sgl) :: vlab2            ! square of LAB recoil velocity
  real(sgl) :: vlabx            ! LAB recoil velocity projection on x-axis // to CM  velocity
  real(sgl) :: vlaby            ! LAB recoil velocity projection on y-axis
!
! EJECTILE TREATMENT
!
! The ejectile LAB angles and energies corresponding to the CM points are deduced
!
  if (type == 0.) then
    Eejlab11 = Eejcm1
    Eejlab12 = Eejcm1
    Eejlab21 = Eejcm2
    Eejlab22 = Eejcm2
    cosejlab11 = cosejcm1
    cosejlab12 = cosejcm2
    cosejlab21 = cosejcm1
    cosejlab22 = cosejcm2
    sinejlab11 = sinejcm1
    sinejlab12 = sinejcm2
    sinejlab21 = sinejcm1
    sinejlab22 = sinejcm2
  else
    coeff = 0.5 * ejectmass
    vlabx = vcm + vejcm1 * cosejcm1
    vlaby = vejcm1 * sinejcm1
    vlab2 = vlabx **2 + vlaby **2
    vlab = max(sqrt(vlab2), 1.e-20)
    cosejlab11 = vlabx / vlab
    sinejlab11 = vlaby / vlab
    Eejlab11 = coeff * vlab2
    vlabx = vcm + vejcm1 * cosejcm2
    vlaby = vejcm1 * sinejcm2
    vlab2 = vlabx **2 + vlaby **2
    vlab = max(sqrt(vlab2), 1.e-20)
    cosejlab12 = vlabx / vlab
    sinejlab12 = vlaby / vlab
    Eejlab12 = coeff * vlab2
    vlabx = vcm + vejcm2 * cosejcm1
    vlaby = vejcm2 * sinejcm1
    vlab2 = vlabx **2 + vlaby **2
    vlab = max(sqrt(vlab2), 1.e-20)
    cosejlab21 = vlabx / vlab
    sinejlab21 = vlaby / vlab
    Eejlab21 = coeff * vlab2
    vlabx = vcm + vejcm2 * cosejcm2
    vlaby = vejcm2 * sinejcm2
    vlab2 = vlabx **2 + vlaby **2
    vlab = max(sqrt(vlab2), 1.e-20)
    cosejlab22 = vlabx / vlab
    sinejlab22 = vlaby / vlab
    Eejlab22 = coeff * vlab2
  endif
!
! RECOIL TREATMENT
!
! The Recoil LAB angles and energies corresponding to the CM points are deduced
!
  coeff = 0.5 * recoilmass
  vlabx = vcm - vreccm1 * cosejcm1
  vlaby = - vreccm1 * sinejcm1
  vlab2 = vlabx **2 + vlaby **2
  vlab = max(sqrt(vlab2), 1.e-20)
  cosreclab11 = vlabx / vlab
  sinreclab11 = - sinejcm1
  Ereclab11 = coeff * vlab2
  vlabx = vcm - vreccm1 * cosejcm2
  vlaby = - vreccm1 * sinejcm2
  vlab2 = vlabx **2 + vlaby **2
  vlab = max(sqrt(vlab2), 1.e-20)
  cosreclab12 = vlabx / vlab
  sinreclab12 = - sinejcm2
  Ereclab12 = coeff * vlab2
  vlabx = vcm - vreccm2 * cosejcm1
  vlaby = - vreccm2 * sinejcm1
  vlab2 = vlabx **2 + vlaby **2
  vlab = max(sqrt(vlab2), 1.e-20)
  cosreclab21 = vlabx / vlab
  sinreclab21 = - sinejcm1
  Ereclab21 = coeff * vlab2
  vlabx = vcm - vreccm2 * cosejcm2
  vlaby = - vreccm2 * sinejcm2
  vlab2 = vlabx **2 + vlaby **2
  vlab = max(sqrt(vlab2), 1.e-20)
  cosreclab22 = vlabx / vlab
  sinreclab22 = - sinejcm2
  Ereclab22 = coeff * vlab2
  return
end subroutine cm2lab
! Copyright A.J. Koning 2021
