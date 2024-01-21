subroutine matrix(A, n)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Matrix element for exciton model
!
! Author    : Arjan Koning and Marieke Duijvestijn
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
! Variables for preequilibrium
!   flag2comp      ! flag for two - component pre - equilibrium model
!   M2constant     ! constant for matrix element in exciton model
!   M2limit        ! constant for asymptotic value for matrix element
!   M2shift        ! constant for energy shift for matrix element
!   preeqadjust    ! logical for energy - dependent pre - eq adjustment
!   preeqmode      ! designator for pre - equilibrium model
!   Rnunu          ! ratio for two - component matrix element
!   Rnupi          ! ratio for two - component matrix element
!   Rpinu          ! ratio for two - component matrix element
!   Rpipi          ! ratio for two - component matrix element
! Variables for main input
!   k0             ! index of incident particle
! Constants
!   parA           ! mass number of particle
! Variables for preequilibrium
!   Ecomp          ! total energy of composite system
! Variables for exciton model
!   M2             ! square of matrix element
!   M2nunu         ! square of neutron - neutron matrix element
!   M2nupi         ! square of neutron - proton matrix element
!   M2pinu         ! square of proton - neutron matrix element
!   M2pipi         ! square of proton - proton matrix element
!   Wompfac        ! adjustable constant for OMP based transition rates
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: key       ! keyword
  integer           :: A         ! mass number of target nucleus
  integer           :: aproj     ! mass number of particle
  integer           :: n         ! exciton number
  real(sgl)         :: factor    ! multiplication factor
  real(sgl)         :: M2c       ! constant for matrix element in exciton model (here used  for MSD model)
!
! ****************** Energy dependent matrix element *******************
!
! adjust     : subroutine for energy-dependent parameter adjustment
!
! We use a parameterization for the matrix element that is reliable between 7 and 200 MeV, see A.J. Koning and M.C. Duijvestijn,
! Nucl. Phys. A744, 15 (2004).
! For the one component exciton model the constant for the matrix element needs to be corrected for two-component effects.
! Empirically, we find that we need to multiply the matrix element by a ratio of 0.50.
! We generally use the exact numerical solution (preeqmode 2) for the transition rates.
! If the analytical solution is used (preeqmode 1), we need to correct the squared matrix element by 20%.
! Finally, the expression is generalized for complex particle emission according to C. Kalbach, Phys. Rev. C00, 004600 (2005).
!
  aproj = max(parA(k0), 1)
  if (preeqadjust) then
    key = 'm2constant'
    call adjust(Ecomp, key, 0, 0, 0, 0, factor)
    M2c = factor * M2constant
  else
    M2c = M2constant
  endif
  M2 = M2c / (A **3) * aproj * (M2limit * 7.48 + 4.62e5 / ((Ecomp / (n * aproj) + M2shift * 10.7) **3))
  if (preeqmode == 1) M2 = 1.20 * M2
  if (flag2comp) then
    M2pipi = Rpipi * M2
    M2nunu = Rnunu * M2
    M2pinu = Rpinu * M2
    M2nupi = Rnupi * M2
  else
    M2 = M2 * 0.50
  endif
!
! ************* Constants for optical model transition rates ***********
!
  if (preeqmode == 3) then
    Wompfac(1) = M2c * 0.55 / (1. + 2. * Rpinu)
    Wompfac(2) = M2c * 0.55 * 2. * Rpinu / (1. + 2. * Rpinu)
    Wompfac(0) = 0.5 * (Wompfac(1) + Wompfac(2))
  endif
  return
end subroutine matrix
! Copyright A.J. Koning 2021
