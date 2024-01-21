subroutine bdef(amass, zchar, epscloc, cou, sym, b)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Computes the binding energy of a deformed nucleus with eccentricity
!
! Author    : Marieke Duijvestijn
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               sgl                ! single precision kind
!
! *** Declaration of local data
!
  real(sgl) :: a3w     ! mass**1/3
  real(sgl) :: amass   ! mass number of residual nucleus
  real(sgl) :: b       ! parameter for slope
  real(sgl) :: cou     ! Coulomb term
  real(sgl) :: epscloc ! Brosa value
  real(sgl) :: fcoul   ! Coulomb self energy factor
  real(sgl) :: fsurf   ! function for form factor for the surface energy
  real(sgl) :: hsc     ! Brosa term
  real(sgl) :: sym     ! symmetry term
  real(sgl) :: zchar   ! charge number of residual nucleus
  external fcoul, fsurf
!
! **********************************************************************
!
  a3w = amass **(1. / 3.)
  cou = 0.7053 / a3w * fcoul(epscloc) - 1.153 / amass
  hsc = 15.4941 * amass - 17.9439 * a3w **2 * fsurf(epscloc)
  sym = 1.7826 / amass **2 * hsc
  b = - hsc + sym * (amass - zchar - zchar) **2 + cou * zchar **2
  return
end subroutine bdef
! Copyright A.J. Koning 2021
