function nubarWahl(Ein)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Wahl systematics for nubar
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
! Definition of single and double precision variables
!   sgl      ! single precision kind
! All global variables
!   Cnubar1    ! adjustable parameter for nubar constant value
!   Cnubar2    ! adjustable parameter for nubar energy slope
! Variables for main input
!   Ainit    ! mass number of initial compound nucleus
!   Ninit    ! neutron number of initial compound nucleus
!   Zinit    ! charge number of initial compound nucleus
! Variables for masses
!   S        ! separation energy
!
! *** Declaration of local data
!
  implicit none
  real(sgl) :: Bn        !
  real(sgl) :: Ein       ! incident energy
  real(sgl) :: Fn        !
  real(sgl) :: Fz        !
  real(sgl) :: nubarWahl !
  real(sgl) :: term1     ! help variable
  real(sgl) :: term2     ! help variable
  real(sgl) :: Th        !
!
! ************************ Wahl systematics ****************************
!
! A.C. Wahl, "Systematics of fission-product yields", LA-13926 (2002), Eq. (21)
!
  Bn = S(0, 0, 1)
  if (mod(Zinit, 2) == 0) then
    Fz = 1.
  else
    Fz = - 1.
  endif
  if (mod(Ninit, 2) == 0) then
    Fn = 1.
  else
    Fn = - 1.
  endif
  term1 = 2.286 + 0.147 * (Zinit - 92) + 0.054 * (Ainit - 236) + 0.040 * (2. - Fz - Fn)
  if (Ein > 0.) then
    Th = 11.47 - 0.166 * Zinit * Zinit / real(Ainit) + 0.093 * (2. - Fz - Fn) - Bn
    term2 = (0.145 - 0.0043 * (Ainit - 236)) * (Ein - Th)
  else
    term2 = 0.
  endif
  nubarWahl = Cnubar1 * term1 + Cnubar2 * term2
  return
end function nubarWahl
! Copyright A.J. Koning 2021
