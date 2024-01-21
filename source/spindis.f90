function spindis(sc, Rspin)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Wigner spin distribution
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               sgl                ! single precision kind
!
! *** Declaration of local data
!
  real(sgl) :: sc      ! spin cutoff factor
  real(sgl) :: Rspin   ! residual spin
  real(sgl) :: sigma22 ! 2 * spin cutoff factor
  real(sgl) :: spindis ! Wigner spin distribution
!
! *********************** Wigner formula ******************************
!
  sigma22 = 2.*sc
  spindis = (2. * Rspin + 1.) / sigma22 * exp( - (Rspin + 0.5) **2 / sigma22)
  return
end function spindis
! Copyright A.J. Koning 2021
