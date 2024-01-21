function thill(ehw, bhw, whw)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Hill-Wheeler penetrability
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
!   sgl    ! single precision kind
! Constants
!   twopi  ! 2 * pi
!
! *** Declaration of local data
!
  implicit none
  real(sgl) :: bhw   ! fission barrier height
  real(sgl) :: ehw   ! E(compound nucleus) - E(transition state)
  real(sgl) :: expo  ! help variable
  real(sgl) :: thill ! Hill-Wheeler penetrability
  real(sgl) :: whw   ! fission barrier width
!
! ************************* Hill-Wheeler formula ***********************
!
  thill = 0.
  expo = twopi * (ehw - bhw) / whw
  if (expo >  - 80.) thill = 1. / (1. + exp( - expo))
  return
end function thill
! Copyright A.J. Koning 2021
