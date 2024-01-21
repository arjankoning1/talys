function evap(rn)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : evap=0 determines the number of evaporated neutrons
!
! Author    : Marieke Duijvestijn
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
! Variables for Brosa model
!   amm    ! parameter for neck rupture
!   ess    ! parameter for neck rupture
!   zee    ! parameter for neck rupture
!
! *** Declaration of local data
!
  implicit none
  real(sgl) :: an    ! neutron level density parameter
  real(sgl) :: an1   ! help variable
  real(sgl) :: b     ! parameter for slope
  real(sgl) :: bn    ! help variable
  real(sgl) :: bn1   ! help variable
  real(sgl) :: dum   ! dummy value
  real(sgl) :: dumm  ! help variable
  real(sgl) :: evap  ! Brosa evaporation function
  real(sgl) :: rn    ! help variable
!
! **********************************************************************
!
! evap: Brosa evaporation function
!
  call bdef(amm, zee, 0., dum, dumm, b)
  an = amm - rn
  call bdef(an, zee, 0., dum, dumm, bn)
  an1 = an - 1.
  call bdef(an1, zee, 0., dum, dumm, bn1)
  evap = ess - 0.5 * (bn + bn1) + b - 0.621 * rn * sqrt(rn + 1.)
  return
end function evap
! Copyright A.J. Koning 2021
