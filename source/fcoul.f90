function fcoul(epscloc)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Form factor for the coulomb self energy
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
  real(sgl) :: epscloc            ! Brosa value
  real(sgl) :: fcoul              ! Coulomb self energy factor
!
! **********************************************************************
!
  if (epscloc) 10, 20, 30
 10   fcoul = (1. + epscloc **2) **(1. / 3.) / epscloc * atan(epscloc)
  return
 20   fcoul = 1.
  return
 30   fcoul = (1. - epscloc **2) **(1. / 3.) / (2. * epscloc) * alog((1. + epscloc) / (1. - epscloc))
  return
end function fcoul
! Copyright A.J. Koning 2021
