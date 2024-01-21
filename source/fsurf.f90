function fsurf(epscloc)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Form factor for the surface energy
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
  real(sgl) :: dum     ! dummy value
  real(sgl) :: epscloc ! Brosa value
  real(sgl) :: fsurf   ! function for form factor for the surface energy
!
! **********************************************************************
!
! fsurf: function for form factor for the surface energy
!
  if (epscloc) 10, 20, 30
 10   dum = sqrt(1. + epscloc **2)
  fsurf = (alog( - epscloc + dum) - epscloc * dum) / ( - 2. * epscloc * dum **(1. / 3.))
  return
 20   fsurf = 1.
  return
 30   dum = sqrt(1. - epscloc **2)
  fsurf = (asin(epscloc) + epscloc * dum) / (2. * epscloc * dum **(1. / 3.))
  return
end function fsurf
! Copyright A.J. Koning 2021
