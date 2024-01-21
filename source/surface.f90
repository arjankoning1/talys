function surface(type, elab)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Surface effects in exciton model
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
! Variables for main input
!   Ainit    ! mass number of initial compound nucleus
! Constants
!   onethird ! 1 / 3
!
! *** Declaration of local data
!
  implicit none
  integer   :: type    ! particle type
  real(sgl) :: elab    ! incident energy
  real(sgl) :: surface ! well depth for first hole
!
! ************************* Kalbach formula ****************************
!
  if (type == 1) then
    surface = 12. + 26. * elab **4 / (elab **4 + (245. / (Ainit **onethird)) **4)
  else
    surface = 22. + 16. * elab **4 / (elab **4 + (450. / (Ainit **onethird)) **4)
  endif
  return
end function surface
! Copyright A.J. Koning 2021
