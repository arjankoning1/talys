subroutine adjustF(E, nrange, en1, en2, Dr, sr, factor)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Local parameter adjustment
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
!   sgl         ! single precision kind
! All global variables
!   numrange    ! number of energy ranges
!
! *** Declaration of local data
!
  implicit none
  integer   :: nr            ! number of radial grid point
  integer   :: nrange        ! number of energy ranges for local adjustment
  real(sgl) :: D             ! depth of local adjustment
  real(sgl) :: Dr(numrange)  ! depth of local adjustment
  real(sgl) :: E             ! incident energy
  real(sgl) :: elow          ! help variable
  real(sgl) :: emid          ! help variable
  real(sgl) :: en1(numrange) ! start energy of local adjustment
  real(sgl) :: en2(numrange) ! end energy of local adjustment
  real(sgl) :: eup           ! help variable
  real(sgl) :: expo          ! help variable
  real(sgl) :: factor        ! multiplication factor
  real(sgl) :: offset        ! offset to guarantee continuity
  real(sgl) :: sigma         ! help variable
  real(sgl) :: sr(numrange)  ! variance of local adjustment
!
! ************************* Woods-Saxon function ***********************
!
  factor = 1.
  do nr = 1, nrange
    elow = en1(nr)
    eup = en2(nr)
    if (E > elow .and. E < eup) then
      emid = 0.5 * (elow + eup)
      D = 0.01 * Dr(nr)
      sigma = sr(nr)
      if (sigma == 0.) sigma = (eup - emid) / 2.
      expo = (eup - emid) **2 / (2. * sigma **2)
      offset = 0.
      if (expo <= 80.) offset = - D * exp( - expo)
      expo = (E - emid) **2 / (2. * sigma **2)
      factor = 1. + offset
      if (expo <= 80.) factor = 1. + D * exp( - expo) + offset
      return
    endif
  enddo
  return
end subroutine adjustF
! Copyright A.J. Koning 2021
