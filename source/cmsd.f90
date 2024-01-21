function cmsd(Ein)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Normalization factor for MSD
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
!   sgl           ! single precision kind
! Variables for preequilibrium
!   M2constant    ! constant for matrix element in exciton model
!   preeqadjust   ! logical for energy - dependent pre - eq adjustment
! Variables for main input
!   Atarget       ! mass number of target nucleus
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: key       ! keyword
  real(sgl)         :: cmsd      ! normalization factor for MSD
  real(sgl)         :: Ein       ! incident energy
  real(sgl)         :: factor    ! multiplication factor
  real(sgl)         :: M2c       ! constant for matrix element in exciton model (here used  for MSD model)
!
! ***************** Normalization factor for MSD ***********************
!
! adjust     : subroutine for energy-dependent parameter adjustment
!
  if (preeqadjust) then
    key = 'm2constant'
    call adjust(Ein, key, 0, 0, 0, 0, factor)
    M2c = factor * M2constant
  else
    M2c = M2constant
  endif
  cmsd = M2c * 3.3e6 / (Atarget **3) / Ein
  return
end function cmsd
! Copyright A.J. Koning 2021
