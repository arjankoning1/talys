subroutine timer
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of execution time
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               sgl                  ! single precision kind
!
! *** Declaration of local data
!
  integer   :: hour      ! number of hours
  integer   :: hundred   ! number of 1/100th of seconds
  integer   :: minute    ! number of minutes
  integer   :: second    ! number of seconds
  real(sgl) :: etime     ! time function
  real(sgl) :: tarray(2) ! help variable
  real(sgl) :: time      ! time
!
! ****** Get elapsed time in seconds from beginning of execution *******
!
! etime  : time function
!
! The returned time should be "charge time" (e.g., cp+pp+sys). This could be machine dependent.
!
  time = etime(tarray)
  hour = int(time / 3600.)
  minute = int((time - hour * 3600) / 60.)
  second = int(time - hour * 3600 - minute * 60)
  hundred = int(100 * (time - int(time)))
  write(*, '(/" Execution time:", i3, " hours ", i2, " minutes ", i2, ".", i2.2, " seconds ")') hour, minute, second, hundred
  write(*, '(/" The TALYS team congratulates you with this successful calculation.")')
  return
end subroutine timer
! Copyright A.J. Koning 2021
