subroutine ompadjust(E, k)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Local optical model parameter adjustment
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
!   sgl            ! single precision kind
! Variables for OMP
!   avadjust       ! adjustable factor for OMP (default 1.)
!   avdadjust      ! adjustable factor for OMP (default 1.)
!   avsoadjust     ! adjustable factor for OMP (default 1.)
!   awadjust       ! adjustable factor for OMP (default 1.)
!   awdadjust      ! adjustable factor for OMP (default 1.)
!   awsoadjust     ! adjustable factor for OMP (default 1.)
!   d1adjust       ! adjustable factor for OMP (default 1.)
!   d2adjust       ! adjustable factor for OMP (default 1.)
!   d3adjust       ! adjustable factor for OMP (default 1.)
!   flagincadj     ! flag for OMP adjustment on incident channel also
!   ompadjustp     ! flag for local optical model parameter adjustment
!   rcadjust       ! adjustable factor for OMP (default 1.)
!   rvadjust       ! adjustable factor for OMP (default 1.)
!   rvdadjust      ! adjustable factor for OMP (default 1.)
!   rvsoadjust     ! adjustable factor for OMP (default 1.)
!   rwadjust       ! adjustable factor for OMP (default 1.)
!   rwdadjust      ! adjustable factor for OMP (default 1.)
!   rwsoadjust     ! adjustable factor for OMP (default 1.)
!   v1adjust       ! adjustable factor for OMP (default 1.)
!   v2adjust       ! adjustable factor for OMP (default 1.)
!   v3adjust       ! adjustable factor for OMP (default 1.)
!   v4adjust       ! adjustable factor for OMP (default 1.)
!   vso1adjust     ! adjustable factor for OMP (default 1.)
!   vso2adjust     ! adjustable factor for OMP (default 1.)
!   w1adjust       ! adjustable factor for OMP (default 1.)
!   w2adjust       ! adjustable factor for OMP (default 1.)
!   w3adjust       ! adjustable factor for OMP (default 1.)
!   w4adjust       ! adjustable factor for OMP (default 1.)
!   wso1adjust     ! adjustable factor for OMP (default 1.)
!   wso2adjust     ! adjustable factor for OMP (default 1.)
! Variables for ECIS
!   flaginvecis    ! logical for calculating inverse channel OMP
! Variables for optical model
!   Fav            ! adjustable factor for OMP (default 1.)
!   Favd           ! adjustable factor for OMP (default 1.)
!   Favso          ! adjustable factor for OMP (default 1.)
!   Faw            ! adjustable factor for OMP (default 1.)
!   Fawd           ! adjustable factor for OMP (default 1.)
!   Fawso          ! adjustable factor for OMP (default 1.)
!   Fd1            ! adjustable factor for OMP (default 1.)
!   Fd2            ! adjustable factor for OMP (default 1.)
!   Fd3            ! adjustable factor for OMP (default 1.)
!   Frc            ! adjustable factor for OMP (default 1.)
!   Frv            ! adjustable factor for OMP (default 1.)
!   Frvd           ! adjustable factor for OMP (default 1.)
!   Frvso          ! adjustable factor for OMP (default 1.)
!   Frw            ! adjustable factor for OMP (default 1.)
!   Frwd           ! adjustable factor for OMP (default 1.)
!   Frwso          ! adjustable factor for OMP (default 1.)
!   Fv1            ! adjustable factor for OMP (default 1.)
!   Fv2            ! adjustable factor for OMP (default 1.)
!   Fv3            ! adjustable factor for OMP (default 1.)
!   Fv4            ! adjustable factor for OMP (default 1.)
!   Fvso1          ! adjustable factor for OMP (default 1.)
!   Fvso2          ! adjustable factor for OMP (default 1.)
!   Fw1            ! adjustable factor for OMP (default 1.)
!   Fw2            ! adjustable factor for OMP (default 1.)
!   Fw3            ! adjustable factor for OMP (default 1.)
!   Fw4            ! adjustable factor for OMP (default 1.)
!   Fwso1          ! adjustable factor for OMP (default 1.)
!   Fwso2          ! adjustable factor for OMP (default 1.)
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: key       ! keyword
  integer           :: k         ! designator for particle
  real(sgl)         :: E         ! incident energy
  real(sgl)         :: factor    ! multiplication factor
!
! ******************* Create adjustable OMP factors *********************
!
! adjust     : subroutine for energy-dependent parameter adjustment
!
  if (ompadjustp(k) .and. (flagincadj .or. flaginvecis)) then
    key = 'v1adjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Fv1 = factor * v1adjust(k)
    key = 'v2adjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Fv2 = factor * v2adjust(k)
    key = 'v3adjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Fv3 = factor * v3adjust(k)
    key = 'v4adjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Fv4 = factor * v4adjust(k)
    key = 'rvadjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Frv = factor * rvadjust(k)
    key = 'avadjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Fav = factor * avadjust(k)
    key = 'w1adjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Fw1 = factor * w1adjust(k)
    key = 'w2adjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Fw2 = factor * w2adjust(k)
    key = 'w3adjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Fw3 = factor * w3adjust(k)
    key = 'w4adjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Fw4 = factor * w4adjust(k)
    key = 'rwadjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Frw = factor * rwadjust(k)
    key = 'awadjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Faw = factor * awadjust(k)
    key = 'rvdadjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Frvd = factor * rvdadjust(k)
    key = 'avdadjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Favd = factor * avdadjust(k)
    key = 'd1adjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Fd1 = factor * d1adjust(k)
    key = 'd2adjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Fd2 = factor * d2adjust(k)
    key = 'd3adjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Fd3 = factor * d3adjust(k)
    key = 'rwdadjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Frwd = factor * rwdadjust(k)
    key = 'awdadjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Fawd = factor * awdadjust(k)
    key = 'vso1adjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Fvso1 = factor * vso1adjust(k)
    key = 'vso2adjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Fvso2 = factor * vso2adjust(k)
    key = 'rvsoadjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Frvso = factor * rvsoadjust(k)
    key = 'avsoadjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Favso = factor * avsoadjust(k)
    key = 'wso1adjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Fwso1 = factor * wso1adjust(k)
    key = 'wso2adjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Fwso2 = factor * wso2adjust(k)
    key = 'rwsoadjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Frwso = factor * rwsoadjust(k)
    key = 'awsoadjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Fawso = factor * awsoadjust(k)
    key = 'rcadjust'
    call adjust(E, key, 0, 0, k, 0, factor)
    Frc = factor * rcadjust(k)
  else
    Fv1 = 1.
    Fv2 = 1.
    Fv3 = 1.
    Fv4 = 1.
    Frv = 1.
    Fav = 1.
    Fw1 = 1.
    Fw2 = 1.
    Fw3 = 1.
    Fw4 = 1.
    Frw = 1.
    Faw = 1.
    Frvd = 1.
    Favd = 1.
    Fd1 = 1.
    Fd2 = 1.
    Fd3 = 1.
    Frwd = 1.
    Fawd = 1.
    Fvso1 = 1.
    Fvso2 = 1.
    Frvso = 1.
    Favso = 1.
    Fwso1 = 1.
    Fwso2 = 1.
    Frwso = 1.
    Fawso = 1.
    Frc = 1.
  endif
  return
end subroutine ompadjust
! Copyright A.J. Koning 2021
