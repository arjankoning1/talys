function plegendre(l, x)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculation of Legendre polynomial
!
! Author    : Marieke Duijvestijn
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               sgl      ! single precision kind
!
! *** Declaration of local data
!
  implicit none
  integer   :: i         ! level
  integer   :: l         ! multipolarity
  real(sgl) :: pl(0:200) ! legendre polynomial
  real(sgl) :: pl0       ! help variable
  real(sgl) :: pl1       ! help variable
  real(sgl) :: plegendre ! function for calculation of Legendre polynomial
  real(sgl) :: x         ! help variable
!
! ************************ Legendre polynomial *************************
!
  pl0 = 1.
  pl1 = x
  pl(0) = pl0
  pl(1) = pl1
  do i = 2, l
    pl(i) = (x * (2 * i - 1) * pl1 - (i - 1) * pl0) / i
    pl0 = pl1
    pl1 = pl(i)
  enddo
  plegendre = pl(l)
  return
end function plegendre
! Copyright A.J. Koning 2021
