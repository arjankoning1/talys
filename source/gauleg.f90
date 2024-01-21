subroutine gauleg(ngl, tgl, wgl)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculation of Gauss-Legendre arrays
!
! Author    : Stephane Hilaire
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
  integer   :: i        ! counter
  integer   :: it       ! counter for tritons
  integer   :: j        ! counter
  integer   :: k        ! designator for particle
  integer   :: ngl      ! number of points for Gauss-Laguerre integration
  integer   :: ns2      ! half of number of points for Gauss-Legendre integration
  real(sgl) :: oi       ! cosine term
  real(sgl) :: p        ! particle number
  real(sgl) :: p1       ! help variable
  real(sgl) :: p2       ! second parameter
  real(sgl) :: pi       ! pi
  real(sgl) :: tgl(ngl) ! points for Gauss-Laguerre integration
  real(sgl) :: ti       ! help variable
  real(sgl) :: wgl(ngl) ! weights for Gauss-Laguerre integration
!
! ****************** Calculation of Gauss-Legendre arrays **************
!
  pi = 3.14159265358979323
  ns2 = ngl / 2
  do i = 1, ns2
    ti = (4 * i - 1) * pi / (4 * ngl + 2.)
    oi = cos(ti + 1. / (8. * ngl * ngl * tan(ti)))
    do it = 1, 10
      p2 = 1.
      p1 = oi
      do k = 2, ngl
        p = ((2 * k - 1) * oi * p1 - (k - 1) * p2) / float(k)
        p2 = p1
        p1 = p
      enddo
      oi = oi - p * (1. - oi * oi) / (ngl * (p2 - oi * p))
    enddo
    tgl(i) = oi
    wgl(i) = (1. - oi * oi) / (ngl * p2 * ngl * p2)
  enddo
  do j = ns2 + 1, ngl
    tgl(j) = - tgl(j - ns2)
    wgl(j) = wgl(j - ns2)
  enddo
  return
end subroutine gauleg
! Copyright A.J. Koning 2021
