function racah(a, b, c, d, e, f, g, numfac)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculation of Racah coefficients
!
! Author    : Stephane Hilaire (adapted from O. Bersillon)
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               dbl  , & ! double precision kind
               sgl      ! single precision kind
!
! *** Declaration of local data
!
  implicit none
  integer   :: i(16)     ! level
  integer   :: il        ! angular momentum
  integer   :: j         ! counter
  integer   :: j1        ! J value
  integer   :: j2        ! J value
  integer   :: j3        ! J value
  integer   :: j4        ! J value
  integer   :: j5        ! J value
  integer   :: j6        ! J value
  integer   :: j7        ! J value
  integer   :: ja        ! J value
  integer   :: jb        ! J value
  integer   :: jc        ! J value
  integer   :: jd        ! J value
  integer   :: je        ! J value
  integer   :: jf        ! J value
  integer   :: k         ! designator for particle
  integer   :: n         ! exciton number
  integer   :: numfac    ! number of terms for factorial logarithm
  real(sgl) :: a         ! fit variables
  real(sgl) :: b         ! parameter for slope
  real(sgl) :: c         ! curvature of neck
  real(sgl) :: d         ! parameter for energy smoothing
  real(sgl) :: e         ! energy
  real(sgl) :: eps       ! help variable
  real(sgl) :: f         ! E-Ef
  real(sgl) :: g(numfac) ! single-particle level density parameter
  real(sgl) :: o         ! J value
  real(sgl) :: p         ! particle number
  real(sgl) :: q         ! help variable
  real(sgl) :: r         ! ratio of neck contribution
  real(sgl) :: racah     ! function for racah coefficients
  real(sgl) :: t         ! help variable
  real(sgl) :: v         ! real volume potential, radius, diffuseness
  real(sgl) :: w         ! imaginary volume potential, radius, diffuseness
  real(sgl) :: x         ! help variable
  real(sgl) :: y         ! coordinates of the point to test
  real(sgl) :: z         ! charge number
  real(dbl) :: h         ! help variable
  real(dbl) :: s         ! help variable
!
! *********** Calculation of Racah coefficients   w(a,b,c,d;e,f) *******
!
! From John.G.Wills  ORNL-TM-1949 (August 1967) and Comp.Phys.Comm. 2(1971)381
! O.Bersillon    August 1977
!
  eps  = 1.0e-03
  racah = 0.
!
!  Convert arguments to integer and make useful combinations
!
  ja  = int(2. * a + eps)
  jb  = int(2. * b + eps)
  jc  = int(2. * c + eps)
  jd  = int(2. * d + eps)
  je  = int(2. * e + eps)
  jf  = int(2. * f + eps)
  i(1)  = ja + jb - je
  i(2)  = jb + je - ja
  i(3)  = je + ja - jb
  i(4)  = jc + jd - je
  i(5)  = jd + je - jc
  i(6)  = je + jc - jd
  i(7)  = ja + jc - jf
  i(8)  = jc + jf - ja
  i(9)  = jf + ja - jc
  i(10) = jb + jd - jf
  i(11) = jd + jf - jb
  i(12) = jf + jb - jd
  i(13) = ja + jb + je
  i(14) = jc + jd + je
  i(15) = ja + jc + jf
  i(16) = jb + jd + jf
!
! Check triangular inequalities, find no. of terms in sum, divide I's by 2
!
  n = i(16)
  do j = 1, 12
    k = i(j) / 2
    if(i(j)  /=  2 * k) return
    if(k  <  0) return
    if(k  <  n) then
      n = k
    endif
    i(j) = k + 1
  enddo
!
! Find minimum value of summation index
!
  il = 0
  do j = 13, 16
    i(j) = i(j) / 2
    if(il  <  i(j)) then
      il = i(j)
    endif
  enddo
  j1 = il  - i(13) + 1
  j2 = il  - i(14) + 1
  j3 = il  - i(15) + 1
  j4 = il  - i(16) + 1
  j5 = i(13) + i(4)  - il
  j6 = i(15) + i(5)  - il
  j7 = i(16) + i(6)  - il
  h  = - exp(0.5 * (g(i(1)) + g(i(2)) + g(i(3)) - g(i(13) + 2) + g(i(4)) + g(i(5)) &
    + g(i(6)) - g(i(14) + 2) + g(i(7)) + g(i(8)) + g(i(9)) - g(i(15) + 2) + g(i(10)) + &
    g(i(11)) + g(i(12)) - g(i(16) + 2)) + g(il + 2) - g(j1) - g(j2) - g(j3) - g(j4) - g(j5) - g(j6) - g(j7))
  if((j5 - 2 * (j5 / 2))  /=  0) h = - h
  if(n  <  0) return
  if(n  ==  0) then
    racah = real(h)
    return
  else
    s = 1.
    q = n  - 1
    p = il + 2
    r = j1
    o = j2
    v = j3
    w = j4
    x = j5 - 1
    y = j6 - 1
    z = j7 - 1
    do j = 1, n
      t = (p + q) / (r + q) * (x - q) / (o + q) * (y - q) / (v + q) * (z - q) / (w + q)
      s = 1. - s * t
      q = q   - 1.
    enddo
    racah = real(h * s)
  endif
  return
end function racah
! Copyright A.J. Koning 2021
