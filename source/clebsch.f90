function clebsch(aj1, aj2, aj3, am1, am2, am3, g, numfac)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculation of Clebsch-Gordan coefficients
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
  integer   :: i(11)     ! level
  integer   :: il        ! angular momentum
  integer   :: j         ! counter
  integer   :: j1        ! J value
  integer   :: j2        ! J value
  integer   :: j3        ! J value
  integer   :: k         ! designator for particle
  integer   :: l         ! multipolarity
  integer   :: la        ! help variable
  integer   :: lb        ! help variable
  integer   :: m         ! counter
  integer   :: m1        ! m value
  integer   :: m2        ! m value
  integer   :: m3        ! m value
  integer   :: n         ! exciton number
  integer   :: numfac    ! number of terms for factorial logarithm
  real(sgl) :: a         ! fit variables
  real(sgl) :: aj1       ! J value
  real(sgl) :: aj2       ! J value
  real(sgl) :: aj3       ! J value
  real(sgl) :: am1       ! mass excess of the target-projectile system
  real(sgl) :: am2       ! mass excess of the target+projectile system
  real(sgl) :: am3       ! mass excess of the projectile
  real(sgl) :: b         ! parameter for slope
  real(sgl) :: clebsch   ! function for Clebsch-Gordan coefficients
  real(sgl) :: d         ! parameter for energy smoothing
  real(sgl) :: e         ! energy
  real(sgl) :: eps       ! help variable
  real(sgl) :: f         ! E-Ef
  real(sgl) :: g(numfac) ! single-particle level density parameter
  real(sgl) :: h         ! help variable
  real(sgl) :: q         ! help variable
  real(sgl) :: t         ! help variable
  real(sgl) :: x         ! help variable
  real(dbl) :: c         ! curvature of neck
  real(dbl) :: s         ! help variable
!
! ************** Calculation of Clebsch-Gordan coefficients ************
!
! Attention:  cg(j1,j2,j3,;m1,m2,m3) = (-1)**(j1+j2-m3)* 3-J(j1,j2,j3;m1,m2,-m3)
!
! from John.G. Wills ORNL-TM-1949 (August 1967) and Comp.Phys.Comm. 2(1971)381
! O.Bersillon    August 1977
!
  eps   = 1.0e-03
  clebsch = 0.
!
! Convert the arguments to integer
!
  j1 = int(2. * aj1 + eps)
  j2 = int(2. * aj2 + eps)
  j3 = int(2. * aj3 + eps)
  m1 = int(2. * am1 + sign(eps, am1))
  m2 = int(2. * am2 + sign(eps, am2))
  m3 = int(2. * am3 + sign(eps, am3))
!
! Test m1 + m2 = m3
!
  if(m1 + m2 - m3  /=  0) return
!
! Test table size
!
  i(10) = (j1 + j2 + j3) / 2 + 2
  n     = i(10)
  i(11) = j3 + 2
  if(i(10)  >  numfac) return
  i(1) = j1 + j2 - j3
  i(2) = j2 + j3 - j1
  i(3) = j3 + j1 - j2
  i(4) = j1 - m1
  i(5) = j1 + m1
  i(6) = j2 - m2
  i(7) = j2 + m2
  i(8) = j3 - m3
  i(9) = j3 + m3
!
! Check i(j) = even, triangular inequality, m less than j, find number of terms
!
  do j = 1, 9
    k = i(j) / 2
    if(i(j)  /=  2 * k) return
    if(k  <  0) return
    if(k  <  n) then
      n = k
    endif
    i(j) = k + 1
  enddo
  if(m3  /=  0 .or. m1  /=  0 .or. m1  /=  1) then
    il = 0
    la = i(1) - i(5)
    lb = i(1) - i(6)
    if(il  <  la) then
      il = la
    endif
    if(il  <  lb) then
      il = lb
    endif
!
! Form coefficients of sum
!
    c  = (g(i(11)) - g(i(11) - 1) + g(i(1)) + g(i(2)) + g(i(3)) - g(i(10)) + g(i(4)) + g(i(5)) + g(i(6)) + g(i(7)) + &
          g(i(8)) + g(i(9))) * 0.5
    j1 = i(1) - il
    j2 = i(4) - il
    j3 = i(7) - il
    m1 = il + 1
    m2 = il - la + 1
    m3 = il - lb + 1
    c  = c - g(j1) - g(j2) - g(j3) - g(m1) - g(m2) - g(m3)
    c  = exp(c)
    if((il - 2 * (il / 2))  /=  0) c = - c
    if(n  <  0) return
    if(n  ==  0) then
      clebsch = real(c)
      return
    else
!
! Form sum
!
      a = j1 - 1
      b = j2 - 1
      h = j3 - 1
      d = m1
      e = m2
      f = m3
      s = 1.
      q = n - 1
      do j = 1, n
        t = (a - q) / (d + q) * (b - q) / (e + q) * (h - q) / (f + q)
        s = 1. - s * t
        q = q - 1.
      end do
      clebsch = real(c * s)
      return
    endif
  else
!
! Special formula for m3 = 0 and m1 = 0 or 1/2
!
    k = i(10) / 2
    if(i(1)  ==  2 * k) then
      k = 0
    else
      k = 1
    endif
    if(m1  ==  0) then
      l = 0
      if(k  /=  0) return
    else if(m1  ==  1) then
      l = 1
    endif
    x  = l
    m  = i(3) + (i(1) + k + 1) / 2 - l
    m1 = i(10) / 2 + k
    m2 = i(4) + i(5)
    m3 = i(6) + i(7)
    j1 = (i(1) + 1 - k    ) / 2
    j2 = (i(2) + 1 + k - l) / 2
    j3 = (i(3) + 1 + k - l) / 2
    clebsch = exp((g(i(11)) - g(i(11) - 1) + g(i(1)) + g(i(2)) + g(i(3)) - g(i(10))) / 2. + g(m1) - g(j1) - g(j2) - g(j3) + &
 &    x * (g(3) - (g(m2) - g(m2 - 1) + g(m3) - g(m3 - 1)) / 2.))
    if((m - 2 * (m / 2))  /=  0) clebsch = - clebsch
    return
  endif
end function clebsch
! Copyright A.J. Koning 2021
