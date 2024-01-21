subroutine barsierk(iz, ia, il, bfis, egs, lbar0)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Fission barrier heights, rotating gs energy and lbar0
!
! Author    : Marieke Duijvestijn
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
! Variables for fission parameters
!   barcof      ! parameter values for barrier heights at l = 0
!   egscof      ! parameter values for rotating ground state energy
!   l20cof      ! parameter values for l20 belonging to 20% barrier height
!   l80cof      ! parameter values for l80 belonging to 80% barrier height
!   lmxcof      ! parameter values for lmax, l - value where barrier disappears
!
! *** Declaration of local data
!
  implicit none
  integer   :: ia        ! mass number from abundance table
  integer   :: il        ! angular momentum
  integer   :: iloop     ! loop counter
  integer   :: iz        ! charge number of residual nucleus
  integer   :: jloop     ! loop counter
  integer   :: kloop     ! loop counter
  real(sgl) :: a         ! fit variables
  real(sgl) :: a1        ! Myers-Swiatecki parameter
  real(sgl) :: a2        ! Myers-Swiatecki parameter
  real(sgl) :: amass     ! mass number of residual nucleus
  real(sgl) :: amax      ! maximal and minimal mass number defining range of the fit
  real(sgl) :: amax2     ! maximum A value
  real(sgl) :: amin      ! maximal and minimal mass number defining range of the fit
  real(sgl) :: amin2     ! minimum A value
  real(sgl) :: bfis      ! barrier height
  real(sgl) :: egs       ! rotating ground state energy
  real(sgl) :: l         ! multipolarity
  real(sgl) :: l20       ! l-value for which bfis is 20% (80% resp.) of bfis(l=0)
  real(sgl) :: l80       ! l-value for which bfis is 20% (80% resp.) of bfis(l=0)
  real(sgl) :: lbar0     ! l-value for which bfis becomes zero
  real(sgl) :: ll        ! angular momentum
  real(sgl) :: p1        ! help variable
  real(sgl) :: p2        ! second parameter
  real(sgl) :: plegendre ! function for calculation of Legendre polynomial
  real(sgl) :: q1        ! help variable
  real(sgl) :: q2        ! help variable
  real(sgl) :: q3        ! help variable
  real(sgl) :: r         ! ratio of neck contribution
  real(sgl) :: x         ! help variable
  real(sgl) :: y         ! coordinates of the point to test
  real(sgl) :: z         ! charge number
  real(sgl) :: zchar     ! charge number of residual nucleus
!
! ******************** Rotating Finite Range Model *********************
!
! plegendre: function for calculation of Legendre polynomial
!
! Barrier height for l=0
!
  bfis = 0.
  egs = 0.
  lbar0 = 0.
  if (iz < 19 .or. iz > 111) then
    write(*,'(/, 10x, "*  *  *  barfit called with  z  less than 19 or greater than 111.  bfis is set to 0.0  *  *  *")')
    return
  endif
  if (iz > 102 .and. il > 0) then
    write(*,'(/, 10x, "*  *  *  barfit called with  z  greater than 102 and  l  not equal to zero.  bfis is set to 0.0  *  *  *")')
    return
  endif
  zchar = real(iz)
  amass = real(ia)
  ll = real(il)
  amin = 1.2 * zchar + 0.01 * zchar * zchar
  amax = 5.8 * zchar - 0.024 * zchar * zchar
  if (amass < amin .or. amass > amax) then
    write(*,'(/, 10x, "*  *  *  barfit called with a = ", i3, ", outside the allowed values for z = ", i3, " *  *  *")') &
 &    int(amass), iz
    return
  endif
  a = amass / 400.
  z = zchar / 100.
  l = ll / 100.
  do iloop = 1, 7
    do jloop = 1, 7
      bfis = bfis + barcof(jloop, iloop) * plegendre(jloop - 1, z) * plegendre(iloop - 1, a)
    enddo
  enddo
  if (il < 1) return
!
! L-values corresponding to fission barrier height which is 20% (80%) of L=0 fission barrier
!
  l80 = 0.
  l20 = 0.
  amin2 = 1.4 * zchar + 0.009 * zchar * zchar
  amax2 = 20. + 3.0 * zchar
  if ((amass < amin2 - 5..or.amass > amax2 + 10.) .and. il > 0) then
    write(*,'(/, 10x, "*  *  *  barfit called with a = ", i3, ", outside the allowed values for z = ", i3, " *  *  *")') &
 &    int(amass), iz
    return
  endif
  do iloop = 1, 4
    do jloop = 1, 5
      l80 = l80 + l80cof(jloop, iloop) * plegendre(jloop - 1, z) * plegendre(iloop - 1, a)
      l20 = l20 + l20cof(jloop, iloop) * plegendre(jloop - 1, z) * plegendre(iloop - 1, a)
    enddo
  enddo
!
! L-value for which the fission barrier vanishes
!
  do iloop = 1, 4
    do jloop = 1, 6
      lbar0 = lbar0 + lmxcof(jloop, iloop) * plegendre(jloop - 1, z) * plegendre(iloop - 1, a)
    enddo
  enddo
!
! L-dependent fission barrier
!
  x = l20 / lbar0
  y = l80 / lbar0
  if (ll > l20) then
    p1 = ( - (20. * x **5) + 25. * x **4 - 4.) * (y - 1.) **2 * y * y
    p2 = ( - (20. * y **5) + 25. * y **4 - 1.) * (x - 1.) **2 * x * x
    q1 = 0.2 / ((y - x) * ((1. - x) * (1. - y) * x * y) **2)
    q2 = q1 * (p1 * y - p2 * x)
    q3 = - (q1 * (p1 * (2. * y + 1.) - p2 * (2. * x + 1.)))
    r = ll / lbar0
    a1 = 4. * r **5 - 5. * r **4 + 1.
    a2 = q2 * (2. * r + 1.)
    bfis = bfis * (a1 + (r - 1.) * (a2 + q3 * r) * r * r * (r - 1.))
  else
    q1 = 0.2 / (l20 **2 * l80 **2 * (l20 - l80))
    q2 = q1 * (4. * l80 **3 - l20 **3)
    q3 = - (q1 * (4. * l80 **2 - l20 **2))
    bfis = bfis * (1. + q2 * ll **2 + q3 * ll **3)
  endif
  if (bfis <= 0..or.ll > lbar0) bfis = 0.
!
! Rotating ground state energy
!
  if(ll > lbar0)return
  do iloop = 1, 4
    do jloop = 1, 6
      do kloop = 1, 5
        egs = egs + egscof(kloop, jloop, iloop) * plegendre(jloop - 1, z) * plegendre(iloop - 1, a) * plegendre(2 * kloop - 2, l)
      enddo
    enddo
  enddo
  if (egs < 0.) egs = 0.
!
! warning messages for attempted use outside validity boundaries
!
  return
end subroutine barsierk
