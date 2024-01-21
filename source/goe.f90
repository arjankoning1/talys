subroutine goe(tjl, numtjl, na, nb, st, nweip, nweis, nweit, tav, x02, x102, &
    x202, x2i, x12i, x22i, fpst1, fpst2, s1, s2, s3, s4, s5, res, numtr, ielas)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : GOE triple integral width fluctuation correction
!
! Author    : Stephane Hilaire
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
  integer   :: i                         ! level
  integer   :: ielas                     ! designator for elastic channel
  integer   :: j                         ! counter
  integer   :: k                         ! designator for particle
  integer   :: na                        ! help variable
  integer   :: nb                        ! help variable
  integer   :: numtjl                    ! number of transmission coefficients
  integer   :: numtr                     ! number of transmission coefficients
  integer   :: nweip                     ! variable for GOE triple integral calculation
  integer   :: nweis                     ! variable for GOE triple integral calculation
  integer   :: nweit                     ! variable for GOE triple integral calculation
  real(sgl) :: dab                       ! help variable
  real(sgl) :: fpst1(nweip, nweis, nweit)! variable for final GOE calculation
  real(sgl) :: fpst2(nweip, nweis, nweit)! variable for final GOE calculation
  real(sgl) :: func1                     ! function for GOE
  real(sgl) :: res                       ! width fluctuation factor
  real(sgl) :: res02                     ! variable for final GOE calculation
  real(sgl) :: res2i                     ! variable for final GOE calculation
  real(sgl) :: s1                        ! variable for final GOE calculation
  real(sgl) :: s2                        ! variable for final GOE calculation
  real(sgl) :: s3                        ! variable for final GOE calculation
  real(sgl) :: s4                        ! variable for final GOE calculation
  real(sgl) :: s5                        ! variable for final GOE calculation
  real(sgl) :: x02(nweip, nweis, nweit)  ! variable for final GOE calculation
  real(sgl) :: x102(nweip, nweis, nweit) ! variable for final GOE calculation
  real(sgl) :: x12i(nweip, nweis, nweit) ! variable for final GOE calculation
  real(sgl) :: x202(nweip, nweis, nweit) ! variable for final GOE calculation
  real(sgl) :: x22i(nweip, nweis, nweit) ! variable for final GOE calculation
  real(sgl) :: x2i(nweip, nweis, nweit)  ! variable for final GOE calculation
  real(sgl) :: y1                        ! variable for final GOE calculation
  real(sgl) :: y2                        ! variable for final GOE calculation
  real(sgl) :: y3                        ! variable for final GOE calculation
  real(sgl) :: y4                        ! variable for final GOE calculation
  real(sgl) :: y5                        ! variable for final GOE calculation
  real(sgl) :: y6                        ! variable for final GOE calculation
  real(sgl) :: y7                        ! variable for final GOE calculation
  real(sgl) :: y8                        ! variable for final GOE calculation
  real(dbl) :: c1s1                      ! term for GOE
  real(dbl) :: c1s2                      ! term for GOE
  real(dbl) :: c1s3                      ! term for GOE
  real(dbl) :: c1s4                      ! term for GOE
  real(dbl) :: c2s1                      ! term for GOE
  real(dbl) :: c2s2                      ! term for GOE
  real(dbl) :: c2s3                      ! term for GOE
  real(dbl) :: c2s4                      ! term for GOE
  real(dbl) :: c2s5                      ! term for GOE
  real(dbl) :: res1                      ! help variable
  real(dbl) :: rhob                      ! help variable
  real(dbl) :: s12                       ! help variable
  real(dbl) :: s13                       ! help variable
  real(dbl) :: s14                       ! help variable
  real(dbl) :: s15                       ! help variable
  real(dbl) :: st                        ! denominator of compound nucleus formula
  real(dbl) :: ta                        ! transmission coefficient
  real(dbl) :: ta2                       ! transmission coefficient ** 2
  real(dbl) :: ta3                       ! transmission coefficient ** 3
  real(dbl) :: ta4                       ! transmission coefficient ** 4
  real(dbl) :: ta5                       ! transmission coefficient ** 5
  real(dbl) :: taptb                     ! term for GOE
  real(dbl) :: taptb2                    ! term for GOE
  real(dbl) :: taptb3                    ! term for GOE
  real(dbl) :: taptb4                    ! term for GOE
  real(dbl) :: taptb5                    ! term for GOE
  real(dbl) :: tav(numtr)                ! average transmission coefficients
  real(dbl) :: tb                        ! transmission coefficient
  real(dbl) :: tb2                       ! term for GOE
  real(dbl) :: tb3                       ! term for GOE
  real(dbl) :: tb4                       ! term for GOE
  real(dbl) :: tb5                       ! term for GOE
  real(dbl) :: tb6                       ! term for GOE
  real(dbl) :: tjl(0:5, numtr)           ! transmission coefficients
  real(dbl) :: tt                        ! term with transmission coefficients
!
! *************************** GOE calculation **************************
!
! Initialization
!
  res1 = 0.
  dab = 0.
  if (ielas == 1) dab = 1.
!
! ****** Numerical calculation of triple integral for few channels *****
!
  if (st < 20.) then
    ta = tav(na)
    tb = tav(nb)
    rhob = tjl(1, nb)
    if (nb == numtjl + 1) tb = dble(0.)
!
! Loop over p
!
    do i = 1, nweip
!
! Loop over s
!
      do j = 1, nweis
!
! Loop over t
!
        do k = 1, nweit
!
! func1         : function for GOE
!
          y1 = fpst1(i, j, k)
          y2 = fpst2(i, j, k)
          y3 = x102(i, j, k)
          y4 = x202(i, j, k)
          y5 = x02(i, j, k)
          y6 = x12i(i, j, k)
          y7 = x22i(i, j, k)
          y8 = x2i(i, j, k)
          res02 = y1 * func1(y3, y4, y5, ta, tb, dab)
          res2i = y2 * func1(y6, y7, y8, ta, tb, dab)
          res1 = res1 + real(ta * (res02 + res2i) * rhob)
        enddo
      enddo
    enddo
  else
!
! ***** Numerical calculation of triple integral for many channels *****
!
! Compute sum of tc,tc**2,tc**3 ....
!
!
    s12 = s1 * s1
    s13 = s12 * s1
    s14 = s13 * s1
    s15 = s14 * s1
    ta = tjl(1, na)
    ta2 = ta * ta
    ta3 = ta2 * ta
    ta4 = ta3 * ta
    ta5 = ta4 * ta
    tb = tjl(1, nb)
    tb2 = tjl(2, nb)
    tb3 = tjl(3, nb)
    tb4 = tjl(4, nb)
    tb5 = tjl(5, nb)
    tb6 = tb * tjl(5, nb) / max(tjl(0, nb), 1.d0)
    taptb = ta * tb + tb2
    taptb2 = ta2 * tb + ta * tb2 + tb3
    taptb3 = ta3 * tb + ta2 * tb2 + ta * tb3 + tb4
    taptb4 = ta4 * tb + ta3 * tb2 + ta2 * tb3 + ta * tb4 + tb5
    taptb5 = ta5 * tb + ta4 * tb2 + ta3 * tb3 + ta2 * tb4 + ta * tb5 + tb6
    c1s1 = ( - 2. - 4. * ta) / s1
    c1s2 = (6. + 3. * s2 + 12. * ta + 34. * ta2) / s12
    c1s3 = ( - 32. - 12. * s2 - 16. * s3 - 64. * ta - 40. * ta * s2 - 136. * ta2 - 304. * ta3) / s13
    c1s4 = (240. + 80. * s2 + 25 * s2 * s2 + 80 * s3 + 100. * s4 + 480. * ta + 200. * ta * s2 + &
      240. * ta * s3 + 864. * ta2 + 524. * ta2 * s2 + 1520. * ta3 + 3508. * ta4) / s14
    c2s1 = - taptb / s1
    c2s2 = (s2 * tb - 2. * taptb + 4. * taptb2) / s12
    c2s3 = (s2 * (2. * tb - 5. * taptb) - 4. * s3 * tb + 4. * taptb + 8. * taptb2 - 20. * taptb3) / s13
    c2s4 = (s2 * ( - 4. * tb - 16. * taptb + 36. * taptb2) + s3 * ( - 8. * tb + 24. * taptb) + &
      20. * s4 * tb - 12. * taptb + 5. * s2 * s2 * tb - 20. * taptb2 - 68. * taptb3 + 148. * taptb4) / s14
    c2s5 = (s2 * (12. * tb + 48. * taptb + 156. * taptb2 - 336. * taptb3) - 60. * s2 * s3 * &
      tb - 148. * s5 * tb + s3 * (20. * tb + 108. * taptb - 228. * taptb2) + s4 * (68. * tb - &
      168. * taptb) + s2 * s2 * (16. * tb - 41. * taptb) + 64. * taptb + 88. * taptb2 + &
      204. * taptb3 + 608. * taptb4 - 1348. * taptb5) / s15
    res1 = dab * 2. * ta2 * (1. - ta) / s12 * (1. + c1s1 + c1s2 + c1s3 + c1s4) + &
      (1. + dab) * ta / s1 * (tb + c2s1 + c2s2 + c2s3 + c2s4 + c2s5)
  endif
  tt = tjl(1, na) * tjl(1, nb) / tjl(0, na)
  if (tt /= 0.) then
    res = real(res1 * st / tt)
  else
    res = 0.
  endif
  return
end subroutine goe
! Copyright A.J. Koning 2021
