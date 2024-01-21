subroutine goeprepare(tjl, numtjl, st, nweip, nweis, nweit, p, s, &
    t, wp, ws, wt, tav, x02, x102, x202, x2i, x12i, x22i, fpst1, fpst2, s1, s2, s3, s4, s5, numtr, Ninc)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Preparation of GOE triple integral width fluctuation
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
  integer   :: i                          ! counter
  integer   :: j                          ! counter
  integer   :: k                          ! designator for particle
  integer   :: Ninc                       ! number of incident energies
  integer   :: numtjl                     ! number of transmission coefficients
  integer   :: numtr                      ! number of transmission coefficients
  integer   :: nweip                      ! variable for GOE triple integral calculation
  integer   :: nweis                      ! variable for GOE triple integral calculation
  integer   :: nweit                      ! variable for GOE triple integral calculation
  real(sgl) :: d                          ! parameter for energy smoothing
  real(sgl) :: ds2                        ! help variable
  real(sgl) :: e                          ! energy
  real(sgl) :: ex02                       ! help variable
  real(sgl) :: ex2i                       ! help variable
  real(sgl) :: fpst1(nweip, nweis, nweit) ! variable for final GOE calculation
  real(sgl) :: fpst2(nweip, nweis, nweit) ! variable for final GOE calculation
  real(sgl) :: p(nweip)                   ! particle number
  real(sgl) :: p1                         ! help variable
  real(sgl) :: p2                         ! second parameter
  real(sgl) :: pm2ps21                    ! term for GOE
  real(sgl) :: pm2s22                     ! term for GOE
  real(sgl) :: pms222                     ! term for GOE
  real(sgl) :: pp                         ! Duflo-Zuker parameter
  real(sgl) :: ps21                       ! term for GOE
  real(sgl) :: ps2t21                     ! term for GOE
  real(sgl) :: pums212                    ! term for GOE
  real(sgl) :: s(nweis)                   ! help variable
  real(sgl) :: s1                         ! variable for final GOE calculation
  real(sgl) :: s2                         ! variable for final GOE calculation
  real(sgl) :: s21                        ! variable for final GOE calculation
  real(sgl) :: s22                        ! variable for final GOE calculation
  real(sgl) :: s2t21                      ! variable for final GOE calculation
  real(sgl) :: s2t22                      ! variable for final GOE calculation
  real(sgl) :: s3                         ! variable for final GOE calculation
  real(sgl) :: s4                         ! variable for final GOE calculation
  real(sgl) :: s5                         ! variable for final GOE calculation
  real(sgl) :: t(nweit)                   ! help variable
  real(sgl) :: t1                         ! help variable
  real(sgl) :: t2                         ! help variable
  real(sgl) :: t21                        ! variable for final GOE calculation
  real(sgl) :: t22                        ! variable for final GOE calculation
  real(sgl) :: um2s21                     ! term for GOE
  real(sgl) :: ume                        ! term for GOE
  real(sgl) :: umes2                      ! term for GOE
  real(sgl) :: umps21                     ! term for GOE
  real(sgl) :: ums21                      ! term for GOE
  real(sgl) :: ums22                      ! term for GOE
  real(sgl) :: umt21                      ! term for GOE
  real(sgl) :: umt22                      ! term for GOE
  real(sgl) :: upes2                      ! term for GOE
  real(sgl) :: uppm2ps2                   ! term for GOE
  real(sgl) :: uppm2s2                    ! term for GOE
  real(sgl) :: wp(nweip)                  ! optical model parameter
  real(sgl) :: wp1                        ! term for GOE
  real(sgl) :: wpspp2                     ! term for GOE
  real(sgl) :: wpws1                      ! term for GOE
  real(sgl) :: wpws2                      ! term for GOE
  real(sgl) :: ws(nweis)                  ! surface component of the imaginary potential
  real(sgl) :: ws1                        ! term for GOE
  real(sgl) :: ws2                        ! term for GOE
  real(sgl) :: wt(nweit)                  ! (inverse) weight
  real(sgl) :: wt2                        ! term for GOE
  real(sgl) :: x02(nweip, nweis, nweit)   ! variable for final GOE calculation
  real(sgl) :: x102(nweip, nweis, nweit)  ! variable for final GOE calculation
  real(sgl) :: x12i(nweip, nweis, nweit)  ! variable for final GOE calculation
  real(sgl) :: x1rat                      ! variable for final GOE calculation
  real(sgl) :: x202(nweip, nweis, nweit)  ! variable for final GOE calculation
  real(sgl) :: x22i(nweip, nweis, nweit)  ! variable for final GOE calculation
  real(sgl) :: x2i(nweip, nweis, nweit)   ! variable for final GOE calculation
  real(sgl) :: x2rat                      ! variable for final GOE calculation
  real(dbl) :: prodm                      ! product function for GOE
  real(dbl) :: prodp                      ! product function for GOE
  real(dbl) :: st                         ! denominator of compound nucleus formula
  real(dbl) :: tav(numtr)                 ! average transmission coefficients
  real(dbl) :: tgamma                     ! gamma transmission coefficient
  real(dbl) :: tjl(0:5, numtr)            ! transmission coefficients
!
! *************************** GOE calculation **************************
!
! Initialisation and average transmission coefficients
!
  tgamma = 0.
  do i = 1, numtr
    tav(i) = 0.
  enddo
  do i = 1, numtjl
    if (tjl(0, i) > 0.) tav(i) = tjl(1, i) / tjl(0, i)
  enddo
  tgamma = tjl(1, numtjl + 1)
!
! ****** Numerical calculation of triple integral for few channels *****
!
  if (st < 20.) then
!
! Loop over p
!
    do i = 1, nweip
      p1 = p(i) + 1.
      pp = 0.5 * p1
      p2 = (3. - pp) / pp
      d = sqrt(0.5 * p2)
      ds2 = 0.5 * d
      wp1 = 2. * wp(i)
      wpspp2 = 3. * wp(i) / (pp * pp)
!
! Special case for capture
!
      if (tgamma == 0.) then
        ex02 = 1.
        ex2i = 1.
      else
        ex02 = real(exp( - p1 * tgamma * 0.5))
        ex2i = real(exp( - p2 * tgamma * 0.5))
      endif
!
! Loop over s
!
      do j = 1, nweis
        s1 = 0.25 * sqrt(2.) * (s(j) + 1.)
        s2 = ds2 * s(j) + ds2
        s21 = real(s1 * s1)
        s22 = real(s2 * s2)
        ums21 = 1. - s21
        ums22 = 1. - s22
        ps21 = p1 * s21
        umps21 = 1. - ps21
        pm2ps21 = p1 - ps21 - ps21
        um2s21 = ums21 - s21
        pums212 = p1 * ums21 * ums21
        uppm2ps2 = pm2ps21 + 1.
        pm2s22 = p2 - s22 - s22
        uppm2s2 = 1. + pm2s22
        pms222 = (p2 - s22) * (p2 - s22)
        ws2 = d * ws(j)
        ws1 = 0.5 * sqrt(2.) * ws(j)
        if (s2 > 1.) then
          e = sqrt(1. - 1. / s22)
        else
          e = 0.
        endif
        ume = 1. - e
        umes2 = 0.5 * ume
        upes2 = 0.5 + 0.5 * e
        wpws1 = wp1 * ws1
        wpws2 = wpspp2 * ws2
!
! Loop over t
!
! prodm,prodp   : product function for GOE
!
        do k = 1, nweit
          t1 = 0.5 * t(k) + 0.5
          t2 = umes2 * t(k) + upes2
          wt2 = ume * wt(k)
          t21 = t1 * t1
          t22 = t2 * t2
          s2t22 = s22 * t22
          s2t21 = s21 * t21
          umt21 = 1. - t21
          umt22 = 1. - t22
          ps2t21 = ps21 * t21
          x02(i, j, k) = ps21 - ps2t21
          x102(i, j, k) = ps2t21
          x202(i, j, k) = pm2ps21 + ps2t21
          x2i(i, j, k) = s22 - s2t22
          x12i(i, j, k) = s2t22
          x22i(i, j, k) = pm2s22 + s2t22
          x1rat = real(prodm(x02(i, j, k), tjl, numtjl, numtr, Ninc) / &
            prodp(x102(i, j, k), tjl, numtjl, numtr, Ninc) / prodp(x202(i, j, k), tjl, numtjl, numtr, Ninc))
          fpst1(i, j, k) = wpws1 * wt(k) * (ps2t21 + umps21) * um2s21 * umt21 * ex02 * x1rat **10 / pums212 / &
            sqrt((um2s21 + s2t21) * (1. + ps2t21) * (ps2t21 + uppm2ps2))
          x2rat = real(prodm(x2i(i, j, k), tjl, numtjl, numtr, Ninc) / &
            prodp(x12i(i, j, k), tjl, numtjl, numtr, Ninc) / prodp(x22i(i, j, k), tjl, numtjl, numtr, Ninc))
          fpst2(i, j, k) = wpws2 * wt2 * umt22 * (ums22 + s2t22) * pm2s22 * ex2i * x2rat **10 / pms222 / &
            sqrt((1. + s2t22) * (pm2s22 + s2t22) * (uppm2s2 + s2t22))
        enddo
      enddo
    enddo
  else
!
! ***** Numerical calculation of triple integral for many channels *****
!
    s1 = 0.
    s2 = 0.
    s3 = 0.
    s4 = 0.
    s5 = 0.
    do i = Ninc + 1, numtjl + 1
      s1 = s1 + real(tjl(1, i))
      s2 = s2 + real(tjl(2, i))
      s3 = s3 + real(tjl(3, i))
      s4 = s4 + real(tjl(4, i))
      s5 = s5 + real(tjl(5, i))
    enddo
  endif
  return
end subroutine goeprepare
! Copyright A.J. Koning 2021
