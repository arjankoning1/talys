subroutine molprepare(tjl, numtjl, st, npmold, xmo, wmo, tav, vnu, product, numtr, Ninc, WFCfactor)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Preparation of Moldauer width fluctuation correction
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
  integer   :: i               ! counter
  integer   :: im              ! counter
  integer   :: npmold          ! variables for Gauss-Legendre integration
  integer   :: Ninc            ! number of incident energies
  integer   :: WFCfactor       ! enhancement factor for WFC: 1: Original, 2: Ernebjerg and Herman 3: Kawano
  integer   :: numtjl          ! number of transmission coefficients
  integer   :: numtr           ! number of transmission coefficients
  real(sgl) :: alpha           ! help variable
  real(sgl) :: beta            ! help variable
  real(sgl) :: gamma           ! help variable
  real(sgl) :: capt            ! help variable
  real(sgl) :: eps             ! help variable
  real(sgl) :: expo            ! help variable
  real(sgl) :: factor          ! multiplication factor
  real(sgl) :: fxmold          ! help variable
  real(sgl) :: fxmsqrt         ! help variable
  real(sgl) :: prod            ! help variable
  real(sgl) :: product(npmold) ! product used in final Moldauer calculation
  real(sgl) :: vnu(numtr)      ! number of degrees of freedom
  real(sgl) :: wmo(npmold)     ! help variable
  real(sgl) :: x               ! help variable
  real(sgl) :: x1              ! coordinates of intersection points inside the bin
  real(sgl) :: xmo(npmold)     ! variables for Gauss-Legendre integration
  real(sgl) :: yy              ! help variable
  real(sgl) :: fT              ! help variable
  real(sgl) :: gT              ! help variable
  real(sgl) :: denom           ! help variable
  real(sgl) :: Ta              ! help variable
  real(sgl) :: f               ! help variable
  real(sgl) :: alpha1          ! help variable
  real(sgl) :: beta1           ! help variable
  real(sgl) :: gamma1          ! help variable
  real(sgl) :: delta1          ! help variable
  real(dbl) :: st              ! denominator of compound nucleus formula
  real(dbl) :: tav(numtr)      ! average transmission coefficients
  real(dbl) :: tgamma          ! gamma transmission coefficient
  real(dbl) :: tjl(0:5, numtr) ! transmission coefficients
!
! **************** Preparation of Moldauer integral ********************
!
! Initialisation and average transmission coefficients
!
  tgamma = 0.
  do i = 1, numtr
    tav(i) = 0.
    vnu(i) = 1.
  enddo
  do i = 1, numtjl
    if (tjl(0, i) > 0.) tav(i) = tjl(1, i) / tjl(0, i)
  enddo
  tgamma = tjl(1, numtjl + 1)
!
! Calculation of number of degrees of freedom
!
  if (WFCfactor == 1) then
    do i = 1, numtjl
      vnu(i) = min(1.78 + real(((tav(i) **1.212) - 0.78) * exp( - 0.228 * st)), 2.)
    enddo
  else
    alpha = 0.177
    beta = 20.337
    gamma = 3.148
    do i = 1,numtjl
      Ta = tav(i)
      fT = alpha/(1.-Ta**beta)
      gT = 1.+gamma*Ta*(1.-Ta)
      denom = 1.+fT*(st**gT)
      vnu(i) = min(2.-1./denom,2.)
      vnu(i) = max(vnu(i), 1.)
    enddo
  endif
  if (WFCfactor.eq.3) then
    do i = 1,numtjl
      Ta = tav(i)
      alpha1 = 0.0287892*Ta+0.245856
      beta1 = 1.+2.5*Ta*(1.-Ta)*exp(-2.*st)
      gamma1 = Ta*Ta-(st-2.*Ta)**2
      if (gamma1 > 0. .and. st < 2*Ta) then
        delta1 = sqrt(gamma1)/Ta
      else
        delta1 = 1.
      endif
      f = alpha1*(st+Ta)/(1.-Ta)*beta1*delta1
      denom = 1.+f
      vnu(i) = 2.-1./denom
    enddo
  endif
  vnu(numtjl + 1) = 1.
!
! Loop over integration points
!
  do im = 1, npmold
    x = xmo(im)
!
! Special case for capture
!
    expo = real(tgamma * x / st)
    if (expo > 80.) then
      capt = 0.
    else
      capt = exp( - expo)
    endif
!
! Calculation of product over tjl
!
    factor = 0.
    do i = Ninc + 1, numtjl
      eps = real(2. * tav(i) * x / (st * vnu(i)))
      yy = real( - vnu(i) * 0.5 * tjl(0, i))
      if (eps <= 1.e-30) cycle
      if (eps < 1.e-5) then
        x1 = eps - 0.5 * eps * eps + (eps **3) / 3. - 0.25 * (eps **4) + 0.2 * (eps **5)
      else
        x1 = log(1. + eps)
      endif
      factor = factor + yy * x1
    enddo
    prod = exp(factor)
    fxmsqrt = wmo(im) * exp(0.5 * x)
    fxmold = fxmsqrt * prod * fxmsqrt
    product(im) = fxmold * capt
  enddo
  return
end subroutine molprepare
! Copyright A.J. Koning 2021
