subroutine hrtwprepare(tjl, numtjl, st, tav, sv, v, w, numtr, Ninc, WFCfactor)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Preparation of HRTW width fluctuation correction
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
  integer   :: ni              ! counter
  integer   :: niter           ! number of iterations
  integer   :: Ninc            ! number of incident energies
  integer   :: WFCfactor       ! enhancement factor for WFC: 1: Original, 2: Ernebjerg and Herman
  integer   :: numtjl          ! number of transmission coefficients
  integer   :: numtr           ! number of transmission coefficients
  real(sgl) :: alpha           ! help variable
  real(sgl) :: beta            ! help variable
  real(sgl) :: gamma           ! help variable
  real(sgl) :: sv              ! variable for width fluctuation
  real(sgl) :: v(numtr)        ! real volume potential, radius, diffuseness
  real(sgl) :: w(numtr)        ! imaginary volume potential, radius, diffuseness
  real(dbl) :: f               ! E-Ef
  real(dbl) :: factor          ! multiplication factor
  real(dbl) :: st              ! denominator of compound nucleus formula
  real(dbl) :: st2             ! help variable
  real(sgl) :: fT              ! help variable
  real(sgl) :: gT              ! help variable
  real(sgl) :: denom           ! help variable
  real(sgl) :: Ta              ! help variable
  real(dbl) :: t               ! help variable
  real(dbl) :: va              ! help variable
  real(dbl) :: tav(numtr)      ! average transmission coefficients
  real(dbl) :: tjl(0:5, numtr) ! transmission coefficients
!
! ************************** HRTW calculation **************************
!
! We use the model of HRTW (1975), revised in 1980.
!
! Initialization and average transmission coefficients
!
  niter = 60
  do i = 1, numtr
    tav(i) = 0.
  enddo
  do i = 1, numtjl + 1
    if (tjl(0, i) > 0.) tav(i) = tjl(1, i) / tjl(0, i)
  enddo
!
! Calculation of sum over t**2
!
  st2 = 0.
  do i = Ninc + 1, numtjl + 1
    st2 = st2 + tjl(2, i)
  enddo
!
! Initialisation for HRTW calculation
!
  if (WFCfactor == 1) then
    factor = real(4. * st2 / (st * st + 3. * st2))
    do i = 1, numtjl + 1
      t = tav(i)
      if (t < st) then
        f = factor * (1. + t / st)
        w(i) = real(1. + 2. / (1. + t **f) + 87. * (((t - st2 / st) / st) **2) * (t / st) **5)
      else
        w(i) = 3.
      endif
    enddo
  else
    alpha = 0.139
    beta = 15.247
    gamma = 4.081
    do i = 1, numtjl + 1
      Ta = tav(i)
      fT = alpha/(1.-Ta**beta)
      gT = 1.+gamma*Ta*(1.-Ta)
      denom = 1.+fT*(st**gT)
      va = 2.-1./denom
      w(i) = 1.+2./va
    enddo
  endif
  do i = 1, numtjl + 1
    v(i) = real(tav(i) / (1. + tav(i) / st * (w(i) - 1.)))
  enddo
!
! Loop over iterations
!
  do ni = 1, niter
    sv = 0.
!
! Sum over v
!
    do i = Ninc + 1, numtjl + 1
      sv = sv + v(i) * tjl(0, i)
    enddo
!
! Determination of new v's
!
    do i = 1, numtjl + 1
      v(i) = real(tav(i) / (1. + v(i) / sv * (w(i) - 1.)))
    enddo
  enddo
  return
end subroutine hrtwprepare
! Copyright A.J. Koning 2021
