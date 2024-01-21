function fmin(extern, k, psh, pd, pa, pe, ncall, eps)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Searches the minimal value fmin of the function extern.
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
  integer, parameter :: mmax=15             ! maximum number of points
  integer, parameter :: nmaxp1=16           ! maximum number of points
  real(sgl)          :: alpha               ! amplitude of the density-dependent term of the DDM3Y interaction
  real(sgl)          :: betap               ! Brosa parameter
  real(sgl)          :: eps                 ! help variable
  real(sgl)          :: extern              ! external function
  real(sgl)          :: f(nmaxp1)           ! E-Ef
  real(sgl)          :: fex                 ! help variable
  real(sgl)          :: fmax                ! maximum function value
  real(sgl)          :: fmax2               ! maximum function value
  real(sgl)          :: fmin                ! minimum function value
  real(sgl)          :: fn                  ! help variable
  real(sgl)          :: fr                  ! help variable
  real(sgl)          :: gam                 ! Brosa parameter
  real(sgl)          :: pa(k)               ! Brosa parameter
  real(sgl)          :: pc(mmax)            ! Brosa parameter
  real(sgl)          :: pd(k)               ! Brosa parameter
  real(sgl)          :: pe(k)               ! Brosa parameter
  real(sgl)          :: pex(mmax)           ! Brosa parameter
  real(sgl)          :: pmax                ! Brosa parameter
  real(sgl)          :: pr(mmax)            ! Brosa parameter
  real(sgl)          :: psh(k)              ! Brosa parameter
  real(sgl)          :: ssh(mmax, nmaxp1)   ! Brosa parameter
  real(sgl)          :: std                 ! Brosa parameter
  real(sgl)          :: sum                 ! help variable
  real(sgl)          :: sum2                ! help variable
!
! **********************************************************************
!
! extern: external function
! fmin  : minimum function value
! fmax  : maximum function value
! fmax2 : maximum function value
!
  alpha = 1.
  betap = 0.5
  gam = 2.
  if (k <= mmax) goto 4
  write ( * , 5) k, mmax
 5    format( / 1x, "from fmin: introduced number of parameters", i4, "    larger than allowed number", i3 / )
  stop
 4    icall = 0
  do i = 1, k
    if (psh(i) < pa(i)) psh(i) = pa(i)
    if (psh(i) > pe(i)) psh(i) = pe(i)
  enddo
  fn = float(k)
  n1 = k + 1

! construct initial simplex.
!
  do  jk = 1, n1
    do  i = 1, k
      ssh(i, jk) = psh(i)
    enddo
  enddo
  do  i = 1, k
    ssh(i, i + 1) = ssh(i, i + 1) + pd(i)
  enddo
  do jk = 1, n1
    if (jk == 1) goto 14
    do i = 1, k
      psh(i) = ssh(i, jk)
    enddo
 14     icall = icall + 1
    if (icall > ncall) goto 80
    f(jk) = extern(psh, k)
  enddo
!
! super loop.
! find largest, second largest, and smallest function value.
!
 17   fmax = f(1)
  maxf = 1
  do  i = 2, n1
    if (f(i) < fmax)  cycle
    maxf = i
    fmax = f(i)
  enddo
  fmin = f(1)
  minf = 1
  do  i = 2, n1
    if (f(i) > fmin)  cycle
    minf = i
    fmin = f(i)
  enddo
  fmax2 = f(minf)
  do  i = 1, n1
    if (i == maxf)  cycle
    if (f(i) < fmax2)  cycle
    fmax2 = f(i)
  enddo
!
! determine centroid.
!
  do  i = 1, k
    pc(i) = 0.
  enddo
  do  jk = 1, n1
    if (jk == maxf)  cycle
    do  i = 1, k
      pc(i) = pc(i) + ssh(i, jk)
    enddo
  enddo
  do  i = 1, k
    pc(i) = pc(i) / fn
  enddo
!
!   check accuracy.
!
  if (eps < 0.) goto 74
  std = 0.
  do i = 1, k
    pmax = abs(pc(i) - ssh(i, maxf))
  enddo
  if (pmax > std) std = pmax
  goto 80
 74   sum = 0.
  do jk = 1, n1
    sum = sum + f(jk)
  enddo
  sum = sum / float(n1)
  sum2 = 0.
  do jk = 1, n1
    sum2 = sum2 + (f(jk) - sum) **2
  enddo
  std = sqrt(sum2 / fn)
 80   if (icall < ncall .and. std > abs(eps)) goto 27
  do i = 1, k
    psh(i) = ssh(i, minf)
  enddo
  return
!
! reflect.
!
 27   do i = 1, k
    pr(i) = (1. + alpha) * pc(i) - alpha * ssh(i, maxf)
 enddo
  do i = 1, k
    if (pr(i) < pa(i) .or. pr(i) > pe(i)) goto 30
  enddo
  icall = icall + 1
  if (icall > ncall) goto 80
  fr = extern(pr, k)
  goto 31
 30   icall = icall + 1
  if (icall > ncall) goto 80
  fr = abs(fmax) * (1. + abs((pr(i) - pc(i)) / pd(i)))
 31   if (fr < fmin)  goto 40
  if (fr > fmax2)  goto 50
  goto 44
!
! expand.
!
 40   do i = 1, k
    pex(i) = gam * pr(i) + (1. - gam) * pc(i)
 enddo
  do i = 1, k
    if (pex(i) < pa(i) .or. pex(i) > pe(i)) goto 33
  enddo
  icall = icall + 1
  if (icall > ncall) goto 80
  fex = extern(pex, k)
  goto 34
 33   icall = icall + 1
  if (icall > ncall) goto 80
  fex = abs(fmax) * (1. + abs((pex(i) - pc(i)) / pd(i)))
 34   if (fex > fmin)  goto 44
 45   do i = 1, k
    ssh(i, maxf) = pex(i)
 enddo
  f(maxf) = fex
  goto 17
 44   do i = 1, k
    ssh(i, maxf) = pr(i)
 enddo
  f(maxf) = fr
  goto 17
 50   if(fr > fmax)  goto 54
  do  i = 1, k
    ssh(i, maxf) = pr(i)
  enddo
  f(maxf) = fr
!
! contract.
!
 54   do i = 1, k
    pex(i) = betap * ssh(i, maxf) + (1. - betap) * pc(i)
 enddo
  icall = icall + 1
  if (icall > ncall) goto 80
  fex = extern(pex, k)
  if (fex < fmax)  goto 45
!
! halve.
!
  do  jk = 1, n1
    if (jk == minf)  cycle
    do  i = 1, k
      ssh(i, jk) = .5 * (ssh(i, jk) + ssh(i, minf))
      psh(i) = ssh(i, jk)
    enddo
    icall = icall + 1
    if (icall > ncall) goto 80
    f(jk) = extern(psh, k)
  enddo
  goto 17
end function fmin
! Copyright A.J. Koning 2021

