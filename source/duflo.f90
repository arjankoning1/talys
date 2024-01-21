subroutine duflo(nn, nz, exc)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Purpose of this module
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
!
!
  call mass10(nn, nz, e)
  exc = nz * 7.28903 + nn * 8.07138 - e
  return
end subroutine duflo
!
subroutine mass10(nx, nz, E)              ! Duflo - Zuker fevrier 1996
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Purpose of this module
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               sgl                        ! single precision kind
!
! *** Declaration of local data
!
  real(sgl) :: b(10)                      !
  real(sgl) :: dei(2)                     !
  real(sgl) :: dx(2)                      !
  real(sgl) :: dyda(10)                   !
  real(sgl) :: oei(2)                     !
  real(sgl) :: onp(0:20, 2, 2)            !
  real(sgl) :: op(2)                      !
  real(sgl) :: os(2)                      !
  real(sgl) :: pp(2)                      !
  real(sgl) :: qx(2)                      !
  real(sgl) :: y(2)                       !
  integer   :: n2(2)                      !
  integer   :: nn(2)                      !
  integer   :: noc(28, 2)                 !
  b =  (/0.7043, 17.7418, 16.2562, 37.5562, 53.9017, 0.4711, 2.1307, 0.0210, 40.5356, 6.0632 /)
!*********
!
  nn(1) = nx
  nn(2) = nz
  a = nx + nz
  t = abs(nx - nz)
  r = a **(1. / 3.)
  rc = r * (1. - .25 * (t / a) **2) !      Charge radius
  ra = (rc * rc) / r
!--------
  z2 = nz * (nz - 1)
  dyda(1) = ( - z2 + .76 * z2 **(2. / 3.)) / rc ! Coulomb energy
!********                          ! beginning of main loop
  do ndef = 1, 2                                !      ndef = 1  spherical
  ju = 0                                        !      ndef = 2  deformed
  y(ndef) = 0.
  if(ndef == 2) ju = 4 !      nucleons associated to deform.
  do kk = 2, 10
    dyda(kk) = 0.
  enddo
!--------                          ! beginning of loop over N and Z
  do j = 1, 2
    do l = 1, 28
      noc(l, j) = 0
    enddo
    do l = 1, 2
      do k = 0, 20
        onp(k, l, j) = 0.
      enddo
    enddo
    n2(j) = 2 * (nn(j) / 2) !      (for pairing calculation)
    ncum = 0
    i = 0
!--------
  20    i = i + 1 !    sub - shells (ssh) j and r filling
    i2 = (i / 2) * 2
    if(i2 /= i)then
      id = i + 1 !             for ssh j
    else
      id = i * (i - 2) / 4 !             for ssc r
    endif
    ncum = ncum + id
    if(ncum < nn(j))then
      noc(i, j) = id !     nb of nucleons in each ssh
      goto 20
    endif
!--------
    imax = i + 1     !     imax = last subshell nb
    ip = (i - 1) / 2 !     HO number (p)
    ipm = i / 2
    pp(j) = ip
    moc = nn(j) - ncum + id
    noc(i, j) = moc - ju !     nb of nucleons in last ssh
    noc(i + 1, j) = ju
    if(i2 /= i)then                    !     ssh j
      oei(j) = moc + ip * (ip - 1)     !     nb of nucleons in last EI shell
      dei(j) = ip * (ip + 1) + 2       !       size of the EI shell
    else                               !     ssh r
      oei(j) = moc - ju                !     nb of nucleons in last EI shell
      dei(j) = (ip + 1) * (ip + 2) + 2 !       size of the EI shell
    endif
    qx(j) = oei(j) * (dei(j) - oei(j) - ju) / dei(j) ! n * (D - n) / D        S3(j)
    dx(j) = qx(j) * (2 * oei(j) - dei(j))            ! n * (D - n) * (2n - D) / D  Q
    if(ndef == 2)qx(j) = qx(j) / sqrt(dei(j))        ! scaling for deformed
!--------
    do i = 1, imax                                   ! Amplitudes
      ip = (i - 1) / 2
      fact = sqrt((ip + 1.) * (ip + 2.))
      onp(ip, 1, j) = onp(ip, 1, j) + noc(i, j) / fact !    for FM term
      vm = - 1.
      if((2 * (i / 2)) /= i)vm = .5 * ip !    for spin - orbit term
      onp(ip, 2, j) = onp(ip, 2, j) + noc(i, j) * vm
    enddo
!--------
    op(j) = 0.
    os(j) = 0.
    do ip = 0, ipm !       FM and SO terms
      pi = ip
      den = ((pi + 1) * (pi + 2)) **(3. / 2.)
      op(j) = op(j) + onp(ip, 1, j) ! FM
      os(j) = os(j) + onp(ip, 2, j) * (1. + onp(ip, 1, j)) * (pi * pi / den) &
                 + onp(ip, 2, j) * (1. - onp(ip, 1, j)) * ((4 * pi - 5) / den) ! SO
    enddo
    op(j) = op(j) * op(j)
  enddo
!--------                          ! end of loop over  N and Z
  dyda(2) = op(1) + op(2)           !   Master term (FM): volume
  dyda(3) = - dyda(2) / ra          !                     surface
  dyda(2) = dyda(2) + os(1) + os(2) !   FM + SO
  dyda(4) = - t * (t + 2) / (r * r) !   isospin term : volume
  dyda(5) = - dyda(4) / ra          !                : surface
  if(ndef == 1)then                 ! sph.
    dyda(6) = dx(1) + dx(2)         !   S3  volume
    dyda(7) = - dyda(6) / ra        !       surface
    px = sqrt(pp(1)) + sqrt(pp(2))
    dyda(8) = qx(1) * qx(2) * (2 **px) !   QQ sph.
  else                                 ! def.
    dyda(9) = qx(1) * qx(2)            !   QQ deform.
  endif
  dyda(5) = t * (1 - t) / (a * ra **3) + dyda(5) !   "Wigner term"
!--------                                 !   PAIRING
  if(n2(1) /= nn(1) .and. n2(2) /= nn(2))dyda(10) = t / a
  if(nx > nz)then
    if(n2(1) == nn(1) .and. n2(2) /= nn(2))dyda(10) = 1 - t / a
    if(n2(1) /= nn(1) .and. n2(2) == nn(2))dyda(10) = 1
  else
    if(n2(1) == nn(1) .and. n2(2) /= nn(2))dyda(10) = 1
    if(n2(1) /= nn(1) .and. n2(2) == nn(2))dyda(10) = 1 - t / a
  endif
  if(n2(2) == nn(2) .and. n2(1) == nn(1))dyda(10) = 2 - t / a
!--------
  do mss = 2, 10
    dyda(mss) = dyda(mss) / ra
  enddo
  do mss = 1, 10
    y(ndef) = y(ndef) + dyda(mss) * b(mss)
  enddo
!--------                            ! end of main loop
  enddo
  de = y(2) - y(1)
  E = y(2)                         ! Binding Energy for def. nuclides
  if(de <= 0..or.nz <= 50)E = y(1) !                spherical nuclides
  return
end subroutine mass10
! Copyright A.J. Koning 2021
