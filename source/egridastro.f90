subroutine egridastro
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculate default incident energy grid for astrophysical
!
! Author    : Stephane Goriely
!
! 2021-12-30: Original code
! 2022-04-07: Changed definition for maximum temperature
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl         ! single precision kind
!   dbl         ! double precision kind
! All global variables
!   numenin     ! number of incident energies
! Variables for astrophysics
!   astroT9     ! temperature, in 10^9 K, for Maxwellian average
!   nTmax       ! effective number of temperatures for Maxwellian
! Variables for input energies
!   eninc       ! incident energy in MeV
!   enincmax    ! maximum incident energy
!   enincmin    ! minimum incident energy
!   Ninc        ! number of incident energies
! Variables for main input
!   k0          ! index of incident particle
!   Ztarget     ! charge number of target nucleus
! Variables for nuclides
!   T9          ! Temperature grid in 10 **9 K
! Constants
!   kT          ! energy kT expressed in MeV corresponding to a temperature T9 = 1
!   parN        ! neutron number of particle
!   parZ        ! charge number of particle
! Variables for masses
!   redumass    ! reduced mass
!   S           ! separation energy
!   specmass    ! specific mass for residual nucleus
!
! *** Declaration of local data
!
  implicit none
  integer   :: neg     ! number of energies
  integer   :: neg1    ! number of energies
  integer   :: neg2    ! number of energies
  integer   :: nen     ! energy counter
  integer   :: Nix     ! neutron number index for residual nucleus
  integer   :: Zix     ! charge number index for residual nucleus
  real(sgl) :: dTgrid  ! temperature increment
  real(sgl) :: T9max   ! Max Temperature in 10**9 K
  real(sgl) :: T9min   ! Min Temperature in 10**9 K
  real(sgl) :: Temp    ! nuclear temperature
  real(sgl) :: Teps    ! temperatures of basic grid in 10**9 K
  real(dbl) :: acm     ! inverse of reduced mass
  real(dbl) :: am      ! reduced mass
  real(dbl) :: b1      ! help variable
  real(dbl) :: b2      ! beta2
  real(dbl) :: b3      ! help variable
  real(dbl) :: de      ! single-nucleon exchange term J00
  real(dbl) :: de1     ! energy increment
  real(dbl) :: de2     ! energy increment
  real(dbl) :: deg     ! constant for Gamow energy
  real(dbl) :: deg1    ! constant for Gamow energy
  real(dbl) :: deg2    ! constant for Gamow energy
  real(dbl) :: del     ! difference
  real(dbl) :: e       ! energy
  real(dbl) :: eg0     ! Gamow energies at T9=1
  real(dbl) :: eg1     ! Gamow energies at T9(min) and T9(max)
  real(dbl) :: eg2     ! Gamow energies at T9(min) and T9(max)
  real(dbl) :: emax    ! maximal emission energy within bin decay
  real(dbl) :: emin    ! minimal emission energy
  real(dbl) :: Q0      ! Q-value
  real(dbl) :: Qa      ! Q-value
  real(dbl) :: qmax    ! maximum Q-value
  real(dbl) :: qmin    ! maximum Q-value
  real(dbl) :: Qn      ! Q-value
  real(dbl) :: Qp      ! Q-value
!
! ******** Set temperature grid for astrophysical calculations *********
!
  ET9 = 0.
  T9 = 0.
  T9(1) = 0.0001
  T9min = T9(1)
  if (astroT9.ne.0.) then
    T9max = astroT9
  else
    T9max = 10.
  endif
  Temp = 0.0001
  dTgrid = 0.0004
  nen = 1
  if (nTmax == 1) then
    T9(nen) = astroT9
    goto 100
  endif
   10 Temp = Temp + dTgrid
  Teps = Temp + 1.e-4
  nen = nen + 1
  if (nen > nTmax) goto 100
  T9(nen) = Temp
  if (Teps > 0.0005) dTgrid = 0.0005
  if (Teps > 0.001) dTgrid = 0.004
  if (Teps > 0.005) dTgrid = 0.005
  if (Teps > 0.01) dTgrid = 0.04
  if (Teps > 0.05) dTgrid = 0.05
  if (Teps > 0.05) dTgrid = 0.05
  if (Teps > 0.30) dTgrid = 0.1
  if (Teps > 1.) dTgrid = 0.50
  if (Teps > 4.) dTgrid = 1.
  if (Teps > 10.) goto 100
  goto 10
!
! ************************ incident energy grid ************************
!
100 do nen = 1, nTmax
    ET9(nen) = kT * T9(nen)
  enddo
  emin = 1.d-12
  if (k0 /= 1) emin = 1.d-3
  emax = 50.
  if (nTmax == 1 .and. k0 == 1) emax = min(emax, dble(50. * kT * astroT9))
  neg = 100
  neg1 = 10
  neg2 = 10
  Zix = parZ(k0)
  Nix = parN(k0)
  acm = 1. / specmass(Zix, Nix, k0)
  am = redumass(Zix, Nix, k0)
  Q0 = S(0, 0, k0)
  Qn = S(0, 0, 1)
  Qp = S(0, 0, 2)
  Qa = S(0, 0, 6)
  qmin = min(Qn, Qp, Qa)
  qmax = max(Qn, Qp, Qa)
  if (k0 == 0) then
    eg1 = max(0.d0, qmin - 4. * kT * T9max)
    eg2 = qmax + 4. * kT * T9max
    eg2 = max(eg2, 20.d0)
    b1 = Qn
    b2 = Qp
    b3 = Qa
  elseif (k0 == 1) then
    eg1 = max(dble(kT * T9min / 4.), 0.001d0)
    eg2 = kT * T9max * 4.
    b1 = - 999.
    b2 = Qp - Q0
    b3 = Qa - Q0
  elseif (k0 > 1) then
    eg0 = 0.122 * am **(1. / 3.) * (Ztarget * parZ(k0)) **(2. / 3.)
    eg1 = 0.122 * am **(1. / 3.) * (Ztarget * parZ(k0) * T9min) **(2. / 3.)
    deg1 = 0.237 * am **(1. / 6.) * (Ztarget * parZ(k0)) **(1. / 3.) * T9min **(5. / 6.)
    eg1 = eg1 - 2. * deg1
    eg2 = 0.122 * am **(1. / 3.) * (Ztarget * parZ(k0) * T9max) **(2. / 3.)
    deg2 = 0.237 * am **(1. / 6.) * (Ztarget * parZ(k0)) **(1. / 3.) * T9max **(5. / 6.)
    eg2 = eg2 + 2. * deg2
    b1 = Qn - Q0
    b2 = Qp - Q0
    if (k0 == 2 .and. k0 /= 6) b2 = Qa - Q0
    b3 = - 999.
    eg2 = max(eg2, b1, b2)
  endif
  if (eg1 < 0.) eg1 = emin
  nen = 0
  e = emin
  if (k0 > 1) then
    de1 = (log10(eg1) - log10(emin)) / float(neg1)
    deg = (log10(eg2) - log10(eg1)) / float(neg)
    de2 = (log10(emax) - log10(eg2)) / float(neg2)
  else
    de1 = (eg1 - emin) / float(neg1)
    deg = (eg2 - eg1) / float(neg)
    de2 = (emax - eg2) / float(neg2)
  endif
  110 continue
  if (e > emax) goto 200
  if (e < eg1) then
    de = de1
  elseif (e < eg2) then
    de = deg
  else
    de = de2
  endif
  if (k0 == 1) then
     de = 0.000002d00
     if (e >= 0.00001d00) de = 0.000010d00
     if (e >= 0.0001d00) de = 0.00010d00
     if (e >= 0.001d00) de = 0.0010d00
     if (e >= 0.010d00) de = 0.0040d00
     if (e >= 0.100d00) de = 0.0500d00
     if (e >= 1.000d00) de = 0.2500d00
     if (e >= 5.000d00) de = 0.5d00
     if (e >= 10.00d00) de = 2.d00
     if (e >= 20.00d00) de = 5.d00
  endif
  if (k0 <= 1) then
    if (abs(e - b1) < 2. * de .or. abs(e - b2) < 2. * de .or. abs(e - b3) < 2. * de) de = de / 10.
    e = e + de
  else
    del = 10. **(log10(e) + de) - e
    if (abs(e - b1) < 2. * del .or. abs(e - b2) < 2. * del .or. &
    abs(e - b3) < 2. * del .or. (e > eg0 .and. e < eg2)) de = de / 5.
    e = 10. **(log10(e) + de)
  endif
  nen = nen + 1
  if (nen > numenin + 2) then
    write(*, *) 'Astro-warning: too many energy points'
    nen = numenin + 2
    goto 200
  endif
  eninc(nen) = e * acm
  goto 110
  200 Ninc = nen
!
! The minimum and maximum value of all the incident energies is determined.
!
  enincmin = eninc(1)
  enincmax = eninc(Ninc)
  return
end subroutine egridastro
! Copyright A.J. Koning 2021
