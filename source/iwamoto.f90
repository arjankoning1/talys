subroutine iwamoto(Z, A, Sn, Ein, Tmadjust, Fsadjust, ETfns, NETfns, pfns, maxpfns, pfnscm, Eav)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Iwamoto systematics for PFNS
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: Ngrid=20000         ! number of energy grid points
  integer            :: NETfns              ! number of energy points
  integer            :: A                   ! mass number of compound nucleus
  integer            :: ifrag               ! counter
  integer            :: nen                 ! energy counter
  integer            :: nen2                ! energy counter
  integer            :: Z                   ! charge number of compound nucleus
  integer            :: AH                  ! mass number of heavy fragment
  integer            :: AL                  ! mass number of light fragment
  real               :: ald                 ! level density parameter
  real               :: dETfns(NETfns)      ! delta energy
  real               :: dE                  ! delta energy
  real               :: e                   ! energy
  real               :: Ef                  ! energy of fragment
  real               :: Eav                 ! average energy
  real               :: EfH                 ! energy of heavy fragment
  real               :: EfL                 ! energy of light fragment
  real               :: Ein                 ! incident energy
  real               :: ETfns(NETfns)       ! energy
  real               :: Ex                  ! excitation energy
  real               :: Emax                ! maximal emission energy for particle channel
  real               :: Fs                  ! fraction of scission neutrons
  real               :: labfns              ! FNS in lab
  real               :: phiTav(Ngrid)       ! phi(T) function
  real               :: Egrid(0:Ngrid)      ! outgoing energy grid
  real               :: sigma(Ngrid)        ! help variable
  real               :: Tfns(NETfns)        ! total fission neutron spectrum
  real               :: pfns0               ! prompt fission neutron spectrum
  real               :: sfns0               ! scission neutron spectrum
  real               :: pfns(NETfns)        ! prompt fission neutron spectrum
  real               :: pfnsfrag(2,NETfns)  ! prompt fission neutron spectrum per fragment
  real               :: prompt              ! prompt part of FNS
  real               :: promptCM            ! prompt part of FNS in CM
  real               :: scission(NETfns)    ! scission part of FNS
  real               :: maxwell(NETfns)     ! Maxwell spectrum
  real               :: pfnscm(NETfns)      ! PFNS in CM
  real               :: maxpfns(NETfns)     ! Maxwell ratio
  real               :: maxfactor           ! factor for Maxwellian
  real               :: pi                  ! pi
  real               :: Sn                  ! neutron separation energy
  real               :: summax              ! help variable
  real               :: sumTfns             ! integrated TFNS
  real               :: EsumTfns            ! integrated TFNS * energy
  real               :: sumpfns             ! integrated PFNS
  real               :: sumpfnscm           ! integrated PFNS
  real               :: sumscission         ! integrated scission
  real               :: sumpfnsfrag(2)      ! integrated PFNS per fragment
  real               :: TKE                 ! total kinetic energy
  real               :: Tm0                 ! temperature
  real               :: Tm                  ! temperature
  real               :: ZA13                ! help variable
  real               :: Tmadjust            ! adjustable parameter for maximum temperature
  real               :: Fsadjust            ! adjustable parameter for scission neutrons
!
! ************************ Iwamoto systematics *************************
!
! PFNS systematics of O. Iwamoto: Journ. Nuc. Sci. Techn 45, 910 (2008).
!
! TKE, average light and heavy fragment mass
!
  ZA13 = Z / (A **(1. / 3.))
  TKE = 0.1189 * Z * ZA13 + 7.3
  AH = 145
  AL = A - AH
  EfH = AL / real(AH) * TKE / A
  EfL = AH / real(AL) * TKE / A
!
! Maximum temperature
!
  Tm0 = Tmadjust * (0.45 * ZA13 - 5.82)
  Ex = max(Ein + Sn, 0.)
  ald = A / 10.
  Tm = sqrt(Tm0 * Tm0 + Ex / ald)
  Emax = ETfns(NETfns)
!
! Scission neutrons
!
  Fs = Fsadjust * (0.11 / (1. + exp((15. - ZA13) / 0.1)))
!
! Inverse reaction cross sections
!
  call inverse_xs(Egrid, sigma, Ngrid, Emax)
!
! Prompt fission neutron spectrum in CM
!
  call cmfns(Egrid, sigma, Ngrid, Tm, phiTav, Emax)
!
! Prompt fission and scission neutron spectrum in Lab
!
  EsumTfns = 0.
  sumTfns = 0.
  sumpfns = 0.
  sumpfnscm = 0.
  sumpfnsfrag = 0.
  sumscission = 0.
  dETfns = 0.
  do nen = 2, NETfns - 1
    dETfns(nen) = 0.5 * (ETfns(nen+1) - ETfns(nen-1))
  enddo
  dETfns(1) = ETfns(1)
  dETfns(NETfns) = 0.5 * (ETfns(NETfns) - ETfns(NETfns-1))
  Tfns = 0.
  scission = 0.
  pfns = 0.
  pfnscm = 0.
  do nen = 1, NETfns
    E = ETfns(nen)
    dE = dETfns(nen)
!
! Scission
!
    nen2 = max(int(E*Ngrid/Emax), 1)
    sfns0 = phiTav(nen2)
    scission(nen) = Fs * sfns0
!
! Loop over light and heavy fragment
!
    do ifrag = 1, 2
      if (ifrag == 1) then
        Ef = EfL
      else
        Ef = EfH
      endif
      pfns0 = labfns(E, Ef, phiTav, Egrid, Ngrid)
      prompt = (1. - Fs) * pfns0
      promptcm = (1. - Fs) * phiTav(nen2)
      pfnsfrag(ifrag, nen) = prompt
      pfns(nen) = pfns(nen) + 0.5 * prompt
      pfnscm(nen) = pfnscm(nen) + 0.5 * promptcm
      sumpfnsfrag(ifrag) = sumpfnsfrag(ifrag) + pfnsfrag(ifrag, nen) * dE
    enddo
    Tfns(nen) = pfns(nen) + scission(nen)
    sumscission = sumscission + scission(nen) * dE
    sumpfns = sumpfns + pfns(nen) * dE
    sumpfnscm = sumpfnscm + pfnscm(nen) * dE
    EsumTfns = EsumTfns + Tfns(nen) * E * dE
    sumTfns = sumTfns + Tfns(nen) * dE
  enddo
  Eav = EsumTfns / sumTfns
!
! Maxwellian ratio
!
  summax = 0.
  pi = 3.14159265358979323
  maxfactor = 2. / sqrt (pi * Eav**3)
  if (Eav > 0.) then
    do nen = 1, NETfns
      E = ETfns(nen)
      dE = dETfns(nen)
      maxwell(nen) = maxfactor * sqrt(E) * exp(-E/(2./3.*Eav))
      summax = summax + maxwell(nen) * dE
    enddo
    if (summax > 0.) then
      do nen = 1, NETfns
        maxwell(nen) = maxwell(nen) / summax
        if (maxwell(nen) > 0.) maxpfns(nen) = Tfns(nen) / maxwell(nen)
      enddo
    endif
  endif
!
! Output
!
  open (unit = 2, file = 'Tfns.tot', status = 'unknown')
  write(2,'("# Fission neutron spectrum in Lab - Iwamoto model ")')
  write(2,'("# Z(CN)   :",i4)') Z
  write(2,'("# A(CN)   :",i4)') A
  write(2,'("# S(N)    :",es12.5, " MeV")') Sn
  write(2,'("# Ein     :",es12.5, " MeV")') Ein
  write(2,'("# Ex      :",es12.5, " MeV")') Ex
  write(2,'("# Tm0     :",es12.5)') Tm0
  write(2,'("# Tm      :",es12.5)') Tm
  write(2,'("# ald     :",es12.5, " MeV^-1")') ald
  write(2,'("# AH      :",i4)') AH
  write(2,'("# AL      :",i4)') AL
  write(2,'("# EfH     :",es12.5, " MeV")') EfH
  write(2,'("# EfL     :",es12.5, " MeV")') EfL
  write(2,'("# Tmadjust:",es12.5)') Tmadjust
  write(2,'("# Fsadjust:",es12.5)') Fsadjust
  write(2,'("# N       :",i4)') NETfns
  write(2,'("#   E-out        TFNS           PFNS          L frag         H frag        Scission     Maxwell ratio    PFNS (CM)")')
  do nen = 1, NETfns
    E = ETfns(nen)
    write(2,'(f10.5,7es15.6)') E, Tfns(nen), pfns(nen), pfnsfrag(1, nen), pfnsfrag(2, nen), scission(nen), maxpfns(nen), pfnscm(nen)
  enddo
  write(2, '("# Eav       :",es12.5," MeV")') Eav
  write(2, '("# FNS int.  :  sumTfns        sumpfns      sumpfnsfrag(1) sumpfnsfrag(2)  sumscission      summax       sumpfnsCM")')
  write(2, '("#         ",7es15.6)') sumTfns, sumpfns, sumpfnsfrag(1), sumpfnsfrag(2), sumscission, summax, sumpfnscm
  close(2)
  return
end subroutine iwamoto
! Copyright A.J. Koning 2021
subroutine inverse_xs(Egrid, sigma, Ngrid, Emax)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Inverse reaction cross section
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  implicit none
  integer :: Ngrid             ! number of energy grid points
  integer :: nen               ! energy counter
  real    :: pi                ! pi
  real    :: hbar              ! Planck's constant / 2.pi in MeV.s
  real    :: clight            ! speed of light in vacuum in m/s
  real    :: hbarc             ! hbar.c in MeV.fm
  real    :: amu               ! atomic mass unit in MeV
  real    :: parmass           ! mass of particle in a.m.u.
  real    :: A100              ! average FF
  real    :: r0                ! effective radius parameter
  real    :: R                 ! radius
  real    :: sigma0            ! constant for inverse c.s.
  real    :: S0                ! neutron strength function
  real    :: sqrt1eV           ! parameter for correct dimension
  real    :: alpha             ! constant
  real    :: Emax              ! maximal emission energy
  real    :: dEgrid            ! energy bin
  real    :: e                 ! energy
  real    :: Egrid(0:Ngrid)    ! outgoing energy grid
  real    :: sigma(Ngrid)      ! inverse cross section
!
! Constants for inverse reaction cross section (Iwamoto)
!
  pi = 3.14159265358979323
  hbar = 6.5821220e-22
  clight = 2.99792458e8
  hbarc = hbar * clight * 1.e15
  amu = 931.49386
  parmass = 1.008664904 * amu
  A100 = 100.
  r0 = 1.2
  R = r0 * (A100 **(1. / 3.))
  sigma0 = pi * R * R
  S0 = 1.e-4
  sqrt1eV = 1.e-3
  alpha = pi * hbarc * hbarc * S0 / (parmass * r0 * r0 * A100 **(2. / 3.) * sqrt1eV)
  dEgrid = Emax / Ngrid
  open (unit = 1, file = 'reaction_pfns.tot', status = 'unknown')
  write(1 ,'("# Neutron reaction cross sections from Iwamoto systematics")')
  write(1,'("# N       :",i6)') Ngrid
  Egrid(0) = 0.
  do nen = 1, Ngrid
    Egrid(nen) = dEgrid * nen
    e = Egrid(nen)
    sigma(nen) = sigma0 * (1. + alpha / sqrt(e)) * 10.
    write(1, * ) e, sigma(nen)
  enddo
  close(1)
  return
end subroutine inverse_xs
! Copyright A.J. Koning 2021
subroutine cmfns(Egrid, sigma, Ngrid, Tm, phiTav, Emax)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : FNS in CM system
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! ****************** Declarations and common blocks ********************
!
  implicit none
  integer, parameter :: nT=100 ! number of temperatures
  integer :: nen               ! energy counter
  integer :: Ngrid             ! number of energy grid points
  integer :: iT                ! temperature counter
  real    :: Egrid(0:Ngrid)    ! outgoing energy
  real    :: dEgrid            ! width of energy bin
  real    :: Emax              ! maximal emission energy for particle channel
  real    :: dT                ! width of temperature bin
  real    :: PT                ! temperature distribution
  real    :: T                 ! temperature
  real    :: Tm                ! maximum temperature
  real    :: e                 ! energy
  real    :: sum               ! help variable
  real    :: sumphi            ! help variable
  real    :: expo              ! exponent
  real    :: kT                ! normalization constant
  real    :: phi0(Ngrid, nT)   ! FNS spectrum
  real    :: phiTav(Ngrid)     ! temperature-averaged FNS spectrum
  real    :: sigma(Ngrid)      ! inverse cross ection
!
! Prompt fission neutron spectrum in CM
!
  open (unit = 2, file = 'phiTav.tot', status = 'unknown')
  write(2, '("# Prompt fission neutron spectrum in CM - Iwamoto model")')
  dT = Tm / nT
  sumphi = 0.
  dEgrid = Emax / Ngrid
  do iT = 1, nT
    T = (iT - 0.5) * dT
    PT = 2. * T / (Tm * Tm)
    sum = 0.
    do nen = 1, Ngrid
      e = Egrid(nen)
      expo = exp( - e / T)
      sum = sum + sigma(nen) * e * expo
    enddo
    kT = sum * dEgrid
    if (kT > 0.) then
      do nen = 1, Ngrid
        e = Egrid(nen)
        expo = exp( - e / T)
        phi0(nen, iT) = sigma(nen) * e * expo / kT
      enddo
    endif
  enddo
  write(2,'("# N       :",i6)') Ngrid
  do nen = 1, Ngrid
    sum = 0.
    do iT = 1, nT
      T = (iT - 0.5) * dT
      PT = 2. * T / (Tm **2)
      sum = sum + phi0(nen, iT) * PT
    enddo
    phiTav(nen) = sum * dT
    sumphi = sumphi + phiTav(nen) * dEgrid
    write(2, * ) Egrid(nen), phiTav(nen)
  enddo
  write(2, '("#Integral: ",es12.5)') sumphi
  close(2)
  return
end subroutine cmfns
! Copyright A.J. Koning 2021
function labfns(E, Ef, phi, Egrid, Ngrid)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : FNS in Lab system
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! ****************** Declarations and common blocks ********************
!
  implicit none
  integer :: nen            ! energy counter
  integer :: Ngrid          ! number of energy grid points
  real    :: dE             ! width of energy bin
  real    :: E              ! incident energy
  real    :: Ef             ! energy of FF
  real    :: Eg             ! outgoing energy
  real    :: Egrid(0:Ngrid) ! outgoing energy
  real    :: Emin           ! minimum energy for integration
  real    :: Eplus          ! maximum energy for integration
  real    :: factor         ! multiplication factor
  real    :: labfns         ! FNS in LAB system
  real    :: phi(Ngrid)     ! FNS in CM system
  real    :: sum            ! help variable
  real    :: term           ! help variable
!
! ************************ Lab fission neutron spectrum ****************
!
  Emin = (sqrt(E) - sqrt(Ef)) **2
  Eplus = (sqrt(E) + sqrt(Ef)) **2
  factor = 1. / (4. * sqrt(Ef))
  sum = 0.
  do nen = 1, Ngrid
    Eg = Egrid(nen)
    if (Eg > Emin .and. Eg <= Eplus) then
      dE = Eg - max(Egrid(nen - 1), Emin)
      term = phi(nen) / sqrt(Eg)
      sum = sum + term * dE
    endif
  enddo
  labfns = factor * sum
  return
end function labfns
! Copyright A.J. Koning 2021
