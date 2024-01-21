subroutine breakup
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Contribution of breakup reactions
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl          ! single precision kind
! All global variables
!   numen        ! maximum number of outgoing energies
!   numpar       ! number of particles
! Variables for preequilibrium
!   Cbreak       ! adjustable parameter for break - up reactions
! Variables for main input
!   Atarget      ! mass number of target nucleus
!   k0           ! index of incident particle
!   Ztarget      ! charge number of target nucleus
! Variables for energy grid
!   deltaE       ! energy bin around outgoing energies
! Variables for energy grid
!   ebegin       ! first energy point of energy grid
!   egrid        ! outgoing energy grid
!   Einc         ! incident energy in MeV
! Variables for energies
!   eend         ! last energy point of energy grid
!   eninccm      ! center - of - mass incident energy in MeV
! Variables for nuclides
!   parskip      ! logical to skip outgoing particle
!   Q            ! Q - value
!   ZZ           ! charge number of residual nucleus
! Constants
!   amu          ! atomic mass unit in MeV
!   onethird     ! 1 / 3
!   parA         ! mass number of particle
!   parmass      ! mass of particle in a.m.u.
!   parZ         ! charge number of particle
!   sqrttwopi    ! sqrt(2. * pi)
! Variables for preequilibrium
!   Ca           ! effective Coulomb barrier
!   Deff         ! effective target - projectile evaporation
!   Ecent        ! centroid energy for emission spectrum
!   xspreeqbu    ! preequilibrium cross section per particle type and outgoing energy for break-up
!   Sab          ! separation energy for projectile
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: key                     ! keyword
  integer           :: nen                     ! energy counter
  integer           :: type                    ! particle type
  integer           :: type2                   ! particle type
  real(sgl)         :: arg                     ! help variable
  real(sgl)         :: Cb                      ! effective Coulomb barrier
  real(sgl)         :: Ehalf                   ! help variable
  real(sgl)         :: Emax                    ! maximal emission energy for particle channel
  real(sgl)         :: Eout                    ! outgoing energy
  real(sgl)         :: Epk                     ! peak energy
  real(sgl)         :: F                       ! full width at half maximum
  real(sgl)         :: fac1                    ! help variable
  real(sgl)         :: fac2                    ! help variable
  real(sgl)         :: factor                  ! multiplication factor
  real(sgl)         :: Feff                    ! effective full width at half maximum
  real(sgl)         :: H                       ! help variable
  real(sgl)         :: Ksig                    ! asymptotic potential
  real(sgl)         :: Nab(numpar, numpar)     ! break-up normalization constants
  real(sgl)         :: Napb(numpar, numpar)    ! break-up normalization constants
  real(sgl)         :: PBU(0:numen)            ! breakup distribution
  real(sgl)         :: PBUint                  ! integral over breakup distribution
  real(sgl)         :: r0                      ! effective radius parameter
  real(sgl)         :: step                    ! step function
  real(sgl)         :: TE                      ! barrier penetrability
  real(sgl)         :: Tpren                   ! barrier-penetrability factor
  real(sgl)         :: wi                      ! half width
  real(sgl)         :: width                   ! Full width at half maximum of the breakup nucleon  energy distribution
  real(sgl)         :: widthsig                ! width of barrier
  real(sgl)         :: wmin                    ! half width
  real(sgl)         :: wplus                   ! half width
  real(sgl)         :: xsbreakup               ! break-up cross section
!
! ************************** Kalbach model *****************************
!
! Break-up model by Kalbach: Phys. Rev. C95, 014606 (2017).
!
! We first determine the possible break-up channels and the other particle resulting from the break-up.
!
  do type = 1, 6
    do type2 = 1, 6
      Nab(type, type2) = 1.
      Napb(type, type2) = 1.
    enddo
  enddo
  Nab(3, 1) = 3.6
  Nab(3, 2) = 3.6
!
! Extra adjustment for (d,n) and (d,p) cross sections, TENDL-2021
!  Later withdrawn due to bad (d,n), (d,2n) etc fits
!
! Nab(3, 1) = Nab(3,1) * max(3. - Atarget/125., 1.)
! Nab(3, 2) = Nab(3,2) * max(6. - Atarget/50., 1.)
!
  Nab(4, 1) = 4.1
  Nab(4, 2) = 2.0
  Nab(4, 3) = 1.3
  Nab(5, 1) = 2.0
  Nab(5, 2) = 4.1
  Nab(5, 3) = 1.3
  Nab(6, 1) = 0.61
  Nab(6, 2) = 0.61
  Nab(6, 3) = 0.23
  Nab(6, 4) = 0.19
  Nab(6, 5) = 0.34
  Napb(6, 2) = 1.2
  Napb(6, 3) = 1.2
  Napb(5, 2) = 1.8
  Napb(5, 3) = 1.8
  Ksig = 112.
  widthsig = 13.
  do type = 1, 6
    if (parskip(type)) cycle
!
! (d,n) and (d,p)
!
    if (k0 == 3) then
      if (type /= 1 .and. type /= 2) cycle
      if (type == 1) type2 = 2
      if (type == 2) type2 = 1
    endif
!
! (t,d) and (t,p)
!
    if (k0 == 4) then
      if (type > 3) cycle
      if (type == 1) type2 = 3
      if (type == 2) type2 = 1
      if (type == 3) type2 = 1
    endif
!
! (h,d) and (h,p)
!
    if (k0 == 5) then
      if (type > 3) cycle
      if (type == 1) type2 = 2
      if (type == 2) type2 = 3
      if (type == 3) type2 = 2
    endif
!
! (a,n), (a,p), (a,d), (a,t) and (a,h)
!
    if (k0 == 6) then
      if (type == 6) cycle
      if (type == 1) type2 = 5
      if (type == 2) type2 = 4
      if (type == 3) type2 = 3
      if (type == 4) type2 = 2
      if (type == 5) type2 = 1
    endif
!
! Calculation of terms independent of emission energy.
!
! Centroid energy
!
    r0 = 1.2 + 5. / (1. + exp(Einc / 30.))
    Deff = r0 * (Atarget **onethird) + 1.2
    Ca = 1.44 * parZ(k0) * Ztarget / Deff
    Cb = 1.44 * parZ(type) * ZZ(0, 0, type) / Deff
    Ecent = parA(type) / real(parA(k0)) * (Einc - Ca) + Cb
!
! Full width at half maximum
!
! step     : step function
!
    arg = parA(k0) - parA(type) - 1.5
    if (arg < 0.) then
      step = 0.
    else
      step = 1.
    endif
    Sab = (parmass(type) + parmass(type2) - parmass(k0)) * amu
    F = 62. * (1. - 1. / exp(Einc / 173.)) * (1. - Atarget / (155. * Sab * Sab)) - 3. * step
!
! Effective full width at half maximum (asymmetric peaks)
!
    H = 0.5 * F
    Emax = eninccm + Q(type)
    Feff = H + min(H, 0.6 * (Emax - Ecent))
    width = max(Feff / 2.35, 0.1 * Einc)
    fac1 = 1. / (width * sqrttwopi)
    wplus = max(0., min(H, 0.6 * (Emax - Ecent))) / 2.35
    wmin = max(0., H - max(0., 0.6 * (Ecent - Emax))) / 2.35
    if (Emax >= Ecent - 1.67 * H .and. Emax <= Ecent) then
      Epk = Emax
    else
      Epk = Ecent
    endif
!
! Breakup cross section
!
! adjust   : subroutine for energy-dependent parameter adjustment
!
    Ehalf = 34. * (parA(k0) - parA(type)) **0.84
    Tpren = 1. / (1. + exp((Ehalf - Einc) / widthsig))
    key = 'cbreak'
    call adjust(Einc, key, 0, 0, type, 0, factor)
!
! Break-up term that depends on emission energy.
!
! First calculate normalization integral and spectrum function
!
    PBUint = 0.
    do nen = ebegin(type), eend(type)
      Eout = egrid(nen)
      if (Eout <= Epk) then
        wi = wmin
      else
        wi = wplus
      endif
      if (wi > 0.) then
        fac2 = 2. * (2. * wi) **2
        if (type > 1) then
          TE = 1. / (1. + exp(3. * (Cb - Eout) / Cb))
        else
          TE = 1.
        endif
        PBU(nen) = fac1 * exp( - (Eout - Epk) **2 / fac2) * TE
        PBUint = PBUint + PBU(nen) * deltaE(nen)
      endif
    enddo
    if (PBUint > 0.) then
      xsbreakup = factor * Cbreak(type) * (Nab(k0, type) * Deff * Deff + Napb(k0, type) * Deff) * exp(Einc / Ksig) * Tpren
      do nen = ebegin(type), eend(type)
        xspreeqbu(type, nen) = xsbreakup * PBU(nen)
      enddo
    endif
  enddo
  return
end subroutine breakup
! Copyright A.J. Koning 2021
