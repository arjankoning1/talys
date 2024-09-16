subroutine foldalpha(Zix, Nix, E)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Purpose of this module
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
!   sgl           ! single precision kind
! All global variables
!   numjlm        ! maximum number of radial points
! Variables for masses
!   beta2         ! deformation parameter
! Variables for OMP
!   adepthcor     ! adjustable parameter for depth of DF alpha potential
!   alphaomp      ! alpha optical model
!   aradialcor    ! adjustable parameter for shape of DF alpha potential
!   avadjust      ! adjustable factor for OMP (default 1.)
!   awdadjust     ! adjustable factor for OMP (default 1.)
!   d1adjust      ! adjustable factor for OMP (default 1.)
!   ompadjustp    ! flag for local optical model parameter adjustment
!   rvadjust      ! adjustable factor for OMP (default 1.)
!   rwdadjust     ! adjustable factor for OMP (default 1.)
!   w1adjust      ! adjustable factor for OMP (default 1.)
! Variables for nuclides
!   AA            ! mass number of residual nucleus
!   ZZ            ! charge number of residual nucleus
! Constants
!   amu           ! atomic mass unit in MeV
!   parmass       ! mass of particle in a.m.u.
!   pi            ! pi
! Variables for levels
!   edis          ! energy of level
! Variables for optical model
!   rc            ! Coulomb radius
! Variables for JLM
!   potjlm        ! JLM potential depth values
!   radjlm        ! radial points for JLM potential
!   rhojlmn       ! density for neutrons
!   rhojlmp       ! density for protons
! Variables for masses
!   expmexc       ! experimental mass excess
!   thmexc        ! theoretical mass excess
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: key               ! keyword
  real(sgl)         :: a                 ! fit variables
  real(sgl)         :: a13               ! mass**1/3
  real(sgl)         :: a2                ! Myers-Swiatecki parameter
  real(sgl)         :: a3                ! mass**3
  real(sgl)         :: aa1               ! help variable
  real(sgl)         :: aaw               ! help variable
  real(sgl)         :: alphav(10)        ! optical model part  alphav = 1 : reduced radius r0=r/a13 of real pot
  real(sgl)         :: am1               ! mass excess of the target-projectile system
  real(sgl)         :: am2               ! mass excess of the target+projectile system
  real(sgl)         :: am3               ! mass excess of the projectile
  real(sgl)         :: as                ! diffuseness of the functional form of the imaginary volume integral
  real(sgl)         :: aws               ! diffuseness of imaginary surface potential
  real(sgl)         :: aww               ! diffuseness of imaginary volume potential
  real(sgl)         :: bb1               ! help variable
  real(sgl)         :: c                 ! curvature of neck
  real(sgl)         :: dvolj             ! volume integral
  real(sgl)         :: dws               ! dispersive contribution to the real potential arising from the surface imaginary pot
  real(sgl)         :: dwv               ! dispersive contribution to the real potential arising from the volume imaginary pot
  real(sgl)         :: e                 ! energy
  real(sgl)         :: e2exp             ! lowest inelastic threshold
  real(sgl)         :: ee                ! energy
  real(sgl)         :: efermia           ! Fermi energy
  real(sgl)         :: eref              ! reference energy
  real(sgl)         :: es                ! energy threshold of the functional form of the imaginary volume integral
  real(sgl)         :: expj0             ! saturation value of the imaginary volume integral
  real(sgl)         :: factor1           ! help variable
  real(sgl)         :: factor2           ! help variable
  real(sgl)         :: hh                ! radius difference
  real(sgl)         :: radmom(numjlm)    ! radial grid for the potential, identical as the one used for JLM
  real(sgl)         :: raw               ! radius
  real(sgl)         :: rb                ! maximum radius value
  real(sgl)         :: rhomom(numjlm)    ! total density (neutron + proton) at a given radius radmom
  real(sgl)         :: rr                ! running variable in integration over the radius
  real(sgl)         :: rva(1000)         ! radial coordinate for the double folding potential
  real(sgl)         :: rws               ! radius of imaginary surface potential
  real(sgl)         :: rww               ! radius of imaginary volume potential
  real(sgl)         :: te                ! energy-dependent term of the imaginary volume integral
  real(sgl)         :: v5d               ! OMP component
  real(sgl)         :: va(1000)          ! double folding potential
  real(sgl)         :: vdis              ! dispersive contribution to the real potential arising from the surface imaginary pot
  real(sgl)         :: vdiv              ! dispersive contribution to the real potential arising from the volume imaginary pot
  real(sgl)         :: vopr(numjlm)      ! Final real part of the alpha optical potential
  real(sgl)         :: wopr(numjlm)      ! Final imaginary part of the alpha optical potential
  real(sgl)         :: ws                ! surface component of the imaginary potential
  real(sgl)         :: ww                ! weight
  real(sgl)         :: wws               ! depth of imaginary surface potential
  real(sgl)         :: www               ! depth of imaginary volume potential
  real(sgl)         :: z                 ! charge number
  integer           :: i                 ! level
  integer           :: k                 ! designator for particle
  integer           :: khi               ! help variable
  integer           :: kk                ! counter
  integer           :: klo               ! help variable
  integer           :: Nix               ! neutron number index for residual nucleus
  integer           :: nradrho           ! number of grid point used in radmom and rhomom
  integer           :: nu                ! number of radial grid point
  integer           :: Zix               ! charge number index for residual nucleus
!
! ************************ Alpha Optical potential *************************
!
  a = AA(Zix, Nix, 0)
  z = ZZ(Zix, Nix, 0)
  ee = E
  do i = 1, numjlm
    radmom(i) = radjlm(Zix, Nix, i)
    rhomom(i) = rhojlmn(Zix, Nix, i, 1) + rhojlmp(Zix, Nix, i, 1)
  enddo
  e2exp = edis(Zix, Nix, 1)
  a2 = a **2
  a3 = a **3
  a13 = a **(1. / 3.)
!
!--------------------------------------------------------------------------------------
! alphav =  1: reduced radius r0=r/a13 of real pot
! as         : diffuseness of the functional form of the imaginary volume integral
! es         : energy threshold of the functional form of the imaginary volume integral
!--------------------------------------------------------------------------------------
!
  v5d = 0.
  goto (10, 20, 30, 200) alphaomp - 2
   10 continue
!--------------------------------------------------------------------------------------
! Case alphaomp=3
! Imaginary potential OMP I of Demetriou et al. (2002)
! Woods-Saxon E-dependent volume imaginary potential
!--------------------------------------------------------------------------------------
  alphav(1) = 1.25
  raw = 0.85281 + 0.02202 * a - 2.14551e-4 * a2 + 7.92942e-7 * a3 - 9.94658e-10 * a2 * a2
  aaw = - 0.13526 + 0.02029 * a - 1.98441e-4 * a2 + 7.35104e-7 * a3 - 9.15272e-10 * a2 * a2
  if (aaw <= 1.e-2) aaw = 1.e-2
  as = 7.65867 - 7.5669 * e2exp + 2.50486 * e2exp **2
  es = 0.1024 * a + 1.1307 * as
  te = 1. / (1. + exp( - (ee - es) / as))
  expj0 = 77.
  if (a <= 90.) expj0 = 135. - 0.644 * a
  ww = 3. / pi / raw **3 * expj0 * te
  alphav(3) = raw
  alphav(4) = aaw
  alphav(8) = ww
  alphav(5) = 0.
  dwv = 0.
  dws = 0.
  goto 50
!
   20 continue
!--------------------------------------------------------------------------------------
! Case alphaomp=4
! Imaginary potential OMP II of Demetriou et al. (2002)
! Woods-Saxon E-dependent volume+surface imaginary potential
!--------------------------------------------------------------------------------------
  alphav(1) = 1.25
  raw = 1.47385 - 0.00134615 * a
  aaw = 0.29
  alphav(5) = 0.9
  as = 7.65867 - 7.5669 * e2exp + 2.50486 * e2exp **2
  es = 0.1024 * a + 1.1307 * as
  te = 1. / (1. + exp( - (ee - es) / as))
  expj0 = 77.
  if(a <= 90.) expj0 = 135. - 0.644 * a
  ww = expj0 * te / (raw **3 / 3. + 7.603 * alphav(5) * aaw * raw **2 / a13) / pi
! ww=expj0*te/((1.-alphav(5))*raw**3/3.+7.603*alphav(5)*aaw*raw**2/a13)/pi
  alphav(3) = raw
  alphav(4) = aaw
  alphav(8) = ww
  dwv = 0.
  dws = 0.
  goto 50
!
   30 continue
!--------------------------------------------------------------------------------------------------
! Case alphaomp=5
! Imaginary potential OMP III of Demetriou et al. (2002)
! Woods-Saxon E-dependent volume+surface imaginary potential plus Dispersion contribution
!  using prescription and functions of Capote et al. J. Phys. G 27 (2001) B15
!
!--------------------------------------------------------------------------------------------------
!
  if (expmexc(Zix - 2, Nix - 2) /= 0..and.expmexc(Zix + 2, Nix + 2) /= 0.) then
    am1 = expmexc(Zix + 2, Nix + 2)
    am2 = expmexc(Zix - 2, Nix - 2)
  else
    am1 = thmexc(Zix + 2, Nix + 2)
    am2 = thmexc(Zix - 2, Nix - 2)
  endif
  am3 = (parmass(6) - 4.) * amu
  efermia = - 0.5 * (am1 - am2 + 2. * am3)
!
  alphav(1) = 1.25
  raw = 1.47385 - 0.00134615 * a
  aaw = 0.30
  alphav(5) = 0.35
  v5d = alphav(5)
  c = 0.005
  if(ee < 13.) c = - 0.165 * ee + 2.15
  alphav(5) = alphav(5) * exp( - c * abs(ee - efermia))
!
  as = 7.65867 - 7.5669 * e2exp + 2.50486 * e2exp **2
  es = 0.0854 * a + 1.1307 * as
  te = 1. / (1. + exp( - (ee - es) / as))
  expj0 = 75.
  if(a <= 96.) expj0 = 135. - 0.644 * a
! volume integral
  ww = expj0 * te / ((1. - alphav(5)) * raw **3 / 3. + 7.603 * alphav(5) * aaw * raw **2 / a13) / pi
  alphav(3) = raw
  alphav(4) = aaw
  alphav(8) = ww
!
   50 continue

! final radius, diffuseness and depth of the imaginary WS-type potential
! They are multiplied by the OMP adjustment keywords of TALYS
! (default 1.)
!   volume term
  rww = rvadjust(6) * alphav(3) * a13
  aww = avadjust(6) * alphav(4)
  if (alphaomp == 4) then
    www = w1adjust(6) * alphav(8)
  else
    www = w1adjust(6) * alphav(8) * (1. - alphav(5))
  endif
!   surface term
  rws = rwdadjust(6) * 1.09 * rww
  aws = awdadjust(6) * 1.6 * aww
  wws = d1adjust(6) * alphav(8) * alphav(5)
!
! dispersive contributions: dwv = volume; dws = surface; dvolj = volume integral
! new method for dispersive relation following Mahaux, Ngo and Satchler, NPA 449 (1986) 354
!
  if (alphaomp == 5) then
    eref = 150.
    call mahaux(a, ee, eref, expj0, as, es, v5d, raw, aaw, dwv, dws, dvolj, efermia)
  endif
!--------------------------------------------------------------------------------------------------
! determination of the real folding potential through the product of the fourier transforms
!--------------------------------------------------------------------------------------------------
  rb = radjlm(Zix, Nix, numjlm)
  nu = 1000
  nradrho = numjlm

  call afold(z, a, ee, rb, nu, rva, va, radmom, rhomom, nradrho)
!
! Depth and shape are multiplied by the OMP adjustment keywords of TALYS
! (default 1.)
!
! adjust    : subroutine for energy-dependent parameter adjustment
!
  if (ompadjustp(6)) then
    key = 'aradialcor'
    call adjust(E, key, Zix, Nix, 0, 0, factor1)
    key = 'adepthcor'
    call adjust(E, key, Zix, Nix, 0, 0, factor2)
  else
    factor1 = 1.
    factor2 = 1.
  endif
  do i = 1, nu
!   rva(i)=factor1*aradialcor*rva(i)
!
!sg Correction of the radius dependence of the real part after an analysis of the
!   of the (a,g) and (a,n) data of deformed nuclei : 27/4/2018 (Brussels)
!   increase of rva by 3% for deformed nuclei but only below typically 18 MeV

    rva(i) = factor1 * aradialcor * rva(i) * (1. + abs(beta2(Zix, Nix, 0)) / 15. / (1. + exp((E - 18.) / 2.)))
    va(i) = factor2 * adepthcor * va(i)
  enddo
!
! interpolation of the alpha optical potential va(rva) on the given radial grid radjlm
!
  do k = 1, numjlm
    rr = radjlm(Zix, Nix, k)
    klo = 1
    khi = nu
    if (rr <= rva(klo)) then
       khi = klo + 1
       goto 110
    endif
    if (rr >= rva(khi)) then
       klo = khi - 1
       goto 110
    endif
  120   if (khi - klo > 1) then
      kk = (khi + klo) / 2
      if (rva(kk) > rr) then
        khi = kk
      else
      klo = kk
      endif
      goto 120
    endif
  110   hh = rva(khi) - rva(klo)
    aa1 = (rva(khi) - rr) / hh
    bb1 = (rr - rva(klo)) / hh
    vopr(k) = aa1 * va(klo) + bb1 * va(khi)
! dispersive contributions to real potential for alphaomp=5  (OMP III)
    if (alphaomp < 5) then
      vdiv = 0.
      vdis = 0.
    else
      if (abs((rr - rww) / aww) < 88.) then
        vdiv = - dwv / (1. + exp((rr - rww) / aww))
      else
        vdiv = 0.
      endif
      if (abs((rr - rws) / aws) < 88.) then
        vdis = - 4. * dws * exp((rr - rws) / aws) / (1. + exp((rr - rws) / aws)) **2
      else
        vdis = 0.
      endif
    endif
    vopr(k) = vopr(k) + vdiv + vdis
! final imaginary potential
    wopr(k) = 0.
    if (abs((rr - rww) / aww) < 88.) wopr(k) = - www / (1. + exp((rr - rww) / aww))
    if (abs((rr - rws) / aws) < 88.) then
      ws = - 4. * wws * exp((rr - rws) / aws) / (1. + exp((rr - rws) / aws)) **2
    else
      ws = 0.
    endif
    wopr(k) = wopr(k) + ws
!-----------------------------------------------------------------
! potjlm   1 : final real central potential
! potjlm   2 : final imaginary central potential
!-----------------------------------------------------------------
    potjlm(Zix, Nix, k, 1) = vopr(k)
    potjlm(Zix, Nix, k, 2) = wopr(k)
    potjlm(Zix, Nix, k, 3) = 0.
    potjlm(Zix, Nix, k, 4) = 0.
    potjlm(Zix, Nix, k, 5) = 0.
    potjlm(Zix, Nix, k, 6) = 0.
  enddo
! calculate the coulomb radius based on elton's formula
  rc = 1.123 + 2.352 * (a **( - .666666)) - 2.07 * (a **( - 1.333333))
  200 continue
  return
end subroutine foldalpha
! Copyright A.J. Koning 2021
