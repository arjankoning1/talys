subroutine neck(Z, A, fmass, fmasscor, fmz, fmzcor, ap, edefo, elt, crel)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Fission fragment mass yields per fission mode based on RNRM
!
! Author    : Marieke Duijvestijn
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
!   numelem       ! number of elements
!   nummass       ! number of masses
! Variables for fission
!   flagffevap    ! flag for calculation of particle evaporati
! Variables for main input
!   Ninit         ! neutron number of initial compound nucleus
!   Zinit         ! charge number of initial compound nucleus
! Constants
!   pi            ! pi
! Variables for Brosa model
!   amm           ! parameter for neck rupture
!   c0            ! curvature of neck
!   di            ! nucleon number density
!   ess           ! parameter for neck rupture
!   r1            ! parameter for neck rupture
!   r2            ! parameter for neck rupture
!   r3            ! parameter for neck rupture
!   rest          ! parameter for neck rupture
!   rp            ! parameter for neck rupture
!   rpt           ! parameter for neck rupture
!   rt            ! parameter for neck rupture
!   totl          ! parameter for neck rupture
!   vtot          ! total volume of the complex
!   z1            ! parameter for neck rupture
!   z2            ! parameter for neck rupture
!   z3            ! parameter for neck rupture
!   zee           ! parameter for neck rupture
!
! *** Declaration of local data
!
  implicit none
  integer   :: A                         ! mass number of target nucleus
  integer   :: i                         ! counter
  integer   :: iaf1                      ! help variable
  integer   :: iaf2                      ! help variable
  integer   :: imax                      ! help variable
  integer   :: irn1                      ! help variable
  integer   :: irn2                      ! help variable
  integer   :: izf1                      ! help variable
  integer   :: izf2                      ! help variable
  integer   :: izloop                    ! counter
  integer   :: izmax                     ! maximum Z value
  integer   :: izstepnum                 ! counter
  integer   :: jimax                     ! help variable
  integer   :: jmx                       ! help variable
  integer   :: k                         ! designator for particle
  integer   :: mcount                    ! counter
  integer   :: Nix                       ! neutron number index for residual nucleus
  integer   :: Z                         ! charge number of target nucleus
  integer   :: Zix                       ! charge number index for residual nucleus
  real(sgl) :: a1                        ! Myers-Swiatecki parameter
  real(sgl) :: a2                        ! Myers-Swiatecki parameter
  real(sgl) :: af1(4000)                 ! help variable
  real(sgl) :: af2(4000)                 ! help variable
  real(sgl) :: ald                       ! level density parameter
  real(sgl) :: aloop                     ! help variable
  real(sgl) :: am1                       ! mass excess of the target-projectile system
  real(sgl) :: am2                       ! mass excess of the target+projectile system
  real(sgl) :: amin                      ! maximal and minimal mass number defining range of the fit
  real(sgl) :: ap                        ! proton level density parameter
  real(sgl) :: astepnum                  ! help variable
  real(sgl) :: astepsize                 ! help variable
  real(sgl) :: at                        ! mass of residual nucleus
  real(sgl) :: atot                      ! mass number
  real(sgl) :: b1                        ! help variable
  real(sgl) :: b2                        ! beta2
  real(sgl) :: bcom                      ! help variable
  real(sgl) :: bind01                    ! binding energy
  real(sgl) :: bind02                    ! binding energy
  real(sgl) :: bind1                     ! binding energy
  real(sgl) :: bind2                     ! binding energy
  real(sgl) :: coul12                    ! Coulomb energy
  real(sgl) :: coulel                    ! Coulomb energy
  real(sgl) :: crel                      ! scaling factor for neck curvature
  real(sgl) :: d                         ! parameter for energy smoothing
  real(sgl) :: de                        ! single-nucleon exchange term J00
  real(sgl) :: delt                      ! help variable
  real(sgl) :: dum                       ! dummy value
  real(sgl) :: dumm                      ! help variable
  real(sgl) :: edefo                     ! deformation energy
  real(sgl) :: elt                       ! help variable
  real(sgl) :: eob                       ! help variable
  real(sgl) :: es                        ! energy threshold of the functional form of the imaginary volume integral
  real(sgl) :: es1                       ! help variable
  real(sgl) :: es2                       ! help variable
  real(sgl) :: expo                      ! help variable
  real(sgl) :: ezdis                     ! energy
  real(sgl) :: ezdisnorm                 ! normalized energy
  real(sgl) :: fimin                     ! fucntion value of fmin
  real(sgl) :: fmass(nummass)            ! fission fragment mass yield
  real(sgl) :: fmasscor(nummass)         ! corrected fission fragment mass yield
  real(sgl) :: fmin                      ! minimum function value
  real(sgl) :: fmz(nummass, numelem)     ! fission fragment isotope yield
  real(sgl) :: fmzcor(nummass, numelem)  ! corrected fission fragment isotope yield
  real(sgl) :: gam                       ! Brosa parameter
  real(sgl) :: ignatyuk                  ! function for energy dependent level density parameter a
  real(sgl) :: pa(7)                     ! Brosa parameter
  real(sgl) :: pd(7)                     ! Brosa parameter
  real(sgl) :: pe(7)                     ! Brosa parameter
  real(sgl) :: psh(7)                    ! Brosa parameter
  real(sgl) :: r0                        ! effective radius parameter
  real(sgl) :: rayl                      ! Brosa constant
  real(sgl) :: rhodi                     ! function that returns the shape of the dinuclear system
  real(sgl) :: rn1(4000, numelem)        ! number of evaporated neutrons from light fragment
  real(sgl) :: rn2(4000, numelem)        ! number of evaporated neutrons from heavy fragment
  real(sgl) :: rnma                      ! help variable
  real(sgl) :: rnmi                      ! help variable
  real(sgl) :: rtbis                     ! function to search for zero crossings of the function
  real(sgl) :: s12                       ! help variable
  real(sgl) :: sform                     ! function for Form factor for the Coulomb interaction energy  between two sphero
  real(sgl) :: sumtmp(4000)              ! help variable
  real(sgl) :: sumw                      ! sum over weights
  real(sgl) :: tmp                       ! temperature
  real(sgl) :: ve1                       ! potential
  real(sgl) :: ve2                       ! potential
  real(sgl) :: vnel                      ! help variable
  real(sgl) :: vr1                       ! function for volume of the projectile-like section
  real(sgl) :: vr2                       ! function for volume of the neck
  real(sgl) :: vr3                       ! help value
  real(sgl) :: wgt(4000)                 ! mass probability distribution
  real(sgl) :: wlog(4000)                ! help variable
  real(sgl) :: x1                        ! coordinates of intersection points inside the bin
  real(sgl) :: x2                        ! coordinates of the 2nd summit of the triangle
  real(sgl) :: xnu                       ! power
  real(sgl) :: zda                       ! charge over mass ratio
  real(sgl) :: zdis(4000, numelem)       ! charge distribution
  real(sgl) :: ze1                       ! help variable
  real(sgl) :: ze2                       ! help variable
  real(sgl) :: zf1(4000, numelem)        ! help variable
  real(sgl) :: zf2(4000, numelem)        ! help variable
  real(sgl) :: zo                        ! help variable
  real(sgl) :: zriss                     ! help variable
  real(sgl) :: zstepsize                 ! help variable
  real(sgl) :: ztot                      ! help variable
  real(sgl) :: zu                        ! help variable
  external fidi, rpoint, evap
!
! Determine mass and charge grid (depending whether evaporation correction is required)
!
  r0 = 1.15
  xnu = 1.0
  rayl = 11.00
  Zix = Zinit - Z
  Nix = Ninit - (A - Z)
  atot = real(A)
  ztot = real(Z)
  if(flagffevap)then
    astepnum = 10.
    astepsize = 0.1
    izstepnum = 30
    zstepsize = 0.1
  else
    astepnum = 1.
    astepsize = 1.
    izstepnum = 3
    zstepsize = 1.
  endif
!
  totl = (1.15 / 1.2249) * 2 * elt
  rayl =  11 * elt / ( 2.4 * 1.2249 * ((atot) **.333333) )
  di = 3. / (4. * pi * r0 **3)
  vtot = atot / di
  zda = ztot / atot
  gam = .9517 * (1. - 1.7826 * ((atot - 2. * ztot) / atot) **2)
  at = atot - ap
  rt = r0 * at **(1. / 3.)
  rp = r0 * ap **(1. / 3.)
!
! fission option
!
  call bdef(atot, ztot, 0., dum, dumm, bcom)
  d = totl - rt - rp
  r2 = totl / rayl
  c0 = crel * 2. * (rt + rp + 2. * (SQRT((rt - r2) * (rp - r2)) - r2)) / d **2
  rpt = (rp / rt) **xnu
  ald = ignatyuk(Zix, Nix, edefo, 0)
  tmp = sqrt(edefo / ald)
!
! Starting values for fmin to compute the shape of the dinuclear system psh(1) for a, psh(2) for z2, psh(3,4) for
! z1,z3, psh(5,6) for r1,r3.
! psh(7) curvature at the smallest cross section of the neck.
!
  psh(1) = .5
  pa(1) = .0
  pe(1) = 1.
  pd(1) = .1
  psh(2) = 0.
  pa(2) = - .2
  pe(2) = .2
  pd(2) = .05
  do i = 3, 4
    psh(i) = .65
    pa(i) = .3
    pe(i) = 1.
    pd(i) = .35
  enddo
  do i = 5, 6
    psh(i) = .8
    pa(i) = .6
    pe(i) = 1.
    pd(i) = .2
  enddo
  psh(7) = 1.
  pa(7) = .5
  pe(7) = 10.
  pd(7) = .5
  r1 = rp
  r3 = rt
  z1 = r1 * psh(3)
  z3 = r3 * psh(4)
!
! fmin needs a few starts to find the right values.
!
  do k = 1, 15
    if (k.ge.2) then
      delt = .2 **((k + 1) / 2)
      do i = 1, 7
        pa(i) = psh(i) - delt
        pe(i) = psh(i) + delt
        pd(i) = delt
      enddo
    endif
    fimin = fmin(fidi, 7, psh, pd, pa, pe, 200, - 1.E-4)
    if (fimin < 1.E-3) exit
  enddo
  d = totl - r1 - r3
!
! Graphical discussion of the rupture shape.
!
  amin = di * vr1(z1)
  imax = int(di * vr2(z1, z2, z3) * astepnum)
!
! In this loop the properties of the different fragmentations are calculated, as tke(a,z), neutron number rn(a,z), and yield wgt(a),
! zdis(a,z).
!
  jmx = 0
  do i = 0, imax
    aloop = amin + i * astepsize
!
! calculate the rupture cut at zriss depending on the mass number.
!
    rest = amin - aloop
    zu = z1
    zo = z3
    zriss = rtbis(rpoint, zu, zo, 2.E-4)
!
! Calculation of the equivalent ellipsoidal shapes and the Coulomb and nuclear proximity repulsion energies
!
    a1 = .5 * (zriss + r1)
    a2 = .5 * (d + r3 - zriss)
    ve1 = vr1(z1) + vr2(z1, z2, zriss)
    ve2 = vr3(z3) + vr2(zriss, z2, z3)
    b1 = sqrt(3. * ve1 / (4. * pi * a1))
    b2 = sqrt(3. * ve2 / (4. * pi * a2))
    de = a1 + a2
    x1 = a1 **2 - b1 **2
    if (x1 >= 0.)then
      x1 = sqrt(x1) / de
    else
      x1 = - sqrt( - x1) / de
    endif
    x2 = a2 **2 - b2 **2
    if(x2 >= 0.)then
      x2 = sqrt(x2) / de
    else
      x2 = - sqrt( - x2) / de
    endif
    s12 = sform(x1, x2)
    coul12 = 1.44 * s12 / de
    x1 = x1 * de / a1
    x2 = x2 * de / a2
    am1 = aloop
    am2 = atot - am1
    jmx = jmx + 1
    mcount = 1
    do izloop = 1, 2 * izstepnum + 1
      ze1 = zda * am1 + (izloop - izstepnum - 1) * zstepsize
      ze2 = ztot - ze1
      if(ze1 < 0. .or. ze2 < 0.) cycle
!
! Excess internal energy of ruptures nucleus es, es1 and es2 excitation energies of separated fragments, rn1 and rn2 number of
! evaporated neutrons from light and heavy fragments.
!
      call bdef(am1, ze1, x1, dum, dumm, bind1)
      call bdef(am2, ze2, x2, dum, dumm, bind2)
      es = edefo
      call bdef(am1, ze1, 0., dum, dumm, bind01)
      call bdef(am2, ze2, 0., dum, dumm, bind02)
      es1 = es / atot * am1
      es2 = es - es1
      es1 = bind1 - bind01 + es1
      es2 = bind2 - bind02 + es2
      wlog(jmx) = es1 + es2
      amm = am1
      zee = ze1
      ess = es1
      rnmi = - 1.
      rnma = wlog(jmx) * .167
      rn1(jmx, izloop) = 0.
      if(rnma > 0.7) rn1(jmx, izloop) = rtbis(evap, rnmi, rnma, 1.E-2)
      if(rn1(jmx, izloop) < 0.)rn1(jmx, izloop) = 0.
      amm = am2
      zee = ze2
      ess = es2
      rnmi = - 1.
      rnma = wlog(jmx) * .167
      rn2(jmx, izloop) = 0.
      if(rnma > 0.7) rn2(jmx, izloop) = rtbis(evap, rnmi, rnma, 1.E-2)
      if(rn2(jmx, izloop) < 0.)rn2(jmx, izloop) = 0.
!
! charge distribution
!
      zdis(jmx, izloop) = 0.
      vnel = 4. * pi * gam * (b1 * b2) **2 / (a1 * b2 * b2 + a2 * b1 * b1) * ( - 1.7817)
      coulel = coul12 * ze1 * ze2
      ezdis = (bind1 + bind2 + coulel + vnel) / tmp
      if(mcount == 1)ezdisnorm = ezdis
      mcount = mcount + 1
      ezdis = ezdis - ezdisnorm
      if (abs(ezdis) <= 80.) zdis(jmx, izloop) = exp( - ezdis)
!
! total kinetic energy tke, mass probability distribution wgt(a).
!
      eob = 2. * pi * gam * (rhodi(zriss) **2 - rhodi(z2) **2)
      wgt(jmx) = 0.
      expo = eob / tmp
      if(expo < 80.) wgt(jmx) = exp( - expo)
      af1(jmx) = am1
      zf1(jmx, izloop) = ze1
      af2(jmx) = am2
      zf2(jmx, izloop) = ze2
    enddo
  enddo
  jimax = jmx
!
! mass and charge distributions, with(out) corrections for evaporated neutrons
!
  sumw = 0.
  izmax = izstepnum * 2 + 1
  do k = 1, jimax
    sumw = sumw + wgt(k)
    sumtmp(k) = 0.
    do i = 1, izmax
      sumtmp(k) = sumtmp(k) + zdis(k, i)
    enddo
  enddo
  do k = 1, jimax
    wgt(k) = wgt(k) / sumw
    do i = 1, izmax
      zdis(k, i) = zdis(k, i) / sumtmp(k)
      iaf1 = max(int(af1(k) + 0.5), 1)
      izf1 = max(int(zf1(k, i) + 0.5), 1)
      iaf2 = max(int(af2(k) + 0.5), 1)
      izf2 = max(int(zf2(k, i) + 0.5), 1)
      fmz(iaf1, izf1) = wgt(k) * zdis(k, i) + fmz(iaf1, izf1)
      fmz(iaf2, izf2) = wgt(k) * zdis(k, i) + fmz(iaf2, izf2)
!
! calculate corrections for neutron evaporation if evapcor=true
!
      if (flagffevap) then
        irn1 = max(int(af1(k) - rn1(k, i) + 0.5), 1)
        irn2 = max(int(af2(k) - rn2(k, i) + 0.5), 1)
        fmzcor(irn1, izf1) = fmzcor(irn1, izf1) + wgt(k) * zdis(k, i)
        fmzcor(irn2, izf2) = fmzcor(irn2, izf2) + wgt(k) * zdis(k, i)
      endif
    enddo
    fmass(int(af1(k) + 0.5)) = fmass(int(af1(k) + 0.5)) + wgt(k)
    fmass(int(af2(k) + 0.5)) = fmass(int(af2(k) + 0.5)) + wgt(k)
  enddo
  do k = 1, nummass
    do i = 1, numelem
      fmasscor(k) = fmasscor(k) + fmzcor(k, i)
    enddo
  enddo
  return
end subroutine neck
! Copyright A.J. Koning 2021
