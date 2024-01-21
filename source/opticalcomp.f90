subroutine opticalcomp(Zix, Nix, k, eopt)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Optical potential for composite particles
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
!   sgl            ! single precision kind
! All global variables
!   numNph         ! maximum number of neutrons away from the initial compound nucleus
!   numrange       ! number of energy ranges
!   numZph         ! maximum number of protons away from the initial compound nucleus
! Variables for OMP
!   alphaomp       ! alpha optical model
!   altomp         ! flag for alternative optical model
!   deuteronomp    ! deuteron optical model
!   ompadjustD     ! depth of local OMP adjustment
!   ompadjustE1    ! start energy of local OMP adjustment
!   ompadjustE2    ! end energy of local OMP adjustment
!   ompadjustF     ! logical for local OMP adjustment
!   ompadjustN     ! number of energy ranges for local OMP adjustment
!   ompadjusts     ! variance of local OMP adjustment
!   optmod         ! file with optical model parameters
! Constants
!   parA           ! mass number of particle
!   parN           ! neutron number of particle
!   parZ           ! charge number of particle
! Variables for optical model
!   av             ! real volume diffuseness
!   avd            ! real surface diffuseness
!   avso           ! real spin - orbit diffuseness
!   aw             ! imaginary volume diffuseness
!   awd            ! imaginary surface diffuseness
!   awso           ! imaginary spin - orbit diffuseness
!   eomp           ! energies on optical model file
!   Eompbeg0       ! upper energy of KD03 OMP
!   Eompbeg1       ! lower energy of alternative OMP
!   Eompend0       ! lower energy of KD03 OMP
!   Eompend1       ! upper energy of alternative
!   Fav            ! adjustable factor for OMP (default 1.)
!   Favd           ! adjustable factor for OMP (default 1.)
!   Favso          ! adjustable factor for OMP (default 1.)
!   Faw            ! adjustable factor for OMP (default 1.)
!   Fawd           ! adjustable factor for OMP (default 1.)
!   Fawso          ! adjustable factor for OMP (default 1.)
!   Fd1            ! adjustable factor for OMP (default 1.)
!   Frc            ! adjustable factor for OMP (default 1.)
!   Frv            ! adjustable factor for OMP (default 1.)
!   Frvd           ! adjustable factor for OMP (default 1.)
!   Frvso          ! adjustable factor for OMP (default 1.)
!   Frw            ! adjustable factor for OMP (default 1.)
!   Frwd           ! adjustable factor for OMP (default 1.)
!   Frwso          ! adjustable factor for OMP (default 1.)
!   Fv1            ! adjustable factor for OMP (default 1.)
!   Fvso1          ! adjustable factor for OMP (default 1.)
!   Fw1            ! adjustable factor for OMP (default 1.)
!   Fwso1          ! adjustable factor for OMP (default 1.)
!   omplines       ! number of lines in optical model file
!   rc             ! Coulomb radius
!   rv             ! real volume radius
!   rvd            ! real surface radius
!   rvso           ! real spin - orbit radius
!   rw             ! imaginary volume radius
!   rwd            ! imaginary surface radius
!   rwso           ! imaginary spin - orbit radius
!   v              ! real volume depth
!   vd             ! real surface depth
!   vomp           ! optical model parameters from file
!   vso            ! real spin - orbit depth
!   w              ! imaginary volume depth
!   wd             ! imaginary surface depth
!   wso            ! imaginary spin - orbit depth
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: optmodfile       ! file with optical model parameters
  integer           :: i                ! counter
  integer           :: ia               ! mass number from abundance table
  integer           :: in               ! counter for neutrons
  integer           :: iz               ! charge number of residual nucleus
  integer           :: k                ! designator for particle
  integer           :: nen              ! energy counter
  integer           :: Nix              ! neutron number index for residual nucleus
  integer           :: nr               ! number of radial grid point
  integer           :: nrange           ! number of energy ranges for local adjustment
  integer           :: omptype          ! type of optical model (spherical or coupled)
  integer           :: Zix              ! charge number index for residual nucleus
  real(sgl)         :: avdkd            ! optical model parameter
  real(sgl)         :: avdn             ! optical model parameter
  real(sgl)         :: avdp             ! optical model parameter
  real(sgl)         :: avkd             ! optical model parameter
  real(sgl)         :: avn              ! optical model parameter
  real(sgl)         :: avp              ! optical model parameter
  real(sgl)         :: avsokd           ! optical model parameter
  real(sgl)         :: avson            ! optical model parameter
  real(sgl)         :: avsop            ! optical model parameter
  real(sgl)         :: awdkd            ! optical model parameter
  real(sgl)         :: awkd             ! optical model parameter
  real(sgl)         :: awn              ! optical model parameter
  real(sgl)         :: awp              ! optical model parameter
  real(sgl)         :: awsokd           ! optical model parameter
  real(sgl)         :: Dr(numrange)     ! depth of local adjustment
  real(sgl)         :: e                ! energy
  real(sgl)         :: Efrac            ! fractional energy
  real(sgl)         :: eint             ! help variable
  real(sgl)         :: elow             ! help variable
  real(sgl)         :: en1(numrange)    ! start energy of local adjustment
  real(sgl)         :: en2(numrange)    ! end energy of local adjustment
  real(sgl)         :: eopt             ! incident energy
  real(sgl)         :: eup              ! help variable
  real(sgl)         :: factor           ! multiplication factor
  real(sgl)         :: rckd             ! optical model parameter
  real(sgl)         :: rvdkd            ! optical model parameter
  real(sgl)         :: rvdn             ! optical model parameter
  real(sgl)         :: rvdp             ! optical model parameter
  real(sgl)         :: rvkd             ! optical model parameter
  real(sgl)         :: rvn              ! optical model parameter
  real(sgl)         :: rvp              ! optical model parameter
  real(sgl)         :: rvsokd           ! optical model parameter
  real(sgl)         :: rvson            ! optical model parameter
  real(sgl)         :: rvsop            ! optical model parameter
  real(sgl)         :: rwdkd            ! optical model parameter
  real(sgl)         :: rwkd             ! optical model parameter
  real(sgl)         :: rwn              ! optical model parameter
  real(sgl)         :: rwp              ! optical model parameter
  real(sgl)         :: rwsokd           ! optical model parameter
  real(sgl)         :: sr(numrange)     ! variance of local adjustment
  real(sgl)         :: vdkd             ! optical model parameter
  real(sgl)         :: vdn              ! optical model parameter
  real(sgl)         :: vdp              ! optical model parameter
  real(sgl)         :: vkd              ! optical model parameter
  real(sgl)         :: vloc(19)         ! interpolated optical model parameters
  real(sgl)         :: vn               ! optical model parameters for neutrons
  real(sgl)         :: vp               ! optical model parameters for protons
  real(sgl)         :: vsokd            ! optical model parameter
  real(sgl)         :: vson             ! Vso and Wso for neutrons
  real(sgl)         :: vsop             ! Vso and Wso for protons
  real(sgl)         :: wdkd             ! optical model parameter
  real(sgl)         :: wdn              ! optical model parameter
  real(sgl)         :: wdp              ! optical model parameter
  real(sgl)         :: wkd              ! optical model parameter
  real(sgl)         :: wn               ! optical model parameter
  real(sgl)         :: wp               ! optical model parameter
  real(sgl)         :: wsokd            ! optical model parameter
  real(sgl)         :: wson             ! Vso and Wso for neutrons
  real(sgl)         :: wsop             ! Vso and Wso for protons
!
! ******************* Optical model input file *************************
!
! 1. In case of an optical model file, we interpolate between the tabulated values
!
  optmodfile = '                                                     '
  if (Zix <= numZph .and. Nix <= numNph) optmodfile = optmod(Zix, Nix, k)
  if (optmodfile(1:1) /= ' ') then
    if (eopt >= eomp(Zix, Nix, k, 1) .and. eopt <= eomp(Zix, Nix, k, omplines(Zix, Nix, k))) then
      call ompadjust(eopt, k)
      do nen = 1, omplines(Zix, Nix, k) - 1
        elow = eomp(Zix, Nix, k, nen)
        eup = eomp(Zix, Nix, k, nen + 1)
        if (elow <= eopt .and. eopt <= eup) then
          eint = (eopt - elow) / (eup - elow)
          do i = 1, 19
            vloc(i) = vomp(Zix, Nix, k, nen, i) + eint * (vomp(Zix, Nix, k, nen + 1, i) - vomp(Zix, Nix, k, nen, i))
          enddo
          v = Fv1 * vloc(1)
          rv = Frv * vloc(2)
          av = Fav * vloc(3)
          w = Fw1 * vloc(4)
          rw = Frw * vloc(5)
          aw = Faw * vloc(6)
          vd = Fd1 * vloc(7)
          rvd = Frvd * vloc(8)
          avd = Favd * vloc(9)
          wd = Fd1 * vloc(10)
          rwd = Frwd * vloc(11)
          awd = Fawd * vloc(12)
          vso = Fvso1 * vloc(13)
          rvso = Frvso * vloc(14)
          avso = Favso * vloc(15)
          wso = Fwso1 * vloc(16)
          rwso = Frwso * vloc(17)
          awso = Fawso * vloc(18)
          rc = Frc * vloc(19)
          goto 200
        endif
      enddo
    endif
  endif
!
! 2. The general energy-dependent form of the optical potential using parameters per nucleus or a global optical model,
!    both from subroutine omppar.
!
! We use the Watanabe method (S. Watanabe, Nucl. Phys. 8, 484 (1958), see also D.G. Madland, Proceedings of a
! Specialists' Meeting on preequilibrium nuclear reactions, Semmering, Austria, February 10-12 1988, p. 103)
! to make a composite particle potential out of the proton and neutron potential. The simplified formula is:
!   V(d,E)=z*V(n,E/a)+n*V(p,E/a)
! with z, n and a, the proton, neutron and mass number of the particle, and similarly for Wd and W.
! VSO and Wso are as in the nucleon potentials.
! We take a similar weighting for the geometry parameters.
!
! 1. Neutron potential.
!
! opticalnp: subroutine for optical potential for neutrons and protons
!
  e = eopt / parA(k)
  call opticalnp(Zix, Nix, 1, e)
  vn = v
  rvn = rv
  avn = av
  wn = w
  rwn = rw
  awn = aw
  vdn = vd
  rvdn = rvd
  avdn = avd
  wdn = wd
!
! 2. Proton potential.
!
  call opticalnp(Zix, Nix, 2, e)
  vp = v
  rvp = rv
  avp = av
  wp = w
  rwp = rw
  awp = aw
  vdp = vd
  rvdp = rvd
  avdp = avd
  wdp = wd
!
! 3. Another 2 calls to opticalnp for eopt=E to determine Vso and Wso.
!
  call opticalnp(Zix, Nix, 1, eopt)
  vson = vso
  rvson = rvso
  avson = avso
  wson = wso
  call opticalnp(Zix, Nix, 2, eopt)
  vsop = vso
  rvsop = rvso
  avsop = avso
  wsop = wso
!
! 4. Final potential depths: construct V, W, Wd, Vso and Wso.
!
  call ompadjust(eopt, k)
  iz = parZ(k)
  in = parN(k)
  ia = parA(k)
  v = Fv1 * (in * vn + iz * vp)
  rv = Frv * (in * rvn + iz * rvp) / ia
  av = Fav * (in * avn + iz * avp) / ia
  w = Fw1 * (in * wn + iz * wp)
  rw = Frw * (in * rwn + iz * rwp) / ia
  aw = Faw * (in * awn + iz * awp) / ia
  vd = Fd1 * (in * vdn + iz * vdp)
  rvd = Frvd * (in * rvdn + iz * rvdp) / ia
  avd = Favd * (in * avdn + iz * avdp) / ia
  wd = Fd1 * (in * wdn + iz * wdp)
  rwd = Frwd * (in * rvdn + iz * rvdp) / ia
  awd = Fawd * (in * avdn + iz * avdp) / ia
  rvso = Frvso * (in * rvson + iz * rvsop) / ia
  avso = Favso * (in * avson + iz * avsop) / ia
  rwso = Frwso * (in * rvson + iz * rvsop) / ia
  awso = Fawso * (in * avson + iz * avsop) / ia
  rc = Frc * rc
  if (k == 3) then
    vso = Fvso1 * (vson + vsop) / 2.
    wso = Fwso1 * (wson + wsop) / 2.
  endif
  if (k == 4 .or. k == 5) then
    vso = Fvso1 * (vson + vsop) / 6.
    wso = Fwso1 * (wson + wsop) / 6.
  endif
  if (k == 6) then
    vso = 0.
    wso = 0.
  endif
!
! Alternative options for complex particle OMP
!
  if (altomp(k)) then
    vkd = v
    rvkd = rv
    avkd = av
    wkd = w
    rwkd = rw
    awkd = aw
    vdkd = vd
    rvdkd = rvd
    avdkd = avd
    wdkd = wd
    rwdkd = rwd
    awdkd = awd
    vsokd = vso
    rvsokd = rvso
    avsokd = avso
    wsokd = wso
    rwsokd = rwso
    awsokd = awso
    rckd = rc
    i = 1
    if (k == 3 .and. deuteronomp >= 2) then
      call opticaldeut(Zix, Nix, eopt)
      i = deuteronomp
    endif
    if (k == 6 .and. (alphaomp == 2 .or. alphaomp >= 6)) then
      call opticalalpha(Zix, Nix, eopt)
      i = alphaomp
    endif
    Efrac = 1.
    if (eopt <= Eompbeg0(k, i)) Efrac = 0.
    if (eopt > Eompbeg0(k, i) .and. eopt <= Eompbeg1(k, i)) Efrac = (eopt - Eompbeg0(k, i)) / (Eompbeg1(k, i) - Eompbeg0(k, i))
    if (eopt > Eompend1(k, i) .and. eopt <= Eompend0(k, i)) Efrac = 1. - (eopt - Eompend1(k, i)) / (Eompend0(k, i) - Eompend1(k, i))
    if (eopt > Eompend0(k, i)) Efrac = 0.
    v = Efrac * v + (1 - Efrac) * vkd
    rv = Efrac * rv + (1 - Efrac) * rvkd
    av = Efrac * av + (1 - Efrac) * avkd
    w = Efrac * w + (1 - Efrac) * wkd
    rw = Efrac * rw + (1 - Efrac) * rwkd
    aw = Efrac * aw + (1 - Efrac) * awkd
    vd = Efrac * vd + (1 - Efrac) * vdkd
    rvd = Efrac * rvd + (1 - Efrac) * rvdkd
    avd = Efrac * avd + (1 - Efrac) * avdkd
    wd = Efrac * wd + (1 - Efrac) * wdkd
    rwd = Efrac * rwd + (1 - Efrac) * rwdkd
    awd = Efrac * awd + (1 - Efrac) * awdkd
    vso = Efrac * vso + (1 - Efrac) * vsokd
    rvso = Efrac * rvso + (1 - Efrac) * rvsokd
    avso = Efrac * avso + (1 - Efrac) * avsokd
    wso = Efrac * wso + (1 - Efrac) * wsokd
    rwso = Efrac * rwso + (1 - Efrac) * rwsokd
    awso = Efrac * awso + (1 - Efrac) * awsokd
    rc = Efrac * rc + (1 - Efrac) * rckd
  endif
!
! Possible additional energy-dependent adjustment of the geometry
!
! adjustF          : subroutine for local parameter adjustment
!
  200 if (ompadjustF(k)) then
    do omptype = 1, 13
      nrange = ompadjustN(k, omptype)
      do nr = 1, nrange
        en1(nr) = ompadjustE1(k, omptype, nr)
        en2(nr) = ompadjustE2(k, omptype, nr)
        Dr(nr) = ompadjustD(k, omptype, nr)
        sr(nr) = ompadjusts(k, omptype, nr)
      enddo
      call adjustF(eopt, nrange, en1, en2, Dr, sr, factor)
      if (omptype == 1) rv = factor * rv
      if (omptype == 2) av = factor * av
      if (omptype == 3) rw = factor * rw
      if (omptype == 4) aw = factor * aw
      if (omptype == 5) rvd = factor * rvd
      if (omptype == 6) avd = factor * avd
      if (omptype == 7) rwd = factor * rwd
      if (omptype == 8) awd = factor * awd
      if (omptype == 9) rvso = factor * rvso
      if (omptype == 10) avso = factor * avso
      if (omptype == 11) rwso = factor * rwso
      if (omptype == 12) awso = factor * awso
      if (omptype == 13) rc = factor * rc
    enddo
  endif
  return
end subroutine opticalcomp
! Copyright A.J. Koning 2021
