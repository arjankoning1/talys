subroutine opticalnp(Zix, Nix, k, eopt)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Optical potential for neutrons and protons
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
!   flagsoukhoinp ! flag for Soukhovitskii OMP for actinides
! Variables for OMP
!   Ejoin          ! joining energy for high energy OMP
!   flagoutkd      ! flag for output of KD03 OMP parameters
!   flagsoukho     ! flag for Soukhovitskii OMP for actinides
!   ompadjustD     ! depth of local OMP adjustment
!   ompadjustE1    ! start energy of local OMP adjustment
!   ompadjustE2    ! end energy of local OMP adjustment
!   ompadjustF     ! logical for local OMP adjustment
!   ompadjustN     ! number of energy ranges for local OMP adjustment
!   ompadjusts     ! variance of local OMP adjustment
!   optmod         ! file with optical model parameters
!   Vinfadjust     ! adj. factor for high energy limit of real centra
! Variables for nuclides
!   AA             ! mass number of residual nucleus
!   ZZ             ! charge number of residual nucleus
! Constants
!   onethird       ! 1 / 3
!   parname        ! name of particle
! Variables for OMP
!   av             ! real volume diffuseness
!   av0            ! diffuseness for real volume OMP
!   avd            ! real surface diffuseness
!   avd0           ! diffuseness for surface OMP
!   avso0          ! diffuseness for real spin - orbit OMP
!   avso           ! real spin - orbit diffuseness
!   aw             ! imaginary volume diffuseness
!   awd            ! imaginary surface diffuseness
!   awso           ! imaginary spin - orbit diffuseness
!   d1             ! parameter for imaginary surface OMP
!   d2             ! parameter for imaginary surface OMP
!   d3             ! parameter for imaginary surface OMP
!   ef             ! Fermi energy
!   eomp           ! energies on optical model file
!   Fav            ! adjustable factor for OMP (default 1.)
!   Favd           ! adjustable factor for OMP (default 1.)
!   Favso          ! adjustable factor for OMP (default 1.)
!   Faw            ! adjustable factor for OMP (default 1.)
!   Fawd           ! adjustable factor for OMP (default 1.)
!   Fawso          ! adjustable factor for OMP (default 1.)
!   Fd1            ! adjustable factor for OMP (default 1.)
!   Fd2            ! adjustable factor for OMP (default 1.)
!   Fd3            ! adjustable factor for OMP (default 1.)
!   Frc            ! adjustable factor for OMP (default 1.)
!   Frv            ! adjustable factor for OMP (default 1.)
!   Frvd           ! adjustable factor for OMP (default 1.)
!   Frvso          ! adjustable factor for OMP (default 1.)
!   Frw            ! adjustable factor for OMP (default 1.)
!   Frwd           ! adjustable factor for OMP (default 1.)
!   Frwso          ! adjustable factor for OMP (default 1.)
!   Fv1            ! adjustable factor for OMP (default 1.)
!   Fv2            ! adjustable factor for OMP (default 1.)
!   Fv3            ! adjustable factor for OMP (default 1.)
!   Fv4            ! adjustable factor for OMP (default 1.)
!   Fvso1          ! adjustable factor for OMP (default 1.)
!   Fvso2          ! adjustable factor for OMP (default 1.)
!   Fw1            ! adjustable factor for OMP (default 1.)
!   Fw2            ! adjustable factor for OMP (default 1.)
!   Fw3            ! adjustable factor for OMP (default 1.)
!   Fw4            ! adjustable factor for OMP (default 1.)
!   Fwso1          ! adjustable factor for OMP (default 1.)
!   Fwso2          ! adjustable factor for OMP (default 1.)
!   ompglobal      ! flag for use of global optical model
!   omplines       ! number of lines in optical model file
!   rc             ! Coulomb radius
!   rc0            ! Coulomb radius
!   rv             ! real volume radius
!   rv0            ! radius for real volume OMP
!   rvd0           ! radius for surface OMP
!   rvd            ! real surface radius
!   rvso           ! real spin - orbit radius
!   rvso0          ! radius for real spin - orbit OMP
!   rw             ! imaginary volume radius
!   rwd            ! imaginary surface radius
!   rwso           ! imaginary spin - orbit radius
!   v              ! real volume depth
!   vd             ! real surface depth
!   vso            ! real spin - orbit depth
!   v1             ! parameter for real volume OMP
!   v2             ! parameter for real volume OMP
!   v3             ! parameter for real volume OMP
!   V0             ! V at zero MeV
!   Vjoin          ! V at joining energy
!   vomp           ! optical model parameters from file
!   vso1           ! parameter for real spin - orbit OMP
!   vso2           ! parameter for real spin - orbit OMP
!   w              ! imaginary volume depth
!   w1             ! parameter for imaginary volume OMP
!   w2             ! parameter for imaginary volume OMP
!   w3             ! parameter for imaginary volume OMP
!   w4             ! parameter for imaginary volume OMP
!   wd             ! imaginary surface depth
!   Wjoin          ! W at joining energy
!   wso            ! imaginary spin - orbit depth
!   wso1           ! parameter for imaginary spin - orbit OMP
!   wso2           ! parameter for imaginary spin - orbit OMP
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: optmodfile       ! file with optical model parameters
  integer           :: A                ! mass number of target nucleus
  integer           :: i                ! counter
  integer           :: k                ! designator for particle
  integer           :: ktype            ! designator for particle
  integer           :: md               ! powers for W and Wd
  integer           :: mw               ! powers for W and Wd
  integer           :: nen              ! energy counter
  integer           :: Nix              ! neutron number index for residual nucleus
  integer           :: nr               ! number of radial grid point
  integer           :: nrange           ! number of energy ranges for local adjustment
  integer           :: omptype          ! type of optical model (spherical or coupled)
  integer           :: Z                ! charge number of target nucleus
  integer           :: Zix              ! charge number index for residual nucleus
  real(sgl)         :: d1loc            ! help variable
  real(sgl)         :: d2loc            ! help variable
  real(sgl)         :: d3loc            ! help variable
  real(sgl)         :: Dr(numrange)     ! depth of local adjustment
  real(sgl)         :: eint             ! help variable
  real(sgl)         :: elow             ! help variable
  real(sgl)         :: en1(numrange)    ! start energy of local adjustment
  real(sgl)         :: en2(numrange)    ! end energy of local adjustment
  real(sgl)         :: eopt             ! incident energy
  real(sgl)         :: eup              ! help variable
  real(sgl)         :: f                ! E-Ef
  real(sgl)         :: factor           ! multiplication factor
  real(sgl)         :: fjoin            ! help variable
  real(sgl)         :: sr(numrange)     ! variance of local adjustment
  real(sgl)         :: V0term           ! high energy V term
  real(sgl)         :: v1loc            ! help variable
  real(sgl)         :: v2loc            ! help variable
  real(sgl)         :: v3loc            ! help variable
  real(sgl)         :: v4loc            ! help variable
  real(sgl)         :: Vc               ! Coulomb term
  real(sgl)         :: Vcoul            ! Coulomb term
  real(sgl)         :: Vinf             ! high energy limit for real central potential
  real(sgl)         :: vloc(19)         ! interpolated optical model parameters
  real(sgl)         :: vso1loc          ! help variable
  real(sgl)         :: vso2loc          ! help variable
  real(sgl)         :: Vterm            ! high energy V term
  real(sgl)         :: w1loc            ! help variable
  real(sgl)         :: w2loc            ! help variable
  real(sgl)         :: w3loc            ! help variable
  real(sgl)         :: w4loc            ! help variable
  real(sgl)         :: wso1loc          ! help variable
  real(sgl)         :: wso2loc          ! help variable
!
! ************************ Calculate parameters ************************
!
! ompadjust    : subroutine for local optical model parameter adjustment
!
! 1. In case of an optical model file, we interpolate between the tabulated values.
!
  if (flagompejec .and. k == k0) then
    ktype = 0
  else
    ktype = k
  endif
  call ompadjust(eopt, ktype)
  optmodfile = '                                                     '
  if (Zix <= numZph .and. Nix <= numNph) optmodfile = optmod(Zix, Nix, k)
  if (optmodfile(1:1) /= ' ' .or. omplines(Zix, Nix, k) > 0) then
    if (eopt >= eomp(Zix, Nix, k, 1) .and. eopt <= eomp(Zix, Nix, k, omplines(Zix, Nix, k))) then
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
! 2. The general energy-dependent form of the optical potential using parameters per nucleus or the global optical model,
!    both from subroutine omppar.
!
! soukhovitskii: subroutine for global optical model parameters for actinides by Soukhovitskii et al.
!
  Z = ZZ(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  if (flagsoukhoinp .or. (flagsoukho .and. Z >= 90 .and. ompglobal(Zix, Nix, k))) then
    call soukhovitskii(k, Z, A, eopt)
  else
    f = eopt - ef(Zix, Nix, k)
    rc = Frc * rc0(Zix, Nix, k)
    v1loc = Fv1 * v1(Zix, Nix, k)
    v2loc = Fv2 * v2(Zix, Nix, k)
    v3loc = Fv3 * v3(Zix, Nix, k)
    v4loc = Fv4 * v4(Zix, Nix, k)
!
! Coulomb term for protons
!
    if (k == 2 .and. ompglobal(Zix, Nix, k)) then
      Vc = 1.73 / rc * Z / (A **onethird)
      vcoul = Vc * v1loc * (v2loc - 2. * v3loc * f + 3. * v4loc * f **2)
    else
      vcoul = 0.
    endif
!
! Extension up to 1 GeV.
!
    if (eopt > Ejoin(k)) fjoin = Ejoin(k) - ef(Zix, Nix, k)
    if (eopt <= Ejoin(k)) then
      v = v1loc * (1. - v2loc * f + v3loc * f **2 - v4loc * f **3) + vcoul
    else
      Vinf = - Vinfadjust(k) * 30.
      v = Vinf
      V0term = V0(k) - Vinf
      if (V0term > 0.) then
        Vterm = (Vjoin(k) - Vinf) / V0term
        if (Vterm > 0.) v = Vinf + V0term * exp(f / fjoin * log(Vterm))
!test     if (k.eq.2.and.ompglobal(Zix,Nix,k))
!    +      vcoul=-Vc*(V0(k)-Vinf)*logterm/fjoin*exp(f/fjoin*logterm)
      endif
    endif
    rv = Frv * rv0(Zix, Nix, k)
    av = Fav * av0(Zix, Nix, k)
    mw = 2
    w1loc = Fw1 * w1(Zix, Nix, k)
    w2loc = Fw2 * w2(Zix, Nix, k)
    w3loc = Fw3 * w3(Zix, Nix, k)
    w4loc = Fw4 * w4(Zix, Nix, k)
    if (eopt <= Ejoin(k)) then
      w = w1loc * f **mw / (f **mw + w2loc **mw)
    else
      w = Wjoin(k) - w3loc * fjoin **4 / (fjoin **4 + w4loc **4) + w3loc * f **4 / (f **4 + w4loc **4)
    endif
    rw = Frw * rv0(Zix, Nix, k)
    aw = Faw * av0(Zix, Nix, k)
    vd = 0.
    rvd = Frvd * rvd0(Zix, Nix, k)
    avd = Favd * avd0(Zix, Nix, k)
    md = 2
    d1loc = Fd1 * d1(Zix, Nix, k)
    d2loc = Fd2 * d2(Zix, Nix, k)
    d3loc = Fd3 * d3(Zix, Nix, k)
    wd = d1loc * f **md * exp( - d2loc * f) / (f **md + d3loc **md)
    rwd = Frwd * rvd0(Zix, Nix, k)
    awd = Fawd * avd0(Zix, Nix, k)
    vso1loc = Fvso1 * vso1(Zix, Nix, k)
    vso2loc = Fvso2 * vso2(Zix, Nix, k)
    vso = vso1loc * exp( - vso2loc * f)
    rvso = Frvso * rvso0(Zix, Nix, k)
    avso = Favso * avso0(Zix, Nix, k)
    wso1loc = Fwso1 * wso1(Zix, Nix, k)
    wso2loc = Fwso2 * wso2(Zix, Nix, k)
    wso = wso1loc * f **2 / (f **2 + wso2loc **2)
    rwso = Frwso * rvso0(Zix, Nix, k)
    awso = Fawso * avso0(Zix, Nix, k)
    if (flagoutkd) then
      write(*, '(" KD03 OMP parameters for ", a8, " E:", f12.5, " Ef:", f12.5)') parname(k), eopt, ef(Zix, Nix, k)
      write(*, '("   rv:", f12.5, "   av:", f12.5, "   v1:", f12.5, "   v2:", f12.5, "   v3:", es12.5, "   v4:", es12.5, &
 &      " Vcoul:", f12.5)') rv,av,v1loc,v2loc,v3loc,v4loc,Vcoul
      write(*, '("   rw:", f12.5, "   aw:", f12.5, "   w1:", f12.5, "   w2:", f12.5, "   w3:", f12.5, "   w4:", f12.5)') &
 &      rw, aw, w1loc, w2loc, w3loc, w4loc
      write(*, '("  rwd:", f12.5, "  awd:", f12.5, "   d1:", f12.5, "   d2:", f12.5, "   d3:", f12.5)') &
 &      rwd, awd, d1loc, d2loc, d3loc
      write(*, '(" rvso:", f12.5, " avso:", f12.5, " vso1:", f12.5, " vso2:", f12.5)') rvso, avso, vso1loc, vso2loc
      write(*, '(" rwso:", f12.5, " awso:", f12.5, " wso1:", f12.5, " wso2:", f12.5)') rwso, awso, wso1loc, wso2loc
    endif
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
end subroutine opticalnp
! Copyright A.J. Koning 2021
