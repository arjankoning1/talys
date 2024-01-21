subroutine opticalalpha(Zix, Nix, eopt)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Other optical potential for alphas
!
! Author    : Arjan Koning, Vlad and Marilena Avrigeanu
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl         ! single precision kind
! Variables for OMP
!   alphaomp    ! alpha optical model
! Variables for nuclides
!   AA          ! mass number of residual nucleus
!   ZZ          ! charge number of residual nucleus
! Constants
!   onethird    ! 1 / 3
! Variables for optical model
!   av          ! real volume diffuseness
!   avd         ! real surface diffuseness
!   avso        ! real spin - orbit diffuseness
!   aw          ! imaginary volume diffuseness
!   awd         ! imaginary surface diffuseness
!   awso        ! imaginary spin - orbit diffuseness
!   Fav         ! adjustable factor for OMP (default 1.)
!   Favd        ! adjustable factor for OMP (default 1.)
!   Favso       ! adjustable factor for OMP (default 1.)
!   Faw         ! adjustable factor for OMP (default 1.)
!   Fawd        ! adjustable factor for OMP (default 1.)
!   Fawso       ! adjustable factor for OMP (default 1.)
!   Fd1         ! adjustable factor for OMP (default 1.)
!   Frc         ! adjustable factor for OMP (default 1.)
!   Frv         ! adjustable factor for OMP (default 1.)
!   Frvd        ! adjustable factor for OMP (default 1.)
!   Frvso       ! adjustable factor for OMP (default 1.)
!   Frw         ! adjustable factor for OMP (default 1.)
!   Frwd        ! adjustable factor for OMP (default 1.)
!   Frwso       ! adjustable factor for OMP (default 1.)
!   Fv1         ! adjustable factor for OMP (default 1.)
!   Fvso1       ! adjustable factor for OMP (default 1.)
!   Fw1         ! adjustable factor for OMP (default 1.)
!   Fwso1       ! adjustable factor for OMP (default 1.)
!   rc          ! Coulomb radius
!   rv          ! real volume radius
!   rvd         ! real surface radius
!   rvso        ! real spin - orbit radius
!   rw          ! imaginary volume radius
!   rwd         ! imaginary surface radius
!   rwso        ! imaginary spin - orbit radius
!   v           ! real volume depth
!   vd          ! real surface depth
!   vso         ! real spin - orbit depth
!   w           ! imaginary volume depth
!   wd          ! imaginary surface depth
!   wso         ! imaginary spin - orbit depth
!
! *** Declaration of local data
!
  implicit none
  integer   :: A     ! mass number of target nucleus
  integer   :: Nix   ! neutron number index for residual nucleus
  integer   :: Z     ! charge number of target nucleus
  integer   :: Zix   ! charge number index for residual nucleus
  real(sgl) :: A13   ! A**1/3
  real(sgl) :: e1opt ! energy for alpha OMP
  real(sgl) :: e2opt ! energy for alpha OMP
  real(sgl) :: e3opt ! energy for alpha OMP
  real(sgl) :: e4opt ! energy for alpha OMP
  real(sgl) :: eopt  ! incident energy
  real(sgl) :: rb    ! maximum radius value
!
! Alternative options for alpha OMP
!
! S. Goriely: inclusion of the alpha OMP of Mc Fadden & Satchler for alphaomp=2, folding model for alphaomp=3,4,5
! alphaomp=6 --> Avrigeanu et al. potential [PRC90,044612(2014)]
!
! Overwrite some of the previous values.
!
  call ompadjust(eopt, 6)
  Z = ZZ(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  A13 = real(A) **onethird
!
! alphaomp=2 --> Global OMP of Mc Fadden & Satchler, Nucl. Phys. A 84, 177 (1966).
!
  if (alphaomp == 2) then
    v = 185.0
    rv = 1.40
    av = 0.52
    w = 25.0
    rw = rv
    aw = av
    vd = 0.
    wd = 0.
    vso = 0.
    wso = 0.
    rc = 1.3
  endif
!
! alphaomp=6 --> Global OMP of Avrigeanu et al., Phys. Rev. C 90, 044612 (2014).
!
  if (alphaomp == 6) then
!
! Norenberg, HIC, N - H, 1980, p.8
!
    rb = 2.66 + 1.36 * A13
    e2opt = (2.59 + 10.4 / A) * Z / rb
    e1opt = - 3.03 - 0.76 * A13 + 1.24 * e2opt
    e3opt = 22.2 + 0.181 * Z / A13
    e4opt = 29.1 - 0.22 * Z / A13
    if (eopt <= e3opt) then
      v = 165.0 + 0.733 * Z / A13 - 2.64 * eopt
    else
      v = 116.5 + 0.337 * Z / A13 - 0.453 * eopt
    endif
    v = max(v, - 100.)
    if (eopt <= 25.0) then
      rv = 1.18 + 0.012 * eopt
    else
      rv = 1.48
    endif
    if (eopt <= e2opt) then
      av = 0.631 + (0.016 - 0.001 * e2opt) * Z / A13
    elseif (eopt <= e4opt) then
      av = 0.631 + (0.016 - 0.001 * eopt) * Z / A13
    else
      av = 0.684 - 0.016 * Z / A13 - (0.0026 - 0.00026 * Z / A13) * eopt
    endif
    av = max(av, 0.1)
    w = amax1(2.73 - 2.88 * A13 + 1.11 * eopt, 0.0)
    rw = 1.34
    aw = 0.50
    vd = 0.
    if (eopt <= e1opt) then
      wd = 4.0
    elseif (eopt <= e2opt) then
      wd = 22.2 + 4.57 * A13 - 7.446 * e2opt + 6.0 * eopt
    else
      wd = 22.2 + 4.57 * A13 - 1.446 * eopt
    endif
    wd = max(wd, 0.)
    if (A <= 152 .or. A >= 190) then
      rwd = 1.52
    else
      rwd = amax1(1.74 - 0.01 * eopt, 1.52)
    endif
    awd = 0.729 - 0.074 * A13
    vso = 0.
    wso = 0.
    rc = 1.3
  endif
!
! alphaomp=7 --> Global OMP of Nolte et al., Phys. Rev. C 36, 1312 (1987).
!
  if (alphaomp == 7) then
    v = 101.1 + 6.051 * Z / A13 - 0.248 * eopt
    rv = 1.245
    av = 0.817 - 0.0085 * A13
    w = 26.82 - 1.706 * A13 + 0.006 * eopt
    rw = 1.57
    aw = 0.692 - 0.02 * A13
    vd = 0.
    wd = 0.
    vso = 0.
    wso = 0.
    rc = 1.3
  endif
!
! alphaomp=8 --> Global OMP of Avrigeanu et al., Phys. Rev. C 49, 2136 (1994).
!
  if (alphaomp == 8) then
    v = 101.1 + 6.051 * Z / A13 - 0.248 * eopt
    rv = 1.245
    av = 0.817 - 0.0085 * A13
    if (eopt <= 73.) then
      w = 12.64 - 1.706 * A13 + 0.20 * eopt
    else
      w = 26.82 - 1.706 * A13 + 0.006 * eopt
    endif
    rw = 1.57
    aw = 0.692 - 0.02 * A13
    vd = 0.
    wd = 0.
    vso = 0.
    wso = 0.
    rc = 1.3
  endif
!
! Possible adjustment of parameters
!
  v = Fv1 * v
  rv = Frv * rv
  av = Fav * av
  w = Fw1 * w
  rw = Frw * rw
  aw = Faw * aw
  vd = Fd1 * vd
  rvd = Frvd * rvd
  avd = Favd * avd
  wd = Fd1 * wd
  rwd = Frwd * rwd
  awd = Fawd * awd
  vso = Fvso1 * vso
  rvso = Frvso * rvso
  avso = Favso * avso
  wso = Fwso1 * wso
  rwso = Frwso * rwso
  awso = Fawso * awso
  rc = Frc * rc
  return
end subroutine opticalalpha
! Copyright A.J. Koning 2021
