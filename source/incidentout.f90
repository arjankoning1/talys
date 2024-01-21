subroutine incidentout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Reaction output for incident channel
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
! Variables for OMP
!   Rprime       ! potential scattering radius
!   Sstrength    ! s, p, d, etc - wave strength function
! Variables for basic reaction
!   flagang      ! flag for output of angular distributions
! Variables for numerics
!   nangle       ! number of angles
! Variables for main input
!   Atarget      ! mass number of target nucleus
!   k0           ! index of incident particle
!   Ltarget      ! excited level of target
! Variables for gamma rays
!   fiso         ! correction factor for isospin forbidden transition
!   gammax       ! number of l - values for gamma multipolarity
! Variables for energy grid
!   angle        ! angle in degrees
!   Einc         ! incident energy in MeV
! Variables for incident channel
!   directad     ! direct angular distribution
!   lmaxinc      ! maximal l - value for transm. coeff. for incident channel
!   Tjlinc       ! transm. coeff. as a function of spin and l for inc. channel
!   Tlinc        ! transm. coeff. as a function of l for incident channel
!   xselasinc    ! total elastic cross section (neutrons only) for inc. channel
!   xsreacinc    ! reaction cross section for incident channel
!   xstotinc     ! total cross section (neutrons only) for incident channel
! Constants
!   parname      ! name of particle
!
! *** Declaration of local data
!
  implicit none
  integer :: iang              ! running variable for angle
  integer :: l                 ! multipolarity
  integer :: type              ! particle type
!
! *********** Total cross sections for incident channel ****************
!
  write(*, '(/" Optical model results"/)')
  if (k0 == 1) then
    write(*, '(" Total cross section   :", es11.4, " mb")') xstotinc
  endif
  write(*, '(" Reaction cross section:", es11.4, " mb")') xsreacinc
  if (k0 == 1) then
    write(*, '(" Elastic cross section :", es11.4, " mb")') xselasinc
  endif
!
! For low energy neutrons we give the resonance parameters.
!
! Sstrength: s,p,d,etc-wave strength function
!
  if (k0 == 1 .and. Einc <= 0.1) then
    write(*, '(/" S-wave and P-wave strength functions and potential scattering radius"/)')
    write(*, '("      A      Value"/)')
    write(*, '(" S0:", i4, f8.4, " .e-4")') Atarget, Sstrength(0)*1.e4
    write(*, '(" S1:", i4, f8.4, " .e-4")') Atarget, Sstrength(1)*1.e4
    write(*, '(" R :", i4, f8.4, " fm")') Atarget, Rprime
  endif
  write(*, '(/" Isospin factors to reduce emission "/)')
  do type = 0, 6
    write(*, '(1x, a8, 1x, f8.5)') parname(type), fiso(type)
  enddo
!
! *********** Transmission coefficients for incident channel ***********
!
  if (lmaxinc ==  - 1) return
  write(*, '(/" Transmission coefficients for incident ", a8, " at ", f8.3, " MeV"/)') parname(k0), Einc
!
! 1. Spin 1/2 particles: Neutrons, protons, tritons and Helium-3
!
  if (k0 == 1 .or. k0 == 2 .or. k0 == 4 .or. k0 == 5) then
    write(*, '("  L    T(L-1/2,L)   T(L+1/2,L)    Tav(L)"/)')
    do l = 0, lmaxinc
      write(*, '(1x, i3, 3es13.5)') l, Tjlinc(-1, l), Tjlinc(1, l), Tlinc(l)
    enddo
  endif
!
! 2. Spin 1 particles: Deuterons
!
  if (k0 == 3) then
    write(*, '("  L     T(L-1,L)      T(L,L)      T(L+1,L)     Tav(L)"/)')
    do l = 0, lmaxinc
      write(*, '(1x, i3, 4es13.5)') l, Tjlinc(-1, l), Tjlinc(0, l), Tjlinc(1, l), Tlinc(l)
    enddo
  endif
!
! 3. Spin 0 particles: Alpha-particles
!
  if (k0 == 6) then
    write(*, '("  L     T(L)"/)')
    do l = 0, lmaxinc
      write(*, '(1x, i3, es13.5)') l, Tjlinc(0, l)
    enddo
  endif
!
! 4. Photons
!
  if (k0 == 0) then
    write(*, '("  L     T(L)"/)')
    do l = 1, gammax
      write(*, '(1x, i3, es13.5)') l, Tjlinc(0, l)
    enddo
  endif
!
! *********** Shape elastic scattering angular distribution ************
!
  if (flagang .and. k0 > 0) then
    write(*, '(/" Shape elastic scattering angular distribution"/)')
    write(*, '(" Angle    Cross section"/)')
    do iang = 0, nangle
      write(*, '(1x, f5.1, es16.5)') angle(iang), directad(k0, Ltarget, iang)
    enddo
  endif
  return
end subroutine incidentout
! Copyright A.J. Koning 2021
