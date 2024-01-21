subroutine mom(Zix, Nix, pZ, e)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Microscopic JLM OMP
!
! Author    : Eric Bauge and Arjan Koning
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
!   dbl           ! double precision kind
! All global variables
!   numjlm        ! maximum number of radial points
! Variables for OMP
!   jlmmode       ! option for JLM imaginary potential normalization
!   lv1adjust     ! adjustable parameter for JLM OMP
!   lvadjust      ! adjustable parameter for JLM OMP
!   lvsoadjust    ! adjustable parameter for JLM OMP
!   lw1adjust     ! adjustable parameter for JLM OMP
!   lwadjust      ! adjustable parameter for JLM OMP
!   lwsoadjust    ! adjustable parameter for JLM OMP
!   ompadjustp    ! flag for local optical model parameter adjustment
! Variables for nuclides
!   AA            ! mass number of residual nucleus
!   ZZ            ! charge number of residual nucleus
! Variables for optical model
!   rc            ! Coulomb radius
! Variables for JLM
!   normjlm       ! JLM potential normalization factors
!   potjlm        ! JLM potential depth values
!   rhojlmn       ! density for neutrons
!   rhojlmp       ! density for protons
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: key                   ! keyword
  integer           :: A                     ! mass number of target nucleus
  integer           :: i                     ! level
  integer           :: j                     ! counter
  integer           :: N                     ! neutron number of residual nucleus
  integer           :: Nix                   ! neutron number index for residual nucleus
  integer           :: Z                     ! charge number of target nucleus
  integer           :: Zix                   ! charge number index for residual nucleus
  real(dbl)         :: e                     ! energy
  real(dbl)         :: lv                    ! real potential depth normalization factor
  real(dbl)         :: lv1                   ! real isovector potential depth normalization factor
  real(dbl)         :: lvso                  ! real spin orbit potential depth normalization factor
  real(dbl)         :: lw                    ! imag potential depth normalization factor
  real(dbl)         :: lw1                   ! imag isovector potential depth normalization factor
  real(dbl)         :: lwso                  ! imaginary spin orbit potential depth normalization factor
  real(dbl)         :: pZ                    ! product of charges
  real(sgl)         :: alam                  ! adjustable parameter
  real(sgl)         :: factor                ! multiplication factor
  real(sgl)         :: rcjlm                 ! coulomb radius for JLM
  real(sgl)         :: rhomomn(numjlm, 6)    ! total neutron density at a given radius radmom
  real(sgl)         :: rhomomp(numjlm, 6)    ! total proton density at a given radius radmom
  real(sgl)         :: vpot(numjlm, 6)       ! optical potential
!
! *********************** Parameterization *****************************
!
! Parameters from prc 63, 024607 2001
! adjust    : subroutine for energy-dependent parameter adjustment
!
  A = AA(Zix, Nix, 0)
  Z = ZZ(Zix, Nix, 0)
  N = A - Z
  lv = .0008 * log(e * 1000.) + .00018 * (log(e * 1000.)) **2. + .951
  lw = (1.24 - (1.0 / (1. + exp(((e) - 4.5) / 2.9)))) * (1. + .06 * exp( - ((e - 14.) / 3.7) **2 )) &
       * (1. - 0.09 * exp( - ((e - 80.) / 78.) **2))
  if(e > 80.)lw = lw * (1 + (e - 80.) / (400.))
  lv1 = 1.5 - (0.65 / (1. + exp((e - 1.3) / 3.0)))
!
! Possible modification of JLMB potential with jlmmode
!  cf Goriely & Delaroche: PLB653, 158 (2007)
!
  alam = 0.44
  if (jlmmode == 1) alam = 1.10 * exp( - 0.4 * e **0.25)
  if (jlmmode >= 2) alam = 1.375 * exp( - 0.2 * e **0.5)
  if (jlmmode == 3) lw = lw * 2.
!
  lw1 = (1.1 + (alam / (1. + (exp(((e) - 40.) / 50.9)) **4 ))) * (1. - .065 * exp( - ((40. - e) / 13.) **2)) &
       * (1. - .083 * exp( - ((200. - e) / 80.) **2))
!
  lvso = 40. + exp( - e * 0.013) * 130.
  lwso = - 0.2 * (e - 20)
  if (ompadjustp(1)) then
    key = 'lvadjust'
    call adjust(real(e), key, 0, 0, 0, 0, factor)
    lv = factor * lvadjust * lv
    key = 'lwadjust'
    call adjust(real(e), key, 0, 0, 0, 0, factor)
    lw = factor * lwadjust * lw
    key = 'lv1adjust'
    call adjust(real(e), key, 0, 0, 0, 0, factor)
    lv1 = factor * lv1adjust * lv1
    key = 'lw1adjust'
    call adjust(real(e), key, 0, 0, 0, 0, factor)
    lw1 = factor * lw1adjust * lw1
    key = 'lvsoadjust'
    call adjust(real(e), key, 0, 0, 0, 0, factor)
    lvso = factor * lvsoadjust * lvso
    key = 'lwsoadjust'
    call adjust(real(e), key, 0, 0, 0, 0, factor)
    lwso = factor * lwsoadjust * lwso
  endif
!
! Calculate and write the potentials in ecis format
!
  do i = 1, numjlm
    do j = 1, 6
      rhomomn(i, j) = rhojlmn(Zix, Nix, i, j)
      rhomomp(i, j) = rhojlmp(Zix, Nix, i, j)
    enddo
  enddo
  call momjlmecis(Z, N, pZ, e, lv1, lw1, 1.25d0, 1.35d0, rhomomn, rhomomp, vpot, rcjlm)
  do i = 1, numjlm
    do j = 1, 6
      potjlm(Zix, Nix, i, j) = min(vpot(i, j), 1.e30)
      potjlm(Zix, Nix, i, j) = max(potjlm(Zix, Nix, i, j), - 1.e30)
    enddo
  enddo
  normjlm(Zix, Nix, 1) = lv
  normjlm(Zix, Nix, 2) = lw
  normjlm(Zix, Nix, 3) = lv1
  normjlm(Zix, Nix, 4) = lw1
  normjlm(Zix, Nix, 5) = lvso
  normjlm(Zix, Nix, 6) = lwso
  rc = rcjlm
  return
end subroutine mom
! Copyright A.J. Koning 2021
