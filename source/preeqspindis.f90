subroutine preeqspindis
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Pre-equilibrium spin distribution
!
! Author    : Arjan Koning
!
! 2023-08-15: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl           ! single precision kind
! Variables for main input
!   Atarget       ! mass number of target nucleus
! Variables for energies
!   flagmulpre    ! flag for multiple pre - equilibrium calculation
! Constants
!   sqrttwopi     ! sqrt(2. * pi)
!   twothird      ! 2 / 3
! Variables for preequilibrium initialization
!   maxJph        ! maximal spin for particle - hole states
!   maxpar        ! maximal particle number
!   RnJ           ! spin distribution for particle - hole stat
!   RnJsum        ! (2J + 1) * sum over spin distributions
!   Rspincutpreeq ! adjustable constant (global) for preequilibrium spin cutoff factor
!
! *** Declaration of local data
!
  implicit none
  integer   :: J        ! spin
  integer   :: n        ! exciton number
  real(sgl) :: denom    ! help variable
  real(sgl) :: expo     ! help variable
  real(sgl) :: sigma2ph ! spin cutoff factor for particle-hole states
  real(dbl) :: newspin! new spin cutoff parameter
  real(dbl) :: spin_wigner ! Wigner distribution
  real(dbl) :: f(3)   ! help variable
  real(dbl) :: eta    ! help variable
!
! ************* Spin distribution for particle-hole states *************
!
! Spin cutoff factor for p-h density: Gruppelaar, Group Meeting on level densities, Brookhaven (1983), p. 143.
!
  maxJph = 30
  RnJ = 0.
  do n = 1, maxexc
    RnJsum(n) = 0.
    sigma2ph = Rspincutpreeq * 0.24 * n * Atarget **twothird
    do J = 0, maxJph
      denom = 2. * sqrttwopi * sigma2ph **1.5
      expo = - (J + 0.5) **2 / (2. * sigma2ph)
      RnJ(n, J) = (2. * J + 1) / denom * exp(expo)
      RnJsum(n) = RnJsum(n) + (2. * J + 1) * RnJ(n, J)
    enddo
  enddo
!
! Marc Dupuis: use a spin cut-off value inferred from JLM/QRPA calculations.
! Microscopic description of target spin distribution after inelastic scattering to the continuum
! Marc  Dupuis, Toshihiko  Kawano, Maelle  Kerveno, Stephane  Hilaire
! EPJ Web of Conf. 284 03003 (2023)
! DOI: 10.1051/epjconf/202328403003
!
  if (pespinmodel == 4) then
    do n = 1, maxexc
      RnJsum(n) = 0.
      f = 0.d0 
      eta = sqrt(Rspincutpreeq) * newspin(Atarget,real(Einc,8),n,1,f)
      do J = 0,maxJph
        RnJ(n, J) =  spin_wigner(real(J,8),eta) / (2. * J + 1)
        RnJsum(n) = RnJsum(n) + (2. * J + 1) * RnJ(n, J)
      enddo
    enddo
  endif
  return
end subroutine preeqspindis
! Copyright A.J. Koning 2023
