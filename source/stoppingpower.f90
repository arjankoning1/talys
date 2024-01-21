subroutine stoppingpower(E, dEdx)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculate stopping power
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
! Variables for medical isotope production
!   rhotarget    ! target material density
! Variables for main input
!   Atarget      ! mass number of target nucleus
!   k0           ! index of incident particle
!   Ztarget      ! charge number of target nucleus
! Constants
!   amu          ! atomic mass unit in MeV
!   emass        ! electron mass in MeV / c^2
!   parmass      ! mass of particle in a.m.u.
!   parZ         ! charge number of particle
!
! *** Declaration of local data
!
  implicit none
  real(sgl) :: beta  ! density-dependenceterm of the DDM3Y interaction
  real(sgl) :: dEdx  ! stopping power
  real(sgl) :: denom ! help variable
  real(sgl) :: E     ! incident energy
  real(sgl) :: enum  ! enumerator of Lorentzian
  real(sgl) :: eta   ! beta*gam
  real(sgl) :: gam   ! Brosa parameter
  real(sgl) :: Imean ! mean excitation potential in MeV
  real(sgl) :: pmass ! rest mass of projectile in MeV (m_0*c^2)
  real(sgl) :: term1 ! help variable
  real(sgl) :: term2 ! help variable
  real(sgl) :: Wmax  ! maximum energy transfer in knock-on collision in MeV
!
! ********************* Velocity of projectile *************************
!
  pmass = parmass(k0)*amu
  enum = E * (E + 2. * pmass)
  denom = (E + pmass) **2
  beta = sqrt(enum / denom)
!
! ********************* Bethe-Bloch formula ****************************
!
! The formula for the stopping power can be found on page 24 of W.R. Leo, "Techniques for nuclear and particle physics experiments",
! Springer Verlag, Berlin, 1994
!
  Imean = Ztarget * (9.76 + 58.8 * Ztarget **( - 1.19)) * 1.e-6
  gam = 1. / sqrt(1. - beta **2)
  eta = beta * gam
  Wmax = 2. * emass * (eta **2)
  term1 = 0.1535 * rhotarget * Ztarget / real(Atarget) * (parZ(k0) **2) / (beta **2)
  term2 = log(2. * emass * (eta **2) * Wmax / (Imean **2)) - 2. * beta **2
  dEdx = term1 * term2
  return
end subroutine stoppingpower
! Copyright A.J. Koning 2021
