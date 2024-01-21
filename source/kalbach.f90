function kalbach(type, Ein, Eout, ang)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Kalbach systematics
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
!   sgl         ! single precision kind
! All global variables
!   numpar      ! number of particles
! Variables for main input
!   k0          ! index of incident particle
! Variables for nuclides
!   Smyers      ! Myers - Swiatecki separation energy
! Constants
!   fourpi      ! 4. * pi
!
! *** Declaration of local data
!
  implicit none
  integer   :: type          ! particle type
  real(sgl) :: akal          ! Kalbach 'a' parameter
  real(sgl) :: ang           ! angle
  real(sgl) :: C1            ! constant
  real(sgl) :: C2            ! constant
  real(sgl) :: C3            ! constant
  real(sgl) :: Cmin(numpar)  ! constant of Kalbach systematics
  real(sgl) :: Cmout(numpar) ! constant of Kalbach systematics
  real(sgl) :: ea            ! help variable
  real(sgl) :: eb            ! help variable
  real(sgl) :: Ein           ! incident energy
  real(sgl) :: Ek1           ! energy
  real(sgl) :: Ek3           ! energy
  real(sgl) :: Eout          ! outgoing energy
  real(sgl) :: Et1           ! constant of Kalbach systematics
  real(sgl) :: Et3           ! constant of Kalbach systematics
  real(sgl) :: kalbach       ! Kalbach function
  real(sgl) :: X1            ! energy
  real(sgl) :: X3            ! energy
!
! ************************ Kalbach systematics *************************
!
! kalbach   : Kalbach function
!
! Systematics of Kalbach: Phys. Rev. C37, 2350, (1987)
!
! Since we separate the pre-equilibrium (xspreeqad) and compound (xscompad) angular distribution in the output we only need to take
! the forward peaked component of the Kalbach formula to calculate the pre-equilibrium angular distribution.
!
  Cmin = (/ 1., 1., 1., 1., 1., 0. /)
  Cmout = (/ 0.5, 1., 1., 1., 1., 2. /)
  Et1 = 130.
  Et3 = 41.
  C1 = 0.04
  C2 = 1.8e-6
  C3 = 6.7e-7
!
! Isotropic distribution for photons
!
  if (k0 == 0 .or. type == 0) then
    kalbach = 1. / fourpi
  else
    ea = Ein + Smyers(k0)
    Ek1 = min(ea, Et1)
    Ek3 = min(ea, Et3)
    eb = Eout + Smyers(type)
    X1 = Ek1 * eb / ea
    X3 = Ek3 * eb / ea
    akal = C1 * X1 + C2 * (X1 **3) + C3 * Cmin(k0) * Cmout(type) * (X3 **4)
    kalbach = 0.
    if (abs(akal) <= 80.) kalbach = 1. / fourpi * akal / sinh(akal) * exp(akal * cos(ang))
  endif
  return
end function kalbach
! Copyright A.J. Koning 2021
