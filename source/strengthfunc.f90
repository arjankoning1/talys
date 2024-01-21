subroutine strengthfunc(Tinc, l, j)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : (l,j) neutron strength function for URR calculations
!
! Author    : Gilles Noguere
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
!   numl          ! number of l values
! Variables for energy grid
!   Einc          ! incident energy in MeV
! Variables for energies
!   wavenum       ! wave number
! Variables for nuclides
!   tarmass       ! mass of target nucleus
! Constants
!   onethird      ! 1 / 3
!   twopi         ! 2 * pi
! Variables for URR
!   strengthlj    ! (l, j) neutron strength function
!
! *** Declaration of local data
!
  implicit none
  integer   :: i            ! counter
  integer   :: j            ! counter
  integer   :: l            ! multipolarity
  real(sgl) :: ac           ! channel radius in ENDF-6 convention
  real(sgl) :: den          ! help variable
  real(sgl) :: num          ! help variable
  real(sgl) :: P(0:numl)    ! pairing energy
  real(sgl) :: ri           ! help variable
  real(sgl) :: Shift(0:numl)! hard sphere shift factor for s-wave
  real(sgl) :: Tinc         ! transmission coefficients as a function of j and l  for the incident channel
  real(sgl) :: vl           ! penetrability factor
!
! ******************** Initialization parameters ***********************
!
! strengthlj: (l,j) neutron strength function
!
  ac = 1.23*tarmass**onethird+0.8
  P(0) = wavenum * ac
  Shift(0) = 0.
!
! *********************  Penetrability factor **************************
!
  do i = 1, l
    ri = real(i)
    num = P(i - 1) * P(0) **2
    den = P(i - 1) **2 + (ri - Shift(i - 1)) **2
    P(i) = num / den
    num = (ri - Shift(i - 1)) * P(0) **2
    den = P(i - 1) **2 + (ri - Shift(i - 1)) **2
    Shift(i) = num / den - ri
  enddo
!
! *********************** Strength function ****************************
!
  vl = P(l) / (P(0) * 10.)
  strengthlj(l, j) = Tinc / (twopi * vl * sqrt(Einc))
  return
end subroutine strengthfunc
! Copyright A.J. Koning 2021
