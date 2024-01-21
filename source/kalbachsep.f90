subroutine kalbachsep
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Separation energies for Kalbach systematics
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
! Variables for nuclides
!   Smyers      ! Myers - Swiatecki separation energy
! Variables for main input
!   Ainit       ! mass number of initial compound nucleus
!   Zinit       ! charge number of initial compound nucleus
! Variables for nuclides
!   AA          ! mass number of residual nucleus
!   ZZ          ! charge number of residual nucleus
! Constants
!   onethird    ! 1 / 3
!   parsym      ! symbol of particle
!   twothird    ! 2 / 3
!
! *** Declaration of local data
!
  implicit none
  integer   :: A               ! mass number of target nucleus
  integer   :: type            ! particle type
  integer   :: Z               ! charge number of target nucleus
  real(sgl) :: Ab              ! help variable
  real(sgl) :: Ac              ! help variable
  real(sgl) :: Ib              ! energy required for breakup
  real(sgl) :: Nb              ! number of solutions
  real(sgl) :: Nc              ! help variable
  real(sgl) :: Zb              ! charge of ejectile
  real(sgl) :: Zc              ! help variable
!
! ************** Spherical Myers-Swiatecki parameters ******************
!
! We use the 'old' Myers-Swiatecki parameters as used by Kalbach in Phys. Rev. C37, 2350, (1988).
!
  do type = 1, 6
    Z = ZZ(0, 0, type)
    A = AA(0, 0, type)
    Ac = real(Ainit)
    Ab = real(A)
    Zc = real(Zinit)
    Zb = real(Z)
    Nc = Ac - Zc
    Nb = Ab - Zb
    Ib = 0.
    if (parsym(type) == 'd') Ib = 2.225
    if (parsym(type) == 't') Ib = 8.482
    if (parsym(type) == 'h') Ib = 7.718
    if (parsym(type) == 'a') Ib = 28.296
    Smyers(type) = 15.68 * (Ac - Ab) - 28.07 * ((Nc - Zc) **2 / Ac - (Nb - Zb) **2 / Ab) - &
 &    18.56 * (Ac **twothird - Ab **twothird) + 33.22 * ((Nc - Zc) **2 / Ac **(2. * twothird) - (Nb - Zb) **2 / Ab ** &
 &    (2. * twothird)) - 0.717 * (Zc **2 / (Ac **onethird) - &
 &    Zb **2 / (Ab **onethird)) + 1.211 * (Zc **2 / Ac - Zb **2 / Ab) - Ib
  enddo
  return
end subroutine kalbachsep
! Copyright A.J. Koning 2021
