function preeqpair(Zix, Nix, n, E, pmodel)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Pairing effects in exciton model
!
! Author    : Marieke Duijvestijn and Arjan Koning
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
! Variables for preequilibrium
!   flag2comp    ! flag for two - component pre - equilibrium model
!   g            ! single - particle level density parameter
!   gn           ! single - particle neutron level density parameter
!   gp           ! single - particle proton level density parameter
! Variables for level density
!   pair         ! pairing energy
!
! *** Declaration of local data
!
  implicit none
  integer   :: n         ! exciton number
  integer   :: Nix       ! neutron number index for residual nucleus
  integer   :: pmodel    ! Fu formula
  integer   :: Zix       ! charge number index for residual nucleus
  real(sgl) :: E         ! incident energy
  real(sgl) :: gs        ! single-particle level density parameter
  real(sgl) :: gsn       ! single-particle neutron level density parameter
  real(sgl) :: gsp       ! single-particle proton level density parameter
  real(sgl) :: ncrit     ! average quasi-particle number at the Tc
  real(sgl) :: pair0     ! ground-state pairing gap
  real(sgl) :: pairCN    ! pairing energy
  real(sgl) :: pairex    ! excited state pairing gap
  real(sgl) :: preeqpair ! pre-equilibrium pairing energy
  real(sgl) :: Tc        ! function for Coulomb dip
!
! **************************** Fu formula ******************************
!
  preeqpair = 0.
  pairCN = pair(Zix, Nix)
  if (pairCN <= 0.) return
!
! pmodel 1: Fu formula
!
  if (pmodel == 1) then
    if (flag2comp) then
      gsp = gp(Zix, Nix)
      gsn = gn(Zix, Nix)
      gs = gsp + gsn
    else
      gs = g(Zix, Nix)
    endif
    pair0 = sqrt(pairCN / (0.25 * gs))
    Tc = 2. * pair0 / 3.5
    ncrit = 2. * gs * Tc * log(2.)
    if (E / pairCN >= (0.716 + 2.44 * (n / ncrit) **2.17)) then
      pairex = pair0 * (0.996 - 1.76 * (n / ncrit) **1.6 / ((E / pairCN) **0.68))
    else
      pairex = 0.
    endif
    preeqpair = pairCN - 0.25 * gs * (pairex **2.)
  else
!
! pmodel 2: compound nucleus value
!
    preeqpair = pairCN
  endif
  return
end function preeqpair
! Copyright A.J. Koning 2021
