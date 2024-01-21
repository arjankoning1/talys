function omega(Zix, Nix, p, h, gs, Eex, rJ)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Particle-hole state density
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
!   sgl       ! single precision kind
! Variables for preequilibrium initialization
!   Efermi    ! depth of Fermi well
!   RnJ       ! spin distribution for particle - hole stat
!
! *** Declaration of local data
!
  implicit none
  logical   :: surfwell ! flag for surface effects in finite well
  integer   :: h        ! help variable
  integer   :: Ji       ! integer of spin
  integer   :: Nix      ! neutron number index for residual nucleus
  integer   :: p        ! particle number
  integer   :: Zix      ! charge number index for residual nucleus
  real(sgl) :: Eex      ! excitation energy
  real(sgl) :: gs       ! single-particle level density parameter
  real(sgl) :: omega    ! particle-hole state density
  real(sgl) :: phdens   ! function for particle-hole state density
  real(sgl) :: rJ       ! help variable
!
! ******************** Particle-hole state density *********************
!
  Ji = int(rJ)
  surfwell = .false.
  omega = (2. * rJ + 1.) * RnJ(2, Ji) * phdens(Zix, Nix, p, h, gs, Eex, Efermi, surfwell)
  return
end function omega
! Copyright A.J. Koning 2021
