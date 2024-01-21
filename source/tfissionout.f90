subroutine tfissionout(Zcomp, Ncomp, nex)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of fission transmission coefficients
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
! Variables for compound nucleus from target
!   Exinc    ! excitation energy of entrance bin
! Variables for excitation energy grid
!   maxJ     ! maximal J - value
! Variables for fission transmission coefficients
!   denfis    ! fission level density
!   gamfis    ! fission width
!   taufis    ! fission lifetime
!   tfis      ! fission transmission coefficient for Hill - Wheeler magnitude
! Variables for nuclides
!   AA       ! mass number of residual nucleus
!   NN       ! neutron number of residual nucleus
!   ZZ       ! charge number of residual nucleus
! Constants
!   nuc      ! symbol of nucleus
!
! *** Declaration of local data
!
  implicit none
  integer :: A     ! mass number of target nucleus
  integer :: J     ! spin of level
  integer :: J2    ! 2 * J
  integer :: N     ! neutron number of residual nucleus
  integer :: Ncomp ! neutron number index for compound nucleus
  integer :: nex   ! excitation energy bin of compound nucleus
  integer :: odd   ! odd (1) or even (0) nucleus
  integer :: Z     ! charge number of target nucleus
  integer :: Zcomp ! proton number index for compound nucleus
!
! *************** Output of fission transmission coefficients **********
!
  Z = ZZ(Zcomp, Ncomp, 0)
  N = NN(Zcomp, Ncomp, 0)
  A = AA(Zcomp, Ncomp, 0)
  write(*, '(/" Fission transmission coefficients for Z=", i3, &
 &  " N=", i3, " (", i3, a2, ") and an excitation energy of ", f8.3, " MeV"/)') Z, N, A, nuc(Z), Exinc
  write(*, '("   J      T(J,-)      T(J,+)    Gamma(J,-)  Gamma(J,+)   tau(J,-)    tau(J,+)  density(J,-) density(J,+)"/)')
  odd = mod(A, 2)
  do J = 0, maxJ(Zcomp, Ncomp, nex)
    J2 = 2 * J + odd
    write(*, '(1x, f4.1, 2x, 8es12.5)') 0.5*J2, tfis(J, -1), tfis(J, 1), gamfis(J, -1), gamfis(J, 1), taufis(J, -1), &
 &  taufis(J, 1), denfis(J, -1), denfis(J, 1)
  enddo
  write( * , * )
  return
end subroutine tfissionout
! Copyright A.J. Koning 2021
