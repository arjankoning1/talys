subroutine residualBU
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Residual production cross sections for break-up
!
! Author    : Marilena Avrigeanu
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl        ! single precision kind
! All global variables
!   numlev     ! maximum number of discrete levels
! Variables for levels
!   tau            ! lifetime of state in seconds
! Variables for numerics
!   maxN            ! maximal number of neutrons away from initial compound nucleus
!   maxZ            ! maximal number of protons away from initial compound nucleus
! Variables for preequilibrium
!   xsBFnuc    ! inelastic breakup enhancement brought by breakup neutrons
! Variables for incident channel
!   maxA            ! maximal number of nucleons away from initial compound nucleus
!   xspopex        ! population cross section summed over spin and parity
!   xspopnuc       ! population cross section per nucleus
! Variables for level density
!   Nlast          ! last discrete level
!
! *** Declaration of local data
!
  implicit none
  integer   :: Acomp            ! mass number index for compound nucleus
  integer   :: Ncomp            ! neutron number index for compound nucleus
  integer   :: nex              ! excitation energy bin of compound nucleus
  integer   :: Zcomp            ! proton number index for compound nucleus
  real(sgl) :: branch(0:numlev) ! branching ratio to a given excited state
!
! ************************ Cross sections ******************************
!
  do Acomp = 0, maxA
    do Zcomp = 0, maxZ
      Ncomp = Acomp - Zcomp
      if (Ncomp < 0 .or. Ncomp > maxN) cycle
      if (xspopnuc(Zcomp, Ncomp) == 0.) cycle
      do nex = 0, Nlast(Zcomp, Ncomp, 0)
        branch(nex) = 0.
        if (nex == 0 .or. tau(Zcomp, Ncomp, nex) /= 0.) then
          branch(nex) = xspopex(Zcomp, Ncomp, nex) / xspopnuc(Zcomp, Ncomp)
        endif
      enddo
      xspopnuc(Zcomp, Ncomp) = xspopnuc(Zcomp, Ncomp) + xsBFnuc(Zcomp, Ncomp)
      do nex = 0, Nlast(Zcomp, Ncomp, 0)
        if (nex == 0 .or. tau(Zcomp, Ncomp, nex) /= 0.) then
          xspopex(Zcomp, Ncomp, nex) = branch(nex) * xspopnuc(Zcomp, Ncomp)
        endif
      enddo
    enddo
  enddo
  return
end subroutine residualBU
! Copyright A.J. Koning 2021
