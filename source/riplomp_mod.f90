subroutine riplomp_mod(NE, Zt, At, Eripl, index1, index2, index3)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Interface for RIPL OMP
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               sgl                  ! single precision kind
  use om_retrieve
!
! *** Declaration of local data
!
  integer   :: NE        ! number of incident energies
  real(sgl) :: Eripl(NE) !
  integer   :: Zt        !
  integer   :: At        !
  integer   :: index1    !
  integer   :: index2    !
  integer   :: index3    !
!
! ************************ Interface to RIPL OMP ***********************
!
  Number_Energies = NE
  Ztarget = Zt
  Atarget = At
  Energies = 0.
  do i = 1, NE
    Energies(i) = Eripl(i)
  enddo
  RIPL_Index = index1
  Calc_Type = index2
  Def_Index = index3
  call retrieve
end subroutine riplomp_mod
! Copyright A.J. Koning 2021
