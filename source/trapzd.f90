subroutine trapzd(Zix, Nix, abeg, bend, snew, nk)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Integration
!
! Author    : Marieke Duijvestijn
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               dbl  , &               ! double precision kind
               sgl                    ! single precision kind
!
! *** Declaration of local data
!
  integer   :: it          ! counter for tritons
  integer   :: jk          ! help variable
  integer   :: Nix         ! neutron number index for residual nucleus
  integer   :: nk          ! help variable
  integer   :: Zix         ! charge number index for residual nucleus
  real(sgl) :: abeg        ! start of A loop
  real(sgl) :: bend        ! end of loop
  real(sgl) :: del         ! difference
  real(sgl) :: x           ! help variable
  real(dbl) :: funcfismode ! function for transmission coefficient per fission mode
  real(dbl) :: snew        ! new trial value
  real(dbl) :: sum         ! help variable
  real(dbl) :: tnm         ! help variable
!
! **********************************************************************
!
  snew = 0.0
  if (nk == 1) then
    snew = 0.5 * (bend - abeg) * (funcfismode(Zix, Nix, abeg) + funcfismode(Zix, Nix, bend))
  else
    it = 2 **(nk - 2)
    tnm = it
    del = (bend - abeg) / tnm
    x = abeg + 0.5 * del
    sum = 0.
    do jk = 1, it
      sum = sum + funcfismode(Zix, Nix, x)
      x = x + del
    enddo
    snew = 0.5 * (snew + (bend - abeg) * sum / tnm)
  endif
  return
end subroutine trapzd
! Copyright A.J. Koning 2021
