subroutine trans(Zix, Nix, transm)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Transmission coefficients per fission mode
!
! Author    : Marieke Duijvestijn
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl    ! single precision kind
!   dbl    ! double precision kind
! Variables for mass distribution
!   excfis ! excitation energy at fission
!
! *** Declaration of local data
!
  implicit none
  integer   :: i           ! counter
  integer   :: icount      ! counter
  integer   :: Nix         ! neutron number index for residual nucleus
  integer   :: Zix         ! charge number index for residual nucleus
  real(sgl) :: xa          ! intersection of the segment with the vertical axis linking  (xl,yl) with (xl,yu)
  real(sgl) :: xb          ! intersection of the segment with the vertical axis linking  (xu,yl) with (xu,yu)
  real(dbl) :: funcfismode ! function for transmission coefficient per fission mode
  real(dbl) :: snew        ! new trial value
  real(dbl) :: sold        ! new trial value
  real(dbl) :: transm      ! function for transmission coefficient per fission mode
!
! **********************************************************************
!
! funcfismode: function for transmission coefficient per fission mode
! transm     : function for transmission coefficient per fission mode
!
  xa = 0.0
  xb = excfis
  if (excfis == 0) then
    transm = funcfismode(Zix, Nix, 0.)
  else
    icount = 0
    i = 0
    snew = 10.0
    sold = 1.0
    do
      if (abs(snew / sold) <= 0.99 .or. abs(snew / sold) >= 1.01) then
        sold = snew
        i = i + 1
        if (i > 10) exit
        call trapzd(Zix, Nix, xa, xb, snew, i)
      else
        icount = icount + 1
        if (icount == 2) exit
        sold = snew
        i = i + 1
        if (i > 10) exit
        call trapzd(Zix, Nix, xa, xb, snew, i)
      endif
    enddo
    transm = snew
  endif
end subroutine trans
! Copyright A.J. Koning 2021
