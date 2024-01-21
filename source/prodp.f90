function prodp(sx, tjl, numtjl, numtr, Ninc)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculation of the product of (1+t*x)**(1/20)
!
! Author    : Stephane Hilaire
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               dbl  , & ! double precision kind
               sgl      ! single precision kind
!
! *** Declaration of local data
!
  implicit none
  integer   :: i               ! counter
  integer   :: Ninc            ! number of incident energies
  integer   :: numtjl          ! number of transmission coefficients
  integer   :: numtr           ! number of transmission coefficients
  real(sgl) :: sx              ! help variable
  real(dbl) :: prodp           ! product function for GOE
  real(dbl) :: tav             ! average transmission coefficients
  real(dbl) :: tjl(0:5, numtr) ! transmission coefficients
!
! **************** Calculation of the product of 1+t*x *****************
!
! prodp : product function for GOE
!
  prodp = 1.
  do i = Ninc + 1, numtjl
    if (tjl(0, i) == 0.) then
      tav = 0.
    else
      tav = tjl(1, i) / tjl(0, i)
    endif
    prodp = prodp * real((1. + tav * sx) **(tjl(0, i) * 0.05))
  enddo
  return
end function prodp
! Copyright A.J. Koning 2021
