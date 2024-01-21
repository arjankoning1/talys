function prodm(sx, tjl, numtjl, numtr, Ninc)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculation of the product of (1-t*x)**(1/10)
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
  real(dbl) :: prodm           ! product function for GOE
  real(dbl) :: tav             ! average transmission coefficients
  real(dbl) :: term            ! help variable
  real(dbl) :: tjl(0:5, numtr) ! transmission coefficients
!
! **************** Calculation of the product of 1-t*x *****************
!
! prodm : product function for GOE
!
  prodm = 1.
  do i = Ninc + 1, numtjl
    if (tjl(0, i) == 0.) then
      tav = 0.
    else
      tav = tjl(1, i) / tjl(0, i)
    endif
    term = 1. - tav * sx
    if (term <= 0.) cycle
    prodm = prodm * real(term **(tjl(0, i) * 0.1))
  enddo
  return
end function prodm
! Copyright A.J. Koning 2021
