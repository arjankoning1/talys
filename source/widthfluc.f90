subroutine widthfluc(ielas)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Width fluctuation corrections
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
! All global variables
!   numhill      ! maximum number of Hill - Wheeler points
!   numtrans     ! number of transmission coefficients
! Variables for compound reactions
!   wmode        ! designator for width fluctuation model
! Variables for compound nucleus from target
!   tnumi        ! counter for width fluctuation calculation
!   tnumo        ! counter for width fluctuation calculation
!   Wab          ! width fluctuation factor
! Variables to prepare information for initial compound nucleus
!   denomhf      ! denominator for compound nucleus formula
!   tnum         ! counter for width fluctuation calculation
!   transjl      ! array for width fluctuation calculation
! Variables for initial compound nucleus
!   ngoep        ! number of points for Gauss - Legendre integration
!   ngoes        ! number of points for Gauss - Legendre integration
!   ngoet        ! number of points for Gauss - Legendre integration
!   nmold        ! number of points for Gauss - Laguerre integration
!   xmold        ! variables for Gauss - Laguerre integration
! Variables for width fluctuation
!   agoe1        ! variable for GOE triple integral calculation
!   agoe2        ! variable for GOE triple integral calculation
!   agoe3        ! variable for GOE triple integral calculation
!   agoe4        ! variable for GOE triple integral calculation
!   agoe5        ! variable for GOE triple integral calculation
!   agoe6        ! variable for GOE triple integral calculation
!   agoe7        ! variable for GOE triple integral calculation
!   agoe8        ! variable for GOE triple integral calculation
!   freedom      ! number of degrees of freedom
!   prodwidth    ! product of widths
!   sgoe1        ! variable for GOE triple integral calculation
!   sgoe2        ! variable for GOE triple integral calculation
!   sgoe3        ! variable for GOE triple integral calculation
!   sgoe4        ! variable for GOE triple integral calculation
!   sgoe5        ! variable for GOE triple integral calculation
!   sumhrtw      ! variable for HRTW calculation
!   tjlav        ! array for width fluctuation calculation
!   vhrtw        ! variable for HRTW calculation
!   whrtw        ! variable for HRTW calculation
!
! *** Declaration of local data
!
  implicit none
  integer :: ielas              ! designator for elastic channel
!
! ******************* Choice of width fluctuation model ****************
!
! moldauer : subroutine for Moldauer width fluctuation correction
! hrtw     : subroutine for HRTW width fluctuation correction
! goe      : subroutine for GOE triple integral width fluctuation correction
!
! Width fluctuation models:
! wmode= 0: Hauser-Feshbach (no width fluctuations)
! wmode= 1: Moldauer
! wmode= 2: HRTW
! wmode= 3: GOE
!
  if (tnumo > numtrans-numhill-2) return
  if (wmode == 1) call moldauer(tnum, tnumi, tnumo, denomhf, nmold, xmold, tjlav, freedom, prodwidth, Wab, numtrans, ielas)
  if (wmode == 2) call hrtw(transjl, tnumi, tnumo, denomhf, sumhrtw, vhrtw, whrtw, Wab, numtrans, ielas)
  if (wmode == 3) call goe(transjl, tnum, tnumi, tnumo, denomhf, ngoep, ngoes, ngoet, tjlav, agoe1, agoe2, agoe3, agoe4, agoe5, &
 &  agoe6, agoe7, agoe8, sgoe1, sgoe2, sgoe3, sgoe4, sgoe5, Wab, numtrans, ielas)
  return
end subroutine widthfluc
! Copyright A.J. Koning 2021
