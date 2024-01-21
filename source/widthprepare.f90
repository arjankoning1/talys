subroutine widthprepare
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Preparation of width fluctuation corrections
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
!   numtrans     ! number of transmission coefficients
! Variables for compound reactions
!   wmode        ! designator for width fluctuation model
! Variables to prepare information for initial compound nucleus
!   denomhf      ! denominator for compound nucleus formula
!   tnum         ! counter for width fluctuation calculation
!   tNinc        ! counter for width fluctuation calculation
!   transjl      ! array for width fluctuation calculation
! Variables for initial compound nucleus
!   ngoep        ! number of points for Gauss - Legendre integration
!   ngoes        ! number of points for Gauss - Legendre integration
!   ngoet        ! number of points for Gauss - Legendre integration
!   nmold        ! number of points for Gauss - Laguerre integration
!   wgoep        ! variables for Gauss - Laguerre integration
!   wgoes        ! variables for Gauss - Laguerre integration
!   wgoet        ! variables for Gauss - Laguerre integration
!   wmold        ! variables for Gauss - Laguerre integration
!   xgoep        ! variables for Gauss - Legendre integration
!   xgoes        ! variables for Gauss - Legendre integration
!   xgoet        ! variables for Gauss - Legendre integration
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
! ******************* Choice of width fluctuation model ****************
!
! molprepare  : subroutine for preparation of Moldauer width fluctuation correction
! hrtwprepare : subroutine for preparation of HRTW width fluctuation correction
! goeprepare  : subroutine for preparation of GOE triple integral width fluctuation correction
!
! Width fluctuation models:
! wmode= 0: Hauser-Feshbach (no width fluctuations)
! wmode= 1: Moldauer
! wmode= 2: HRTW
! wmode= 3: GOE
!
! First, all width fluctuation variables that only depend on J and P and not on the other angular momentum quantum numbers
! are calculated.
!
  if (wmode == 1) call molprepare(transjl, tnum, denomhf, nmold, xmold, wmold, tjlav, freedom, prodwidth, numtrans, tNinc, &
 &  WFCfactor)
  if (wmode == 2) call hrtwprepare(transjl, tnum, denomhf, tjlav, sumhrtw, vhrtw, whrtw, numtrans, tNinc, WFCfactor)
  if (wmode == 3) call goeprepare(transjl, tnum, denomhf, ngoep, ngoes, ngoet, xgoep, xgoes, xgoet, wgoep, wgoes, wgoet, tjlav, &
 &  agoe1, agoe2, agoe3, agoe4, agoe5, agoe6, agoe7, agoe8, sgoe1, sgoe2, sgoe3, sgoe4, sgoe5, numtrans, tNinc)
  return
end subroutine widthprepare
! Copyright A.J. Koning 2022
