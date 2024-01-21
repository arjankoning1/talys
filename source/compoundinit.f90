subroutine compoundinit
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialization of compound model parameters
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
!   numfact     ! number of terms for factorial logarithm
! Variables for basic reaction
!   flagang     ! flag for output of angular distributions
! Variables for compound reactions
!   ewfc        ! off - set incident energy for width fluctuation
!   wmode       ! designator for width fluctuation model
! Variables for input energies
!   enincmin    ! minimum incident energy
! Variables for initial compund nucleus
!   logfact     ! factorial logarithm
!   ngoep       ! number of points for Gauss - Legendre integration
!   ngoes       ! number of points for Gauss - Legendre integration
!   ngoet       ! number of points for Gauss - Legendre integration
!   nmold       ! number of points for Gauss - Laguerre integration
!   wgoep       ! variables for Gauss - Laguerre integration
!   wgoes       ! variables for Gauss - Laguerre integration
!   wgoet       ! variables for Gauss - Laguerre integration
!   wmold       ! variables for Gauss - Laguerre integration
!   wpower      ! power used for rho * (t **wpower)
!   xgoep       ! variables for Gauss - Legendre integration
!   xgoes       ! variables for Gauss - Legendre integration
!   xgoet       ! variables for Gauss - Legendre integration
!   xmold       ! variables for Gauss - Laguerre integration
!
! *** Declaration of local data
!
  implicit none
  integer :: k                 ! counter
!
! *********** Initialization for width fluctuation calculation *********
!
! gaulag         : subroutine for calculation of Gauss-Laguerre arrays
! gauleg         : subroutine for calculation of Gauss-Legendre arrays
!
! Generate weight and nodes for Gauss-Legendre/Laguerre integration
!
  if (enincmin <= ewfc) then
    wpower = 1
!
! 1. Moldauer
!
    if (wmode == 1) then
      nmold = 32
      call gaulag(nmold, xmold, wmold)
    endif
!
! 2. HRTW
!
    if (wmode == 2) wpower = 2
!
! 3. GOE
!
    if (wmode == 3) then
      wpower = 5
      ngoep = 50
      ngoes = 50
      ngoet = 50
      call gauleg(ngoep, xgoep, wgoep)
      call gauleg(ngoes, xgoes, wgoes)
      call gauleg(ngoet, xgoet, wgoet)
    endif
  endif
!
! ***** Initialization for compound nucleus angular distributions ******
!
  if (flagang) then
    logfact(1) = 0.
    logfact(2) = 0.
    do k = 3, numfact
      logfact(k) = logfact(k - 1) + log(float(k - 1))
    enddo
  endif
  return
end subroutine compoundinit
! Copyright A.J. Koning 2021
