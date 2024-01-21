subroutine isoprod
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculate isotope production
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
! prodinitial    : subroutine for initialization of isotope production info
! prodres        : subroutine for residual production cross sections
! prodrates      : subroutine to calculate reaction rates
! prodyield      : subroutine to calculate production yields
! prodout        : subroutine for output of isotope production
!
  call prodinitial
  call prodres
  call prodrates
  call prodyield
  call prodout
  return
end subroutine isoprod
! Copyright A.J. Koning 2021
