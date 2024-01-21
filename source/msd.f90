subroutine msd
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Multi-step direct model
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
! msdinit    : subroutine for initialization of MSD model parameters
! dwbaecis   : subroutine for ECIS calculations of DWBA for MSD
! msdcalc    : subroutine for MSD calculation
! msdtotal   : subroutine for total multi-step direct cross sections
! msdout     : subroutine for output of multi-step direct cross sections
! msdplusmsc : subroutine for total quantum-mechanical pre-equilibrium cross sections
!
! The MSD model is implemented for neutrons and protons only.
!
  call msdinit
  call dwbaecis
  call msdcalc
  call msdtotal
  call msdout
  call msdplusmsc
  return
end subroutine msd
! Copyright A.J. Koning 2021
