subroutine astroinit
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialization of arrays for astro calculations
!
! Author    : Stephane Hilaire and Stephane Goriely
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
! Variables for astro initialization
!   xsastro           ! cross section for astrophysical calculatio
!   xsastroex         !  cross section for astrophysical calculati
!   xsastrofis        ! astrophysical fission cross section
! Variables for astro initialization
!   macsastro         ! Maxwellian-averaged thermonuclear reaction
!   macsastroex       ! thermonuclear reaction cross section to a
!   macsastrofis      ! thermonuclear reaction cross section for f
!   macsastroracap    ! Maxwellian-averaged thermonuc. reac. c.s.
!   maxNastro         ! maximal number of neutrons away from initi
!   maxZastro         ! maximal number of protons away from initia
!   partf             ! integrated partition function
!   rateastro         ! thermonuclear reaction rate
!   rateastroex       ! thermonuclear reaction rate to a given exc
!   rateastrofis      ! thermonuclear reaction rate factor for fis
!   rateastroracap    ! thermonuclear reaction rate for direct cap
!
! *** Initialization
!
  maxZastro = numZastro
  maxNastro = numNastro
  rateastro = 0.
  macsastro = 0.
  rateastroex = 0.
  macsastroex = 0.
  partf = 0.
  rateastrofis = 0.
  rateastroracap = 0.
  macsastrofis = 0.
  macsastroracap = 0.
  xsastro = 0.
  xsastroex = 0.
  xsastrofis = 0.
  return
end subroutine astroinit
! Copyright A.J. Koning 2021
