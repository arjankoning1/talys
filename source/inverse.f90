subroutine inverse(Zcomp, Ncomp)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculation of total, reaction and elastic cross sections
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
! Variables for output
!   flaginverse     ! flag for output of transmission coefficients and inverse cross sections
! Variables for direct reactions
!   flageciscalc    ! flag for new ECIS calculation for outgoing particles
! Variables for inverse channel data
!   csfile          ! file with inverse reaction cross sections
!   transfile       ! file with transmission coefficients
! Variables for nuclides
!   AA              ! mass number of residual nucleus
!   invexist        ! logical to state necessity of new inverse cross section calc.
!   ZZ              ! charge number of residual nucleus
! Variables for ECIS
!   flagecisinp     ! flag for existence of ecis input file
!
! *** Declaration of local data
!
  implicit none
  logical :: lexist ! logical to determine existence
  integer :: A      ! mass number of target nucleus
  integer :: Ncomp  ! neutron number index for compound nucleus
  integer :: Z      ! charge number of target nucleus
  integer :: Zcomp  ! proton number index for compound nucleus
!
! ************************** ECIS calculation **************************
!
! All transmission coefficients calculated by ECIS will be written on a file trZZZAAA, where ZZZ and AAA are the charge and
! mass number in (i3.3) format. The reaction cross sections will be written to csZZZAAA.
!
  transfile = 'tr000000     '
  csfile = 'cs000000     '
  Z = ZZ(Zcomp, Ncomp, 0)
  A = AA(Zcomp, Ncomp, 0)
  if (Z < 10) then
    write(transfile(5:5), '(i1)') Z
    write(csfile(5:5), '(i1)') Z
  endif
  if (Z >= 10 .and. Z < 100) then
    write(transfile(4:5), '(i2)') Z
    write(csfile(4:5), '(i2)') Z
  endif
  if (Z >= 100) then
    write(transfile(3:5), '(i3)') Z
    write(csfile(3:5), '(i3)') Z
  endif
  if (A < 10) then
    write(transfile(8:8), '(i1)') A
    write(csfile(8:8), '(i1)') A
  endif
  if (A >= 10 .and. A < 100) then
    write(transfile(7:8), '(i2)') A
    write(csfile(7:8), '(i2)') A
  endif
  if (A >= 100) then
    write(transfile(6:8), '(i3)') A
    write(csfile(6:8), '(i3)') A
  endif
!
! Calculate transmission coefficients and inverse reaction cross sections.
!
! inverseecis: subroutine for ECIS calculation for outgoing particles and energy grid
!
  if ( .not. invexist(Zcomp, Ncomp)) call inverseecis(Zcomp, Ncomp)
!
! Modification 5/5/11 by Kevin Kelley
! If the user has specified 'eciscalc n' and a cs/tr file pair is not present, go ahead and call the inverseecis subroutine
! to produce them. This is a little more user friendly than terminating 3/4 of the way through a calculation.
! This may save a lot of time in cases with 'optmodall y' when many cross sections and transmission coefficients are calculated,
! and which could then be transferred from directory to directory.
!
  flagecisinp = .true.
  if ( .not. flageciscalc) then
    inquire(file = csfile, exist = lexist)
    if (lexist) inquire(file = transfile, exist = lexist)
    if ( .not. lexist) then
      flageciscalc = .true.
      call inverseecis(Zcomp, Ncomp)
      flageciscalc = .false.
    endif
    invexist(Zcomp, Ncomp) = .true.
  endif
! end modification 5/5/11 by Kevin Kelley
!
! Read transmission coefficients and inverse reaction cross sections.
!
! inverseread: subroutine to read ECIS results for outgoing particles and energy grid
! inversenorm: subroutine for normalization of reaction cross sections and transmission coefficients
! inverseout : subroutine for reaction output for outgoing channels
!
  if (flagecisinp .and. invexist(Zcomp, Ncomp)) call inverseread(Zcomp, Ncomp)
  call inversenorm(Zcomp, Ncomp)
  if (flaginverse) call inverseout(Zcomp, Ncomp)
  return
end subroutine inverse
! Copyright A.J. Koning 2021
