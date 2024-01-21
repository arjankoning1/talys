subroutine msdout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of multi-step direct cross sections
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
!   flagddx         ! flag for output of double - differential cross sections
! Variables for numerics
!   nanglecont      ! number of angles for continuum
! Variables for main input
!   k0              ! index of incident particle
! Variables for energy grid
!   anglecont       ! angle in degrees for continuum
!   ebegin          ! first energy point of energy grid
!   egrid           ! outgoing energy grid
! Variables for energies
!   eend            ! last energy point of energy grid
! Variables for nuclides
!   parskip         ! logical to skip outgoing particle
! Constants
!   parname         ! name of particle
! Constants
!   parsym          ! symbol of particle
! Variables for MSD
!   maxmsd          ! number of MSD steps
!   msdall          ! total multi - step direct cross section
!   msdstep         ! continuum n - step direct cross section
!   msdstepad       ! continuum n - step direct angular distribution
!   msdstepint      ! n - step direct cross section integrated over energy
!   msdstepintad    ! n - step direct angular distribution integrated over energy
!   msdsum          ! multi - step direct cross section summed over steps and inte
!   msdtot          ! multi - step direct cross section summed over steps
!   msdtotad        ! multi - step direct angular distribution summed over steps
!   msdtotintad     ! multi - step direct angular distribution summed over steps a
!
! *** Declaration of local data
!
  implicit none
  integer :: iang ! running variable for angle
  integer :: nen  ! energy counter
  integer :: ns   ! counter
  integer :: type ! particle type
!
! ************** Angle-integrated multi-step cross sections ************
!
  write(*, '(/" ++++++++++ MULTI-STEP DIRECT MODEL ++++++++++")')
  write(*, '(/" 1. Total multi-step direct cross sections"/)')
  write(*, '(" Step ", 2(5x, a8)/)') (parname(type), type = 1, 2)
  do ns = 1, maxmsd
    write(*, '(1x, i3, 3x, 2(1x, es12.5))') ns, (msdstepint(type, ns), type = 1, 2)
  enddo
  write(*, '(/" Total ", 2(1x, es12.5))') (msdsum(type), type = 1, 2)
  write(*, '(/" Total MSD cross section:", f12.5/)') msdall
  write(*, '(" 2. Multi-step direct spectra")')
  do type = 1, 2
    if (parskip(type)) cycle
    write(*, '(/, " Angle-integrated (", a1, ",", a1, ") MSD spectra"/)') parsym(k0), parsym(type)
    write(*, '("   Energy    Total", 4(5x, i2, "-step")/)') (ns, ns = 1, 4)
    do nen = ebegin(type), eend(type)
      write(*, '(1x, f8.3, 5es12.5)') egrid(nen), msdtot(type, nen), (msdstep(type, ns, nen), ns = 1, 4)
    enddo
  enddo
!
! ******************* Multi-step angular distributions *****************
!
  if ( .not. flagddx) return
  write(*, '(/" 3. Multi-step direct angular distributions")')
  do type = 1, 2
    if (parskip(type)) cycle
    do nen = ebegin(type), eend(type)
      write(*, '(/, " (", a1, ",", a1, ") MSD angular distribution for E-out= ", f8.3, /)') parsym(k0), parsym(type), egrid(nen)
      write(*, '(" Angle    Total", 4(5x, i2, "-step")/)') (ns, ns = 1, 4)
      do iang = 0, nanglecont
        write(*, '(1x, f5.1, 5es12.5)') anglecont(iang), msdtotad(type, nen, iang), (msdstepad(type, ns, nen, iang), ns = 1, 4)
      enddo
    enddo
    write(*, '(/, " (", a1, ",", a1, ") MSD angular distribution integrated over energy"/)') parsym(k0), parsym(type)
    write(*, '(" Angle    Total", 4(5x, i2, "-step")/)') (ns, ns = 1, 4)
    do iang = 0, nanglecont
      write(*, '(1x, f5.1, 5es12.5)') anglecont(iang), msdtotintad(type, iang), (msdstepintad(type, ns, iang), ns = 1, 4)
    enddo
  enddo
  return
end subroutine msdout
! Copyright A.J. Koning 2021
