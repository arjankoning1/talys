subroutine abundance
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Natural abundances
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
  use A1_error_handling_mod
!
! Definition of single and double precision variables
!   sgl          ! single precision kind
! Variables for main input
!   Starget    ! symbol of target nucleus
!   Ztarget    ! charge number of target nucleus
! Variables for abundance
!   abun         ! Natural abundance
!   isonum       ! number of isotopes
!   isotope      ! isotope number of residual nucleus
! Variables for reading input lines
!   inline    ! input line
!   nlines    ! number of input lines
! Variables for files
!   path         ! directory containing files to be read
! Constants
!   natstring    ! string extension for file names
!   nuc          ! symbol of nucleus
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist    ! logical to determine existence
  character(len=1)  :: ch        ! character
  character(len=7)  :: abchar    ! help variable
  character(len=132):: abfile    ! isotopic abundance file
  integer           :: i         ! level
  integer           :: i2        ! value
  integer           :: ia        ! mass number from abundance table
  integer           :: istat     ! error code
  real(sgl)         :: ab        ! isotopic abundance
  real(sgl)         :: abtot     ! summed abundances for normalization
!
! ****************************** Abundances ****************************
!
! Note that for non-natural elements we take the longest-lived isotope as default.
!
! 1. Isotopic abundances from user file
!
  abun = 0.
  isotope = 0
  do i = 1, nlines
    if (inline(i)(1:10) == 'abundance ') then
      do i2 = 11, 132
        ch = inline(i)(i2:i2)
        if (ch /= ' ') then
          read(inline(i)(i2:132), * , iostat = istat) abfile
          if (istat /= 0) call read_error(abfile, istat)
          goto 100
        endif
      enddo
    endif
  enddo
!
! 2. Isotopic abundances from abundance directory
!
  abchar = trim(nuc(Ztarget))//'.abun'
  abfile = trim(path)//'abundance/'//abchar
  inquire (file = abfile, exist = lexist)
  if ( .not. lexist) then
    write(*, '(" TALYS-error: No natural isotopes for this element, the mass keyword must be different from 0")')
    stop
  endif
!
! Read abundances from file
!
  100 open (unit = 2, file = abfile, status = 'old')
  i = 1
  do
    read(2, '(4x, i4, f11.6)', iostat = istat) ia, ab
    if (istat == -1) exit
    if (istat > 0) call read_error(abfile, istat)
    isotope(i) = ia
    abun(i) = 0.01 * ab
    i = i + 1
  enddo
  close (unit = 2)
  isonum = i - 1
!
! Normalize abundances to 1.
!
  abtot = 0.
  do i = 1, isonum
    abtot = abtot + abun(i)
  enddo
  do i = 1, isonum
    abun(i) = abun(i) / abtot
  enddo
  write(*, '(/" Calculation for multi-isotope case"/)')
  write(*, '("  Isotope Abundance"/)')
  do i = 1, isonum
    write(*, '(2x, i3, a2, f11.6)') isotope(i), Starget, abun(i)
  enddo
!
! Create file name extensions
!
  do i = 1, isonum
    write(natstring(i)(1:4), '(".", i3.3)') isotope(i)
  enddo
  return
end subroutine abundance
! Copyright A.J. Koning 2021
