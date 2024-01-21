subroutine thermalxs
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Cross sections at thermal energies
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
!   sgl             ! single precision kind
! All global variables
!   numisom       ! number of isomers
! Variables for compound reactions
!   xsalphatherm    ! thermal (n, a) cross section
!   xscaptherm      ! thermal capture cross section
!   xsptherm        ! thermal (n, p) cross section
! Variables for level density
!   alev            ! level density parameter
!   pair            ! pairing energy
! Variables for nuclides
!   AA              ! mass number of residual nucleus
!   ZZ              ! charge number of residual nucleus
! Constants
!   nuc             ! symbol of nucleus
! Variables for files
!   path            ! directory containing files to be read
! Variables for masses
!   S               ! separation energy
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist     ! logical to determine existence
  character(len=8)  :: thchar     ! help variable
  character(len=132):: thfile     ! thermal cross section file
  integer           :: A          ! mass number of target nucleus
  integer           :: i          ! counter
  integer           :: ia         ! mass number from abundance table
  integer           :: isoR       ! isomer of residual
  integer           :: isoT       ! isomer of target
  integer           :: istat      ! error code
  integer           :: Z          ! charge number of target nucleus
  real(sgl)         :: ald        ! level density parameter
  real(sgl)         :: Spair      ! help variable
  real(sgl)         :: sumalpha   ! help variable
  real(sgl)         :: sumcap     ! help variable
  real(sgl)         :: sump       ! help variable
  real(sgl)         :: xs         ! help variable
  real(sgl)         :: xsalpha    ! (n,a) cross section
  real(sgl)         :: xsp        ! (n,p) cross section
!
! ********** Resonance spacings and total radiative widths *************
!
! Read experimental values from thermal cross section file.
! Values from the table can always be overruled by values given in the input file.
!
! 1. Inquire whether file is present
!
  Z = ZZ(0, 0, 1)
  A = AA(0, 0, 1)
  thchar = trim(nuc(Z))//'.ther'
  thfile = trim(path)//'thermal/'//thchar
  inquire (file = thfile, exist = lexist)
  if (lexist) then
    open (unit = 2, file = thfile, status = 'old', iostat = istat)
    if (istat /= 0) call read_error(thfile, istat)
!
! 2. Search for the isotope under consideration and read information
!
    xs = 0.
    xsp = 0.
    xsalpha = 0.
    do
      read(2, '(4x, 3i4, 3(e9.2, 9x))', iostat = istat) ia, isoT, isoR, xs, xsp, xsalpha
      if (istat == -1) exit
      if (istat /= 0) call read_error(thfile, istat)
      if (A == ia) then
        if (xs /= 0.) xscaptherm(isoR) = xs
        if (xsp /= 0.) xsptherm(isoR) = xsp
        if (xsalpha /= 0.) xsalphatherm(isoR) = xsalpha
      endif
    enddo
    close (unit = 2)
    sumcap = 0.
    sump = 0.
    sumalpha = 0.
    do i = 0, numisom
      sumcap = sumcap + xscaptherm(i)
      sump = sump + xsptherm(i)
      sumalpha = sumalpha + xsalphatherm(i)
    enddo
    if (xscaptherm(-1) == 0.) xscaptherm(-1) = sumcap
    if (xsptherm(-1) == 0.) xsptherm(-1) = sump
    if (xsalphatherm(-1) == 0.) xsalphatherm(-1) = sumalpha
  endif
!
! 2. Systematics
!
! Kopecky's value for (n,gamma) cross section at thermal energy.
! J. Kopecky, M.G. Delfini, H.A.J. van der Kamp and D. Nierop
!
  if (xscaptherm(-1) == 0.) then
    ald = alev(0, 0)
    Spair = S(0, 0, 1) - pair(0, 0)
    Spair = max(Spair, 1.)
    xscaptherm(-1) = 1.5e-3 * (ald * Spair) **3.5
  endif
  return
end subroutine thermalxs
! Copyright A.J. Koning 2021
