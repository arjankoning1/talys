subroutine resonancepar(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : S-wave resonance parameters
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
! Variables for level density
!   D0              ! s - wave resonance spacing in eV
! Variables for gamma rays
!   gamgam          ! total radiative width in eV
!   gamgamadjust    ! adjustable factor for radiative parameter
!   gnorm           ! gamma normalization factor
! Variables for nuclides
!   AA              ! mass number of residual nucleus
!   ZZ              ! charge number of residual nucleus
! Constants
!   nuc             ! symbol of nucleus
! Variables for files
!   path            ! directory containing files to be read
!  Variables for gamma-ray strength functions
!   gamkopecky      ! radiative width in eV by spline fit of Kopecky
! Variables for resonance parameters
!   dD0             ! uncertainty in D0
!   dgamgam         ! uncertainty in gamgam
!   Eavres          ! average resonance energy
!   Nrr             ! number of resonances
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist      ! logical to determine existence
  character(len=6)  :: reschar     ! help variable
  character(len=132):: resfile     ! file with residual production cross sections
  integer           :: A           ! mass number of target nucleus
  integer           :: ia          ! mass number from abundance table
  integer           :: istat       ! logical for file access
  integer           :: Nix         ! neutron number index for residual nucleus
  integer           :: Nrrf        ! number of resonances
  integer           :: Z           ! charge number of target nucleus
  integer           :: Zix         ! charge number index for residual nucleus
  real(sgl)         :: D0f         ! help variable
  real(sgl)         :: dD0f        ! help variable
  real(sgl)         :: dgamgamf    ! uncertainty in gamgam
  real(sgl)         :: gamgamf     ! experimental total radiative width in eV
!
! ********** Resonance spacings and total radiative widths *************
!
! Read experimental values from resonance parameter file.
! Resonance parameters from the table can always be overruled by a value given in the input file.
!
! 1. Inquire whether file is present
!
  Eavres = 0.01
  Z = ZZ(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  reschar = trim(nuc(Z))//'.res'
  resfile = trim(path)//'resonances/'//reschar
  inquire (file = resfile, exist = lexist)
  if (lexist) then
    open (unit = 2, file = resfile, status = 'old', iostat = istat)
    if (istat /= 0) call read_error(resfile, istat)
!
! 2. Search for the isotope under consideration and read information
!
    do
      read(2, '(4x, i4, 2e9.2, 10x, 2f9.5, i4)', iostat = istat) ia, D0f, dD0f, gamgamf, dgamgamf, Nrrf
      if (istat == -1) exit
      if (A == ia) then
        if (dD0f /= 0..and.D0(Zix, Nix) == 0.) dD0(Zix, Nix) = dD0f * 1000.
        if (D0f /= 0..and.D0(Zix, Nix) == 0.) D0(Zix, Nix) = D0f * 1000.
        if (dgamgamf /= 0..and.gamgam(Zix, Nix) == 0.) dgamgam(Zix, Nix) = dgamgamf
        if (gamgamf /= 0..and.gamgam(Zix, Nix) == 0.) gamgam(Zix, Nix) = gamgamf
        if (Nrrf /= 0) Nrr(Zix, Nix) = Nrrf
        if (Zix == 0 .and. Nix == 0) then
          if (Nrrf > 0 .and. D0(Zix, Nix) > 0.) Eavres = 0.5 * ((Nrrf - 1) * D0(Zix, Nix)) * 1.e-6
        endif
        exit
      endif
    enddo
    close (unit = 2)
  endif
  gamgam(Zix, Nix) = gamgamadjust(Zix, Nix) * gamgam(Zix, Nix)
  return
end subroutine resonancepar
! Copyright A.J. Koning 2021
