subroutine resonancepar(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : S- and P-wave resonance parameters
!
! Author    : Arjan Koning
!
! 2025-12-13: Original code
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
  character(len=4)  :: D0ext       ! help variable
  character(len=6)  :: reschar     ! help variable
  character(len=132):: resfile     ! file with residual production cross sections
  integer           :: A           ! mass number of target nucleus
  integer           :: L           ! orbital angular momentum
  integer           :: iz          ! charge number from table
  integer           :: ia          ! mass number from table
  integer           :: istat       ! logical for file access
  integer           :: Nix         ! neutron number index for residual nucleus
  integer           :: Nrrf        ! number of resonances
  integer           :: Z           ! charge number of target nucleus
  integer           :: Zix         ! charge number index for residual nucleus
  real(sgl)         :: Df          ! help variable
  real(sgl)         :: dDf         ! help variable
  real(sgl)         :: Sf          ! help variable
  real(sgl)         :: dSf         ! help variable
  real(sgl)         :: Rf          ! help variable
  real(sgl)         :: dRf         ! help variable
  real(sgl)         :: dgamgamf    ! uncertainty in gamgam
  real(sgl)         :: gamgamf     ! experimental total radiative width in eV
  real(sgl)         :: D0glob
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
      read(2, '(4x, 2i4, 8es15.6, i4)', iostat = istat) ia, L, Df, dDf, Sf, dSf, gamgamf, dgamgamf, Rf, dRf, Nrrf
      if (istat == -1) exit
      if (A == ia) then
        if (L == 0) then
          if (dSf /= 0. .and. S0(Zix, Nix) == 0.) dS0(Zix, Nix) = dSf
          if (Sf /= 0. .and. S0(Zix, Nix) == 0.) S0(Zix, Nix) = Sf
          if (dDf /= 0. .and. D0(Zix, Nix) == 0.) dD0(Zix, Nix) = dDf * 1000.
          if (Df /= 0. .and. D0(Zix, Nix) == 0.) D0(Zix, Nix) = Df * 1000.
          if (dgamgamf /= 0. .and. gamgam(Zix, Nix) == 0.) dgamgam(Zix, Nix) = dgamgamf
          if (gamgamf /= 0. .and. gamgam(Zix, Nix) == 0.) gamgam(Zix, Nix) = gamgamf
          if (dRf /= 0. .and. Rscat(Zix, Nix) == 0.) dRscat(Zix, Nix) = dRf
          if (Rf /= 0. .and. Rscat(Zix, Nix) == 0.) Rscat(Zix, Nix) = Rf
          if (Nrrf /= 0) Nrr(Zix, Nix) = Nrrf
          if (Zix == 0 .and. Nix == 0) then
            if (Nrrf > 0 .and. D0(Zix, Nix) > 0.) Eavres = 0.5 * ((Nrrf - 1) * D0(Zix, Nix)) * 1.e-6
          endif
        else
          if (dDf /= 0..and.D1r(Zix, Nix) == 0.) dD1r(Zix, Nix) = dDf * 1000.
          if (Df /= 0..and.D1r(Zix, Nix) == 0.) D1r(Zix, Nix) = Df * 1000.
          if (dSf /= 0. .and. S1(Zix, Nix) == 0.) dS1(Zix, Nix) = dSf
          if (Sf /= 0. .and. S1(Zix, Nix) == 0.) S1(Zix, Nix) = Sf
        endif
      endif
    enddo
    close (unit = 2)
  endif
  gamgam(Zix, Nix) = gamgamadjust(Zix, Nix) * gamgam(Zix, Nix)
!
! 3. Read global D0 values from theory
!
  if (ldmodel(Zix, Nix) == 1 .and. .not.flagcol(Zix, Nix)) D0ext='ld1n'
  if (ldmodel(Zix, Nix) == 1 .and. flagcol(Zix, Nix)) D0ext='ld1y'
  if (ldmodel(Zix, Nix) == 2 .and. .not.flagcol(Zix, Nix)) D0ext='ld2n'
  if (ldmodel(Zix, Nix) == 2 .and. flagcol(Zix, Nix)) D0ext='ld2y'
  if (ldmodel(Zix, Nix) == 3 .and. .not.flagcol(Zix, Nix)) D0ext='ld3n'
  if (ldmodel(Zix, Nix) == 3 .and. flagcol(Zix, Nix)) D0ext='ld3y'
  if (ldmodel(Zix, Nix) == 4) D0ext='ld4'
  if (ldmodel(Zix, Nix) == 5) D0ext='ld5'
  if (ldmodel(Zix, Nix) == 6) D0ext='ld6'
  if (ldmodel(Zix, Nix) == 7) D0ext='ld7'
  resfile = trim(path)//'resonances/D0global.'//D0ext
  inquire (file = resfile, exist = lexist)
  if (lexist) then
    open (unit = 2, file = resfile, status = 'old', iostat = istat)
    if (istat /= 0) call read_error(resfile, istat)
    do
      read(2, *, iostat = istat) iz, ia, D0glob
      if (istat == -1) exit
      if (Z == iz .and. A == ia) then
        D0global(Zix, Nix) = D0glob
        dD0global(Zix, Nix) = D0glob
        exit
      endif
    enddo
  endif
  return
end subroutine resonancepar
! Copyright A.J. Koning 2025
