function densitytot(Zix, Nix, Eex, ibar, ldmod)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Total level density
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
! Definition of single and double precision variables
!   sgl           ! single precision kind
!   dbl           ! double precision kind
! Variables for level density
!   ctable        ! constant to adjust tabulated level densities
!   ldadjust      ! logical for energy - dependent level density adjustment
!   ptable        ! constant to adjust tabulated level densities
! Variables for level density
!   delta         ! energy shift
!   edens         ! energy grid for tabulated level densities
!   Edensmax      ! maximum energy on level density table
!   ldexist       ! flag for existence of level density table
!   ldtottable    ! total level density per parity from table
!   nendens       ! number of energies for level density grid
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: key           ! keyword
  integer           :: ibar          ! fission barrier
  integer           :: ldmod         ! level density model
  integer           :: nex2          ! counter
  integer           :: Nix           ! neutron number index for residual nucleus
  integer           :: Zix           ! charge number index for residual nucleus
  real(sgl)         :: ald           ! level density parameter
  real(sgl)         :: ct            ! constant to adjust tabulated level densities
  real(sgl)         :: Eex           ! excitation energy
  real(sgl)         :: Eshift        ! shifted excitation energy
  real(sgl)         :: expo          ! help variable
  real(sgl)         :: factor        ! multiplication factor
  real(sgl)         :: ignatyuk      ! function for energy dependent level density parameter a
  real(sgl)         :: Kcoll         ! total collective enhancement
  real(sgl)         :: Krot          ! rotational enhancement factor
  real(sgl)         :: Kvib          ! vibrational enhancement factor
  real(sgl)         :: P             ! pairing energy
  real(sgl)         :: pt            ! constant to adjust tabulated level densities
  real(dbl)         :: bsfgmodel     ! level density
  real(dbl)         :: cctable       ! constant to adjust tabulated level densities
  real(dbl)         :: dens          ! total level density
  real(dbl)         :: densitytot    ! total level density
  real(dbl)         :: eb            ! help variable
  real(dbl)         :: ee            ! energy
  real(dbl)         :: gilcam        ! Gilbert-Cameron level density formula
  real(dbl)         :: ldb           ! level density
  real(dbl)         :: lde           ! level density
  real(dbl)         :: lldb          ! log of level density
  real(dbl)         :: llde          ! log of level density
  real(dbl)         :: superfluid    ! Superfluid model level density formula
!
! *********************** Total level density **************************
!
! adjust     : subroutine for energy-dependent parameter adjustment
!
  densitytot = 0.
  do
    if (Eex < 0.) exit
!
! Possible adjustment of final level densities
!
    if (ldadjust(Zix, Nix)) then
      key = 'ptable'
      call adjust(Eex, key, Zix, Nix, 0, ibar, factor)
      pt = ptable(Zix, Nix, ibar) + factor - 1.
      key = 'ctable'
      call adjust(Eex, key, Zix, Nix, 0, ibar, factor)
      ct = ctable(Zix, Nix, ibar) + factor - 1.
    else
      pt = ptable(Zix, Nix, ibar)
      ct = ctable(Zix, Nix, ibar)
    endif
    Eshift = Eex - pt
    if (Eshift <= 0.) exit
    if (ldmod <= 3 .or. .not. ldexist(Zix, Nix, ibar)) then
!
! ignatyuk  : function for energy dependent level density parameter a
! colenhance: subroutine for collective enhancement
!
! Kcoll will be determined.
!
      ald = ignatyuk(Zix, Nix, Eex, ibar)
      call colenhance(Zix, Nix, Eex, ald, ibar, Krot, Kvib, Kcoll)
      P = delta(Zix, Nix, ibar)
!
! 1. Gilbert and Cameron
!
      if (ldmod == 1 .or. (ldmod >= 4 .and. .not. ldexist(Zix, Nix, ibar))) dens = gilcam(Zix, Nix, ald, Eex, P, ibar)
!
! 2. Back-shifted Fermi gas
!
      if (ldmod == 2) dens = bsfgmodel(Zix, Nix, ald, Eex, P, ibar)
!
! 3. Superfluid model
!
      if (ldmod == 3) dens = superfluid(Zix, Nix, ald, Eex, P, ibar)
      dens = Kcoll * dens
    else
!
! 4. Tabulated level densities
!
! locate       : subroutine to find value in ordered table
!
      if (Eshift <= Edensmax(Zix, Nix)) then
        call locate(edens, 0, nendens(Zix, Nix), Eshift, nex2)
        eb = edens(nex2)
        ee = edens(nex2 + 1)
        ldb = ldtottable(Zix, Nix, nex2, ibar)
        lde = ldtottable(Zix, Nix, nex2 + 1, ibar)
      else
        eb = edens(nendens(Zix, Nix) - 1)
        ee = edens(nendens(Zix, Nix))
        ldb = ldtottable(Zix, Nix, nendens(Zix, Nix) - 1, ibar)
        lde = ldtottable(Zix, Nix, nendens(Zix, Nix), ibar)
      endif
      if (ldb > 1..and.lde > 1.) then
        lldb = log(ldb)
        llde = log(lde)
        dens = exp(lldb + (Eshift - eb) / (ee - eb) * (llde - lldb))
      else
        dens = ldb + (Eshift - eb) / (ee - eb) * (lde - ldb)
      endif
    endif
    expo = min(ct * sqrt(Eshift), 80.)
    cctable = exp(dble(expo))
    densitytot = cctable * dens
    exit
  enddo
  densitytot = max(densitytot, 1.d-30)
  return
end function densitytot
! Copyright A.J. Koning 2021
