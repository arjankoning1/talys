function density(Zix, Nix, Eex, Rspin, parity, ibar, ldmod)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Level density
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
!   sgl         ! single precision kind
!   dbl         ! double precision kind
! Variables for level density
!   ctable      ! constant to adjust tabulated level densities
!   ldadjust    ! logical for energy - dependent level density adjustment
!   ptable      ! constant to adjust tabulated level densities
! Constants
!   pardis      ! parity distribution
! Variables for level density
!   edens       ! energy grid for tabulated level densities
!   Edensmax    ! maximum energy on level density table
!   ldexist     ! flag for existence of level density table
!   ldtable     ! level density from table
!   nendens     ! number of energies for level density grid
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: key           ! keyword
  integer           :: ibar          ! fission barrier
  integer           :: jj            ! counter
  integer           :: ldmod         ! level density model
  integer           :: nex2          ! counter
  integer           :: Nix           ! neutron number index for residual nucleus
  integer           :: parity        ! parity
  integer           :: Zix           ! charge number index for residual nucleus
  real(sgl)         :: ald           ! level density parameter
  real(sgl)         :: ct            ! constant to adjust tabulated level densities
  real(sgl)         :: Eex           ! excitation energy
  real(sgl)         :: Eshift        ! shifted excitation energy
  real(sgl)         :: expo          ! help variable
  real(sgl)         :: factor        ! multiplication factor
  real(sgl)         :: ignatyuk      ! function for energy dependent level density parameter a
  real(sgl)         :: pt            ! constant to adjust tabulated level densities
  real(sgl)         :: Rspin         ! residual spin
  real(sgl)         :: spindis       ! Wigner spin distribution
  real(sgl)         :: spincut       ! spin cutoff factor
  real(sgl)         :: sc            ! spin cutoff factor
  real(dbl)         :: cctable       ! constant to adjust tabulated level densities
  real(dbl)         :: density       ! level density
  real(dbl)         :: densitytot    ! total level density
  real(dbl)         :: eb            ! help variable
  real(dbl)         :: ee            ! energy
  real(dbl)         :: ldb           ! level density
  real(dbl)         :: lde           ! level density
  real(dbl)         :: ldtab         ! tabulated level density
!
! ************************** Level density *****************************
!
! ignatyuk  : function for energy dependent level density parameter a
!
! 1. Gilbert and Cameron
! 2. Back-shifted Fermi gas
! 3. Superfluid model
!
  density = 0.
  if (Eex < 0.) goto 100
  if (ldmod <= 3 .or. .not. ldexist(Zix, Nix, ibar)) then
    ald = ignatyuk(Zix, Nix, Eex, ibar)
    sc = spincut(Zix, Nix, ald, Eex, ibar, 0)
    density = densitytot(Zix, Nix, Eex, ibar, ldmod) * pardis * spindis(sc, Rspin)
  else
!
! 4. Tabulated level densities
!
! adjust       : subroutine for energy-dependent parameter adjustment
! locate       : subroutine to find value in ordered table
!
    jj = min(numJ - 1, int(Rspin))
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
!
! Interpolate from tables
!
    if (Eshift <= 0.) goto 100
    if (Eshift <= Edensmax(Zix, Nix)) then
      call locate(edens, 0, nendens(Zix, Nix), Eshift, nex2)
      eb = edens(nex2)
      ee = edens(nex2 + 1)
      ldb = ldtable(Zix, Nix, nex2, jj, parity, ibar)
      lde = ldtable(Zix, Nix, nex2 + 1, jj, parity, ibar)
    else
      eb = edens(nendens(Zix, Nix) - 1)
      ee = edens(nendens(Zix, Nix))
      ldb = ldtable(Zix, Nix, nendens(Zix, Nix) - 1, jj, parity, ibar)
      lde = ldtable(Zix, Nix, nendens(Zix, Nix), jj, parity, ibar)
    endif
    if (ldb <= 1..or.lde <= 1.) then
      ldtab = ldb + (Eshift - eb) / (ee - eb) * (lde - ldb)
    else
      ldtab = exp(log(ldb) + (Eshift - eb) / (ee - eb) * (log(lde) - log(ldb)))
    endif
    expo = min(ct * sqrt(Eshift), 80.)
    cctable = exp(dble(expo))
    density = cctable * ldtab
  endif
  100 density = max(density, 1.d-30)
  return
end function density
! Copyright A.J. Koning 2021
