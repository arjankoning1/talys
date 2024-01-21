function densitytotP(Zix, Nix, Eex, parity, ibar, ldmod)
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
! Constants
!   pardis           ! parity distribution
! Variables for level density
!   ctable        ! constant to adjust tabulated level densities
!   ldadjust      ! logical for energy - dependent level density adjustment
!   ptable        ! constant to adjust tabulated level densities
! Variables for level density
!   edens         ! energy grid for tabulated level densities
!   Edensmax      ! maximum energy on level density table
!   ldexist       ! flag for existence of level density table
!   ldtottableP   ! total level density per parity from table
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
  integer           :: parity        ! parity
  integer           :: Zix           ! charge number index for residual nucleus
  real(sgl)         :: ct            ! constant to adjust tabulated level densities
  real(sgl)         :: Eex           ! excitation energy
  real(sgl)         :: Eshift        ! shifted excitation energy
  real(sgl)         :: expo          ! help variable
  real(sgl)         :: factor        ! multiplication factor
  real(sgl)         :: pt            ! constant to adjust tabulated level densities
  real(dbl)         :: cctable       ! constant to adjust tabulated level densities
  real(dbl)         :: dens          ! total level density
  real(dbl)         :: densitytot    ! total level density
  real(dbl)         :: densitytotP   ! total level density per parity
  real(dbl)         :: eb            ! help variable
  real(dbl)         :: ee            ! energy
  real(dbl)         :: ldb           ! level density
  real(dbl)         :: lde           ! level density
  real(dbl)         :: lldb          ! log of level density
  real(dbl)         :: llde          ! log of level density
!
! *********************** Total level density **************************
!
! adjust     : subroutine for energy-dependent parameter adjustment
!
  densitytotP = 0.
  if (Eex >= 0.) then
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
    if (Eshift > 0.) then
      if (ldmod <= 3 .or. .not. ldexist(Zix, Nix, ibar)) then
        dens = densitytot(Zix, Nix, Eex, ibar, ldmod) * pardis
      else
!
! Tabulated level densities
!
! locate       : subroutine to find value in ordered table
!
        if (Eshift <= Edensmax(Zix, Nix)) then
          call locate(edens, 0, nendens(Zix, Nix), Eshift, nex2)
          eb = edens(nex2)
          ee = edens(nex2 + 1)
          ldb = ldtottableP(Zix, Nix, nex2, parity, ibar)
          lde = ldtottableP(Zix, Nix, nex2 + 1, parity, ibar)
        else
          eb = edens(nendens(Zix, Nix) - 1)
          ee = edens(nendens(Zix, Nix))
          ldb = ldtottableP(Zix, Nix, nendens(Zix, Nix) - 1, parity, ibar)
          lde = ldtottableP(Zix, Nix, nendens(Zix, Nix), parity, ibar)
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
      densitytotP = cctable * dens
    endif
  endif
  densitytotP = max(densitytotP, 1.d-30)
  return
end function densitytotP
! Copyright A.J. Koning 2021
