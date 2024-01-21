function spincut(Zix, Nix, ald, Eex, ibar, ipop)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Spin cutoff factor
!
! Author    : Arjan Koning and Stephane Hilaire
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl             ! single precision kind
! Variables for level density
!   alimit          ! asymptotic level density parameter
!   Exmatch         ! matching point for Ex
!   flagcolldamp    ! flag for damping of coll. effects in eff. level density (without explicit coll. enh.)
!   ldadjust        ! logical for energy - dependent level density adjustment
!   ldmodel         ! level density model
!   pair            ! pairing energy
!   Pshift          ! adjustable pairing shift
!   Rspincut        ! adjustable constant (global) for spin cutoff factor
!   s2adjust        ! adjustable constant (Z, A, barrier - dependent) for spin
!   spincutmodel    ! model for spin cutoff factor for ground state
! Variables for masses
!   beta2           ! deformation parameter
! Variables for deformation parameters
!   Irigid0         ! undeformed rigid body value of moment of inertia
! Variables for level density
!   aldcrit         ! critical level density parameter
!   delta           ! energy shift
!   Ediscrete       ! energy of middle of discrete level region
!   scutoffdisc     ! spin cutoff factor for discrete level region
!   Tcrit           ! critical temperature
!   Ucrit           ! critical U
! Variables for masses
!   S               ! separation energy
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: key          ! keyword
  integer           :: ibar         ! fission barrier
  integer           :: ipop         ! 0: normal level density 1: initial population
  integer           :: ldmod        ! level density model
  integer           :: Nix          ! neutron number index for residual nucleus
  integer           :: Zix          ! charge number index for residual nucleus
  real(sgl)         :: ald          ! level density parameter
  real(sgl)         :: aldm         ! level density parameter at matching energy
  real(sgl)         :: Ed           ! energy of middle of discrete level region
  real(sgl)         :: Eex          ! excitation energy
  real(sgl)         :: Em           ! intermediate energy of local adjustment
  real(sgl)         :: factor       ! multiplication factor
  real(sgl)         :: ignatyuk     ! function for energy dependent level density parameter a
  real(sgl)         :: Rs           ! sensitivity factor
  real(sgl)         :: s2           ! variable for final GOE calculation
  real(sgl)         :: s2d          ! help variable
  real(sgl)         :: s2m          ! help variable
  real(sgl)         :: scutconst    ! constant for spin cutoff factor
  real(sgl)         :: spincut      ! spin cutoff factor
  real(sgl)         :: U            ! excitation energy minus pairing energy
  real(sgl)         :: Umatch       ! Exmatch - pairing
!
! *********************** Spin cutoff parameter ************************
!
! Models for spin cutoff factor:
!
! spincutmodel 1: s.c. = I0 * (a/alimit) * sqrt (U/a) =
!   I0/alimit * sqrt (a.U)
! spincutmodel 2: s.c. = I0 * sqrt (U/a)
!   where I0 = 0.01389 * A^(5/3)
!
! We interpolate between the discrete level spin cutoff factor and the value at the matching energy.
!
! 1. Below the matching energy
!
! adjust      : subroutine for energy-dependent parameter adjustment cutoff parameter
! ignatyuk    : function for energy dependent level density parameter a
!
  if (ipop == 1) then
    key = 'rspincutff'
    call adjust(Eex, key, 0, 0, 0, 0, factor)
    Rs = factor * Rspincutff
  else
    key = 'rspincut'
    call adjust(Eex, key, 0, 0, 0, 0, factor)
    Rs = factor * Rspincut
  endif
  if (ldadjust(Zix, Nix)) then
    key = 's2adjust'
    call adjust(Eex, key, Zix, Nix, 0, ibar, factor)
    s2 = factor * s2adjust(Zix, Nix, ibar)
  else
    s2 = s2adjust(Zix, Nix, ibar)
  endif
  ldmod = ldmodel(Zix, Nix)
  spincut = 1.
  if (spincutmodel == 1) then
    scutconst = Rs * s2 * Irigid0(Zix, Nix) / alimit(Zix, Nix)
  else
    scutconst = Rs * s2 * Irigid0(Zix, Nix)
  endif
  Em = Exmatch(Zix, Nix, ibar)
  if (ldmod == 2 .or. ldmod >= 4) Em = S(Zix, Nix, 1)
  if (ldmod == 3) Em = Ucrit(Zix, Nix, ibar) - pair(Zix, Nix) - Pshift(Zix, Nix, ibar)
  if (Eex <= Em) then
    aldm = ignatyuk(Zix, Nix, Em, ibar)
    Umatch = Em - pair(Zix, Nix) - Pshift(Zix, Nix, ibar)
    if (Umatch > 0.) then
      if (spincutmodel == 1) then
        if (ldmod == 3) then
          s2m = scutconst * aldcrit(Zix, Nix, ibar) * Tcrit(Zix, Nix)
        else
          s2m = scutconst * sqrt(aldm * Umatch)
        endif
      else
        if (ldmod == 3) then
          s2m = scutconst * Tcrit(Zix, Nix)
        else
          s2m = scutconst * sqrt(Umatch / aldm)
        endif
      endif
    else
      s2m = scutoffdisc(Zix, Nix, ibar)
    endif
    if (flagcolldamp .and. ibar /= 0) s2m = (1. + beta2(Zix, Nix, ibar) / 3.) * s2m
    s2d = scutoffdisc(Zix, Nix, ibar)
    Ed = Ediscrete(Zix, Nix, ibar)
    spincut = s2d
    if (Em /= Ed .and. Eex > Ed) spincut = s2d+(Eex - Ed) / (Em - Ed) * (s2m - s2d)
  else
!
! 2. Above the matching energy
!
    U = Eex - delta(Zix, Nix, ibar)
    if (U > 0.) then
      if (spincutmodel == 1) then
        spincut = scutconst * sqrt(ald * U)
      else
        spincut = scutconst * sqrt(U / ald)
      endif
    else
      spincut = scutoffdisc(Zix, Nix, ibar)
    endif
  endif
  if (flagcolldamp .and. ibar /= 0) spincut = (1. + beta2(Zix, Nix, ibar) / 3.) * spincut
  spincut = max(scutoffdisc(Zix, Nix, ibar), spincut)
  return
end function spincut
! Copyright A.J. Koning 2021
