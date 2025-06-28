subroutine t1barrier(Zcomp, Ncomp, J2, parity, ibar, trfis, rhof, Eex, iloop)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Fission transmission coefficient for one barrier
!
! Author    : Stephane Hilaire and Marieke Duijvestijn
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
!   dbl             ! double precision kind
! All global variables
!   numhill         ! maximum number of Hill - Wheeler points
! Variables for fission
!   fbaradjust      ! adjustable factor for fission parameters
!   fbarrier        ! height of fission barrier
!   fisadjust       ! logical for energy - dependent fission adjustment
!   fismodelx       ! fission model
!   fwidth          ! width of fission barrier
!   fwidthadjust    ! adjustable factor for fission parameters
! Variables for level density
!   deltaW          ! shell correction in nuclear mass
! Variables for energy grid, level densities and transmission coefficients
!   eintfis         ! excitation energy for fission
!   nbintfis        ! number of bins
!   rhofis          ! integrated level density corresponding to tfisA
! Variables for fission transmission coefficients
!   rhofisA         ! integrated level density corresponding to tfisA
!   tfisA           ! transmission coefficient for Hill - Wheeler magnitude
! Variables for nuclides
!   primary         ! flag to designate primary (binary) reaction
! Variables for fission parameters
!   efistrrot       ! energy of rotational transition states
!   fecont          ! start of continuum energy
!   jfistrrot       ! spin of rotational transition states
!   nfisbar         ! number of fission barrier parameters
!   nfistrrot       ! number of rotational transition states for barrier
!   pfistrrot       ! parity of rotational transition states
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: key         ! keyword
  integer           :: i           ! level
  integer           :: ibar        ! fission barrier
  integer           :: ihill       ! counter for Hill-Wheeler magnitude
  integer           :: iloop       ! loop counter
  integer           :: itr         ! counter
  integer           :: J           ! spin of level
  integer           :: J2          ! 2 * J
  integer           :: j2trans     ! help variable
  integer           :: Ncomp       ! neutron number index for compound nucleus
  integer           :: parity      ! parity
  integer           :: pitrans     ! help variable
  integer           :: Zcomp       ! proton number index for compound nucleus
  real(sgl)         :: bfis        ! barrier height
  real(sgl)         :: dE1         ! help variable
  real(sgl)         :: dE2         ! help variable
  real(sgl)         :: Eeff        ! help variable
  real(sgl)         :: Eex         ! excitation energy
  real(sgl)         :: elow        ! help variable
  real(sgl)         :: emid        ! help variable
  real(sgl)         :: etrans      ! help variable
  real(sgl)         :: eup         ! help variable
  real(sgl)         :: factor1     ! help variable
  real(sgl)         :: factor2     ! help variable
  real(sgl)         :: fbar        ! fission barrier
  real(sgl)         :: thill       ! Hill-Wheeler penetrability
  real(sgl)         :: twkbint     ! WKB penetrability
  real(sgl)         :: wfis        ! help variable
  real(dbl)         :: r1log       ! help variable
  real(dbl)         :: r2log       ! help variable
  real(dbl)         :: r3log       ! help variable
  real(dbl)         :: rho         ! integrated level density
  real(dbl)         :: rho1        ! help variable
  real(dbl)         :: rho2        ! help variable
  real(dbl)         :: rho3        ! help variable
  real(dbl)         :: rhof        ! value for total level density
  real(dbl)         :: rhotr       ! level density x transmission coefficient
  real(dbl)         :: trfis       ! transmission coefficient
  real(dbl)         :: trfisone    ! help variable
!
! ********** Fission transmission coefficient for one barrier **********
!
  J = J2 / 2
  trfis = 0.
  rhof = 0.
!
! Correct LDM barrier height with ground state shell correction
!
! adjust    : subroutine for energy-dependent parameter adjustment (default 1.)
!
  if (fisadjust(Zcomp, Ncomp)) then
    key = 'fisbar'
    call adjust(Eex, key, Zcomp, Ncomp, 0, ibar, factor1)
    key = 'fisbaradjust'
    call adjust(Eex, key, Zcomp, Ncomp, 0, ibar, factor2)
    fbar = factor1 * factor2 * fbarrier(Zcomp, Ncomp, ibar) * fbaradjust(Zcomp, Ncomp, ibar)
  else
    fbar = fbarrier(Zcomp, Ncomp, ibar)
  endif
  if ((fismodelx(Zcomp, Ncomp) >= 3) .or. ((fismodelx(Zcomp, Ncomp) < 3) .and. (nfisbar(Zcomp, Ncomp) == 1))) then
    bfis = fbar - (deltaW(Zcomp, Ncomp, 0) - deltaW(Zcomp, Ncomp, 1))
  else
    bfis = fbar
  endif
  if (fisadjust(Zcomp, Ncomp)) then
    key = 'fishw'
    call adjust(Eex, key, Zcomp, Ncomp, 0, ibar, factor1)
    key = 'fishwadjust'
    call adjust(Eex, key, Zcomp, Ncomp, 0, ibar, factor2)
    wfis = factor1 * factor2 * fwidth(Zcomp, Ncomp, ibar) * fwidthadjust(Zcomp, Ncomp, ibar)
  else
    wfis = fwidth(Zcomp, Ncomp, ibar)
  endif
!
! 1. Discrete states
!
  do itr = 1, nfistrrot(Zcomp, Ncomp, ibar)
    etrans = efistrrot(Zcomp, Ncomp, ibar, itr)
    if (Eex < etrans) cycle
    j2trans = int(2. * jfistrrot(Zcomp, Ncomp, ibar, itr))
    pitrans = pfistrrot(Zcomp, Ncomp, ibar, itr)
    Eeff = Eex - etrans
    if ((J2 == j2trans) .and. (parity == pitrans)) then
      if (fismodelx(Zcomp, Ncomp) >= 5) then
        trfisone = twkbint(Eeff, ibar, Zcomp, Ncomp)
      else
        trfisone = thill(Eeff, bfis, wfis)
      endif
      if ((ibar == 1) .and. primary .and. iloop == 2) then
        ihill = min(int(numhill * trfisone) + 1, numhill)
        tfisA(J, parity, ihill) = tfisA(J, parity, ihill) + trfisone
        tfisA(J, parity, 0) = tfisA(J, parity, 0) + trfisone
        rhofisA(J, parity, ihill) = rhofisA(J, parity, ihill) + 1.
      endif
      trfis = trfis + trfisone
      rhof = rhof + 1.
    endif
  enddo
!
! 2. Continuum
!
  if (Eex >= fecont(Zcomp, Ncomp, ibar)) then
    do i = 1, nbintfis(ibar) - 2, 2
      elow = eintfis(i, ibar)
      if (elow > Eex) cycle
      emid = min(eintfis(i + 1, ibar), Eex)
      eup = min(eintfis(i + 2, ibar), Eex)
      dE1 = emid - elow
      dE2 = eup - emid
      rho1 = rhofis(i, J, parity, ibar) * (1. + 1.d-10) + 1.d-10
      rho2 = rhofis(i + 1, J, parity, ibar) + 1.d-10
      rho3 = rhofis(i + 2, J, parity, ibar) * (1. + 1.d-10) + 1.d-10
      r1log = log(rho1)
      r2log = log(rho2)
      r3log = log(rho3)
      if (r2log /= r1log .and. r2log /= r3log) then
        rho = (rho1 - rho2) / (r1log - r2log) * dE1 + (rho2 - rho3) / (r2log - r3log) * dE2
      else
        rho = rho2 * (dE1 + dE2)
      endif
      Eeff = Eex - emid
      if (fismodelx(Zcomp, Ncomp) >= 5) then
        trfisone = twkbint(Eeff, ibar, Zcomp, Ncomp)
      else
        trfisone = thill(Eeff, bfis, wfis)
      endif
      rhotr = rho * trfisone
      trfis = trfis + rhotr
      rhof = rhof + rho
      if ((ibar == 1) .and. primary .and. iloop == 2) then
        ihill = min(int(numhill * trfisone) + 1, numhill)
        tfisA(J, parity, ihill) = tfisA(J, parity, ihill) + rhotr
        tfisA(J, parity, 0) = tfisA(J, parity, 0) + rhotr
        rhofisA(J, parity, ihill) = rhofisA(J, parity, ihill) + rho
      endif
    enddo
    if (iloop == 2) tfisA(J, parity, 0) = max(tfisA(J, parity, 0), 1.d-30)
  endif
  return
end subroutine t1barrier
! Copyright A.J. Koning 2021
