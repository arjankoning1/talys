subroutine tdirbarrier(Zcomp, Ncomp, J2, parity, ibar, ibar2, trfis, rhof, Eex)
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
!   sgl          ! single precision kind
!   dbl          ! double precision kind
! Variables for energy grid, level densities and transmission coefficients
!   eintfis         ! excitation energy for fission
!   nbintfis        ! number of bins
!   rhofis          ! integrated level density corresponding to tfisA
! Variables for fission parameters
!   efistrrot       ! energy of rotational transition states
!   fecont          ! start of continuum energy
!   jfistrrot       ! spin of rotational transition states
!   nfistrrot       ! number of rotational transition states for barrier
!   pfistrrot       ! parity of rotational transition states
! Variables for fission
!   flagfispartdamp ! flag for fission partial damping
!
! *** Declaration of local data
!
  implicit none
!
! +---------------------------------------------------------------------
! | Author: Stephane Hilaire and Marieke Duijvestijn
! | Date  : October 31, 2019
! | Task  : Fission transmission coefficient for one barrier
! +---------------------------------------------------------------------
!
! ****************** Declarations and common blocks ********************
!
  integer   :: i            ! counter
  integer   :: ibar         ! fission barrier
  integer   :: ibar2        !
  integer   :: ibar3        !
  integer   :: itr          ! counter
  integer   :: J            ! spin of level
  integer   :: J2           ! 2 * J
  integer   :: j2trans      ! help variable
  integer   :: Ncomp        ! neutron number index for compound nucleus
  integer   :: parity       ! parity
  integer   :: pitrans      ! help variable
  integer   :: Zcomp        ! proton number index for compound nucleus
  real(sgl) :: dE1          ! help variable
  real(sgl) :: dE2          ! help variable
  real(sgl) :: Eeff         ! help variable
  real(sgl) :: Eex          ! excitation energy
  real(sgl) :: elow         ! help variable
  real(sgl) :: emid         ! help variable
  real(sgl) :: etrans       ! help variable
  real(sgl) :: eup          ! help variable
  real(sgl) :: twkbint      ! WKB penetrability
  real(sgl) :: twkbphaseint !
  real(dbl) :: r1log        ! help variable
  real(dbl) :: r2log        ! help variable
  real(dbl) :: r3log        ! help variable
  real(dbl) :: rho          ! integrated level density
  real(dbl) :: rho1         ! help variable
  real(dbl) :: rho2         ! help variable
  real(dbl) :: rho3         ! help variable
  real(dbl) :: rhof         ! value for total level density
  real(dbl) :: trfis        ! transmission coefficient
  real(dbl) :: trfisone     ! help variable
  real(dbl) :: trfisonetwo  !
  real(dbl) :: trfisthree   !
  real(dbl) :: trfistwo     !
  external twkbphaseint
!
! ********** Fission transmission coefficient for one barrier **********
!
  J = J2 / 2
  trfis = 0.
  rhof = 0.
!
! Correct LDM barrier height with ground state shell correction
!
  if ( .not. flagfispartdamp) return
!
! 1. Discrete states
!
  if (abs(ibar - ibar2) == 2) then
    ibar3 = (ibar2 + ibar) / 2
  else
    ibar3 = ibar
  endif
  do itr = 1, nfistrrot(Zcomp, Ncomp, ibar3)
    etrans = efistrrot(Zcomp, Ncomp, ibar3, itr)
    if (Eex < etrans) cycle
    j2trans = int(2. * jfistrrot(Zcomp, Ncomp, ibar3, itr))
    pitrans = pfistrrot(Zcomp, Ncomp, ibar3, itr)
    Eeff = Eex - etrans
    if ((J2 == j2trans) .and. (parity == pitrans)) then
      if (abs(ibar - ibar2) == 1) then
        trfisone = twkbint(Eeff, ibar, Zcomp, Ncomp)
        trfistwo = twkbint(Eeff, ibar2, Zcomp, Ncomp)
        if (twkbphaseint(Eeff, ibar, Zcomp, Ncomp) > 0) then
          trfis = trfis + trfisone * trfistwo / (1 + (1. - trfisone) * (1. - trfistwo))
        else
          trfis = trfis + trfisone * trfistwo
        endif
      elseif (abs(ibar - ibar2) == 2) then
        trfisone = twkbint(Eeff, ibar, Zcomp, Ncomp)
          trfistwo = twkbint(Eeff, ibar2, Zcomp, Ncomp)
        ibar3 = (ibar2 + ibar) / 2
          trfisthree = twkbint(Eeff, ibar3, Zcomp, Ncomp)
        if (twkbphaseint(Eeff, ibar, Zcomp, Ncomp) > 0) then
          trfisonetwo = trfisone * trfisthree / (1 + (1. - trfisone) * (1. - trfisthree))
        else
          trfisonetwo = trfisone * trfisthree
        endif
        if (twkbphaseint(Eeff, ibar3, Zcomp, Ncomp) > 0) then
            trfis = trfis + trfisonetwo * trfistwo / (1 + (1. - trfisonetwo) * (1. - trfistwo))
        else
          trfis = trfis + trfisonetwo * trfistwo
        endif
      endif
      rhof = rhof + 1.
    endif
  enddo
!
! 2. Continuum
!
  if (abs(ibar - ibar2) == 2) then
    ibar3 = (ibar2 + ibar) / 2
  else
    ibar3 = ibar
  endif
  if (Eex >= fecont(Zcomp, Ncomp, ibar3)) then
    do i = 1, nbintfis(ibar3) - 2, 2
      elow = eintfis(i, ibar3)
      if (elow > Eex) cycle
      emid = min(eintfis(i + 1, ibar3), Eex)
      eup = min(eintfis(i + 2, ibar3), Eex)
      dE1 = emid - elow
      dE2 = eup - emid
      if (abs(ibar - ibar2) == 1) then
        rho1 = rhofis(i, J, parity, ibar3) * (1. + 1.d-10)
        rho2 = rhofis(i + 1, J, parity, ibar3)
        rho3 = rhofis(i + 2, J, parity, ibar3) * (1. + 1.d-10)
      else
        rho1 = min(rhofis(i, J, parity, ibar) * (1. + 1.d-10), rhofis(i, J, parity, ibar2) * (1. + 1.d-10))
        rho2 = min(rhofis(i + 1, J, parity, ibar), rhofis(i + 1, J, parity, ibar2))
        rho3 = min(rhofis(i + 2, J, parity, ibar) * (1. + 1.d-10), rhofis(i + 2, J, parity, ibar2) * (1. + 1.d-10))
      endif
      r1log = log(rho1)
      r2log = log(rho2)
      r3log = log(rho3)
      if (r2log /= r1log .and. r2log /= r3log) then
        rho = (rho1 - rho2) / (r1log - r2log) * dE1 + (rho2 - rho3) / (r2log - r3log) * dE2
      else
        rho = rho2 * (dE1 + dE2)
      endif
      Eeff = Eex - emid
      if (abs(ibar - ibar2) == 1) then
        trfisone = twkbint(Eeff, ibar, Zcomp, Ncomp)
        trfistwo = twkbint(Eeff, ibar2, Zcomp, Ncomp)
        if ( twkbphaseint(Eeff, ibar, Zcomp, Ncomp)  >  0 ) then
          trfis = trfis + rho * trfisone * trfistwo / (1 + (1. - trfisone) * (1. - trfistwo))
        else
          trfis = trfis + rho * trfisone * trfistwo
        endif
      elseif (abs(ibar - ibar2) == 2) then
        trfisone = twkbint(Eeff, ibar, Zcomp, Ncomp)
        trfistwo = twkbint(Eeff, ibar2, Zcomp, Ncomp)
        ibar3 = (ibar2 + ibar) / 2
          trfisthree = twkbint(Eeff, ibar3, Zcomp, Ncomp)
        if (twkbphaseint(Eeff, ibar, Zcomp, Ncomp) > 0) then
          trfisonetwo = rho * trfisone * trfisthree / (1 + (1. - trfisone) * (1. - trfisthree))
        else
          trfisonetwo = rho * trfisone * trfisthree
        endif
        if (twkbphaseint(Eeff, ibar3, Zcomp, Ncomp) > 0) then
          trfis = trfis + rho * trfisonetwo * trfistwo / (1 + (1. - trfisonetwo) * (1. - trfistwo))
        else
          trfis = trfis + rho * trfisonetwo * trfistwo
        endif
      endif
      rhof = rhof + rho
    enddo
  endif
  return
end subroutine tdirbarrier
! Copyright A.J. Koning 2021
