function aldmatch(Zix, Nix, Eex, ibar)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Effective level density parameter
!
! Author    : Stephane Hilaire, Marieke Duijvestijn and Arjan Koning
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
! Variables for level density
!   pair         ! pairing energy
! Variables for nuclides
!   AA           ! mass number of residual nucleus
! Constants
!   sqrttwopi    ! sqrt(2. * pi)
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: jmax=40    ! maximum j-value
  integer            :: A          ! mass number of target nucleus
  integer            :: ibar       ! fission barrier
  integer            :: j          ! counter
  integer            :: Nix        ! neutron number index for residual nucleus
  integer            :: Zix        ! charge number index for residual nucleus
  real(sgl)          :: ald1       ! boundary values for search
  real(sgl)          :: ald2       ! boundary values for search
  real(sgl)          :: aldacc     ! accuracy
  real(sgl)          :: aldmatch   ! function to determine effective level density parameter
  real(sgl)          :: aldmid     ! help variable
  real(sgl)          :: aldref     ! level density parameter
  real(sgl)          :: dald       ! difference in level density parameter
  real(sgl)          :: Eex        ! excitation energy
  real(sgl)          :: ignatyuk   ! function for energy dependent level density parameter a
  real(sgl)          :: Kcoll      ! total collective enhancement
  real(sgl)          :: Krot       ! rotational enhancement factor
  real(sgl)          :: Kvib       ! vibrational enhancement factor
  real(sgl)          :: rj         ! help variable
  real(sgl)          :: sc         ! help variable
  real(sgl)          :: rjbegin    ! help variable
  real(sgl)          :: sigma      ! help variable
  real(sgl)          :: sigmamid   ! help variable
  real(sgl)          :: spincut    ! spin cutoff factor
  real(sgl)          :: spindis    ! Wigner spin distribution
  real(dbl)          :: factor     ! multiplication factor
  real(dbl)          :: fdiff      ! difference in level density
  real(dbl)          :: fermi      ! function for Fermi gas level density formula
  real(dbl)          :: fmid       ! help variable
  real(dbl)          :: rhoref     ! help variable
  real(dbl)          :: rhosum     ! help variable
!
! ************ Search for effective level density parameter ************
!
! The effective level density parameter is obtained in three steps:
!
! 1. Create total level density in Fermi gas region
!
! ignatyuk     : function for energy dependent level density parameter a
! fermi        : function for Fermi gas level density formula
!
  A = AA(Zix, Nix, 0)
  rjbegin = 0.5 * mod(A, 2)
  rj = rjbegin - 1.
  rhosum = 0.
  aldref = ignatyuk(Zix, Nix, Eex, ibar)
  do
    rj = rj + 1.
    sc=spincut(Zix, Nix, aldref, Eex, ibar, 0)
    factor = (2. * rj + 1.) * fermi(Zix, Nix, aldref, Eex, pair(Zix, Nix), ibar) * spindis(sc, rj)
    rhosum = rhosum + factor
    if (factor < 0.00001) exit
  enddo
!
! 2. Apply a rotational enhancement to the total level density
!
! colenhance: subroutine for collective enhancement
!
  call colenhance(Zix, Nix, Eex, aldref, ibar, Krot, Kvib, Kcoll)
  rhoref = Kcoll * rhosum
!
! 3. Determine effective level density parameter by equating the rotational enhanced level density by a new effective
!    total level density.
!
! aldmatch    : function to determine effective level density parameter
!
  aldacc = 0.001
  ald1 = 0.5 * aldref
  ald2 = 2.0 * aldref
  sigma = sqrt(spincut(Zix, Nix, ald1, Eex, ibar, 0))
  fdiff = fermi(Zix, Nix, ald1, Eex, pair(Zix, Nix), ibar) * sigma * sqrttwopi - rhoref
  if (fdiff < 0.) then
    aldmatch = ald1
    dald = ald2 - ald1
  else
    aldmatch = ald2
    dald = ald1 - ald2
  endif
  do j = 1, jmax
    dald = dald * 0.5
    aldmid = aldmatch + dald
    sigmamid = sqrt(spincut(Zix, Nix, aldmid, Eex, ibar, 0))
    fmid = fermi(Zix, Nix, aldmid, Eex, pair(Zix, Nix), ibar) * sigmamid * sqrttwopi - rhoref
    if (fmid <= 0.) aldmatch = aldmid
    if (abs(dald) < aldacc .or. fmid == 0.) return
  enddo
  return
end function aldmatch
! Copyright A.J. Koning 2021
