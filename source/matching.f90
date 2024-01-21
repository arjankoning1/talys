subroutine matching(Zix, Nix, Exm, ibar)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Determine matching between temperature and Fermi-gas region
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
!   sgl       ! single precision kind
! Variables for level density
!   pair      ! pairing energy
! Variables for nuclides
!   AA        ! mass number of residual nucleus
! Variables for level density
!   Exmemp    ! empirical estimate for matching point for Ex
!
! *** Declaration of local data
!
  implicit none
  integer   :: A         ! mass number of target nucleus
  integer   :: ibar      ! fission barrier
  integer   :: Nb        ! number of solutions
  integer   :: Nix       ! neutron number index for residual nucleus
  integer   :: nseg      ! number of segments
  integer   :: Zix       ! charge number index for residual nucleus
  real(sgl) :: ald       ! level density parameter
  real(sgl) :: Exm       ! maximal attainable energy
  real(sgl) :: Exm1      ! matching energy
  real(sgl) :: Exm2      ! matching energy
  real(sgl) :: ignatyuk  ! function for energy dependent level density parameter a
  real(sgl) :: match     ! function to search for zero crossings of the function
  real(sgl) :: rtbis     ! function to search for zero crossings of the function
  real(sgl) :: x1        ! coordinates of intersection points inside the bin
  real(sgl) :: x2        ! coordinates of the 2nd summit of the triangle
  real(sgl) :: xacc      ! help variable
  real(sgl) :: xb1(2)    ! help variable
  real(sgl) :: xb2(2)    ! help variable
  external match
!
! ************************ Search for zeroes ***************************
!
! ignatyuk: function for energy dependent level density parameter a
! match   : function to search for zero crossings of the function
! zbrak   : function to bracket the function
! rtbis   : function to search for zero crossings of the function
!
  A = AA(Zix, Nix, 0)
  ald = ignatyuk(Zix, Nix, Exmemp, ibar)
  xacc = 0.0001
!
! Set possible region for solution
!
  x1 = max(2.25 / ald + pair(Zix, Nix), 0.) + 0.11
  x2 = 19. + 300. / A
  nseg = 100
  Nb = 2
  call zbrak(match, x1, x2, nseg, xb1, xb2, Nb)
!
! Look for 0,1 or 2 solutions
!
  if (Nb == 0) Exm = 0.
  if (Nb == 1) Exm = rtbis(match, xb1(1), xb2(1), xacc)
!
! If there are 2 solutions, we choose the one closest to the empirical expression.
!
  if (Nb == 2) then
    Exm1 = rtbis(match, xb1(1), xb2(1), xacc)
    Exm2 = rtbis(match, xb1(2), xb2(2), xacc)
    if (abs(Exm1 - Exmemp) > abs(Exm2 - Exmemp)) then
      Exm = Exm2
    else
      Exm = Exm1
    endif
  endif
!
! If the solution is unphysical, we choose the empirical expression.
!
  if (Exm < x1 .or. Exm > x2) Exm = 0.
  return
end subroutine matching
! Copyright A.J. Koning 2021
