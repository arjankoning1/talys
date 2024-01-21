function gammaxs(Zcomp, Ncomp, Egamma)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Gamma ray cross sections
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
! Variables for gamma rays
!   gammax      ! number of l - values for gamma multipolarity
! Variables for energy grid
!   Einc        ! incident energy in MeV
!  Variables for gamma-ray strength functions
!   kgr         ! constant for gamma - ray strength functio
!
! *** Declaration of local data
!
  implicit none
  integer   :: irad          ! variable to indicate M(=0) or E(=1) radiation
  integer   :: l             ! multipolarity
  integer   :: Ncomp         ! neutron number index for compound nucleus
  integer   :: Zcomp         ! proton number index for compound nucleus
  real(sgl) :: Egamma        ! gamma energy
  real(sgl) :: fstrength     ! gamma ray strength function
  real(sgl) :: gammaxs       ! function for gamma ray cross sections
  real(sgl) :: quasideuteron ! Quasi-deuteron function of Chadwick and Oblozinsky
  real(sgl) :: xsgdr         ! photo-absorption cross section from GDR part
  real(sgl) :: xsqd          ! photo-absorption cross section from QD part
!
! ************** Calculate photo-absorption cross section **************
!
! fstrength    : gamma ray strength function
! kgr          : constant for gamma-ray strength function
! quasideuteron: function for quasi-deuteron cross section
!
! 1. GDR part
!
  xsgdr = 0.
  do irad = 0, 1
    do l = 1, gammax
      xsgdr = xsgdr + fstrength(Zcomp, Ncomp, Einc, Egamma, irad, l) / kgr(l) * Egamma
    enddo
  enddo
!
! 2. QD part
!
  xsqd = quasideuteron(Egamma)
!
! Total absorption cross section
!
  gammaxs = xsgdr + xsqd
  return
end function gammaxs
! Copyright A.J. Koning 2021
