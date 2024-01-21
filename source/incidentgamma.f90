subroutine incidentgamma
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Incident photons
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
! Variables for incident channel
!   lmaxinc     ! maximal l - value for transm. coeff. for incident channel
!   Tjlinc      ! transm. coeff. as a function of spin and l for inc. channel
!   xsreacinc   ! reaction cross section for incident channel
! Constants
!   twopi       ! 2 * pi
!
! *** Declaration of local data
!
  implicit none
  integer   :: irad          ! variable to indicate M(=0) or E(=1) radiation
  integer   :: l             ! multipolarity
  real(sgl) :: factor        ! multiplication factor
  real(sgl) :: fstrength     ! gamma ray strength function
  real(sgl) :: gammaxs       ! function for gamma ray cross sections
  real(sgl) :: quasideuteron ! Quasi-deuteron function of Chadwick and Oblozinsky
  real(sgl) :: Tgamma        ! gamma transmission coefficient
  real(sgl) :: xsgdr         ! photo-absorption cross section from GDR part
  real(sgl) :: xsqd          ! photo-absorption cross section from QD part
!
! **** Photo-absorption cross section and transmission coefficients ****
!
! gammaxs  : function for gamma ray cross sections
! fstrength: gamma ray strength function
!
! Note that we use the first index of Tjlinc for the radiation type, instead of the particle spin index.
!
  lmaxinc = gammax
  xsreacinc = gammaxs(0, 0, Einc)
  xsqd = quasideuteron(Einc)
  xsgdr = xsreacinc - xsqd
  if (xsgdr > 0.) then
    factor = xsreacinc / xsgdr
  else
    factor = 1.
  endif
  do irad = 0, 1
    do l = 1, gammax
      Tgamma = twopi * (Einc **(2 * l + 1)) * fstrength(0, 0, Einc, Einc, irad, l)
      Tjlinc(irad, l) = Tgamma * factor
    enddo
  enddo
  return
end subroutine incidentgamma
! Copyright A.J. Koning 2021
