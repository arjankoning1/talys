subroutine incidentnorm
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Normalization of reaction cross sections and transmission
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
! All global variables
!   numl          ! number of l values
! Variables for direct reactions
!   flagsys       ! flag for reaction cross section from systematics
! Variables for main input
!   Atarget       ! mass number of target nucleus
!   k0            ! index of incident particle
!   Ztarget       ! charge number of target nucleus
! Variables for energy grid
!   Einc          ! incident energy in MeV
! Variables for incident channel
!   Tjlinc        ! transm. coeff. as a function of spin and l for inc. channel
!   Tlinc         ! transm. coeff. as a function of l for incident channel
!   xselasinc     ! total elastic cross section (neutrons only) for inc. channel
!   xsoptinc      ! optical model reaction cross section for incident channel
!   xsreacinc     ! reaction cross section for incident channel
! Variables for inverse channel data
!   threshnorm    ! normalization factor at trheshold
! Constants
!   parA          ! mass number of particle
!   parZ          ! charge number of particle
!
! *** Declaration of local data
!
  implicit none
  integer   :: ispin    ! spin index
  integer   :: l        ! multipolarity
  real(sgl) :: enuc     ! incident energy in MeV per nucleon
  real(sgl) :: norm     ! normalization factor
  real(sgl) :: tripathi ! function for semi-empirical reaction cross section of
  real(sgl) :: xs       ! help variable
!
! ************ Normalization with semi-empirical results ***************
!
! tripathi  : function for semi-empirical reaction cross section of Tripathi et al.
!
! The normalization is only performed if the option for semi-empirical reaction cross sections is enabled.
!
  if ( .not. flagsys(k0)) return
  if (xsoptinc == 0.) return
  enuc = Einc / parA(k0)
  xs = tripathi(parZ(k0), parA(k0), Ztarget, Atarget, enuc)
  if (xs == 0.) then
    norm = threshnorm(k0)
    xs = xsoptinc * threshnorm(k0)
  else
    norm = xs / xsoptinc
  endif
  xsreacinc = xs
  if (k0 == 1) xselasinc = xselasinc + xsoptinc - xs
  do l = 0, numl
    Tlinc(l) = Tlinc(l) * norm
    do ispin = - 1, 1
      Tjlinc(ispin, l) = Tjlinc(ispin, l) * norm
    enddo
  enddo
  return
end subroutine incidentnorm
! Copyright A.J. Koning 2021
