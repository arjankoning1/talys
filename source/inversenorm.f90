subroutine inversenorm(Zcomp, Ncomp)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Normalization and extrapolation of reaction cross sections
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
!   numen         ! maximum number of outgoing energies
!   numl          ! number of l values
! Variables for direct reactions
!   flagsys       ! flag for reaction cross section from systematics
! Variables for energy grid
!   ebegin        ! first energy point of energy grid
!   eendmax       ! last energy point of energy grid
!   egrid         ! outgoing energy grid
! Variables for inverse channel data
!   threshnorm    ! normalization factor at threshold
!   Tjl           ! transmission coefficient per particle, energy, spin and l - value
!   Tl            ! transmission coefficients per particle, energy and l - value
!   xselas        ! total elastic cross section (shape + compound)
!   xsopt         ! optical model reaction cross section
!   xsreac        ! reaction cross section
! Variables for nuclides
!   AA            ! mass number of residual nucleus
!   parskip       ! logical to skip outgoing particle
!   ZZ            ! charge number of residual nucleus
! Constants
!   parA          ! mass number of particle
!   parZ          ! charge number of particle
!
! *** Declaration of local data
!
  implicit none
  integer   :: A                            ! mass number of target nucleus
  integer   :: ispin                        ! spin index
  integer   :: l                            ! multipolarity
  integer   :: Ncomp                        ! neutron number index for compound nucleus
  integer   :: nen                          ! energy counter
  integer   :: type                         ! particle type
  integer   :: Z                            ! charge number of target nucleus
  integer   :: Zcomp                        ! proton number index for compound nucleus
  real(sgl) :: enuc                         ! incident energy in MeV per nucleon
  real(sgl) :: norm                         ! normalization factor
  real(sgl) :: tripathi                     ! function for semi-empirical reaction cross section of
  real(sgl) :: xs                           ! help variable
  real(sgl) :: xsprev                       ! help variable
  real(sgl) :: xstripathi(0:numen)          ! help variable
!
! ************ Normalization with semi-empirical results ***************
!
! tripathi  : function for semi-empirical reaction cross section of Tripathi et al.
!
! The normalization is only performed if the option for semi-empirical reaction cross sections is enabled.
! The semi-empirical results have a too sharp cutoff at low energies.
! Therefore, for the lowest energies the optical model results are renormalized with the ratio at the threshold.
!
  do type = 1, 6
    if (parskip(type)) cycle
    if ( .not. flagsys(type)) cycle
    Z = ZZ(Zcomp, Ncomp, type)
    A = AA(Zcomp, Ncomp, type)
    xsprev = 0.
    do nen = ebegin(type), eendmax(type)
      if (xsopt(type, nen) == 0.) cycle
      enuc = egrid(nen) / parA(type)
      xs = tripathi(parZ(type), parA(type), Z, A, enuc)
      xstripathi(nen) = xs
      if (xsprev == 0 .and. xs /= 0.) threshnorm(type) = xs / xsopt(type, nen)
      xsprev = xs
    enddo
    do nen = ebegin(type), eendmax(type)
      if (xsopt(type, nen) == 0.) cycle
      xs = xstripathi(nen)
      if (xs == 0.) then
        norm = threshnorm(type)
        xs = xsopt(type, nen) * threshnorm(type)
      else
        norm = xs / xsopt(type, nen)
      endif
      xsreac(type, nen) = xs
      if (type == 1) xselas(type, nen) = xselas(type, nen) + xsopt(type, nen) - xs
      do l = 0, numl
        Tl(type, nen, l) = Tl(type, nen, l) * norm
        do ispin = - 1, 1
          Tjl(type, nen, ispin, l) = Tjl(type, nen, ispin, l) * norm
        enddo
      enddo
    enddo
  enddo
  return
end subroutine inversenorm
! Copyright A.J. Koning 2021
