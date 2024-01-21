subroutine tgamma(Zcomp, Ncomp)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Photon transmission coefficients
!
! Author    : Stephane Hilaire and Arjan Koning
!
! 2021-12-30: Original code
! 2022-04-04: Removed gamma normalization
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl             ! single precision kind
! Variables for gamma rays
!   gammax          ! number of l - values for gamma multipolarity
! Variables for energy grid
!   ebegin          ! first energy point of energy grid
!   egrid           ! outgoing energy grid
!   Einc            ! incident energy in MeV
! Variables for energies
!   eend            ! last energy point of energy grid
! Variables for inverse channel data
!   Tjl             ! transmission coefficient per particle, energy, spin and l - value
!   xsreac          ! reaction cross section
!  Variables for gamma-ray strength functions
!   lmax            ! maximal l - value for transmission coefficients
! Constants
!   twopi           ! 2 * pi
!
! *** Declaration of local data
!
  implicit none
  integer           :: irad         ! variable to indicate M(=0) or E(=1) radiation
  integer           :: l            ! multipolarity
  integer           :: Ncomp        ! neutron number index for compound nucleus
  integer           :: nen          ! energy counter
  integer           :: Zcomp        ! proton number index for compound nucleus
  real(sgl)         :: Egamma       ! gamma energy
  real(sgl)         :: fstrength    ! gamma ray strength function
  real(sgl)         :: gammaxs      ! function for gamma ray cross sections
!
! ************** Normalization of transmission coefficients ***********
!
! fstrength   : gamma ray strength function
!
  do nen = ebegin(0), eend(0)
    Egamma = egrid(nen)
    lmax(0, nen) = gammax
    do l = 1, gammax
      do irad = 0, 1
        Tjl(0, nen, irad, l) = twopi * (Egamma **(2 * l + 1)) * fstrength(Zcomp, Ncomp, Einc, Egamma, irad, l)
    enddo
  enddo
!
! Photo-absorption cross section
!
! gammaxs: function for gamma ray cross sections
!
    xsreac(0, nen) = gammaxs(Zcomp, Ncomp, Egamma)
  enddo
  return
end subroutine tgamma
! Copyright A.J. Koning 2021
