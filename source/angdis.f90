subroutine angdis
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculation of angular distributions for discrete states
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
!   sgl            ! single precision kind
! All global variables
!   numang         ! maximum number of angles
!   numJ           ! maximum J - value
! Variables for basic reaction
!   flagrecoil     ! flag for calculation of recoils
! Variables for numerics
!   nangle         ! number of angles
! Variables for main input
!   k0             ! index of incident particle
!   Ltarget        ! excited level of target
! Variables for energy grid
!   angle          ! angle in degrees
! Variables for energies
!   flagcompang    ! flag for compound angular distribution calculation
! Variables for binary reactions
!   xscompdisc     ! compound cross section for discrete state
!   xscompel       ! compound elastic cross section
!   xselastot      ! total elastic cross section (shape + compound)
!   xspopex0       ! binary population cross section
! Variables for incident channel
!   cleg           ! compound nucleus Legendre coefficient
!   directad       ! direct angular distribution
!   dleg           ! direct reaction Legendre coefficient
! Variables for angular distributions
!   cleg0          ! Legendre coefficient normalized to the first one
!   compad         ! compound angular distribution
!   discad         ! discrete state angular distribution
!   tleg           ! total Legendre coefficient
!   tlegnor        ! total Legendre coefficient normalized to 1
! Variables to normalize compound nucleus cross section
!   J2end          ! end of J summation
! Variables for nuclides
!   parskip        ! logical to skip outgoing particle
! Constants
!   deg2rad        ! conversion factor for degrees to radians
!   fourpi         ! 4. * pi
!   parN           ! neutron number of particle
!   parZ           ! charge number of particle
! Variables for level density
!   Nlast          ! last discrete level
!
! *** Declaration of local data
!
  implicit none
  integer   :: iang                    ! running variable for angle
  integer   :: LL                      ! counter for l value
  integer   :: nex                     ! excitation energy bin of compound nucleus
  integer   :: type                    ! particle type
  real(sgl) :: ang                     ! angle
  real(sgl) :: leg(0:2*numJ, 0:numang) ! Legendre polynomial
  real(sgl) :: plegendre               ! function for calculation of Legendre polynomial
  real(sgl) :: x                       ! help variable
  real(sgl) :: xs                      ! help variable
!
! ************* Create arrays with Legendre polynomials ****************
!
! plegendre: function for calculation of Legendre polynomial
!
  do iang = 0, nangle
    ang = angle(iang) * deg2rad
    x = cos(ang)
    do LL = 0, J2end, 2
      leg(LL, iang) = plegendre(LL, x)
    enddo
  enddo
!
! ******** Calculation of compound and total angular distribution ******
!
  do type = 0, 6
    if (parskip(type)) cycle
    do nex = 0, Nlast(parZ(type), parN(type), 0)
!
! For some unknown reason, compound inelastic inelastic scattering to 0+ states can have a L=2 Legendre coefficient that is larger
! than the L=0 Legendre coefficient.
! This only happens for high incident energies, typically above 10 MeV, i.e. for cases that compound inelastic scattering to
! individual states is negligible.
! Nevertheless, to avoid complaints by ENDF6 checking codes we put the L=2 Legendre coefficient equal to the L=0 coefficient
! for these rare cases.
!
      if (type /= k0 .or. nex > 0) cleg(type, nex, 2) = min(cleg(type, nex, 2), cleg(type, nex, 0))
      if (flagcompang .and. (type == k0 .or. cleg(type, nex, 0) /= 0.)) then
        do LL = 0, J2end, 2
          do iang = 0, nangle
            compad(type, nex, iang) = compad(type, nex, iang) + (2 * LL + 1) * cleg(type, nex, LL) * leg(LL, iang)
          enddo
        enddo
      else
        if (k0 == type .and. nex == Ltarget) then
          xs = xscompel
        else
          xs = xscompdisc(type, nex)
        endif
        do iang = 0, nangle
          compad(type, nex, iang) = xs / fourpi
        enddo
      endif
      do iang = 0, nangle
        discad(type, nex, iang) = directad(type, nex, iang) + compad(type, nex, iang)
      enddo
    enddo
  enddo
!
! The nuclear term is removed from charged-particle elastic scattering,
! so the nuclear + interference term remains (for data libraries).
!
  if (k0 > 1) then
    do iang = 0, nangle
      if (ruth(iang) > 0.) elasni(iang) = discad(k0, 0, iang) * (1. - 1./ruth(iang))
    enddo
  endif
!
! ************ Total Legendre coefficients and normalization ***********
!
  do type = 0, 6
    if (parskip(type)) cycle
    do nex = 0, Nlast(parZ(type), parN(type), 0)
      do LL = 0, J2end
        tleg(type, nex, LL) = dleg(type, nex, LL) + cleg(type, nex, LL)
        if (xspopex0(type, nex) >= 1.e-20) tlegnor(type, nex, LL) = tleg(type, nex, LL) / xspopex0(type, nex)
      enddo
    enddo
  enddo
  if (k0 == 1) then
    do LL = 0, J2end
      tlegnor(k0, 0, LL) = tleg(k0, 0, LL) / xselastot
    enddo
  endif
!
! Normalization to first Legendre coefficient (for ENDF-6 files)
!
  do type = 0, 6
    if (parskip(type)) cycle
    do nex = 0, Nlast(parZ(type), parN(type), 0)
      if (tlegnor(type, nex, 0) == 0.) cycle
      do LL = 0, J2end
        cleg0(type, nex, LL) = tlegnor(type, nex, LL) / tlegnor(type, nex, 0)
      enddo
    enddo
  enddo
!
! ***************************** Recoils ********************************
!
! angdisrecoil: subroutine for recoil angular distributions for discrete states
!
  if (flagrecoil) call angdisrecoil
  return
end subroutine angdis
! Copyright A.J. Koning 2021
