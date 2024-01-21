subroutine giant
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Giant resonance contribution
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
! 2022-11-14: Current version
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl              ! single precision kind
! All global variables
!   numen            ! maximum number of outgoing energies
!   numlev2          ! maximum number of levels
! Variables for direct reactions
!   elwidth          ! width of elastic peak in MeV
! Variables for output
!   flagddx          ! flag for output of double - differential cross sections
! Variables for numerics
!   nangle           ! number of angles
!   nanglecont       ! number of angles for continuum
! Variables for main input
!   k0               ! index of incident particle
! Variables for energy grid
!   anglecont        ! continuum angle
!   deltaE           ! energy bin around outgoing energies
!   ebegin           ! first energy point of energy grid
!   egrid            ! outgoing energy grid
! Variables for energies
!   eend             ! last energy point of energy grid
!   eninccm          ! center - of - mass incident energy in MeV
!   eoutdis          ! outgoing energy of discrete state reaction
! Variables for incident channel
!   directad         ! direct angular distribution
!   xscollconttot    ! total collective cross section in the continuum
!   xsdirdisc        ! direct cross section for discrete state direct cross section
!   xsgr             ! total smoothed giant resonance cross section
!   xsgrsum          ! sum over giant resonance cross sections
!   xsgrtot          ! total smoothed giant resonance cross section
! Variables for nuclides
!   Nindex           ! neutron number index for residual nucleus
!   Zindex           ! charge number index for residual nucleus
! Variables for deformation parameters
!   betagr           ! deformation parameter for giant resonance
!   Egrcoll          ! energy of giant resonance
!   Ggrcoll          ! width of giant resonance
! Variables for giant resonances
!   collcontad       ! collective angular distribution in the continuum
!   eoutgr           ! emission energy
!   grcollad         ! giant resonance angular distribution
!   xscollcont       ! total collective cross section in the continuum
!   xsgrad           ! smoothed giant resonance angular distribution
!   xsgrcoll         ! giant resonance cross section
!   xsgrstate        ! smoothed giant resonance cross section per state
! Constants
!   sqrttwopi        ! sqrt(2. * pi)
! Variables for deformation parameters
!   deform           ! deformation parameter
! Variables for level density
!   Nlast            ! last discrete level
!
! *** Declaration of local data
!
  implicit none
  integer   :: i                   ! counter
  integer   :: iang                ! running variable for angle
  integer   :: iangd(0:nanglecont) ! running variable for angle
  integer   :: l                   ! multipolarity
  integer   :: J                   ! spin
  integer   :: parity              ! parity
  integer   :: nen                 ! energy counter
  integer   :: Nix                 ! neutron number index for residual nucleus
  integer   :: Zix                 ! charge number index for residual nucleus
  real(sgl) :: diswidth            ! width of discrete level peak
  real(sgl) :: dang                ! delta angle
  real(sgl) :: edist               ! help variable
  real(sgl) :: fac1                ! help variable
  real(sgl) :: fac2                ! help variable
  real(sgl) :: term                ! help variable
  real(sgl) :: gauss(0:numen)      ! Gaussian contribution
  real(sgl) :: grwidth             ! width of giant resonance
  real(sgl) :: sumgauss            ! sum over Gaussians (for normalization)
  real(sgl) :: weight              ! weight of emission energy bin
  real(sgl) :: wscale              ! scaling factor for giant resonance to Gaussian width
!
! ************* Smearing of giant resonances into spectra **************
!
  wscale = 0.42
  do l = 0, 3
    do i = 1, 2
      if (betagr(l, i) == 0.) cycle
      eoutgr(k0, l, i) = eninccm - Egrcoll(l, i)
      grwidth = Ggrcoll(l, i) * wscale
      fac1 = 1. / (grwidth * sqrttwopi)
      fac2 = 1. / (2. * grwidth **2)
      sumgauss = 0.
      do nen = ebegin(k0), eend(k0)
        gauss(nen) = 0.
        edist = abs(eoutgr(k0, l, i) - egrid(nen))
        if (edist > 5. * grwidth) cycle
        gauss(nen) = fac1 * exp( - (edist **2) * fac2)
        sumgauss = sumgauss + gauss(nen)
      enddo
      if (sumgauss /= 0.) then
        do nen = ebegin(k0), eend(k0)
          weight = gauss(nen) / sumgauss
          xsgrstate(k0, l, i, nen) = weight * xsgrcoll(k0, l, i)
          xsgr(k0, nen) = xsgr(k0, nen) + xsgrstate(k0, l, i, nen)
          if (flagddx) then
            do iang = 0, nanglecont
              xsgrad(k0, nen, iang) = xsgrad(k0, nen, iang) + weight * grcollad(k0, l, i, iang)
            enddo
          endif
        enddo
      endif
      xsgrtot(k0) = xsgrtot(k0) + xsgrcoll(k0, l, i)
    enddo
  enddo
  xsgrsum = xsgrtot(k0)
!
! *********** Other collective contributions to the continuum **********
!
  if (xscollconttot(k0) == 0.) return
  xsgrtot(k0) = xsgrtot(k0) + xscollconttot(k0)
  xsgrsum = xsgrsum + xscollconttot(k0)
  if (flagddx) then
    dang = 180. / nangle
    do iang = 0, nanglecont
      iangd(iang) = int(anglecont(iang) / dang)
    enddo
  endif
  Zix = Zindex(0, 0, k0)
  Nix = Nindex(0, 0, k0)
  do i = Nlast(Zix, Nix, 0) + 1, numlev2
    if (deform(Zix, Nix, i) == 0.) cycle
    if (eoutdis(k0, i) <= 0.) cycle
    diswidth = elwidth * (eoutdis(k0, i) / eninccm) **1.5
    fac1 = 1. / (diswidth * sqrttwopi)
    fac2 = 1. / (2. * diswidth **2)
    sumgauss = 0.
    do nen = ebegin(k0), eend(k0)
      gauss(nen) = 0.
      edist = abs(egrid(nen) - eoutdis(k0, i))
      if (edist > 5. * diswidth) cycle
      gauss(nen) = fac1 * exp( - (edist **2) * fac2)
      sumgauss = sumgauss + gauss(nen)
    enddo
    if (sumgauss /= 0.) then
      J = int(jdis(Zix,Nix,i))
      parity = parlev(Zix,Nix,i)
      do nen = ebegin(k0), eend(k0)
        weight = gauss(nen) / sumgauss / deltaE(nen)
        term = weight * xsdirdisc(k0, i)
        xscollcont(k0, nen) = xscollcont(k0, nen) + term
        xscollcontJP(k0, J, parity, nen) = xscollcontJP(k0, J, parity, nen) + term
        if (flagddx) then
          do iang = 0, nanglecont
            collcontad(k0, nen, iang) = collcontad(k0, nen, iang) + weight * directad(k0, i, iangd(iang))
          enddo
        endif
      enddo
    endif
  enddo
!
! Add collective contribution to giant resonance results
!
  do nen = ebegin(k0), eend(k0)
    xsgr(k0, nen) = xsgr(k0, nen) + xscollcont(k0, nen)
    if (flagddx) then
      do iang = 0, nanglecont
        xsgrad(k0, nen, iang) = xsgrad(k0, nen, iang) + collcontad(k0, nen, iang)
      enddo
    endif
  enddo
  return
end subroutine giant
! Copyright A.J. Koning 2021
