subroutine pseudo_resonance
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Compound nucleus pseudo-resonances
!
! Author    : Arjan Koning
!
! 2024-05-08: Original code
! 2026-05-06: Revised normalization of pseudo-resonance grid
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl            ! single precision kind
! All global variables
!   numlev2        ! maximum number of levels
!   numresgrid     ! number of energies on resonance grid
! Variables for input energies
!   Emaxpseudores ! max lab incident energy for pseudo-resonances
!   eninc          ! incident energy in MeV
!   Ninc           ! number of incident energies
! Variables for levels
!   edis           ! energy of level
!   jdis           ! spin of level
!   parlev         ! parity of level
! Variables for Q values
!   S              ! separation energy
! Constants
!   k0             ! index of incident particle
!   parN           ! neutron number of particle
!   parZ           ! charge number of particle
! Variables for masses
!   specmass       ! specific mass for residual nucleus
! Variables for pseudo-resonance treatment
!   Eresgrid       ! energies on resonance grid
!   pseudoresfade  ! fade width to return smoothly to optical model
!   pseudoreswidth ! Lorentzian full width
!   resgrid        ! resonance grid
!
! *** Declaration of local data
!
  implicit none
  integer           :: i                    ! energy counter
  integer           :: igrid                ! resonance-grid index
  integer           :: j                    ! energy counter
  integer           :: k                    ! resonance counter
  integer           :: nen                  ! incident energy counter
  integer           :: Nres                 ! number of resonances
  real(sgl)         :: active               ! active normalization weight
  real(sgl)         :: Ecut                 ! cutoff energy for pseudo-resonances
  real(sgl)         :: fgr(0:numresgrid)    ! unnormalised fluctuation strength
  real(sgl)         :: Eres(0:numlev2)      ! resonance energies
  real(sgl)         :: Ec                   ! compound excitation energy
  real(sgl)         :: Ecm                  ! incident energy in the center-of-mass system
  real(sgl)         :: Er                   ! resonance energy relative to separation threshold
  real(sgl)         :: E                    ! energy difference
  real(sgl)         :: degrid               ! resonance-grid energy step
  real(sgl)         :: Emax                 ! maximum energy of resonance grid
  real(sgl)         :: peak                 ! strength of fluctuations
  real(sgl)         :: weight               ! resonance weight
  real(sgl)         :: fac1                 ! 1 / pi
  real(sgl)         :: hgam                 ! half resonance width
  real(sgl)         :: lorentz              ! Lorentzian strength
  real(sgl)         :: fsum                 ! summed fluctuation strength
  real(sgl)         :: fave                 ! average fluctuation strength
  real(sgl)         :: favegrid(0:numresgrid) ! local average fluctuation strength
  real(sgl)         :: fade                 ! fade width in center-of-mass energy
  real(sgl)         :: gate                 ! damping factor to optical model
  real(sgl)         :: normwidth            ! width for local envelope normalization
!
! ******************** Default nuclear levels **************************
!
! The pseudo-resonance grid is a multiplicative factor applied later in
! incidentnorm. Therefore its average value should be one. The parameter
! peak damps or enhances the fluctuation around unity:
!
!   peak = 1 : full fluctuation strength
!   peak = 0 : no pseudo-resonance structure, resgrid = 1
!
  peak = 1.
!
! Initialise arrays
!
  Eres = 0.
  fgr = 0.
  Eresgrid = 0.
  resgrid = 1.
!
! Define resonance grid
!
  Emax = 10.
  degrid = Emax / real(numresgrid, sgl)
  do i = 0, numresgrid
    Eresgrid(i) = real(i, sgl) * degrid
  enddo
!
! Collect discrete compound-nucleus levels above the incident-channel
! separation energy. At this point in the initialization, Q has not yet
! been assigned, so use S(0, 0, k0) directly.
!
  Nres = 0
  do i = 1, nlevmax2(0, 0)
    Ec = edis(0, 0, i)
    Er = Ec - S(0, 0, k0)
    if (Er > 0.) then
      if (Nres < numlev2) then
        Nres = Nres + 1
        Eres(Nres) = Er
      endif
    endif
  enddo
!
! If no levels are available above threshold, leave resgrid equal to one.
!
  if (Nres == 0) return
  if (Ninc == 0) return
!
! Determine the energy above which the optical-model result should take
! over. The default cutoff is the last stored pseudo-resonance energy;
! this prevents the tail of an incomplete level list from suppressing the
! high-energy optical-model reaction cross section. A positive
! Emaxpseudores input value overrides this automatic cutoff. Input
! energies are laboratory energies, while the resonance grid is in the
! center-of-mass system.
!
  if (Emaxpseudores > 0.) then
    Ecut = Emaxpseudores * specmass(parZ(k0), parN(k0), k0)
  else
    Ecut = Eres(Nres)
  endif
  fade = pseudoresfade * specmass(parZ(k0), parN(k0), k0)
!
! Construct Lorentzian fluctuation strength on the resonance grid.
! The Lorentzian is area-normalised:
!
!   L(E) = 1/pi * hgam / ((E-Er)**2 + hgam**2)
!
! Give each selected discrete level equal strength. A 1/k weight makes
! the first alpha-unbound levels dominate the whole excitation function
! and suppresses the higher resonances that should provide the structure.
!
  fac1 = 1. / pi
  hgam = 0.5 * pseudoreswidth
!
  do i = 0, numresgrid
    do k = 1, Nres
      E = Eresgrid(i) - Eres(k)
      lorentz = fac1 * hgam / (E * E + hgam * hgam)
      weight = 1.
      fgr(i) = fgr(i) + lorentz * weight
    enddo
  enddo
!
! Remove the slow envelope caused by the finite and non-uniform level
! list. Without this local normalization, the pseudo-resonance factor can
! have a long trend, so an arithmetic mean of one still gives a biased
! cross-section average when the optical cross section rises with energy.
!
  normwidth = max(0.5, 10. * pseudoreswidth)
  do i = 0, numresgrid
    fsum = 0.
    active = 0.
    do j = 0, numresgrid
      if (abs(Eresgrid(j) - Eresgrid(i)) <= normwidth) then
        fsum = fsum + fgr(j)
        active = active + 1.
      endif
    enddo
    if (active > 0.) then
      favegrid(i) = fsum / active
    else
      favegrid(i) = 0.
    endif
  enddo
  do i = 0, numresgrid
    if (favegrid(i) > 0.) fgr(i) = fgr(i) / favegrid(i)
  enddo
!
! Normalise the fluctuation strength to average unity on the actual
! incident-energy grid. This is the essential step: incidentnorm later
! multiplies xsreacinc, Tlinc and Tjlinc by resgrid sampled at those
! incident energies. Hence <resgrid> on that same grid should be one,
! otherwise the pseudo-resonance option changes the smooth TALYS average.
!
  fsum = 0.
  active = 0.
  do nen = 1, Ninc
    Ecm = eninc(nen) * specmass(parZ(k0), parN(k0), k0)
    call locate(Eresgrid, 0, numresgrid, Ecm, igrid)
    gate = 1.
    if (Ecm > Ecut - fade) gate = (Ecut - Ecm) / fade
    gate = max(0., min(1., gate))
    if (gate > 0.) then
      fsum = fsum + gate * fgr(igrid)
      active = active + gate
    endif
  enddo
  if (active > 0.) then
    fave = fsum / active
  else
    resgrid = 1.
    return
  endif
!
! Convert the Lorentzian strength into a multiplicative fluctuation
! factor with mean value one:
!
!   resgrid = (1 - peak) + peak * fgr / <fgr>
!
! For peak = 1, this reduces to resgrid = fgr / <fgr>.
!
  if (fave > 0.) then
    do i = 0, numresgrid
      gate = 1.
      if (Eresgrid(i) > Ecut - fade) gate = (Ecut - Eresgrid(i)) / fade
      gate = max(0., min(1., gate))
      resgrid(i) = 1. + gate * peak * (fgr(i) / fave - 1.)
      resgrid(i) = max(resgrid(i), 0.)
    enddo
  else
    resgrid = 1.
  endif
!
  return
end subroutine pseudo_resonance
! Copyright A.J. Koning 2026
