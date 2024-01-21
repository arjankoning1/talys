subroutine exgrid(Zcomp, Ncomp)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Set excitation energy grid
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
!   dbl            ! double precision kind
! All global variables
!   numex          ! maximum number of excitation energies
!   numJ           ! maximum J - value
! Variables for compound reactions
!   flagcomp       ! flag for compound angular distribution calculation
! Variables for basic parameters
!   flagequi       ! flag to use equidistant excitation bins instead of logarithmic
! Variables for main input
!   k0             ! index of incident particle
!   Ltarget        ! excited level of target
! Variables for level density
!   flagparity     ! flag for non - equal parity distribution
!   ldmodel        ! level density model
! Variables for OMP
!   flagomponly    ! flag to execute ONLY an optical model calculation
! Variables for energies
!   Etotal         ! total energy of compound system (target + projectile)
! Variables for energies
!   Ethresh        ! threshold incident energy for residual nucleus
!   nbins          ! number of continuum excitation energy bins
!   Qres           ! Q - value for residual nucleus
! Variables for excitation energy grid
!   deltaEx        ! excitation energy bin for population arrays
!   Ex             ! excitation energy
!   Exmax          ! maximum excitation energy for residual nucleus
!   Exmax0         ! maximum excitation energy for res. nucleus (incl. neg. en.)
!   maxex          ! maximum excitation energy bin for residual nucleus
!   maxJ           ! maximal J - value
!   nexmax         ! maximum excitation energy bin for residual nucleus
!   rhogrid        ! integrated level density
! Variables for nuclides
!   AA             ! mass number of residual nucleus
!   Nindex         ! neutron number index for residual nucleus
!   parskip        ! logical to skip outgoing particle
!   targetE        ! excitation energy of target
!   Zindex         ! charge number index for residual nucleus
! Constants
!   parN           ! neutron number of particle
!   parZ           ! charge number of particle
! Variables for levels
!   edis           ! energy of level
! Variables for level density
!   Nlast          ! last discrete level
! Variables for masses
!   S              ! separation energy
!   specmass       ! specific mass for residual nucleus
!
! *** Declaration of local data
!
  implicit none
  integer   :: A            ! mass number of target nucleus
  integer   :: Aix          ! mass number index for residual nucleus
  integer   :: i            ! counter
  integer   :: ipop         ! help varible
  integer   :: Ir           ! residual spin
  integer   :: ldmod        ! level density model
  integer   :: Ncomp        ! neutron number index for compound nucleus
  integer   :: Ndeep        ! neutron number index for lightest residual nucleus
  integer   :: nex          ! excitation energy bin of compound nucleus
  integer   :: nexbins      ! number of continuum excitation energy bins
  integer   :: Nix          ! neutron number index for residual nucleus
  integer   :: NL           ! last discrete level
  integer   :: Nmother      ! neutron number index for mother nucleus
  integer   :: odd          ! odd (1) or even (0) nucleus
  integer   :: Pbeg         ! begin and end of parity summation
  integer   :: Pprime       ! parity
  integer   :: type         ! particle type
  integer   :: Zcomp        ! proton number index for compound nucleus
  integer   :: Zdeep        ! charge number index for lightest residual nucleus
  integer   :: Zix          ! charge number index for residual nucleus
  integer   :: Zmother      ! charge number index for mother nucleus
  real(sgl) :: ald          ! level density parameter
  real(sgl) :: dEx          ! excitation energy bin for population arrays
  real(sgl) :: eb           ! help variable
  real(sgl) :: Edif         ! difference in energy
  real(sgl) :: ee           ! energy
  real(sgl) :: Eup(0:numex) ! upper energy of bin
  real(sgl) :: Ex1min       ! lower boundary of residual bin
  real(sgl) :: Ex1plus      ! upper boundary of residual bin
  real(sgl) :: Exout        ! excitation energy
  real(sgl) :: Rodd         ! term to determine integer or half-integer spins
  real(sgl) :: Rspin        ! residual spin
  real(sgl) :: spincut      ! spin cutoff factor
  real(dbl) :: density      ! level density
  real(dbl) :: r1log        ! help variable
  real(dbl) :: r2log        ! help variable
  real(dbl) :: r3log        ! help variable
  real(dbl) :: rho1         ! help variable
  real(dbl) :: rho2         ! help variable
  real(dbl) :: rho3         ! help variable
!
! ********************* Set maximum excitation energy ******************
!
! All possible routes to all residual nuclei are followed, to determine the maximum possible excitation energy for each nucleus,
! given the incident energy.
!
  if (flagomponly .and. .not. flagcomp) return
  Zdeep = Zcomp
  Ndeep = Ncomp
  do type = 1, 6
    if (parskip(type)) cycle
    Zix = Zindex(Zcomp, Ncomp, type)
    Nix = Nindex(Zcomp, Ncomp, type)
    Zdeep = max(Zdeep, Zix)
    Ndeep = max(Ndeep, Nix)
  enddo
  do Zix = Zcomp, Zdeep
    do Nix = 0, Ndeep
      if (Zix == Zcomp .and. Nix == Ncomp) cycle
      if (Exmax(Zix, Nix) /= 0.) cycle
      do type = 1, 6
        if (parskip(type)) cycle
        Zmother = Zix - parZ(type)
        if (Zmother < 0) cycle
        Nmother = Nix - parN(type)
        if (Nmother < 0) cycle
        Exmax0(Zix, Nix) = Exmax0(Zmother, Nmother) - S(Zmother, Nmother, type)
        Exmax(Zix, Nix) = max(Exmax0(Zix, Nix), 0.)
      enddo
    enddo
  enddo
!
! ******* Define excitation energy grid for residual nuclei ************
!
! The excitation energies are given by the array Ex.
! The first NL values of Ex correspond to the discrete level excitation energies of the residual nucleus specified by (Zix,Nix).
! The NL+1th value corresponds to the first continuum energy bin.
! The continuum bins have a width of deltaEx.
! The continuum part of the target and initial compound nucleus are divided into nbins equidistant energy bins where nbins was given
! in the input file or set by default.
! The continuum parts of the other residual nuclides are also divided in equidistant bins.
! For the first generation of nuclides (within 4 mass units of the initial compound nucleus), the number of bins is equal
! to that of target.
! For nuclides more than 8 mass units away, half the number of bins is chosen.
! For intermediate nuclides, an interpolated number is adopted.
! In sum, for each nucleus the excitation energy range is completely filled by equidistant bins.
! The bin widths thus gradually change.
!
! The Q-value for residual nuclides is determined, both for the ground state and for possible isomers.
!
  do type = 0, 6
    if (parskip(type)) cycle
    Zix = Zindex(Zcomp, Ncomp, type)
    Nix = Nindex(Zcomp, Ncomp, type)
    NL = Nlast(Zix, Nix, 0)
    if (maxex(Zix, Nix) /= 0) cycle
    do nex = 0, numex
      deltaEx(Zix, Nix, nex) = 0.
    enddo
    if (Qres(Zix, Nix, 0) == 0.) then
      Edif = Exmax0(Zix, Nix) - Etotal
      Qres(Zix, Nix, 0) = S(0, 0, k0) + targetE + Edif
    endif
    if (Ltarget /= 0 .and. Zix == parZ(k0) .and. Nix == parN(k0)) Qres(Zix, Nix, 0) = targetE
    do nex = 0, NL
      if (Ethresh(Zix, Nix, nex) == 0.) then
        Qres(Zix, Nix, nex) = Qres(Zix, Nix, 0) - edis(Zix, Nix, nex)
        Ethresh(Zix, Nix, nex) = - (Qres(Zix, Nix, nex) / specmass(parZ(k0), parN(k0), k0))
        Ethresh(Zix, Nix, nex) = max(Ethresh(Zix, Nix, nex), 0.d0)
      endif
    enddo
!
! The highest possible excitation energy could be a discrete state.
!
    do nex = 0, NL
      if (nex > 0 .and. edis(Zix, Nix, nex) >= Exmax(Zix, Nix)) then
        maxex(Zix, Nix) = nex - 1
        goto 180
      endif
      Ex(Zix, Nix, nex) = edis(Zix, Nix, nex)
      if (nex > 0) then
        deltaEx(Zix, Nix, nex) = 0.5 * (edis(Zix, Nix, min(NL, nex + 1)) - edis(Zix, Nix, nex - 1))
      else
        deltaEx(Zix, Nix, nex) = 0.5 * edis(Zix, Nix, 1)
      endif
    enddo
!
! Division of the continuum into bins.
!
    Aix = Zix + Nix
    if (Aix <= 4) then
      nexbins = nbins
    else
      if (Aix <= 8) then
        nexbins = int(real(1 - 0.1 * (Aix - 4)) * nbins)
      else
        nexbins = nbins / 2
      endif
    endif
    nexbins = max(nexbins, 2)
    eb = Ex(Zix, Nix, NL)
    ee = max(Exmax(Zix, Nix), eb + 0.001)
    if (flagequi .or. eb == 0.) then
      do i = 0, nexbins
        Eup(i) = eb + real(i) / nexbins * (ee - eb)
      enddo
    else
      eb = log(eb)
      ee = log(ee)
      do i = 0, nexbins
        Eup(i) = exp(eb + real(i) / nexbins * (ee - eb))
      enddo
    endif
    do i = 1, nexbins
      nex = NL + i
      Ex(Zix, Nix, nex) = 0.5 * (Eup(i - 1) + Eup(i))
      deltaEx(Zix, Nix, nex) = Eup(i) - Eup(i - 1)
    enddo
    maxex(Zix, Nix) = NL + nexbins
    if (Zix == 0 .and. Nix == 0) Ex(0, 0, maxex(0, 0) + 1) = Etotal
  180   nexmax(type) = maxex(Zix, Nix)
!
! ****** Determine level densities on basic excitation energy grid *****
!
! The calculation of level densities can be done outside many loops of various quantum numbers performed in other subroutines.
! Therefore, we store the level density as function of residual nucleus, excitation energy (nex), spin (Ir) and parity (Pprime)
! in the array rhogrid.
!
    A = AA(Zcomp, Ncomp, type)
    odd = mod(A, 2)
    Rodd = 0.5 * odd
    ald = real(A) / 8.
    ldmod = ldmodel(Zix, Nix)
    do nex = NL + 1, maxex(Zix, Nix)
      dEx = deltaEx(Zix, Nix, nex)
      Exout = Ex(Zix, Nix, nex)
      Ex1min = Exout - 0.5 * dEx
      Ex1plus = Exout + 0.5 * dEx
!
! Here we define the maximum J that can be reached in a given excitation energy bin.
! By default, this maxJ value is set to 3*sigma where sigma is the square root of the spin cut-off parameter of the level density.
!
      if (flagffruns) then
        ipop=1
      else
        ipop=0
      endif
      maxJ(Zix, Nix, nex) = int(4.+ 3. * sqrt(spincut(Zix, Nix, ald, Exout, 0, ipop)))
      maxJ(Zix, Nix, nex) = min(maxJ(Zix, Nix, nex), numJ)
!
! In the compound nucleus subroutines, the particle widths are determined by means of products of level densities and transmission
! coefficients.
! Instead of taking this product exactly at the middle of the excitation energy bins, we get a better numerical result by
! performing a logarithmic average over the bin for the level density, using the middle, top and bottom.
!
! For an equiprobable parity distribution for level densities, the loop over Pprime only needs to be performed once and the result
! for the level density is equal for both parities.
!
      if (flagparity) then
        Pbeg = - 1
      else
        Pbeg = 1
      endif
      do Pprime = Pbeg, 1, 2
        do Ir = 0, maxJ(Zix, Nix, nex)
          Rspin = real(Ir) + Rodd
          rho1 = density(Zix, Nix, Ex1min, Rspin, Pprime, 0, ldmod) * (1. + 1.d-10)
          rho2 = density(Zix, Nix, Exout, Rspin, Pprime, 0, ldmod)
          rho3 = density(Zix, Nix, Ex1plus, Rspin, Pprime, 0, ldmod) * (1. + 1.d-10)
          r1log = log(rho1)
          r2log = log(rho2)
          r3log = log(rho3)
          if (r2log /= r1log .and. r2log /= r3log) then
            rhogrid(Zix, Nix, nex, Ir, Pprime) = 0.5 * dEx * ((rho1 - rho2) / (r1log - r2log) + (rho2 - rho3) / (r2log - r3log))
          else
            rhogrid(Zix, Nix, nex, Ir, Pprime) = dEx * rho2
          endif
          if ( .not. flagparity) rhogrid(Zix, Nix, nex, Ir, - 1) = rhogrid(Zix, Nix, nex, Ir, 1)
        enddo
      enddo
    enddo
  enddo
end subroutine exgrid
! Copyright A.J. Koning 2021
