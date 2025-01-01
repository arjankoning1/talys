subroutine energies
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Energies
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
!   numbins        ! maximum number of continuum excitation energy bins
!   numlev2        ! maximum number of levels
!   numN           ! maximum number of neutrons from initial compound nucleus
!   numZ           ! maximum number of protons from initial compound nucleus
! Variables for numerics
!   nbins0         ! number of continuum excitation energy bins
!   xseps          ! limit for cross sections
! Variables for basic reaction
!   flagang        ! flag for output of angular distributions
!   flagrel        ! flag for relativistic kinematics
! Variables for main input
!   k0             ! index of incident particle
! Variables for compound reactions
!   eurr           ! off - set incident energy for URR calculation
!   ewfc           ! off - set incident energy for width fluctuation
!   flagurr        ! flag for output of unresolved resonance parameters
! Variables for direct reactions
!   eadd           ! on - set incident energy for addition of discrete states
!   eaddel         ! on - set incident energy for addition of elastic peak
!   flagcpang      ! flag for compound angular distribution calculation for
!   flaggiant0     ! flag for collective contribution from giant resonances
! Variables for preequilibrium
!   emulpre        ! on - set incident energy for multiple preequilibrium
!   epreeq         ! on - set incident energy for preequilibrium calculation
! Variables for energy grid
!   Einc           ! incident energy in MeV
! Variables for energy grid
!   ebegin         ! first energy point of energy grid
!   Ebottom        ! bottom of outgoing energy bin
!   egrid          ! outgoing energy grid
!   maxen          ! total number of energies
! Variables for energies
!   eend           ! last energy point of energy grid
!   eendhigh       ! last energy point for energy grid for any particle
!   eninccm        ! center - of - mass incident energy in MeV
!   eoutdis        ! outgoing energy of discrete state reaction
!   Etotal         ! total energy of compound system (target + projectile)
!   flagadd        ! flag for addition of discrete states to spectra flag
!   flagaddel      ! flag for addition of elastic peak to spectra
!   flagcompang    ! flag for compound angular distribution calculation
!   flaggiant      ! flag for collective contribution from giant resonances
!   flagmulpre     ! flag for multiple pre - equilibrium calculation
!   flagpreeq      ! flag for pre - equilibrium calculation
!   flagwidth      ! flag for width fluctuation calculation
!   mulpreZN       ! logical for multiple pre - equilibrium per nucleus
!   nendisc        ! last discrete bin
!   nbins          ! number of continuum excitation energy bins
!   speceps        ! limit for cross section spectra
!   wavenum        ! wave number
! Variables for nuclides
!   Nindex         ! neutron number index for residual nucleus
!   parskip        ! logical to skip outgoing particle
!   targetE        ! excitation energy of target
!   tarmass        ! mass of target nucleus
!   Zindex         ! charge number index for residual nucleus
! Constants
!   amu            ! atomic mass unit in MeV
!   hbarc          ! hbar.c in MeV.fm
!   parmass        ! mass of particle in a.m.u.
!   parN           ! neutron number of particle
!   parZ           ! charge number of particle
! Variables for levels
!   edis           ! energy of level
! Variables for level density
!   Nlast          ! last discrete level
! Variables for masses
!   redumass       ! reduced mass
!   S              ! separation energy
!   specmass       ! specific mass for residual nucleus
!
! *** Declaration of local data
!
  implicit none
  integer   :: i       ! counter
  integer   :: maxbins ! maximum number of bins
  integer   :: minbins ! minimum number of bins
  integer   :: nen     ! energy counter
  integer   :: Nix     ! neutron number index for residual nucleus
  integer   :: NL      ! last discrete level
  integer   :: type    ! particle type
  integer   :: Zix     ! charge number index for residual nucleus
  real(sgl) :: b       ! parameter for slope
  real(sgl) :: Elast   ! help variable
  real(dbl) :: eps     ! help variable
  real(dbl) :: ma      ! help variable
  real(dbl) :: ma1     ! help variable
  real(dbl) :: ma2     ! help variable
!
! ****************************** Energies ******************************
!
! Center-of-mass energy and wave number
!
! 1. Relativistic case: Ecm = sqrt[ (m+M)**2 + 2ME ] -(m+M)
!   = 2ME / { sqrt [ (m+M)**2 + 2ME ] + (m+M) }
!
! For incident photons, we always take the relativistic case.
!
  if (flagrel .or. k0 == 0) then
    ma1 = parmass(k0)
    ma2 = tarmass
    ma = ma1 + ma2
    eps = 2. * ma2 * Einc / amu
    eninccm = real(amu * eps / (sqrt(ma **2 + eps) + ma))
    wavenum = real(amu * ma2 / hbarc * sqrt(Einc / amu * (Einc / amu + 2. * ma1) / ((ma1 + ma2) **2 + 2. * ma2 * Einc / amu)))
  else
!
! 2. Non-relativistic case: Ecm = E*M/(m+M)
!
    eninccm = real(Einc * specmass(parZ(k0), parN(k0), k0))
    wavenum = sqrt(real(2. * amu * redumass(parZ(k0), parN(k0), k0) * eninccm)) / hbarc
  endif
  if (flaginitpop) eninccm = Einc
  Etotal = eninccm + S(0, 0, k0) + targetE
!
! ***************** Set upper limit for energy grid ********************
!
  eendhigh = 0
  do type = 0, 6
    eend(type) = maxen - 1
    if (parskip(type)) cycle
    do nen = 0, maxen
      if (egrid(nen) > Etotal - S(0, 0, type)) then
        eend(type) = nen
        exit
      endif
    enddo
    if (eend(type) > ebegin(type)) eend(type) = max(eend(type), ebegin(type) + 3)
    eendhigh = max(eendhigh, eend(type))
  enddo
!
! **** Set limit for cross section spectra ******
!
  speceps = xseps / eninccm
!
! ********* Set outgoing energies for discrete level scattering ********
!
! locate : subroutine to find value in ordered table
!
! The variable nendisc is used to determine the last discrete level of the population of a residual nucleus.
! This will be needed for the calculation of boundary effects in spectra.
!
  do type = 0, 6
    nendisc(type) = 1
    if (parskip(type)) cycle
    Zix = Zindex(0, 0, type)
    Nix = Nindex(0, 0, type)
    NL = Nlast(Zix, Nix, 0)
    do i = 0, numlev2
      eoutdis(type, i) = Etotal - S(0, 0, type) - edis(Zix, Nix, i)
    enddo
    if (ebegin(type) >= eend(type)) cycle
    Elast = eoutdis(type, NL)
    if (Elast > 0.) then
      call locate(Ebottom, 1, eend(type), Elast, nen)
      nendisc(type) = nen
    endif
  enddo
!
! ************ Set pre-equilibrium and compound nucleus flags **********
!
  if (Einc <= ewfc)  then
    flagwidth = .true.
  else
    flagwidth = .false.
  endif
  if (Einc <= eurr)  then
    flagurr = .true.
  else
    flagurr = .false.
  endif
  if ((k0 == 1 .or. flagcpang) .and. flagang .and. Einc <= 50.)  then
    flagcompang = .true.
  else
    flagcompang = .false.
  endif
  if (Einc < epreeq) then
    flagpreeq = .false.
    flaggiant = .false.
  else
    flagpreeq = .true.
    if (flaggiant0) then
      flaggiant = .true.
    else
      flaggiant = .false.
    endif
  endif
  if (Einc < emulpre) then
    flagmulpre = .false.
  else
    flagmulpre = .true.
  endif
  if (Einc < eadd) then
    flagadd = .false.
  else
    flagadd = .true.
  endif
  if (Einc < eaddel) then
    flagaddel = .false.
  else
    flagaddel = .true.
  endif
  if (flagffruns.or.flagrpruns.or.flaginitpop) then
    flagaddel = .false.
    flagadd = .false.
  endif
  do Zix = 0, numZ
    do Nix = 0, numN
      mulpreZN(Zix, Nix) = .false.
    enddo
  enddo
!
! Progressive bin width for excitation energy bins
!
  if (nbins0 == 0) then
    minbins = 30
    maxbins = numbins
    b = 60.
    nbins = minbins + int((maxbins - minbins) * Einc * Einc / (Einc * Einc + b * b))
  else
    nbins = nbins0
  endif
  return
end subroutine energies
! Copyright A.J. Koning 2021
