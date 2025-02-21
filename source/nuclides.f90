subroutine nuclides
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Properties of nuclides
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
! All global variables
!   numN            ! maximum number of neutrons from initial compound nucleus
!   numZ            ! maximum number of protons from initial compound nucleus
! Variables for compound reactions
!   eurr            ! off - set incident energy for URR calculation
!   ewfc            ! off - set incident energy for width fluctuation
!   flagcomp        ! flag for compound angular distribution calculation
! Variables for direct reactions
!   flaggiant0      ! flag for collective contribution from giant resonances
! Variables for basic reaction
!   flagastro       ! flag for calculation of astrophysics reaction rate
!   flagpartable    ! flag for output of model parameters on separate file
! Variables for medical isotope production
!   flagprod        ! flag for isotope production
! Variables for numerics
!   maxN            ! maximal number of neutrons away from initial compound nucleus
!   maxZ            ! maximal number of protons away from initial compound nucleus
! Variables for input energies
!   enincmin        ! minimum incident energy
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
!   Ltarget         ! excited level of target
!   Ninit           ! neutron number of initial compound nucleus
!   Zinit           ! charge number of initial compound nucleus
!   Ztarget         ! charge number of target nucleus
! Variables for preequilibrium
!   epreeq          ! on - set incident energy for preequilibrium calculation
! Variables for discrete levels
!   flagbestbr      ! flag to use only best set of branching ratios
!   nlev            ! number of levels for nucleus
!   nlevbin         ! number of excited levels for binary nucleus
!   nlevmax         ! maximum number of included discrete levels for target
!   nlevmaxres      ! maximum number of included discrete levels for residual nuclides
! Variables for energy grid
!   Einc            ! incident energy in MeV
! Variables for compound reactions
!   flagurr        ! flag for output of unresolved resonance parameters
! Variables for nuclides
!   AA              ! mass number of residual nucleus
!   coulbar         ! Coulomb barrier
!   invexist        ! logical to state necessity of new inverse cross section calc.
!   Nindex          ! neutron number index for residual nucleus
!   NN              ! neutron number of residual nucleus
!   parinclude      ! logical to include outgoing particle
!   parskip         ! logical to skip outgoing particle
!   primary         ! flag to designate primary (binary) reaction
!   Q               ! Q - value
!   strucexist      ! flag to state whether structure info for this nucleus exists
!   strucwrite      ! flag for output of nuclear structure info
!   targetE         ! excitation energy of target
!   targetP         ! parity of target
!   targetspin      ! spin of target
!   targetspin2     ! 2 * target spin
!   tarmass         ! mass of target nucleus
!   Zindex          ! charge number index for residual nucleus
!   ZZ              ! charge number of residual nucleus
! Constants
!   e2              ! square of elementary charge in MeV.fm
!   onethird        ! 1 / 3
!   parN            ! neutron number of particle
!   parZ            ! charge number of particle
! Variables for levels
!   edis            ! energy of level
!   jdis            ! spin of level
!   parlev          ! parity of level
! Variables for resonance parameters
!   Eavres        ! number of resonances
! Variables for level density
!   Nlast           ! last discrete level
! Variables for masses
!   nucmass         ! mass of nucleus
!   S               ! separation energy
!
! *** Declaration of local data
!
  implicit none
  character(len=3) :: massstring ! mass string
  integer          :: A          ! mass number of target nucleus
  integer          :: Ncomp      ! neutron number index for compound nucleus
  integer          :: Nix        ! neutron number index for residual nucleus
  integer          :: NL         ! last discrete level
  integer          :: odd        ! odd (1) or even (0) nucleus
  integer          :: type       ! particle type
  integer          :: Zcomp      ! proton number index for compound nucleus
  integer          :: Zix        ! charge number index for residual nucleus
!
! ************************ Properties of nuclides **********************
!
! Assignment of Z and N of all possible residual nuclei.
! The initial compound nucleus (created by projectile + target) has the indices (0,0).
! The first index represents the number of protons and the second index the number of neutrons away from the initial compound
! nucleus.
! Example: For the reaction p + 208Pb, the set (0,0) represents 209Bi and the set (1,2) represents 206Pb.
! Many arrays have Zindex and Nindex as their first two indices.
! At any point in the reaction calculation, given Zcomp, Ncomp and the particle type, these variables are directly known
! through the arrays we initialize here.
! ZZ is the charge number of the nucleus that is reached through emission of ejectile type from nucleus (Zcomp,Ncomp).
! Note that ZZ, NN, AA, Zinit and Ninit represent true values of the charge and neutron number and that Zindex, Nindex, Zcomp and
! Ncomp are indices relative to the initial compound nucleus.
!
! Zcomp=0 (primary compound nucleus)
! Ncomp=0 (primary compound nucleus)
! Zindex=1
! Nindex=0
! Zinit=26
! Ninit=31
! ZZ=25
! NN=31
!
  do Zcomp = 0, numZ
    do Ncomp = 0, numN
      do type = 0, 6
        Zindex(Zcomp, Ncomp, type) = Zcomp + parZ(type)
        Nindex(Zcomp, Ncomp, type) = Ncomp + parN(type)
        ZZ(Zcomp, Ncomp, type) = Zinit - Zindex(Zcomp, Ncomp, type)
        NN(Zcomp, Ncomp, type) = Ninit - Nindex(Zcomp, Ncomp, type)
        AA(Zcomp, Ncomp, type) = ZZ(Zcomp, Ncomp, type) + NN(Zcomp, Ncomp, type)
      enddo
    enddo
  enddo
!
! We make sure the lightest possible residual nucleus is heavier than an alpha particle.
! Also, we set that at the beginning no structure information is available at all .
!
  maxZ = min(maxZ, Zinit - 3)
  maxN = min(maxN, Ninit - 3)
  do Zix = 0, numZ
    do Nix = 0, numN
      strucexist(Zix, Nix) = .false.
      strucwrite(Zix, Nix) = .false.
      invexist(Zix, Nix) = .false.
    enddo
  enddo
!
! Set the maximum number of discrete levels for each residual nucleus
!
! Note the order of priority: nlevmax (keyword: maxlevels) overrules the value set for nlevbin (keyword: maxlevelsbin)
! for the inelastic channel.
! nlevbin overrules nlevmaxres (keyword: maxlevelsres) for particular binary channels.
!
  Zcomp = 0
  Ncomp = 0
  do type = 0, 6
    Zix = Zindex(Zcomp, Ncomp, type)
    Nix = Nindex(Zcomp, Ncomp, type)
    if (nlev(Zix, Nix) == 0) nlev(Zix, Nix) = nlevbin(type)
  enddo
  do Zix = 0, numZ
    do Nix = 0, numN
      if (Zix == parZ(k0) .and. Nix == parN(k0)) cycle
      if (nlev(Zix, Nix) == 0) nlev(Zix, Nix) = nlevmaxres
    enddo
  enddo
  nlev(parZ(k0), parN(k0)) = nlevmax
!
! ************ Nuclear structure for first set of nuclides *************
!
! strucinitial: subroutine for initialization of arrays for various structure parameters
! masses      : subroutine for nuclear masses
! separation  : subroutine for separation energies and reduced and specific masses
! branching   : subroutine for best set of branching ratios
! structure   : subroutine for nuclear structure parameters
! weakcoupling: subroutine for weak coupling model
! radwidtheory: subroutine for theoretical calculation of total radiative width
! sumrules    : subroutine for giant resonance sum rules
!
! The nuclear masses and separation energies are read/calculated first, for all nuclides that can possibly be formed in
! multiple reactions.
! For all nuclides that can be formed by the first binary reaction, we determine the nuclear structure information (discrete levels,
! parameters for resonances, photons, fission and level densities).
! The logical strucexist ensures that we only need to do this once for each nucleus.
! (i.e., we do not waste time if the same nucleus is encountered later in the reaction chain).
! All parameters for each nucleus are written to file 'parameters.dat'.
!
  call strucinitial
  call masses
  call separation
  primary = .true.
  if (flagbestbr) call branching
  massstring = '   '
  write(massstring,'(i3)') Atarget
  targetnuclide = trim(Starget) // adjustl(massstring)
  targetnuclide0 = targetnuclide
!
! Open file to write all nuclear model parameters
!
  if (flagpartable) then
    open (unit = 51, file = 'parameters.dat', status = 'unknown')
    write(51, '("## header:")')
    write(51, '("##   title: ", a, "(", a1, ",x) nuclear model parameters")') trim(targetnuclide),parsym(k0)
    write(51, '("##   source: ", a)') trim(source)
    write(51, '("##   user: ", a)') trim(user)
    write(51, '("##   date: ", a)') trim(date)
    write(51, '("##   format: ", a)') trim(oformat)
    write(51, '("## target:")')
    write(51, '("##   Z: ", i0)') Ztarget
    write(51, '("##   A: ", i0)') Atarget
    write(51, '("##   nuclide: ", a)') targetnuclide
  endif
!
! Retrieve nuclear structure information
!
  Einc = enincmin
  do type = 0, 6
    if (parskip(type) .and. type /= 0) cycle
    Zix = Zindex(Zcomp, Ncomp, type)
    Nix = Nindex(Zcomp, Ncomp, type)
    call structure(Zix, Nix)
    strucexist(Zix, Nix) = .true.
    A = AA(Zcomp, Ncomp, type)
    odd = mod(A, 2)
    if (type == k0 .and. odd == 1) call weakcoupling(Zix, Nix, type)
  enddo
  if (parinclude(0) .or. flagcomp) call radwidtheory(Zcomp, Ncomp, Eavres)
  if (flaggiant0) call sumrules
  targetspin = jdis(parZ(k0), parN(k0), Ltarget)
  targetspin2 = int(2. * targetspin)
  targetP = parlev(parZ(k0), parN(k0), Ltarget)
  targetE = edis(parZ(k0), parN(k0), Ltarget)
!
! ******************* Q-values and Coulomb barriers ********************
!
  tarmass = nucmass(parZ(k0), parN(k0))
  Q = 0.
  do type = 0, 6
    if (parskip(type)) cycle
    Q(type) = S(0, 0, k0) - S(0, 0, type)
    coulbar(type) = Ztarget * parZ(type) * e2 / (1.25 * Atarget **onethird)
  enddo
!
! * Off- and on-set energies for preequilibrium and width fluctuations *
!
! Assignment of default values: Width fluctuation corrections and URR are included for incident energies
! up to the separation energy.
!
  if (k0 >= 1 .and. ewfc ==  -1.) ewfc = S(parZ(k0), parN(k0), k0)
  if (k0 == 1 .and. eurr ==  -1. .and. flagurr) eurr = S(parZ(k0), parN(k0), k0)
!
! Pre-equilibrium reactions are included for incident energies above the last discrete level.
!
! kalbachsep: subroutine for separation energy for Kalbach systematics
!
  NL = Nlast(parZ(k0), parN(k0), 0)
  if (epreeq ==  - 1.) epreeq = max(edis(parZ(k0), parN(k0), NL), 1.)
  call kalbachsep
!
! For astrophysical calculations, the incident energy grid is determined
! egridastro: subroutine to calculate default incident energy grid for astrophysical rate
!
  if (flagastro) call egridastro
!
! decaydata : subroutine for decay data
!
  if (flagprod) call decaydata
  return
end subroutine nuclides
! Copyright A.J. Koning 2021
