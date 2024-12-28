subroutine checkvalue
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Check for errors in values
!
! Author    : Arjan Koning
!
! 2024-12-08: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
  use A1_error_handling_mod
!
! Definition of single and double precision variables
!   sgl               ! single precision kind
! All global variables
!   numang            ! maximum number of angles
!   numangcont        ! maximum number of angles for continuum
!   numangrec         ! maximum number of recoil angles
!   numbar            ! number of fission barriers
!   numbins           ! maximum number of continuum excitation energy bins
!   numelem           ! number of elements
!   numenmsd          ! maximum number of energy points for DWBA calculation for MSD
!   numenrec          ! maximum number of recoil energies
!   numgam            ! maximum number of l - values for gamma multipolarity
!   numisom           ! number of isomers
!   numl              ! number of l values
!   numlev            ! maximum number of discrete levels
!   nummass           ! number of masses
!   nummt             ! number of MT numbers
!   numN              ! maximum number of neutrons from initial compound nucleus
!   numneu            ! number of neutrons
!   numNph            ! maximum number of neutrons away from the initial compound nucleus
!   numompadj         ! number of adjustable ranges for OMP
!   numZ              ! maximum number of protons from initial compound nucleus
!   numZph            ! maximum number of protons away from the initial compound nucleus
! Variables for adjustment
!   adjustfile        ! file for local adjustment
!   adjustkey         ! keyword for local adjustment
!   adjustpar         ! local adjustment parameters
!   Nadjust           ! number of adjustable parameters
! Variables for basic reaction
!   flagastro         ! flag for calculation of astrophysics reaction rate
!   flagmassdis       ! flag for calculation of fission fragment mass yields
!   ompenergyfile     ! file with energies for OMP calculation (ENDF files)
! Variables for best files
!   grescue           ! global multipl. factor for incident energy dep. adj. factors
!   rescuefile        ! file with incident energy dependent adjustment factors
! Variables for numerics
!   maxchannel        ! maximal number of outgoing particles in individual channel description
!   maxenrec          ! number of recoil energies
!   maxN              ! maximal number of neutrons away from initial compound nucleus
!   maxZ              ! maximal number of protons away from initial compound nucleus
!   maxNrp            ! maximal number of neutrons away from initial compound nucleus for residual production
!   maxZrp            ! maximal number of protons away from initial compound nucleus for residual production
!   nangle            ! number of angles
!   nanglecont        ! number of angles for continuum
!   nbins0            ! number of continuum excitation energy bins
!   nanglerec         ! number of recoil angles
!   popeps            ! limit for population cross sections
!   segment           ! help array for storing segment intersection points
!   transeps          ! absolute limit for transmission coefficient
!   transpower        ! power for transmission coefficient limit
!   xseps             ! limit for cross sections
! Variables for basic parameters
!   eninclow          ! minimal incident energy for nuclear model calculations
!   isomer            ! definition of isomer in seconds
!   Lisoinp           ! user assignment of target isomer number
! Variables for input energies
!   eninc             ! incident energy in MeV
!   enincmax          ! maximum incident energy
!   enincmin          ! minimum incident energy
!   Estop             ! incident energy above which TALYS stops
! Variables for main input
!   Ainit             ! mass number of initial compound nucleus
!   Atarget           ! mass number of target nucleus
!   flagnatural       ! flag for calculation of natural element
!   k0                ! index of incident particle
!   Ltarget           ! excited level of target
!   Ninit             ! neutron number of initial compound nucleus
!   ptype0            ! type of incident particle
!   Starget           ! symbol of target nucleus
!   Zinit             ! charge number of initial compound nucleus
! Variables for level density
!   D0                ! s - wave resonance spacing in eV
! Variables for output
!   ddxacount         ! counter for double - differential cross section files
!   ddxecount         ! counter for double - differential cross section files
!   ddxmode           ! mode for DDX: 0: None, 1: Angular distributions, 2: Spectra per angle, 3: Both
!   fileddxa          ! designator for double - differential cross sections on separate file
!   fileddxe          ! designator for double - differential cross sections on separate file
!   flagdecay         ! flag for output of decay of each population bin
!   flagpop           ! flag for output of population
! Variables for compound reactions
!   eurr              ! off - set incident energy for URR calculation
!   ewfc              ! off - set incident energy for width fluctuation
!   flagcomp          ! flag for compound angular distribution calculation
!   flageciscomp      ! flag for compound nucleus calculation by ECIS
!   flagres           ! flag for output of low energy resonance cross sections
!   flagurr           ! flag for output of unresolved resonance parameters
!   lurr              ! maximal orbital angular momentum for URR
!   reslib            ! library with resonance parameters
!   Tres              ! temperature for broadening low energy cross sections
!   wmode             ! designator for width fluctuation model
!   xsalphatherm      ! thermal (n, a) cross section
!   xscaptherm        ! thermal capture cross section
!   xsptherm          ! thermal (n, p) cross section
! Variables for direct reactions
!   core              ! even - even core for weakcoupling ( - 1 or 1)
!   eadd              ! on - set incident energy for addition of discrete states
!   eaddel            ! on - set incident energy for addition of elastic peak
!   elwidth           ! width of elastic peak in MeV
!   flaggiant0        ! flag for collective contribution from giant resonances
!   maxband           ! highest vibrational band added to rotational model
!   maxrot            ! number of included excited rotational levels
!   soswitch          ! switch for deformed spin - orbit calculation
! Variables for preequilibrium
!   breakupmodel      ! model for break - up reaction: 1. Kalbach 2. Avrigeanu
!   Cbreak            ! adjustable parameter for break - up reactions
!   Cknock            ! adjustable parameter for knockout reactions
!   Cstrip            ! adjustable parameter for stripping / pick - up reactions
!   Emsdmin           ! minimal outgoing energy for MSD calculation
!   emulpre           ! on - set incident energy for multiple preequilibrium
!   epreeq            ! on - set incident energy for preequilibrium calculation
!   Esurf0            ! well depth for surface interaction
!   flagpecomp        ! flag for Kalbach complex particle emission model
!   g                 ! single - particle level density parameter
!   gadjust           ! adjustable factor for single - particle particle-hole states
!   gn                ! single - particle neutron level density parameter
!   gnadjust          ! adjustable factor for single - particle proton parameter
!   gp                ! single - particle proton level density parameter
!   gpadjust          ! adjustable factor for single - particle neutron parameter
!   Kph               ! constant for single - particle level density par. (g = A / Kph)
!   M2constant        ! constant for matrix element in exciton model
!   M2limit           ! constant for asymptotic value for matrix element
!   M2shift           ! constant for energy shift for matrix element
!   mpreeqmode        ! designator for multiple pre - equilibrium model
!   msdbins           ! number of energy points for DWBA calculation for MSD
!   pairmodel         ! model for preequilibrium pairing energy
!   pespinmodel       ! model for pre - equilibrium or compound spin distribution
!   phmodel           ! particle - hole state density model
!   preeqmode         ! designator for pre - equilibrium model
!   Rgamma            ! adjustable parameter for pre - equilibrium gamma decay
!   Rnunu             ! ratio for two - component matrix element
!   Rnupi             ! ratio for two - component matrix element
!   Rpinu             ! ratio for two - component matrix element
!   Rpipi             ! ratio for two - component matrix element
!   Rspincutpreeq     ! adjustable constant (global) for preequilibrium spin cutoff factor
! Variables for fission
!   axtype            ! type of axiality of barrier
!   bdamp             ! fission partial damping parameter
!   bdampadjust       ! correction for fission partial damping parameter
!   betafiscor        ! adjustable factor for fission path width
!   clas2file         ! file with class 2 transition states
!   Cnubar1           ! adjustable parameter for nubar constant value
!   Cnubar2           ! adjustable parameter for nubar constant value
!   Tmadjust          ! adjustable parameter for PFNS temperature
!   Fsadjust          ! adjustable parameter for PFNS scission fraction
!   fbaradjust        ! adjustable factor for fission parameters
!   fbarrier          ! height of fission barrier
!   flagfission       ! flag for fission
!   flagfisout        ! flag for output of fission information
!   fismodel          ! fission model alternative fission model for default barriers
!   fismodelalt       ! alternative fission model for default barriers
!   fwidth            ! width of fission barrier
!   fwidthadjust      ! adjustable factor for fission parameters
!   fymodel           ! fission yield model, 1: Brosa 2: GEF
!   ffmodel           ! fission fragment model, 1: GEF 2: HF3D (Okumura) 3: SPY 4: Langevin-4D
!   pfnsmodel         ! PFNS  model, 1: Iwamoto 2: from FF decay
!   gefran            ! number of random events for GEF calculation
!   hbtransfile       ! file with head band transition states
!   Rfiseps           ! ratio for limit for fission cross section per nucleus
!   vfiscor           ! adjustable factor for fission path height
!   yieldfile         ! fission yield file
! Variables for astrophysics
!   astroE            ! energy, in MeV, for Maxwellian average
!   astroT9           ! temperature, in 10^9 K, for Maxwellian average
!   nonthermlev       ! non - thermalized level in the calculation of astrophysic
! Variables for medical isotope production
!   Area              ! target area in cm^2
!   Eback             ! lower end of energy range in MeV for isotope
!   Ebeam             ! incident energy in MeV for isotope production
!   flagprod          ! flag for isotope production
!   Ibeam             ! beam current in mA for isotope production
!   radiounit         ! unit for radioactivity: Bq, kBq, MBq, Gbq, mCi, Ci or kCi
!   rhotarget         ! target material density
!   Tcool             ! cooling time per unit cooling time unit (y, d, h, m, s)
!   Tirrad            ! irradiation time per unit irradiation time unit (y, d, h, m, s)
!   unitTcool         ! cooling time unit (y, d, h, m, s)
!   unitTirrad        ! irradiation time unit (y, d, h, m, s)
!   yieldunit         ! unit for isotope yield: num (number), mug, mg, g, or kg
! Variables for gamma rays
!   egr               ! energy of GR
!   egradjust         ! adjustable factor for energy of GR
!   epr               ! energy of PR
!   epradjust         ! adjustable factor for energy of PR
!   etable            ! constant to adjust tabulated strength functions
!   etableadjust      ! correction to adjust tabulated strength functions
!   Exlfile           ! tabulated gamma ray strength function
!   fiso              ! correction factor for isospin forbidden transitions
!   fisom             ! correction factor for isospin forbidden transitions for multiple emission
!   flagracap         ! flag for radiative capture model
!   ftable            ! constant to adjust tabulated strength functions
!   ftableadjust      ! correction to adjust tabulated strength functions
!   gamgam            ! total radiative width in eV
!   gamgamadjust      ! adjustable factor for radiative parameter
!   ggr               ! width of GR
!   ggradjust         ! adjustable factor for width of GR
!   gammax            ! number of l - values for gamma multipolarity
!   gpr               ! width of PR
!   gpradjust         ! adjustable factor for width of PR
!   ldmodelracap      ! level density model for direct radiative capture
!   sgr               ! strength of GR
!   sgradjust         ! adjustable factor for strength of GR
!   spectfacexp       ! experimental spectroscopic factor
!   spectfacth        ! theoretical spectroscopic factor
!   strength          ! E1 strength function model
!   strengthM1        ! model for M1 gamma - ray strength function
!   tpr               ! strength of PR
!   tpradjust         ! adjustable factor for strength of PR
!   upbend            ! properties of the low - energy upbend of given multipolarity
!   wtable            ! constant to adjust tabulated strength functions
!   wtableadjust      ! correction to adjust tabulated strength functions
! Variables for discrete levels
!   disctable         ! table with discrete levels
!   deformfile        ! deformation parameter file
!   levelfile         ! discrete level file
!   nlev              ! number of levels for nucleus
!   nlevbin           ! number of excited levels for binary nucleus
!   nlevmax           ! maximum number of included discrete levels for target
!   nlevmaxres        ! maximum number of included discrete levels for residual nuclides
! Variables for level density
!   aadjust           ! adjustable factor for level density parameter
!   alev              ! level density parameter
!   alimit            ! asymptotic level density parameter
!   alphald           ! alpha - constant for asymptotic level density parameter
!   betald            ! beta - constant for asymptotic level density parameter
!   cfermi            ! width of Fermi distribution for damping
!   cglobal           ! global constant to adjust tabulated level densities
!   ctable            ! constant to adjust tabulated level densities
!   ctableadjust      ! correction to adjust tabulated level densities
!   deltaW            ! shell correction in nuclear mass
!   E0                ! particle constant of temperature formula
!   E0adjust          ! adjustable factor for E0
!   Exmatch           ! matching point for Ex
!   Exmatchadjust     ! adjustable factor for matching energy
!   gammald           ! gamma - constant for asymptotic level density parameter
!   gammashell1       ! gamma - constant for asymptotic level density parameter
!   gammashell2       ! gamma - constant for asymptotic level density parameter
!   Krotconstant      ! normalization constant for rotational enhancement
!   kvibmodel         ! model for vibrational enhancement
!   ldmodel           ! level density model
!   ldmodelCN         ! level density model for compound nucleus
!   Nlow              ! lowest discrete level for temperature matching
!   Ntop              ! highest discrete level for temperature matching
!   pair              ! pairing energy
!   pairconstant      ! constant for pairing energy systematics
!   pglobal           ! global constant to adjust tabulated level densities
!   Pshift            ! adjustable pairing shift
!   Pshiftadjust      ! adjustable correction to pairing shift
!   Pshiftconstant    ! global constant for pairing shift
!   ptable            ! constant to adjust tabulated level densities
!   ptableadjust      ! correction to adjust tabulated level densities
!   Rclass2mom        ! norm. constant for moment of inertia for class 2 states
!   Risomer           ! adjustable correction to level branching ratios
!   Rspincut          ! adjustable constant (global) for spin cutoff factor
!   Rtransmom         ! norm. constant for moment of inertia for transition states
!   s2adjust          ! adjustable constant (Z, A, barrier - dependent) for spin
!   shellmodel        ! model for shell correction energies
!   spincutmodel      ! model for spin cutoff factor for ground state
!   T                 ! temperature
!   Tadjust           ! adjustable factor for temperature
!   Ufermi            ! energy of Fermi distribution for damping
! Variables for masses
!   beta2             ! deformation parameter
!   massdir           ! directory with mass tables
!   massexcess        ! mass excess in MeV as read from user input file
!   massmodel         ! model for theoretical nuclear mass
!   massnucleus       ! mass of nucleus in amu as read from user input file
! Variables for OMP
!   adepthcor         ! adjustable parameter for depth of DF alpha potential
!   alphaomp          ! alpha optical model
!   aradialcor        ! adjustable parameter for shape of DF alpha potential
!   avadjust          ! adjustable factor for OMP (default 1.)
!   avdadjust         ! adjustable factor for OMP (default 1.)
!   avsoadjust        ! adjustable factor for OMP (default 1.)
!   awadjust          ! adjustable factor for OMP (default 1.)
!   awdadjust         ! adjustable factor for OMP (default 1.)
!   awsoadjust        ! adjustable factor for OMP (default 1.)
!   d1adjust          ! adjustable factor for OMP (default 1.)
!   d2adjust          ! adjustable factor for OMP (default 1.)
!   d3adjust          ! adjustable factor for OMP (default 1.)
!   deuteronomp       ! deuteron optical model
!   Ejoin             ! joining energy for high energy OMP
!   jlmmode           ! option for JLM imaginary potential normalization
!   lv1adjust         ! adjustable parameter for JLM OMP
!   lvadjust          ! adjustable parameter for JLM OMP
!   lvsoadjust        ! adjustable parameter for JLM OMP
!   lw1adjust         ! adjustable parameter for JLM OMP
!   lwadjust          ! adjustable parameter for JLM OMP
!   lwsoadjust        ! adjustable parameter for JLM OMP
!   ompadjustD        ! depth of local OMP adjustment
!   ompadjustE1       ! start energy of local OMP adjustment
!   ompadjustE2       ! end energy of local OMP adjustment
!   ompadjustN        ! number of energy ranges for local OMP adjustment
!   ompadjusts        ! variance of local OMP adjustment
!   optmod            ! file with optical model parameters
!   optmodfileN       ! optical model parameter file for neutrons
!   optmodfileP       ! optical model parameter file for protons
!   radialfile        ! radial matter density file
!   radialmodel       ! model for radial matter densities (JLM OMP only)
!   RprimeU           ! potential scattering radius
!   rcadjust          ! adjustable factor for OMP (default 1.)
!   rvadjust          ! adjustable factor for OMP (default 1.)
!   rvdadjust         ! adjustable factor for OMP (default 1.)
!   rvsoadjust        ! adjustable factor for OMP (default 1.)
!   rwadjust          ! adjustable factor for OMP (default 1.)
!   rwdadjust         ! adjustable factor for OMP (default 1.)
!   rwsoadjust        ! adjustable factor for OMP (default 1.)
!   v1adjust          ! adjustable factor for OMP (default 1.)
!   v2adjust          ! adjustable factor for OMP (default 1.)
!   v3adjust          ! adjustable factor for OMP (default 1.)
!   v4adjust          ! adjustable factor for OMP (default 1.)
!   Vinfadjust        ! adj. factor for high energy limit of real central potential
!   vso1adjust        ! adjustable factor for OMP (default 1.)
!   vso2adjust        ! adjustable factor for OMP (default 1.)
!   w1adjust          ! adjustable factor for OMP (default 1.)
!   w2adjust          ! adjustable factor for OMP (default 1.)
!   w3adjust          ! adjustable factor for OMP (default 1.)
!   w4adjust          ! adjustable factor for OMP (default 1.)
!   wso1adjust        ! adjustable factor for OMP (default 1.)
!   wso2adjust        ! adjustable factor for OMP (default 1.)
! Variables for files
!   path              ! directory containing files to be read
! Constants
!   Emaxtalys         ! maximum acceptable energy for TALYS
!   nuc               ! symbol of nucleus
!   parsym            ! symbol of particle
! Constants
!   kT                ! energy kT expressed in MeV corresponding to a temperature T9 = 1
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   range_real_error    ! Test if real variable is out of range
!
! *** Declaration of local data
!
  implicit none
  logical            :: lexist      ! logical to determine existence
  character(len=132) :: massdir0    ! mass directory
  character(len=132) :: massfile    ! mass file
  integer            :: A           ! mass number of target nucleus
  integer            :: i           ! counter
  integer            :: ibar        ! fission barrier
  integer            :: igr         ! giant resonance
  integer            :: irad        ! variable to indicate M(=0) or E(=1) radiation
  integer            :: is          ! isotope counter: -1=total, 0=ground state 1=isomer
  integer            :: k           ! designator for particle
  integer            :: l           ! multipolarity
  integer            :: lval        ! multipolarity
  integer            :: m           ! counter
  integer            :: mt          ! MT number
  integer            :: n           ! exciton number
  integer            :: Nix         ! neutron number index for residual nucleus
  integer            :: nr          ! number of radial grid point
  integer            :: nr2         ! counter
  integer            :: omptype     ! type of optical model (spherical or coupled)
  integer            :: type        ! particle type
  integer            :: Z           ! charge number of target nucleus
  integer            :: Zix         ! charge number index for residual nucleus
  real(sgl)          :: D           ! depth of local adjustment
  real(sgl)          :: Ea          ! start energy of local adjustment
  real(sgl)          :: Ea2         ! start energy of local adjustment
  real(sgl)          :: Eb          ! end energy of local adjustment
  real(sgl)          :: Eb2         ! end energy of local adjustment
  real(sgl)          :: Em          ! intermediate energy of local adjustment
!
! All parameters need to fall within certain ranges.
! These ranges are specified in this subroutine and in the manual.
!
! ******************* Check for wrong input variables ******************
!
! 1. Check of values for four main keywords.
!
  do type = 0, 6
    if (ptype0 == parsym(type)) goto 20
  enddo
  if (ptype0 == '0') goto 20
  write(*, '(" TALYS-error: Wrong symbol for projectile: ", a1)') ptype0
  stop
   20 do i = 3, numelem
    if (Starget == nuc(i)) goto 40
   enddo
  write(*, '(" TALYS-error: Wrong symbol for element: ", a2)') Starget
  stop
40 call range_integer_error('Target mass', Atarget, 6, nummass)
  call range_integer_error('CN Z', Zinit, 3, numelem)
  call range_integer_error('CN N', Ninit, 3, numneu)
  call range_real_error('Incident energy', enincmin, 1.e-11, Emaxtalys, unit = 'MeV')
  call range_real_error('Incident energy', enincmax, 1.e-11, Emaxtalys, unit = 'MeV')
  call range_real_error('Estop', Estop, 1.e-11, Emaxtalys, unit = 'MeV')
!
! 2. Check of values for basic physical and numerical parameters
!
  call range_integer_error('maxZ', maxZ, 0, numZ - 2)
  call range_integer_error('maxN', maxN, 0, numN - 2)
  call range_integer_error('maxZrp', maxZrp, 0, numZ - 2)
  call range_integer_error('maxNrp', maxNrp, 0, numN - 2)
  call range_integer_error('bins', nbins0, 2, numbins, default = 0)
  call range_integer_error('segment', segment, 1, 4)
  if (segment > 1 .and. enincmax > 100.) then
    write(*, '(" TALYS-error: segment = 1 for incident energy of ", f8.3, " MeV")') enincmax
    stop
  endif
  if (segment > 2 .and. enincmax > 40.) then
    write(*, '(" TALYS-error: 1 <= segment = 2 for incident energy of ", f8.3, " MeV")') enincmax
    stop
  endif
  if (segment > 3 .and. enincmax > 20.) then
    write(*, '(" TALYS-error: 1 <= segment = 3 for incident energy of ", f8.3, " MeV")') enincmax
    stop
  endif
  if (segment > 1 .and. flagastro) then
    write(*, '(" TALYS-error: segment = 1 for astrophysical calculations")')
    stop
  endif
  if (flagnatural .and. Ltarget > 0) then
    write(*,'(" TALYS-error: Excited level for target not possible for natural targets")')
    stop
  endif
  if (Lisoinp ==  - 1) then
    call range_integer_error('maxlevelstar', nlevmax, 0, numlev)
    call range_integer_error('maxlevelsres', nlevmaxres, 0, numlev)
    do type = 0, 6
      call range_integer_error('maxlevelsbin', nlevbin(type), 0, numlev, index1 = type, name1 = 'type')
    enddo
  endif
  do Zix = 0, numZ
    do Nix = 0, numN
      Z = Zinit - Zix
      A = Ainit - Zix - Nix
      call range_integer_error('nlevels', nlev(Zix, Nix), 0, numlev, index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('massnucleus', massnucleus(Zix, Nix), real(A) - 1., real(A) + 1., default = 0., &
 &      index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('massexcess', massexcess(Zix, Nix), -600., 600., default = 0., &
 &      index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
    enddo
  enddo
  if (massdir(1:1) /= ' ') then
    massdir0 = massdir
    massdir = trim(path)//'masses/'//trim(massdir0)
    massfile = trim(massdir)//'/'//'Fe.mass'
    inquire (file = massfile, exist = lexist)
    if (.not.lexist) then
      write(*, '(" TALYS-error: Non-existent mass file ", a)') trim(massfile)
      stop
    endif
    massfile = trim(massdir)//'/'//trim(Starget)//'.mass'
    inquire (file = massfile, exist = lexist)
    if (.not.lexist) write(*, '(" TALYS-warning: Non-existent mass file ",a)') trim(massfile)
  endif
  if (Lisoinp ==  - 1) call range_integer_error('Ltarget', Ltarget, 0, numlev)
  call range_integer_error('Liso', Lisoinp, 0, 9, default = -1)
  call range_real_error('isomer', isomer, 0., 1.e38, unit ='s')
  call range_integer_error('core', core, -1, 1)
  if (core ==  0) then
    write(*, '(" TALYS-error: core = -1 or 1")')
    stop
  endif
  call range_integer_error('transpower', transpower, 2, 20)
  call range_real_error('transeps', real(transeps), 0., 1.)
  call range_real_error('xseps', xseps, 0., 1000.)
  call range_real_error('popeps', popeps, 0., 1000.)
  call range_real_error('Rfiseps', Rfiseps, 0., 1.)
  call range_real_error('Elow', eninclow, 1.e-11, 1., default = 0.)
  call range_integer_error('angles', nangle, 1, numang)
  call range_integer_error('anglescont', nanglecont, 1, numangcont)
  call range_integer_error('anglesrec', nanglerec, 1, numangrec)
  call range_integer_error('maxenrec', maxenrec, 1, numenrec)
  call range_integer_error('maxchannel', maxchannel, 1, 8)
  call range_integer_error('massmodel', massmodel, 0, 3)
  call range_integer_error('disctable', disctable, 1, 3)
  call range_real_error('astroT', astroT9, 0.0001, 10., default = 0.)
  call range_real_error('astroE', astroE, 0.00001, 1., default = 0.)
  if (astroE /= 0..and.astroT9 /= 0.) then
    write(*, '(" TALYS-error: Only astroE OR astroT can be given")')
    stop
  endif
  if (astroE /= 0.) astroT9 = astroE / kT
  if (astroT9 /= 0.) astroE = astroT9 * kT
  call range_integer_error('nonthermlev', nonthermlev, 0, numlev, default = -1)
  if (flagprod) then
    if (k0 <= 1) then
      write(*, '(" TALYS-error: isotope production not yet enabled for incident photons or neutrons)")')
      stop
    endif
    if (Ebeam ==  - 1.) then
      write(*, '(" TALYS-error: accelerator energy Ebeam must be given for isotope production (production y)")')
      stop
    endif
    call range_real_error('Ebeam', Ebeam, 0., Emaxtalys, unit = 'MeV')
    if (Eback ==  - 1.) then
      Eback = max(Ebeam - 5., 0.1)
    else
      call range_real_error('Eback', Eback, 0., Emaxtalys, unit = 'MeV')
    endif
    call range_real_error('Ebeam', Ebeam, Eback, Emaxtalys, unit = 'MeV')
    if (Ebeam > enincmax + 1.e-4) then
      write(*, '(" TALYS-error: Ebeam is not in the energy range ", &
 &      "with TALYS results, Ebeam = ", f10.5, " Ein(max) = ", f10.5, ". Rerun with wider energy grid")') Ebeam, enincmax
      stop
    endif
    if (Eback < eninc(1) - 1.e-4) then
      write(*, '(" TALYS-error: Eback is not in the energy range ", &
 &      "with TALYS results, Eback = ", f10.5, " Ein(1) = ", f10.5, ". Rerun with wider energy grid")') Eback, eninc(1)
      stop
    endif
    if (radiounit /= 'bq' .and. radiounit /= 'kbq' .and. radiounit /= 'mbq' .and. radiounit /= 'gbq' .and. &
      radiounit /= 'mci' .and. radiounit /= 'ci' .and. radiounit /= 'kci') then
      write(*, '(" TALYS-error: radiounit should be equal to Bq, kBq, MBq, Gbq, mCi, Ci or kCi")')
      stop
    endif
    if (yieldunit /= 'num' .and. yieldunit /= 'mug' .and. yieldunit /= 'mg' .and. yieldunit /= 'g' .and. yieldunit /= 'kg') &
      then
      write(*, '(" TALYS-error: yieldunit should be equal to num (number), mug (micro-gram), mg, g, or kg")')
      stop
    endif
    call range_real_error('Ibeam', Ibeam, 0., 10000., unit = 'mA')
    call range_real_error('Area', Area, 0., 10000., unit = 'cm^2')
    do k = 1, 5
      call range_integer_error('Tirrad', Tirrad(k), 0, 1000000, index1 = k, name1 = 'k')
      call range_integer_error('Tcool', Tcool(k), 0, 1000000, index1 = k, name1 = 'k')
    enddo
    do k = 1, 5
      if (unitTirrad(k) /= ' ' .and. unitTirrad(k) /= 'y' .and. unitTirrad(k) /= 'd' .and. unitTirrad(k) /= 'h' .and. &
        unitTirrad(k) /= 'm' .and. unitTirrad(k) /= 's') then
        write(*, '(" TALYS-error: wrong unit for Tirrad = ", i9)') Tirrad(k)
        stop
      endif
      if (unitTcool(k) /= ' ' .and. unitTcool(k) /= 'y' .and. unitTcool(k) /= 'd' .and. unitTcool(k) /= 'h' .and. &
        unitTcool(k) /= 'm' .and. unitTcool(k) /= 's') then
        write(*, '(" TALYS-error: wrong unit for Tcool = ", i9)') Tcool(k)
        stop
      endif
    enddo
    call range_real_error('rhotarget', rhotarget, 0., 100., default = -1.)
  endif
  call range_real_error('Tres', Tres, 0., 1.e12)
!
! 3. Check of values of optical model
!
  do Zix = 0, numZph
    do Nix = 0, numNph
      do type = 1, 6
        if (optmod(Zix, Nix, type)(1:1) == ' ') cycle
        inquire (file = optmod(Zix, Nix, type), exist = lexist)
        if ( .not. lexist) then
          write(*, '(" TALYS-error: Non-existent optical model file: ", a)') trim(optmod(Zix, Nix, type))
          stop
        endif
      enddo
    enddo
    if (optmodfileN(Zix)(1:1) /= ' ') then
      inquire (file = optmodfileN(Zix), exist = lexist)
      if ( .not. lexist) then
        write(*, '(" TALYS-error: Non-existent optical model file: ", a)') trim(optmodfileN(Zix))
        stop
      endif
    endif
    if (optmodfileP(Zix)(1:1) /= ' ') then
      inquire (file = optmodfileP(Zix), exist = lexist)
      if ( .not. lexist) then
        write(*, '(" TALYS-error: Non-existent optical model file: ", a)') trim(optmodfileP(Zix))
        stop
      endif
    endif
    if (radialfile(Zix)(1:1) /= ' ') then
      inquire (file = radialfile(Zix), exist = lexist)
      if ( .not. lexist) then
        write(*, '(" TALYS-error: Non-existent radial file: ", a)') trim(radialfile(Zix))
        stop
      endif
    endif
!
! Check other parameter input files
!
    if (levelfile(Zix)(1:1) /= ' ') then
      inquire (file = levelfile(Zix), exist = lexist)
      if ( .not. lexist) then
        write(*, '(" TALYS-error: Non-existent level file: ", a)') trim(levelfile(Zix))
        stop
      endif
    endif
    if (deformfile(Zix)(1:1) /= ' ') then
      inquire (file = deformfile(Zix), exist = lexist)
      if ( .not. lexist) then
        write(*, '(" TALYS-error: Non-existent deformation parameter file: ", a)') trim(deformfile(Zix))
        stop
      endif
    endif
    do Nix = 0, numN
      do irad = 0, 1
        do l = 1, numgam
          if (Exlfile(Zix, Nix, irad, l)(1:1) /= ' ') then
            inquire (file = Exlfile(Zix, Nix, irad, l), exist = lexist)
            if ( .not. lexist) then
              write(*, '(" TALYS-error: Non-existent strength ", &
 &              "function file irad = ", i1, " l = ", i2, " : ", a)') irad, l, trim(Exlfile(Zix, Nix, irad, l))
              stop
            endif
          endif
        enddo
      enddo
      if (densfile(Zix,Nix)(1:1) /= ' ') then
        inquire (file=densfile(Zix,Nix),exist=lexist)
        if (.not.lexist) then
          write(*,'(" TALYS-error: Non-existent level density file: ",a)') trim(densfile(Zix,Nix))
          stop
        endif
        if (ctable(Zix,Nix,0) == 1.e-20) ctable(Zix,Nix,0)=0.
        if (ptable(Zix,Nix,0) == 1.e-20) ptable(Zix,Nix,0)=0.
        if (ldmodel(Zix,Nix) <= 3) then
          if (flagparity) then
            ldmodel(Zix,Nix)=5
          else
            ldmodel(Zix,Nix)=4
          endif
        endif
      endif
      if (hbtransfile(Zix, Nix)(1:1) /= ' ') then
        inquire (file = hbtransfile(Zix, Nix), exist = lexist)
        if ( .not. lexist) then
          write(*, '(" TALYS-error: Non-existent head band transition state file: ", a)') trim(hbtransfile(Zix, Nix))
          stop
        endif
      endif
      if (clas2file(Zix, Nix)(1:1) /= ' ') then
        inquire (file = clas2file(Zix, Nix), exist = lexist)
        if ( .not. lexist) then
          write(*, '(" TALYS-error: Non-existent class 2 transition state file: ", a)') trim(clas2file(Zix, Nix))
          stop
        endif
      endif
    enddo
  enddo
  if (ompenergyfile(1:1) /= ' ') then
    inquire (file = ompenergyfile, exist = lexist)
    if ( .not. lexist) then
      write(*, '(" TALYS-error: Non-existent ompenergyfile: ", a)') trim(ompenergyfile)
      stop
    endif
  endif
  do mt = 1, nummt
    do is = - 1, numisom
      if (rescuefile(mt, is)(1:1) /= ' ') then
        inquire (file = rescuefile(mt, is), exist = lexist)
        if ( .not. lexist) then
          write(*, '(" TALYS-error: Non-existent rescue file: ", a)') trim(rescuefile(mt, is))
          stop
        endif
      endif
      call range_real_error('grescue', grescue(mt, is), 0.001, 1000., index1 = mt, name1 = 'MT', index2 = is, name2 = 'iso')
    enddo
  enddo
  call range_integer_error('alphaomp', alphaomp, 1, 8)
  call range_integer_error('deuteronomp', deuteronomp, 1, 5)
  call range_integer_error('radialmodel', radialmodel, 1, 2)
!
! Check adjustable OMP parameters
!
  do type = 0, 6
    call range_real_error('rvadjust', rvadjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('avadjust', avadjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('v1adjust', v1adjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('v3adjust', v2adjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('v3adjust', v3adjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('v4adjust', v4adjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('rwadjust', rwadjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('awadjust', awadjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('w1adjust', w1adjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('w2adjust', w2adjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('w3adjust', w3adjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('w4adjust', w4adjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('rvdadjust', rvdadjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('avdadjust', avdadjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('d1adjust', d1adjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('d2adjust', d2adjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('d3adjust', d3adjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('rwdadjust', rwdadjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('awdadjust', awdadjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('rvsoadjust', rvsoadjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('avsoadjust', avsoadjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('vso1adjust', vso1adjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('vso2adjust', vso2adjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('rwsoadjust', rwsoadjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('awsoadjust', awsoadjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('wso1adjust', wso1adjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('wso2adjust', wso2adjust(type), 0.1, 10., index1 = type, name1 = 'type')
    call range_real_error('rcadjust', rcadjust(type), 0.1, 10., index1 = type, name1 = 'type')
    do omptype = 1, numompadj
      do nr = 1, ompadjustN(type, omptype)
        call range_real_error('ompadjustE1', ompadjustE1(type, omptype, nr), 0., Emaxtalys, &
 &        index1 = type, name1 = 'type', index2 = omptype, name2 = 'omptype', index3 = nr, name3 = 'nr')
        call range_real_error('ompadjustE2', ompadjustE2(type, omptype, nr), 0., Emaxtalys, &
 &        index1 = type, name1 = 'type', index2 = omptype, name2 = 'omptype', index3 = nr, name3 = 'nr')
        call range_real_error('ompadjustE2', ompadjustE2(type, omptype, nr),  ompadjustE1(type, omptype, nr), Emaxtalys, &
 &        index1 = type, name1 = 'type', index2 = omptype, name2 = 'omptype', index3 = nr, name3 = 'nr')
        do nr2 = 1, ompadjustN(type, omptype)
          if (nr == nr2) cycle
          if (ompadjustE1(type, omptype, nr) > ompadjustE1(type, omptype, nr2) .and. &
            ompadjustE1(type, omptype, nr) < ompadjustE2(type, omptype, nr2)) then
            write(*, '(" TALYS-error: ompadjustE1 and ompadjustE2 overlapping")')
            stop
          endif
        enddo
        call range_real_error('ompadjustD', ompadjustD(type, omptype, nr), -100., 100., &
 &        index1 = type, name1 = 'type', index2 = omptype, name2 = 'omptype', index3 = nr, name3 = 'nr')
        call range_real_error('ompadjusts', ompadjusts(type, omptype, nr), 0., 100., &
 &        index1 = type, name1 = 'type', index2 = omptype, name2 = 'omptype', index3 = nr, name3 = 'nr')
      enddo
    enddo
  enddo
  call range_integer_error('jlmmode', jlmmode, 0, 3)
  call range_real_error('lvadjust', lvadjust, 0.5, 1.5)
  call range_real_error('lwadjust', lwadjust, 0.5, 1.5)
  call range_real_error('lv1adjust', lv1adjust, 0.5, 1.5)
  call range_real_error('lw1adjust', lw1adjust, 0.5, 1.5)
  call range_real_error('lvsoadjust', lvsoadjust, 0.5, 1.5)
  call range_real_error('lwsoadjust', lwsoadjust, 0.5, 1.5)
  call range_real_error('aradialcor', aradialcor, 0.5, 1.5)
  call range_real_error('adepthcor', adepthcor, 0.5, 1.5)
  call range_real_error('soswitch', soswitch, 0.1, 10.)
  do type = 1, 2
    call range_real_error('Ejoin', Ejoin(type), 0., Emaxtalys, index1 = type, name1 = 'type')
    call range_real_error('Vinfadjust', Vinfadjust(type), 0.01, 10., index1 = type, name1 = 'type')
  enddo
  call range_integer_error('pruittset', pruittset, 0, 416)
!
! Check direct reaction parameters
!
  call range_integer_error('maxband', maxband, 0, 10)
  call range_integer_error('maxrot', maxrot, 0, 20)
  if (k0 == 0 .and. flaggiant0) then
    write(*, '(" TALYS-error: No giant resonance sumrules for photonuclear reactions")')
    stop
  endif
!
! 4. Check of values for compound nucleus
!
  call range_real_error('ewfc', ewfc, 0., 20., default = -1.)
  call range_real_error('eurr', eurr, 0., 20., default = -1.)
  call range_integer_error('wmode', wmode, 0, 3)
  call range_integer_error('wfcfactor', wfcfactor, 1, 3)
  call range_integer_error('lurr', lurr, 0, numl)
  if (k0 == 0 .and. ewfc > 0.) then
    write(*, '(" TALYS-error: No width fluctuations for photonuclear reactions")')
    stop
  endif
  if (k0 /= 1 .and. flagres) then
    write(*, '(" TALYS-error: resonance calculation only possible for incident neutrons")')
    stop
  endif
  if (k0 /= 1 .and. (eurr > 0..or.flagurr)) then
    write(*, '(" TALYS-error: URR calculation only possible for incident neutrons")')
    stop
  endif
  if ( .not. flagcomp .and. eurr > 0.) then
    write(*, '(" TALYS-error: URR calculation only possible if compound nucleus model enabled")')
    stop
  endif
  if (k0 == 0..and.flageciscomp) then
    write(*, '(" TALYS-error: No compound calculation by ECIS for incident photons")')
    stop
  endif
  if (enincmax > 20..and.flageciscomp) then
    write(*, '(" TALYS-error: No compound calculation by ECIS for E > 20 MeV")')
    stop
  endif
!
! 5. Check of values for gamma emission
!
  call range_integer_error('gammax', gammax, 1, 6)
  call range_integer_error('strength', strength, 1, 12)
  if (strength == 3 .or. strength == 4) then
    inquire (file = trim(path)//'gamma/hfb/Sn.psf', exist = lexist)
    if ( .not. lexist) then
      write(*, '(" TALYS-error: Microscopic HFB tables are not installed: download the full TALYS package from www.talys.eu")')
      stop
    endif
  endif
  if ((strengthM1 < 1 .or. strengthM1 > 4) .and. strengthM1 /= 8 .and. strengthM1 /= 10 .and. strengthM1 /= 11 .and. &
 &  strengthM1 /= 12) then
    write(*,'(" TALYS-error: strengthM1 = 1, 2, 3, 4, 8, 10, 11 or 12")')
    stop
  endif
  do Zix = 0, numZ
    do Nix = 0, numN
      Z = Zinit - Zix
      A = Ainit - Zix - Nix
      do irad = 0, 1
        do lval = 1, gammax
          call range_real_error('etable', etable(Zix, Nix, irad, lval), -10., 10., &
 &          index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L')
          call range_real_error('ftable', ftable(Zix, Nix, irad, lval), 0.1, 10., &
 &          index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L')
          call range_real_error('wtable', wtable(Zix, Nix, irad, lval), 0.1, 10., &
 &          index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L')
          call range_real_error('etableadjust', etableadjust(Zix, Nix, irad, lval), -10., 10., &
 &          index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L')
          call range_real_error('ftableadjust', ftableadjust(Zix, Nix, irad, lval), 0.1, 10., &
 &          index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L')
          call range_real_error('wtableadjust', wtableadjust(Zix, Nix, irad, lval), 0.1, 10., &
 &          index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L')
          do igr = 1, 2
            call range_real_error('energy of GR', egr(Zix, Nix, irad, lval, igr), 1., 100., default = 0., index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('width of GR', ggr(Zix, Nix, irad, lval, igr), 0.5, 100., default = 0., index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('strength of GR', sgr(Zix, Nix, irad, lval, igr), 0., 10000., index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('energy of PR', epr(Zix, Nix, irad, lval, igr), 1., 100., default = 0., index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('width of PR', gpr(Zix, Nix, irad, lval, igr), 0.1, 100., default = 0., index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('strength of PR', tpr(Zix, Nix, irad, lval, igr), 0., 10000., index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('egradjust', egradjust(Zix, Nix, irad, lval, igr), 0.05, 20., index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('ggradjust', ggradjust(Zix, Nix, irad, lval, igr), 0.05, 20., index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('sgradjust', sgradjust(Zix, Nix, irad, lval, igr), 0.05, 20., index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('epradjust', epradjust(Zix, Nix, irad, lval, igr), 0.05, 20., index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('gpradjust', gpradjust(Zix, Nix, irad, lval, igr), 0.05, 20., index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('spradjust', tpradjust(Zix, Nix, irad, lval, igr), 0.05, 20., index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
          enddo
          call range_real_error('upbendc', upbend(Zix, Nix, irad, lval, 1), 0., 1.e-5, index1 = Z, name1 = 'Z', &
 &          index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = 1, name5 = 'igr')
          call range_real_error('upbende', upbend(Zix, Nix, irad, lval, 2), 0., 10., index1 = Z, name1 = 'Z', &
 &          index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = 2, name5 = 'igr')
          call range_real_error('upbendf', upbend(Zix, Nix, irad, lval, 3), -10., 10., index1 = Z, name1 = 'Z', &
 &          index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = 2, name5 = 'igr')
          call range_real_error('upbendcadjust', upbendadjust(Zix, Nix, irad, lval, 1), 0.05, 20., index1 = Z, name1 = 'Z', &
 &          index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = 1, name5 = 'igr')
          call range_real_error('upbendeadjust', upbendadjust(Zix, Nix, irad, lval, 2), 0.05, 20., index1 = Z, name1 = 'Z', &
 &          index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = 2, name5 = 'igr')
          call range_real_error('upbendfadjust', upbendadjust(Zix, Nix, irad, lval, 3), 0.05, 20., index1 = Z, name1 = 'Z', &
 &          index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = 2, name5 = 'igr')
        enddo
      enddo
      call range_real_error('gamgam', gamgam(Zix, Nix), 0., 10., default = 0., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('D0', D0(Zix, Nix), 1.e-3, 1.e7, default = 0., unit = 'keV', index1 = Z, name1 = 'Z', &
 &      index2 = A, name2 = 'A')
      do type = -1, 6
        call range_real_error('fiso', fiso(type), 0.01, 100., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', &
 &        index3 = type, name3 = 'type', default = -1.)
        call range_real_error('fisom', fisom(type), 0.01, 100., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', &
 &        index3 = type, name3 = 'type', default = -1.)
      enddo
      call range_real_error('gamgamadjust', gamgamadjust(Zix, Nix), 0.01, 20., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
    enddo
  enddo
  call range_real_error('RprimeU', RprimeU, 0., 10.)
  if (flagracap .and. k0 == 0) then
    write(*, '(" TALYS-error: Radiative capture model not possible for incident photons")')
    stop
  endif
  call range_integer_error('ldmodelracap', ldmodelracap, 1, 3)
  do Zix = 0, numZ
    do Nix = 0, numN
      call range_real_error('sfth', spectfacth(Zix, Nix), 0., 10., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      do i = 0, numlev
        call range_real_error('sfexp', spectfacexp(Zix, Nix, i), 0., 10., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', &
 &        index3 = i, name3 = 'level')
      enddo
    enddo
  enddo
!
! 6. Check of values for pre-equilibrium
!
  call range_real_error('epreeq', epreeq, 0., Emaxtalys, default = -1.)
  call range_integer_error('preeqmode', preeqmode, 1, 4)
  call range_integer_error('mpreeqmode', mpreeqmode, 1, 2)
  call range_integer_error('breakupmodel', breakupmodel, 1, 2)
  call range_integer_error('phmodel', phmodel, 1, 2)
  call range_integer_error('pairmodel', pairmodel, 1, 2)
  call range_integer_error('pespinmodel', pespinmodel, 1, 4)
  call range_real_error('emulpre', emulpre, 0., Emaxtalys)
  call range_real_error('M2constant', M2constant, 0., 100.)
  call range_real_error('M2limit', M2limit, 0., 100.)
  call range_real_error('M2shift', M2shift, 0., 100.)
  call range_real_error('Rpipi', Rpipi, 0., 100.)
  call range_real_error('Rnunu', Rnunu, 0., 100.)
  call range_real_error('Rpinu', Rpinu, 0., 100.)
  call range_real_error('Rnupi', Rnupi, 0., 100.)
  call range_real_error('Rgamma', Rgamma, 0., 100.)
  call range_real_error('Esurf', Esurf0, 0., 38., default = -1.)
  call range_integer_error('msdbins', msdbins, 2, numenmsd/2 - 1, default = 0)
  call range_real_error('E-in', Emsdmin, 0., Emaxtalys)
  call range_real_error('elwidth', elwidth, 1.e-6, 100.)
  call range_real_error('xscaptherm', xscaptherm(-1), 1.e-20, 1.e10, default = 0.)
  call range_real_error('xsptherm', xsptherm(-1), 1.e-20, 1.e10, default = 0.)
  call range_real_error('xsalphatherm', xsalphatherm(-1), 1.e-20, 1.e10, default = 0.)
  if (k0 == 0 .and. flagpecomp) then
    write(*, '(" TALYS-error: No pick-up and knock-out mechanism for photonuclear reactions")')
    stop
  endif
  do type = 0, 6
    call range_real_error('Cstrip', Cstrip(type), 0., 100., index1 = type, name1 = 'type')
    call range_real_error('Cknock', Cknock(type), 0., 100., index1 = type, name1 = 'type')
    call range_real_error('Cbreak', Cbreak(type), 0., 100., index1 = type, name1 = 'type')
  enddo
!
! 7. Check of values for level densities
!
  call range_integer_error('spincutmodel', spincutmodel, 1, 2)
  call range_integer_error('shellmodel', shellmodel, 1, 2)
  call range_integer_error('kvibmodel', kvibmodel, 1, 2)
  call range_integer_error('ldmodelcn', ldmodelCN, 1, 8)
  do Zix = 0, numZ
    do Nix = 0, numN
      Z = Zinit - Zix
      A = Ainit - Zix - Nix
      call range_integer_error('ldmodel', ldmodel(Zix, Nix), 1, 8, index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      if (ldmodel(Zix, Nix) >= 4) then
        inquire (file = trim(path) // 'density/ground/hilaire/Sn.tab', exist = lexist)
        if ( .not. lexist) then
          write(*, '(" TALYS-error: Microscopic HFB tables are not installed: download the full TALYS package from www.talys.eu")')
          stop
        endif
      endif
      call range_real_error('a', alev(Zix, Nix), 1., 100., default = 0.,  index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('alimit', alimit(Zix, Nix), 1., 100., default = 0.,  index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('gammald', gammald(Zix, Nix), 0., 1., default = -1.,  index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('risomer', Risomer(Zix, Nix), 0.1, 10., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      do ibar = 0, numbar
        call range_real_error('deltaW', deltaW(Zix, Nix, ibar), -20., 20., default = 0.,  index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_integer_error('Nlow', Nlow(Zix, Nix, ibar), 0, 200, default = -1,  index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_integer_error('Ntop', Ntop(Zix, Nix, ibar), 0, 200, default = -1,  index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_integer_error('Ntop', Ntop(Zix, Nix, ibar), Nlow(Zix, Nix, ibar), 200, default = -1,  index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('E0', E0(Zix, Nix, ibar), -15., 15., default = 0.,  index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('beta2', beta2(Zix, Nix, ibar), -0.5, 1.5, index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('s2adjust', s2adjust(Zix, Nix, ibar), 0.02, 50., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('Krotconstant', Krotconstant(Zix, Nix, ibar), 0.001, 1000., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('Ufermi', cfermi(Zix, Nix, ibar), 0., 1000., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('cfermi', cfermi(Zix, Nix, ibar), 0., 1000., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('T', T(Zix, Nix, ibar), 1.e-3, 10., default = 0.,  index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('Exmatch', Exmatch(Zix, Nix, ibar), 0.05, 20., default = 0.,  index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('Tadjust', Tadjust(Zix, Nix, ibar), 0.05, 20., default = 0.,  index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('E0adjust', E0adjust(Zix, Nix, ibar), 0.02, 50., default = 0.,  index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('Exmatchadjust', Exmatchadjust(Zix, Nix, ibar), 0.2, 2., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('Pshift', Pshift(Zix, Nix, ibar), -10., 10., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('Pshiftadjust', Pshiftadjust(Zix, Nix, ibar), -10., 10., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('ctable', ctable(Zix, Nix, ibar), -10., 10., default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('ptable', ptable(Zix, Nix, ibar), -10., 10., default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('ctableadjust', ctableadjust(Zix, Nix, ibar), -10., 10., default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('ptableadjust', ptableadjust(Zix, Nix, ibar), -10., 10., default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
      enddo
      call range_real_error('aadjust', aadjust(Zix, Nix), 0.1, 10., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('gadjust', gadjust(Zix, Nix), 0.1, 10., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('gnadjust', gnadjust(Zix, Nix), 0.1, 10., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('gpadjust', gpadjust(Zix, Nix), 0.1, 10., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('pair', pair(Zix, Nix), -10., 10., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('g', g(Zix, Nix), 0.1, 100., default = 0., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('gn', gn(Zix, Nix), 0.1, 100., default = 0., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('gp', gp(Zix, Nix), 0.1, 100., default = 0., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('alphald', alphald(Zix, Nix), 0.01, 0.2, index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('betald', betald(Zix, Nix), -0.5, 0.5, index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      if (betald(Zix, Nix) < 0..and. abs(betald(Zix, Nix)) > alphald(Zix, Nix)) then
        write(*, '(" TALYS-error: if betald<0, |betald|<alphald")')
        stop
      endif
      call range_real_error('gammashell1', gammashell1(Zix, Nix), 0., 1., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('Pshiftconstant', Pshiftconstant(Zix, Nix), -5., 5., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
    enddo
  enddo
  call range_real_error('cglobal', cglobal, -10., 10., default = 0.)
  call range_real_error('pglobal', pglobal, -10., 10., default = 0.)
!
! There are many input possibilities for the energy dependent level density parameter of the Ignatyuk formula.
! The required parameters are alev, alimit, gammald and deltaW.
! The Ignatyuk formula implies that they can not all be given at the same time in the input file.
!
  do Zix = 0, numZ
    do Nix = 0, numN
      do ibar = 0, numbar
        if (alev(Zix, Nix) /= 0 .and. deltaW(Zix, Nix, ibar) /= 0..and. &
          alimit(Zix, Nix) /= 0 .and. gammald(Zix, Nix) /=  -1.) then
          write(*, '(" TALYS-error: Level density conflict - a, deltaW, alimit and gammald are ALL given", &
 &          " in the input for Z=", i3, " A=", i3, " fission barrier = ", i3)') Zinit-Zix, Ainit-Zix-Nix, ibar
          stop
        endif
      enddo
    enddo
  enddo
  call range_real_error('gammashell2', gammashell2, 0., 0.2)
  call range_real_error('pairconstant', pairconstant, 0., 30.)
  call range_real_error('Kph', Kph, 1., 100.)
  call range_real_error('Rspincut', Rspincut, 0., 10.)
  call range_real_error('Rspincutpreeq', Rspincutpreeq, 0., 10.)
  call range_real_error('Rspincutff', Rspincut, 0., 20.)
!
! 8. Check of values for fission
!
  if ((flagfission .or. flagfisout) .and. Atarget <= 150) then
    write(*, '(" TALYS-error: Fission not allowed for A <= 150")')
    stop
  endif
  if (flagfission .and. flagmassdis .and. flagnatural) then
    write(*, '(" TALYS-error: Fission yield calculation not possible for natural targets")')
    stop
  endif
  call range_integer_error('fismodel', fismodel, 1, 5)
  call range_integer_error('fismodelalt', fismodelalt, 3, 4)
  if (fismodel /= 5 .and. flagfispartdamp) then
    write(*,'(" TALYS-error: Fission partial damping only allowed for fismodel 5")')
    stop
  endif
  call range_integer_error('fymodel', fymodel, 1, 5)
  call range_integer_error('ffmodel', ffmodel, 0, 4)
  call range_integer_error('pfnsmodel', pfnsmodel, 1, 2)
  call range_integer_error('gefran', gefran, 1000, 1000000)
  call range_real_error('Cnubar1', Cnubar1, 0.1, 10.)
  call range_real_error('Cnubar2', Cnubar2, 0.1, 10.)
  call range_real_error('Tmadjust', Tmadjust, 0.1, 10.)
  call range_real_error('Fsadjust', Fsadjust, 0.1, 10.)
  call range_real_error('Cbarrier', Cbarrier, 0.1, 10.)
  if (yieldfile(1:1) /= ' ') then
    inquire (file = yieldfile, exist = lexist)
    if ( .not. lexist) then
      write(*, '(" TALYS-error: Non-existent yieldfile: ", a)') trim(yieldfile)
      stop
    endif
  endif
  do Zix = 0, numZ
    do Nix = 0, numN
      do ibar = 1, numbar
        call range_integer_error('type of axiality', axtype(Zix, Nix, ibar), 1, 5, default = 0, index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('fission barrier', fbarrier(Zix, Nix, ibar), 0., 100., default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('fisbaradjust', fbaradjust(Zix, Nix, ibar), 0.02, 50., default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('fission width', fwidth(Zix, Nix, ibar), 0.01, 10., default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('fwidthadjust', fwidthadjust(Zix, Nix, ibar), 0.02, 50., default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('bdamp', bdamp(Zix, Nix, ibar), 0., 50., default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('bdampadjust', bdampadjust(Zix, Nix, ibar), 0.01, 100., default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('Rtransmom', Rtransmom(Zix, Nix, ibar), 0.05, 20., default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('Rclass2mom', Rclass2mom(Zix, Nix, ibar), 0.05, 20., default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
      enddo
      call range_real_error('betafiscor', betafiscor(Zix, Nix), 0.05, 20., default = 0., index1 = Z, name1 = 'Z', &
 &      index2 = A, name2 = 'A')
      call range_real_error('vfiscor', vfiscor(Zix, Nix), 0.05, 20., default = 0., index1 = Z, name1 = 'Z', &
 &      index2 = A, name2 = 'A')
      call range_real_error('betafiscoradjust', betafiscoradjust(Zix, Nix), 0.1, 10., default = 0., index1 = Z, name1 = 'Z', &
 &      index2 = A, name2 = 'A')
      call range_real_error('vfiscoradjust', vfiscoradjust(Zix, Nix), 0.1, 10., default = 0., index1 = Z, name1 = 'Z', &
 &      index2 = A, name2 = 'A')
    enddo
  enddo
!
! 9. Check of values for output
!
  call range_real_error('eadd', eadd, 0., Emaxtalys, unit = 'MeV')
  call range_real_error('eaddel', eaddel, 0., Emaxtalys, unit = 'MeV')
  call range_integer_error('ddxmode', ddxmode, 0, 3)
  do type = 0, 6
    do i = 1, ddxecount(type)
      call range_real_error('fileddxe', fileddxe(type, i), 0., enincmax, default = 0., index1 = type, name1 = 'type', &
 &      index2 = i, name2 = 'i')
    enddo
    do i = 1, ddxacount(type)
      call range_real_error('fileddxa', fileddxa(type, i), 0., 180., default = 0., index1 = type, name1 = 'type', &
 &      index2 = i, name2 = 'i')
    enddo
  enddo
  if (flagdecay) flagpop = .true.
!
! 10. Check of values energy-dependent parameter adjustment
!
  do n = 1, Nadjust
    if (adjustfile(i)(1:1) /= ' ') cycle
    Ea = adjustpar(n, 1)
    Eb = adjustpar(n, 2)
    Em = adjustpar(n, 3)
    D = adjustpar(n, 4)
    if ((Ea >= Eb) .or. (Ea >= Em) .or. (Em >= Eb)) then
      write(*, '(" TALYS-error: energy range for adjustment should", &
 &      " be given as follows: Ea Eb Em D, with Ea < Em < Eb for keyword ", a)') trim(adjustkey(n))
      stop
    endif
    call range_real_error(adjustkey(n), D, 0., 10., index1 = n, name1 = 'n')
    do m = 1, Nadjust
      if (m == n) cycle
      if (adjustkey(m) /= adjustkey(n)) cycle
      Ea2 = adjustpar(m, 1)
      Eb2 = adjustpar(m, 2)
      if ((Ea2 > Ea .and. Ea2 < Eb) .or. (Eb2 > Ea .and. Eb2 < Eb)) then
        write(*, '(" TALYS-error: overlapping energy ranges for keyword ", a)') trim(adjustkey(n))
        stop
      endif
    enddo
  enddo
!
! 11. Check for correct name of libraries for resonance parameters
!
  if (flagres) then
    if (trim(reslib) == 'tendl.2023') return
    if (trim(reslib) == 'jeff3.3') return
    if (trim(reslib) == 'endfb8.1') return
    if (trim(reslib) == 'cendl3.2') return
    if (trim(reslib) == 'jendl5.0') return
    write(*, '(" TALYS-error: Wrong library name: ", a)') trim(reslib)
    stop
  endif
  return
end subroutine checkvalue
! Copyright A.J. Koning 2021
