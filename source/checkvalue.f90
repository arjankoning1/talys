subroutine checkvalue
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Check for errors in values
!
! Author    : Arjan Koning
!
! 2026-03-26: Current version
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
!   ecisstep          ! integration step size for ECIS OMP calculation
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
  call loadkeyranges
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
  if (Lisoinp ==  -1) then
    call range_integer_error('maxlevelstar', nlevmax, nlevmaxLIM(1), nlevmaxLIM(2))
    call range_integer_error('maxlevelsres', nlevmaxres, nlevmaxresLIM(1), nlevmaxresLIM(2))
    do type = 0, 6
      call range_integer_error('maxlevelsbin', nlevbin(type), nlevbinLIM(1), nlevbinLIM(2), index1 = type, name1 = 'type')
    enddo
  endif
  call range_real_error('emaxpseudores', emaxpseudores, EmaxpseudoresLIM(1), EmaxpseudoresLIM(2), default = -1.)
  call range_real_error('pseudoreswidth', pseudoreswidth, pseudoreswidthLIM(1), pseudoreswidthLIM(2))
  call range_real_error('pseudoresfade', pseudoresfade, pseudoresfadeLIM(1), pseudoresfadeLIM(2))
  do Zix = 0, numZ
    do Nix = 0, numN
      Z = Zinit - Zix
      A = Ainit - Zix - Nix
      call range_integer_error('nlevels', nlev(Zix, Nix), nlevLIM(1), nlevLim(2), index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('massnucleus', massnucleus(Zix, Nix), massnucleusLIM(1), massnucleusLIM(2), default = 0., &
 &      index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('massexcess', massexcess(Zix, Nix), massexcessLIM(1), massexcessLIM(2), default = 0., &
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
  if (Lisoinp ==  -1) call range_integer_error('Ltarget', Ltarget, LtargetLIM(1), LtargetLIM(2))
  call range_integer_error('Liso', Lisoinp, LisoinpLIM(1), LisoinpLIM(2), default = -1)
  call range_real_error('isomer', isomer, isomerLIM(1), isomerLIM(2), unit ='s')
  call range_integer_error('core', core, coreLIM(1), coreLIM(2))
  if (core ==  0) then
    write(*, '(" TALYS-error: core = -1 or 1")')
    stop
  endif
  call range_integer_error('transpower', transpower, transpowerLIM(1), transpowerLIM(2))
  call range_real_error('transeps', real(transeps), transepsLIM(1), transepsLIM(2))
  call range_real_error('xseps', xseps, xsepsLIM(1), xsepsLIM(2))
  call range_real_error('popeps', popeps, popepsLIM(1), popepsLIM(2))
  call range_real_error('Rfiseps', Rfiseps, RfisepsLIM(1), RfisepsLIM(2))
  call range_real_error('Elow', eninclow, eninclowLIM(1), eninclowLIM(2), default = 0.)
  call range_integer_error('angles', nangle, nangleLIM(1), nangleLIM(2))
  call range_integer_error('anglescont', nanglecont, nanglecontLIM(1), nanglecontLIM(2))
  call range_integer_error('anglesrec', nanglerec, nanglerecLIM(1), nanglerecLIM(2))
  call range_integer_error('maxenrec', maxenrec, maxenrecLIM(1), maxenrecLIM(2))
  call range_integer_error('maxchannel', maxchannel, maxchannelLIM(1), maxchannelLIM(2))
  call range_integer_error('massmodel', massmodel, massmodelLIM(1), massmodelLIM(2))
  call range_integer_error('disctable', disctable, disctableLIM(1), disctableLIM(2))
  call range_real_error('astroT', astroT9, astroT9LIM(1), astroT9LIM(2), default = 0.)
  call range_real_error('astroE', astroE, astroELIM(1), astroELIM(2), default = 0.)
  if (astroE /= 0..and.astroT9 /= 0.) then
    write(*, '(" TALYS-error: Only astroE OR astroT can be given")')
    stop
  endif
  if (astroE /= 0.) astroT9 = astroE / kT
  if (astroT9 /= 0.) astroE = astroT9 * kT
  call range_integer_error('nonthermlev', nonthermlev, nonthermlevLIM(1), nonthermlevLIM(2), default = -1)
  if (flagprod) then
    if (k0 <= 1) then
      write(*, '(" TALYS-error: isotope production not yet enabled for incident photons or neutrons)")')
      stop
    endif
    if (Ebeam ==  -1.) then
      write(*, '(" TALYS-error: accelerator energy Ebeam must be given for isotope production (production y)")')
      stop
    endif
    call range_real_error('Ebeam', Ebeam, EbeamLIM(1), EbeamLIM(2), unit = 'MeV')
    if (Eback ==  -1.) then
      Eback = max(Ebeam - 5., 0.1)
    else
      call range_real_error('Eback', Eback, EbackLIM(1), EbackLIM(2), unit = 'MeV')
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
    call range_real_error('Ibeam', Ibeam, IbeamLIM(1), IbeamLIM(2), unit = 'mA')
    call range_real_error('Area', Area, AreaLIM(1), AreaLIM(2), unit = 'cm^2')
    do k = 1, 5
      call range_integer_error('Tirrad', Tirrad(k), TirradLIM(1), TirradLIM(2), index1 = k, name1 = 'k')
      call range_integer_error('Tcool', Tcool(k), TcoolLIM(1), TcoolLIM(2), index1 = k, name1 = 'k')
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
    call range_real_error('rhotarget', rhotarget, rhotargetLIM(1), rhotargetLIM(2), default = -1.)
  endif
  call range_real_error('Tres', Tres, TresLIM(1), TresLIM(2))
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
      call range_real_error('grescue', grescue(mt, is), grescueLIM(1), grescueLIM(2), index1 = mt, name1 = 'MT', index2 = is, name2 = 'iso')
    enddo
  enddo
  call range_integer_error('alphaomp', alphaomp, alphaompLIM(1), alphaompLIM(2))
  call range_integer_error('deuteronomp', deuteronomp, deuteronompLIM(1), deuteronompLIM(2))
  call range_integer_error('radialmodel', radialmodel, radialmodelLIM(1), radialmodelLIM(2))
!
! Check adjustable OMP parameters
!
  do type = 0, 6
    call range_real_error('rvadjust', rvadjust(type), rvadjustLIM(1), rvadjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('avadjust', avadjust(type), avadjustLIM(1), avadjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('v1adjust', v1adjust(type), v1adjustLIM(1), v2adjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('v2adjust', v2adjust(type), v2adjustLIM(1), v2adjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('v3adjust', v3adjust(type), v3adjustLIM(1), v3adjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('v4adjust', v4adjust(type), v4adjustLIM(1), v4adjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('rwadjust', rwadjust(type), rwadjustLIM(1), rwadjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('awadjust', awadjust(type), awadjustLIM(1), awadjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('w1adjust', w1adjust(type), w1adjustLIM(1), w1adjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('w2adjust', w2adjust(type), w2adjustLIM(1), w2adjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('w3adjust', w3adjust(type), w3adjustLIM(1), w3adjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('w4adjust', w4adjust(type), w4adjustLIM(1), w4adjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('rvdadjust', rvdadjust(type), rvdadjustLIM(1), rvdadjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('avdadjust', avdadjust(type), avdadjustLIM(1), avdadjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('d1adjust', d1adjust(type), d1adjustLIM(1), d1adjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('d2adjust', d2adjust(type), d2adjustLIM(1), d2adjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('d3adjust', d3adjust(type), d3adjustLIM(1), d3adjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('rwdadjust', rwdadjust(type), rwdadjustLIM(1), rwdadjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('awdadjust', awdadjust(type), awdadjustLIM(1), awdadjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('rvsoadjust', rvsoadjust(type), rvsoadjustLIM(1), rvsoadjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('avsoadjust', avsoadjust(type), avsoadjustLIM(1), avsoadjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('vso1adjust', vso1adjust(type), vso1adjustLIM(1), vso1adjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('vso2adjust', vso2adjust(type), vso2adjustLIM(1), vso2adjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('rwsoadjust', rwsoadjust(type), rwsoadjustLIM(1), rwsoadjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('awsoadjust', awsoadjust(type), awsoadjustLIM(1), awsoadjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('wso1adjust', wso1adjust(type), wso1adjustLIM(1), wso1adjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('wso2adjust', wso2adjust(type), wso2adjustLIM(1), wso2adjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('rcadjust', rcadjust(type), rcadjustLIM(1), rcadjustLIM(2), index1 = type, name1 = 'type')
    call range_real_error('ecisstep', ecisstep, ecisstepLIM(1), ecisstepLIM(2), default = 0.)
    do omptype = 1, numompadj
      do nr = 1, ompadjustN(type, omptype)
        call range_real_error('ompadjustE1', ompadjustE1(type, omptype, nr), ompadjustE1LIM(1), ompadjustE1LIM(2), &
 &        index1 = type, name1 = 'type', index2 = omptype, name2 = 'omptype', index3 = nr, name3 = 'nr')
        call range_real_error('ompadjustE2', ompadjustE2(type, omptype, nr), ompadjustE2LIM(1), ompadjustE2LIM(2), &
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
        call range_real_error('ompadjustD', ompadjustD(type, omptype, nr), ompadjustDLIM(1), ompadjustDLIM(2), &
 &        index1 = type, name1 = 'type', index2 = omptype, name2 = 'omptype', index3 = nr, name3 = 'nr')
        call range_real_error('ompadjusts', ompadjusts(type, omptype, nr), ompadjustsLIM(1), ompadjustsLIM(2), &
 &        index1 = type, name1 = 'type', index2 = omptype, name2 = 'omptype', index3 = nr, name3 = 'nr')
      enddo
    enddo
  enddo
  call range_integer_error('jlmmode', jlmmode, jlmmodeLIM(1), jlmmodeLIM(2))
  call range_real_error('lvadjust', lvadjust, lvadjustLIM(1), lvadjustLIM(2))
  call range_real_error('lwadjust', lwadjust, lwadjustLIM(1), lwadjustLIM(2))
  call range_real_error('lv1adjust', lv1adjust, lv1adjustLIM(1), lv1adjustLIM(2))
  call range_real_error('lw1adjust', lw1adjust, lw1adjustLIM(1), lw1adjustLIM(2))
  call range_real_error('lvsoadjust', lvsoadjust, lvsoadjustLIM(1), lvsoadjustLIM(2))
  call range_real_error('lwsoadjust', lwsoadjust, lwsoadjustLIM(1), lwsoadjustLIM(2))
  call range_real_error('aradialcor', aradialcor, aradialcorLIM(1), aradialcorLIM(2))
  call range_real_error('adepthcor', adepthcor, adepthcorLIM(1), adepthcorLIM(2))
  call range_real_error('soswitch', soswitch, soswitchLIM(1), soswitchLIM(2))
  do type = 1, 2
    call range_real_error('Ejoin', Ejoin(type), EjoinLIM(1), EjoinLIM(2), index1 = type, name1 = 'type')
    call range_real_error('Vinfadjust', Vinfadjust(type), VinfadjustLIM(1), VinfadjustLIM(2), index1 = type, name1 = 'type')
  enddo
  call range_integer_error('pruittset', pruittset, pruittsetLIM(1), pruittsetLIM(2))
!
! Check direct reaction parameters
!
  call range_integer_error('maxband', maxband, maxbandLIM(1), maxbandLIM(2))
  call range_integer_error('maxrot', maxrot, maxrotLIM(1), maxrotLIM(2))
  if (k0 == 0 .and. flaggiant0) then
    write(*, '(" TALYS-error: No giant resonance sumrules for photonuclear reactions")')
    stop
  endif
!
! 4. Check of values for compound nucleus
!
  call range_real_error('ewfc', ewfc, ewfcLIM(1), ewfcLIM(2), default = -1.)
  call range_real_error('eurr', eurr, eurrLIM(1), eurrLIM(2), default = -1.)
  call range_integer_error('wmode', wmode, wmodeLIM(1), wmodeLIM(2))
  call range_integer_error('wfcfactor', wfcfactor, 1, 3)
  call range_integer_error('lurr', lurr, lurrLIM(1), lurrLIM(2))
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
  call range_integer_error('gammax', gammax, gammaxLIM(1), gammaxLIM(2))
  call range_integer_error('strength', strength, strengthLIM(1), strengthLIM(2))
  if (strength /= 11 .and. flagstrengthjp) then
    write(*,'(" TALYS-error: strengthjp must be n for strength /= 11")')
    stop
  endif
  if ((strengthM1 < 1 .or. strengthM1 > 4) .and. strengthM1 /= 8 .and. strengthM1 /= 10 .and. &
 &    strengthM1 /= 11 .and. strengthM1 /= 12) then
    write(*,'(" TALYS-error: strengthM1 = 1, 2, 3, 4, 8, 10, 11 or 12")')
    stop
  endif
  if (strengthM1 /= 11 .and. flagstrengthjp) then
    write(*,'(" TALYS-error: strengthjp must be n for strengthM1 /= 11")')
    stop
  endif
  if (.not. flaglegacy) then
    if (strength <= 7) then
      write(*,'(" TALYS-error: strength = 8, 9, 10, 11, 12 or 13 are recommended. If you want to use legacy models put legacy y")')
      stop
    endif
    if (strengthM1 /= 3 .and. strengthM1 /= 8 .and. strengthM1 /= 10 .and. strengthM1 /= 11 .and. strengthM1 /= 12) then
      write(*,'(" TALYS-error: strengthM1 = 3, 8, 10, 11 or 12 are recommended. If you want to use legacy models put legacy y")')
      stop
    endif
  endif
  do Zix = 0, numZ
    do Nix = 0, numN
      Z = Zinit - Zix
      A = Ainit - Zix - Nix
      do irad = 0, 1
        do lval = 1, gammax
          call range_real_error('etable', etable(Zix, Nix, irad, lval), etableLIM(1), etableLIM(2), &
 &          index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L')
          call range_real_error('ftable', ftable(Zix, Nix, irad, lval), ftableLIM(1), ftableLIM(2), &
 &          index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L')
          call range_real_error('wtable', wtable(Zix, Nix, irad, lval), wtableLIM(1), wtableLIM(2), &
 &          index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L')
          call range_real_error('etableadjust', etableadjust(Zix, Nix, irad, lval), etableadjustLIM(1), etableadjustLIM(2), &
 &          index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L')
          call range_real_error('ftableadjust', ftableadjust(Zix, Nix, irad, lval), ftableadjustLIM(1), ftableadjustLIM(2), &
 &          index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L')
          call range_real_error('wtableadjust', wtableadjust(Zix, Nix, irad, lval), wtableadjustLIM(1), wtableadjustLIM(2), &
 &          index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L')
          do igr = 1, 2
            call range_real_error('energy of GR', egr(Zix, Nix, irad, lval, igr), egrLIM(1), egrLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('width of GR', ggr(Zix, Nix, irad, lval, igr), ggrLIM(1), ggrLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('strength of GR', sgr(Zix, Nix, irad, lval, igr), sgrLIM(1), sgrLIM(2), index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('energy of PR', epr(Zix, Nix, irad, lval, igr), eprLIM(1), eprLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('width of PR', gpr(Zix, Nix, irad, lval, igr), gprLIM(1), gprLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('strength of PR', tpr(Zix, Nix, irad, lval, igr), tprLIM(1), tprLIM(2), index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('egradjust', egradjust(Zix, Nix, irad, lval, igr), egradjustLIM(1), egradjustLIM(2), index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('ggradjust', ggradjust(Zix, Nix, irad, lval, igr), ggradjustLIM(1), ggradjustLIM(2), index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('sgradjust', sgradjust(Zix, Nix, irad, lval, igr), sgradjustLIM(1), sgradjustLIM(2), index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('epradjust', epradjust(Zix, Nix, irad, lval, igr), epradjustLIM(1), epradjustLIM(2), index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('gpradjust', gpradjust(Zix, Nix, irad, lval, igr), gpradjustLIM(1), gpradjustLIM(2), index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
            call range_real_error('spradjust', tpradjust(Zix, Nix, irad, lval, igr), tpradjustLIM(1), tpradjustLIM(2), index1 = Z, name1 = 'Z', &
 &            index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = igr, name5 = 'igr')
          enddo
          call range_real_error('upbendc', upbend(Zix, Nix, irad, lval, 1), upbendcLIM(1), upbendcLIM(2), index1 = Z, name1 = 'Z', &
 &          index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = 1, name5 = 'igr')
          call range_real_error('upbende', upbend(Zix, Nix, irad, lval, 2), upbendeLIM(1), upbendeLIM(2), index1 = Z, name1 = 'Z', &
 &          index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = 2, name5 = 'igr')
          call range_real_error('upbendf', upbend(Zix, Nix, irad, lval, 3), upbendfLIM(1), upbendfLIM(2), index1 = Z, name1 = 'Z', &
 &          index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = 2, name5 = 'igr')
          call range_real_error('upbendcadjust', upbendadjust(Zix, Nix, irad, lval, 1), upbendcadjustLIM(1), upbendcadjustLIM(2), index1 = Z, name1 = 'Z', &
 &          index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = 1, name5 = 'igr')
          call range_real_error('upbendeadjust', upbendadjust(Zix, Nix, irad, lval, 2), upbendeadjustLIM(1), upbendcadjustLIM(2), index1 = Z, name1 = 'Z', &
 &          index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = 2, name5 = 'igr')
          call range_real_error('upbendfadjust', upbendadjust(Zix, Nix, irad, lval, 3), upbendfadjustLIM(1), upbendfadjustLIM(2), index1 = Z, name1 = 'Z', &
 &          index2 = A, name2 = 'A', index3 = irad, name3 = 'irad', index4 = lval, name4 = 'L', index5 = 2, name5 = 'igr')
        enddo
      enddo
      call range_real_error('gamgam', gamgam(Zix, Nix), gamgamLIM(1), gamgamLIM(2), default = 0., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('D0', D0(Zix, Nix), D0LIM(1), D0LIM(2), default = 0., unit = 'keV', index1 = Z, name1 = 'Z', &
 &      index2 = A, name2 = 'A')
      do type = -1, 6
        call range_real_error('fiso', fiso(type), fisoLIM(1), fisoLIM(2), index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', &
 &        index3 = type, name3 = 'type', default = -1.)
        call range_real_error('fisom', fisom(type), fisomLIM(1), fisomLIM(2), index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', &
 &        index3 = type, name3 = 'type', default = -1.)
      enddo
      call range_real_error('gamgamadjust', gamgamadjust(Zix, Nix), gamgamadjustLIM(1), gamgamadjustLIM(2), index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
    enddo
  enddo
  call range_real_error('RprimeU', RprimeU, RprimeULIM(1), RprimeULIM(2))
  if (flagracap .and. k0 == 0) then
    write(*, '(" TALYS-error: Radiative capture model not possible for incident photons")')
    stop
  endif
  call range_integer_error('ldmodelracap', ldmodelracap, ldmodelracapLIM(1), ldmodelracapLIM(2))
  call range_real_error('levinger', levinger, levingerLIM(1), levingerLIM(2))
  do Zix = 0, numZ
    do Nix = 0, numN
      call range_real_error('sfth', spectfacth(Zix, Nix), spectfacthLIM(1), spectfacthLIM(2), index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      do i = 0, numlev
        call range_real_error('sfexp', spectfacexp(Zix, Nix, i), spectfacexpLIM(1), spectfacexpLIM(2), index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', &
 &        index3 = i, name3 = 'level')
      enddo
    enddo
  enddo
!
! 6. Check of values for pre-equilibrium
!
  call range_real_error('epreeq', epreeq, epreeqLIM(1), epreeqLIM(2), default = -1.)
  call range_integer_error('preeqmode', preeqmode, preeqmodeLIM(1), preeqmodeLIM(2))
  call range_integer_error('mpreeqmode', mpreeqmode, mpreeqmodeLIM(1), mpreeqmodeLIM(2))
  call range_integer_error('breakupmodel', breakupmodel, breakupmodelLIM(1), breakupmodelLIM(2))
  call range_integer_error('phmodel', phmodel, phmodelLIM(1), phmodelLIM(2))
  call range_integer_error('pairmodel', pairmodel, pairmodelLIM(1), pairmodelLIM(2))
  call range_integer_error('pespinmodel', pespinmodel, pespinmodelLIM(1), pespinmodelLIM(2))
  call range_real_error('emulpre', emulpre, emulpreLIM(1), emulpreLIM(2))
  call range_real_error('M2constant', M2constant, M2constantLIM(1), M2constantLIM(2))
  call range_real_error('M2limit', M2limit, M2limitLIM(1), M2limitLIM(2))
  call range_real_error('M2shift', M2shift, M2shiftLIM(1), M2shiftLIM(2))
  call range_real_error('Rpipi', Rpipi, RpipiLIM(1), RpipiLIM(2))
  call range_real_error('Rnunu', Rnunu, RnunuLIM(1), RnunuLIM(2))
  call range_real_error('Rpinu', Rpinu, RpinuLIM(1), RpinuLIM(2))
  call range_real_error('Rnupi', Rnupi, RnupiLIM(1), RnupiLIM(2))
  call range_real_error('Rgamma', Rgamma, RgammaLIM(1), RgammaLIM(2))
  call range_real_error('Esurf', Esurf0, Esurf0LIM(1), Esurf0LIM(2), default = -1.)
  call range_integer_error('msdbins', msdbins, msdbinsLIM(1), msdbinsLIM(2), default = 0)
  call range_real_error('E-in', Emsdmin, EmsdminLIM(1), EmsdminLIM(2))
  call range_real_error('elwidth', elwidth, elwidthLIM(1), elwidthLIM(2))
  call range_real_error('xscaptherm', xscaptherm(-1), xscapthermLIM(1), xscapthermLIM(2), default = 0.)
  call range_real_error('xsptherm', xsptherm(-1), xspthermLIM(1), xspthermLIM(2), default = 0.)
  call range_real_error('xsalphatherm', xsalphatherm(-1), xsalphathermLIM(1), xsalphathermLIM(2), default = 0.)
  if (k0 == 0 .and. flagpecomp) then
    write(*, '(" TALYS-error: No pick-up and knock-out mechanism for photonuclear reactions")')
    stop
  endif
  do type = 0, 6
    call range_real_error('Cstrip', Cstrip(type), CstripLIM(1), CstripLIM(2), index1 = type, name1 = 'type')
    call range_real_error('Cknock', Cknock(type), CknockLIM(1), CknockLIM(2), index1 = type, name1 = 'type')
    call range_real_error('Cbreak', Cbreak(type), CbreakLIM(1), CbreakLIM(2), index1 = type, name1 = 'type')
  enddo
  call range_real_error('GMRadjustE', GMRadjustE, GMRadjustELIM(1), GMRadjustELIM(2))
  call range_real_error('GQRadjustE', GQRadjustE, GQRadjustELIM(1), GQRadjustELIM(2))
  call range_real_error('LEORadjustE', LEORadjustE, LEORadjustELIM(1), LEORadjustELIM(2))
  call range_real_error('HEORadjustE', HEORadjustE, HEORadjustELIM(1), HEORadjustELIM(2))
  call range_real_error('GMRadjustG', GMRadjustG, GMRadjustGLIM(1), GMRadjustGLIM(2))
  call range_real_error('GQRadjustG', GQRadjustG, GQRadjustGLIM(1), GQRadjustGLIM(2))
  call range_real_error('LEORadjustG', LEORadjustG, LEORadjustGLIM(1), LEORadjustGLIM(2))
  call range_real_error('HEORadjustG', HEORadjustG, HEORadjustGLIM(1), HEORadjustGLIM(2))
  call range_real_error('GMRadjustD', GMRadjustD, GMRadjustDLIM(1), GMRadjustDLIM(2))
  call range_real_error('GQRadjustD', GQRadjustD, GQRadjustDLIM(1), GQRadjustDLIM(2))
  call range_real_error('LEORadjustD', LEORadjustD, LEORadjustDLIM(1), LEORadjustDLIM(2))
  call range_real_error('HEORadjustD', HEORadjustD, HEORadjustDLIM(1), HEORadjustDLIM(2))
!
! 7. Check of values for level densities
!
  call range_integer_error('spincutmodel', spincutmodel, spincutmodelLIM(1), spincutmodelLIM(2))
  call range_integer_error('shellmodel', shellmodel, shellmodelLIM(1), shellmodelLIM(2))
  call range_integer_error('kvibmodel', kvibmodel, kvibmodelLIM(1), kvibmodelLIM(2))
  call range_integer_error('ldmodelcn', ldmodelCN, ldmodelCNLIM(1), ldmodelCNLIM(2))
  if (.not. flaglegacy) then
    if (ldmodelCN == 3 .or. ldmodelCN == 4 .or. ldmodelCN == 6) then
      write(*,'(" TALYS-error: ldmodelCN = 1, 2, 5, or 7 are recommended. If you want to use legacy models put legacy y")')
      stop
    endif
  endif
  do Zix = 0, numZ
    do Nix = 0, numN
      Z = Zinit - Zix
      A = Ainit - Zix - Nix
      call range_integer_error('ldmodel', ldmodel(Zix, Nix), 1, 7, index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      if (.not. flaglegacy) then
        if (ldmodel(Zix, Nix) == 3 .or. ldmodel(Zix, Nix) == 4 .or. ldmodel(Zix, Nix) == 6) then
          write(*,'(" TALYS-error: ldmodel = 1, 2, 5, or 7 are recommended. If you want to use legacy models put legacy y")')
          stop
        endif
      endif
      if (ldmodel(Zix, Nix) >= 4) then
        inquire (file = trim(path) // 'density/ground/hilaire/Sn.tab', exist = lexist)
        if ( .not. lexist) then
          write(*, '(" TALYS-error: Microscopic HFB tables are not installed: download the full TALYS package from www.talys.eu")')
          stop
        endif
      endif
      call range_real_error('a', alev(Zix, Nix), alevLIM(1), alevLIM(2), default = 0.,  index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('alimit', alimit(Zix, Nix), alimitLIM(1), alimitLIM(2), default = 0.,  index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('gammald', gammald(Zix, Nix), gammaldLIM(1), gammaldLIM(2), default = -1.,  index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('risomer', Risomer(Zix, Nix), RisomerLIM(1), RisomerLIM(2), index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      do ibar = 0, numbar
        call range_real_error('deltaW', deltaW(Zix, Nix, ibar), deltaWLIM(1), deltaWLIM(2), default = 0.,  index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_integer_error('Nlow', Nlow(Zix, Nix, ibar), NlowLIM(1), NlowLIM(2), default = -1,  index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_integer_error('Ntop', Ntop(Zix, Nix, ibar), NtopLIM(1), NtopLIM(2), default = -1,  index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_integer_error('Ntop', Ntop(Zix, Nix, ibar), Nlow(Zix, Nix, ibar), 200, default = -1,  index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('E0', E0(Zix, Nix, ibar), E0LIM(1), E0LIM(2), default = 0.,  index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('beta2', beta2(Zix, Nix, ibar), beta2LIM(1), beta2LIM(2), index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('s2adjust', s2adjust(Zix, Nix, ibar), s2adjustLIM(1), s2adjustLIM(2), index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('Krotconstant', Krotconstant(Zix, Nix, ibar), KrotconstantLIM(1), KrotconstantLIM(2), index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('Ufermi', Ufermi(Zix, Nix, ibar), UfermiLIM(1), UfermiLIM(2), index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('cfermi', cfermi(Zix, Nix, ibar), cfermiLIM(1), cfermiLIM(2), index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('T', T(Zix, Nix, ibar), TLIM(1), TLIM(2), default = 0.,  index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('Exmatch', Exmatch(Zix, Nix, ibar), ExmatchLIM(1), ExmatchLIM(2), default = 0.,  index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('Tadjust', Tadjust(Zix, Nix, ibar), TadjustLIM(1), TadjustLIM(2), default = 0.,  index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('E0adjust', E0adjust(Zix, Nix, ibar), E0adjustLIM(1), E0adjustLIM(2), default = 0.,  index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('Exmatchadjust', Exmatchadjust(Zix, Nix, ibar), ExmatchadjustLIM(1), ExmatchadjustLIM(2), index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('Pshift', Pshift(Zix, Nix, ibar), PshiftLIM(1), PshiftLIM(2), index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('Pshiftadjust', Pshiftadjust(Zix, Nix, ibar), PshiftadjustLIM(1), PshiftadjustLIM(2), index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('ctable', ctable(Zix, Nix, ibar), ctableLIM(1), ctableLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('ptable', ptable(Zix, Nix, ibar), ptableLIM(1), ptableLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('ctableadjust', ctableadjust(Zix, Nix, ibar), ctableadjustLIM(1), ctableadjustLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('ptableadjust', ptableadjust(Zix, Nix, ibar), ptableadjustLIM(1), ptableadjustLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
      enddo
      call range_real_error('aadjust', aadjust(Zix, Nix), aadjustLIM(1), aadjustLIM(2), index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('gadjust', gadjust(Zix, Nix), gadjustLIM(1), gadjustLIM(2), index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('gnadjust', gnadjust(Zix, Nix), gnadjustLIM(1), gnadjustLIM(2), index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('gpadjust', gpadjust(Zix, Nix), gpadjustLIM(1), gpadjustLIM(2), index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('pair', pair(Zix, Nix), pairLIM(1), pairLIM(2), index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('g', g(Zix, Nix), gLIM(1), gLIM(2), default = 0., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('gn', gn(Zix, Nix), gnLIM(1), gnLIM(2), default = 0., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('gp', gp(Zix, Nix), gpLIM(1), gpLIM(2), default = 0., index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('alphald', alphald(Zix, Nix), alphaldLIM(1), alphaldLIM(2), index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('betald', betald(Zix, Nix), betaldLIM(1), betaldLIM(2), index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      if (betald(Zix, Nix) < 0..and. abs(betald(Zix, Nix)) > alphald(Zix, Nix)) then
        write(*, '(" TALYS-error: if betald<0, |betald|<alphald")')
        stop
      endif
      call range_real_error('gammashell1', gammashell1(Zix, Nix), gammashell1LIM(1), gammashell1LIM(2), index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
      call range_real_error('Pshiftconstant', Pshiftconstant(Zix, Nix), PshiftconstantLIM(1), PshiftconstantLIM(2), index1 = Z, name1 = 'Z', index2 = A, name2 = 'A')
    enddo
  enddo
  call range_real_error('cglobal', cglobal, cglobalLIM(1), cglobalLIM(2), default = 0.)
  call range_real_error('pglobal', pglobal, pglobalLIM(1), pglobalLIM(2), default = 0.)
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
  call range_real_error('gammashell2', gammashell2, gammashell2LIM(1), gammashell2LIM(2))
  call range_real_error('pairconstant', pairconstant, pairconstantLIM(1), pairconstantLIM(2))
  call range_real_error('Kph', Kph, KphLIM(1), KphLIM(2))
  call range_real_error('Rspincut', Rspincut, RspincutLIM(1), RspincutLIM(2))
  call range_real_error('Rspincutpreeq', Rspincutpreeq, RspincutpreeqLIM(1), RspincutpreeqLIM(2))
  call range_real_error('Rspincutff', Rspincutff, RspincutffLIM(1), RspincutffLIM(2))
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
  call range_integer_error('fismodel', fismodel, fismodelLIM(1), fismodelLIM(2))
  call range_integer_error('fismodelalt', fismodelalt, fismodelaltLIM(1), fismodelaltLIM(2))
  if (fismodel < 5 .and. flagfispartdamp) then
    write(*,'(" TALYS-error: Fission partial damping only allowed for fismodel 5 or 6")')
    stop
  endif
  call range_integer_error('fymodel', fymodel, fymodelLIM(1), fymodelLIM(2))
  call range_integer_error('ffmodel', ffmodel, ffmodelLIM(1), ffmodelLIM(2))
  call range_integer_error('pfnsmodel', pfnsmodel, pfnsmodelLIM(1), pfnsmodelLIM(2))
  call range_integer_error('gefran', gefran, gefranLIM(1), gefranLIM(2))
  call range_real_error('Cnubar1', Cnubar1, Cnubar1LIM(1), Cnubar1LIM(2))
  call range_real_error('Cnubar2', Cnubar2, Cnubar2LIM(1), Cnubar2LIM(2))
  call range_real_error('Tmadjust', Tmadjust, TmadjustLIM(1), TmadjustLIM(2))
  call range_real_error('Fsadjust', Fsadjust, FsadjustLIM(1), FsadjustLIM(2))
  call range_real_error('Cbarrier', Cbarrier, CbarrierLIM(1), CbarrierLIM(2))
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
        call range_integer_error('type of axiality', axtype(Zix, Nix, ibar), axtypeLIM(1), axtypeLIM(2), default = 0, index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('fission barrier', fbarrier(Zix, Nix, ibar), fbarrierLIM(1), fbarrierLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('fisbaradjust', fbaradjust(Zix, Nix, ibar), fbaradjustLIM(1), fbaradjustLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('fission width', fwidth(Zix, Nix, ibar), fwidthLIM(1), fwidthLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('fishwadjust', fwidthadjust(Zix, Nix, ibar), fwidthadjustLIM(1), fwidthadjustLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('bdamp', bdamp(Zix, Nix, ibar), bdampLIM(1), bdampLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('bdampadjust', bdampadjust(Zix, Nix, ibar), bdampadjustLIM(1), bdampadjustLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('Rtransmom', Rtransmom(Zix, Nix, ibar), RtransmomLIM(1), RtransmomLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('Rclass2mom', Rclass2mom(Zix, Nix, ibar), Rclass2momLIM(1), Rclass2momLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
        call range_real_error('class2width', widthc2(Zix, Nix, ibar), widthc2LIM(1), widthc2LIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &        index2 = A, name2 = 'A', index3 = ibar, name3 = 'barrier')
      enddo
      call range_real_error('betafiscor', betafiscor(Zix, Nix), betafiscorLIM(1), betafiscorLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &      index2 = A, name2 = 'A')
      call range_real_error('rmiufiscor', rmiufiscor(Zix, Nix), rmiufiscorLIM(1), rmiufiscorLIM(2), default = -1., index1 = Z, name1 = 'Z', &
 &      index2 = A, name2 = 'A')
      call range_real_error('vfiscor', vfiscor(Zix, Nix), vfiscorLIM(1), vfiscorLIM(2), default = -1., index1 = Z, name1 = 'Z', &
 &      index2 = A, name2 = 'A')
      call range_real_error('betafiscoradjust', betafiscoradjust(Zix, Nix), betafiscoradjustLIM(1), betafiscoradjustLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &      index2 = A, name2 = 'A')
      call range_real_error('vfiscoradjust', vfiscoradjust(Zix, Nix), vfiscoradjustLIM(1), vfiscoradjustLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &      index2 = A, name2 = 'A')
      call range_real_error('rmiufiscoradjust', rmiufiscoradjust(Zix, Nix), rmiufiscoradjustLIM(1), rmiufiscoradjustLIM(2), default = 0., index1 = Z, name1 = 'Z', &
 &      index2 = A, name2 = 'A')
      if (flagsffactor) then
        if (rmiufiscor(Zix, Nix) /= -1. .and. vfiscor(Zix, Nix) /= -1.) then
          write(*,'(" TALYS-error: if sffactor y one can not give both vfiscor and rmiufiscor in the input")')
          stop
        endif
      else
        if (vfiscor(Zix, Nix) == -1.) vfiscor(Zix, Nix) = 1.
        if (rmiufiscor(Zix, Nix) == -1.) rmiufiscor(Zix, Nix) = 1.
      endif
    enddo
  enddo
!
! 9. Check of values for output
!
  call range_real_error('eadd', eadd, eaddLIM(1), eaddLIM(2), unit = 'MeV')
  call range_real_error('eaddel', eaddel, eaddelLIM(1), eaddelLIM(2), unit = 'MeV')
  call range_integer_error('ddxmode', ddxmode, ddxmodeLIM(1), ddxmodeLIM(2))
  do type = 0, 6
    do i = 1, ddxecount(type)
      call range_real_error('fileddxe', fileddxe(type, i), fileddxeLIM(1), fileddxeLIM(2), default = 0., index1 = type, name1 = 'type', &
 &      index2 = i, name2 = 'i')
    enddo
    do i = 1, ddxacount(type)
      call range_real_error('fileddxa', fileddxa(type, i), fileddxaLIM(1), fileddxaLIM(2), default = 0., index1 = type, name1 = 'type', &
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
    call range_real_error(adjustkey(n), D, DLIM(1), DLIM(2), index1 = n, name1 = 'n')
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
    if (trim(reslib) == 'tendl.2025') return
    if (trim(reslib) == 'jeff4.0') return
    if (trim(reslib) == 'endfb8.1') return
    if (trim(reslib) == 'cendl3.2') return
    if (trim(reslib) == 'jendl5.0') return
    write(*, '(" TALYS-error: Wrong library name: ", a)') trim(reslib)
    stop
  endif
  return
end subroutine checkvalue
! Copyright A.J. Koning 2021
