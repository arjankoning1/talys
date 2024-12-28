subroutine inputout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write input parameters
!
! Author    : Arjan Koning
!
! 2024-10-25: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Variables for numerics
!   maxchannel      ! maximal number of outgoing particles in individual channel description
!   maxenrec        ! number of recoil energies
!   maxN            ! maximal number of neutrons away from initial compound nucleus
!   maxNrp          ! maximal number of neutrons away from the initial compound nucleus
!   maxZ            ! maximal number of protons away from initial compound nucleus
!   maxZrp          ! maximal number of protons away from the initial compound nucleus
!   nangle          ! number of angles
!   nanglecont      ! number of angles for continuum
!   nbins0          ! number of continuum excitation energy bins
!   nanglerec       ! number of recoil angles
!   popeps          ! limit for population cross sections
!   segment         ! help array for storing segment intersection points
!   transeps        ! absolute limit for transmission coefficient
!   transpower      ! power for transmission coefficient limit
!   xseps           ! limit for cross sections
! Variables for basic reaction
!   flagang         ! flag for output of angular distributions
!   flagastro       ! flag for calculation of astrophysics reaction rate
!   flagbasic       ! flag for output of basic information and results
!   flagoutall      ! flag for output of all data in main output file
!   flagchannels    ! flag for exclusive channels calculation
!   flagendf        ! flag for information for ENDF - 6 file
!   flagendfdet     ! flag for detailed ENDF - 6 information per channel
!   flagendfecis    ! flag for new ECIS calculation for ENDF - 6 files
!   flagEchannel    ! flag for channel energy for emission spectrum
!   flaglabddx      ! flag for calculation of DDX in LAB system
!   flagmassdis     ! flag for calculation of fission fragment mass yields
!   flagmicro       ! flag for completely microscopic Talys calculation
!   flagfit         ! flag for using fitted nuclear model parameters
!   flagngfit       ! flag for using fitted (n,g) nuclear model parameters
!   flagnffit       ! flag for using fitted (n,f) nuclear model parameters
!   flagnnfit       ! flag for using fitted (n,n'), (n,2n) and (n,p) nuclear model parameters
!   flagnafit       ! flag for using fitted (n,a) nuclear model parameters
!   flagndfit       ! flag for using fitted (n,d) nuclear model parameters
!   flagpnfit       ! flag for using fitted (p,n) nuclear model parameters
!   flaggnfit       ! flag for using fitted (g,n) nuclear model parameters
!   flagdnfit       ! flag for using fitted (d,n) nuclear model parameters
!   flagpartable    ! flag for output of model parameters on separate file
!   flagreaction    ! flag for calculation of nuclear reactions
!   flagrecoil      ! flag for calculation of recoils
!   flagrecoilav    ! flag for average velocity in recoil calculation
!   flagrel         ! flag for relativistic kinematics
!   flagrpevap      ! flag for evaporation of residual products at high incident energies
! Variables for best files
!   flagbest        ! flag to use best set of adjusted parameters
!   flagbestend     ! flag to put best set of parameters at end of input file
! Variables for basic parameters
!   eninclow        ! minimal incident energy for nuclear model calculations
!   flagequi        ! flag to use equidistant excitation bins instead of logarithmic
!   flagequispec    ! flag to use equidistant bins for emission spectra
!   isomer          ! definition of isomer in seconds
!   outtype         ! type of outgoing particles
! Variables for input energies
!   energyfile      ! file with energies for OMP calculation
!   eninc           ! incident energy in MeV
!   enincmax        ! maximum incident energy
!   enincmin        ! minimum incident energy
!   Estop           ! incident energy above which TALYS stops
!   flaginitpop     ! flag for initial population distribution
!   flagpopMeV      ! flag to use initial population per MeV instead of histogram
!   Ninc            ! number of incident energies
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   Ltarget         ! excited level of target
!   ptype0          ! type of incident particle
!   Starget         ! symbol of target nucleus
! Variables for output
!   ddxmode         ! mode for DDX: 0: None, 1: Angular distributions, 2: Spectra per angle, 3: Both
!   flagcompo       ! flag for output of cross section components
!   flagsacs        ! flag for statistical analysis of cross sections
!   flagbinspec     ! flag for output of emission spectrum per excitation bin
!   flagblock       ! flag to block spectra, angle and gamma files
!   flagcheck       ! flag for output of numerical checks
!   flagdecay       ! flag for output of decay of each population bin
!   flagexc         ! flag for output of excitation functions
!   flaggamdis      ! flag for output of discrete gamma - ray intensities
!   flaginverse     ! flag for output of transmission coefficients and inverse cross sections
!   flagmain        ! flag for main output
!   flagpop         ! flag for output of population
!   flagspec        ! flag for output of spectra
! Variables for compound reactions
!   ewfc            ! off - set incident energy for width fluctuation
!   eurr            ! off - set incident energy for URR calculation
!   flagcomp        ! flag for compound angular distribution calculation
!   flageciscomp    ! flag for compound nucleus calculation by ECIS
!   flaggroup       ! flag for output of low energy groupwise cross sections
!   flagfullhf      ! flag for full spin dependence of transmission coefficie
!   flagres         ! flag for output of low energy resonance cross sections
!   flagurr         ! flag for output of unresolved resonance parameters
!   flagurrnjoy     ! normalization of URR parameters with NJOY method
!   lurr            ! maximal orbital angular momentum for URR
!   wmode           ! designator for width fluctuation model
! Variables for direct reactions
!   core            ! even - even core for weakcoupling ( - 1 or 1)
!   eadd            ! on - set incident energy for addition of discrete states
!   eaddel          ! on - set incident energy for addition of elastic peak
!   flagautorot     ! flag for automatic rotational coupled channels
!   flagcoulomb     ! flag for Coulomb excitation calculation with ECIS
!   flagcpang       ! flag for compound angular distribution calculation for
!   flagdirect      ! flag for output of direct reaction results
!   flagdisc        ! flag for output of discrete state cross sections
!   flageciscalc    ! flag for new ECIS calculation for outgoing particles
!   flagecissave    ! flag for saving ECIS input and output files
!   flaggiant0      ! flag for collective contribution from giant resonances
!   flaginccalc     ! flag for new ECIS calculation for incident channel
!   flaglegendre    ! flag for output of Legendre coefficients
!   flagoutecis     ! flag for output of ECIS results
!   flagrot         ! flag for use of rotational optical model per outgoing p
!   flagspher       ! flag to force spherical optical model
!   flagstate       ! flag for optical model potential for each excited state
!   flagsys         ! flag for reaction cross section from systematics
!   flagtransen     ! flag for output of transmission coefficients per energy
!   maxband         ! highest vibrational band added to rotational model
!   maxrot          ! number of included excited rotational levels
! Variables for preequilibrium
!   breakupmodel    ! model for break - up reaction: 1. Kalbach 2. Avrigeanu
!   flagecisdwba    ! flag for new ECIS calculation for DWBA for MSD
!   emulpre         ! on - set incident energy for multiple preequilibrium
!   epreeq          ! on - set incident energy for preequilibrium calculation
!   flag2comp       ! flag for two - component pre - equilibrium model
!   flagpeout       ! flag for output of pre - equilibrium results
!   flaggshell      ! flag for energy dependence of single particle level den
!   flagonestep     ! flag for continuum one - step direct only
!   flagoutdwba     ! flag for output of DWBA cross sections for MSD
!   flagpecomp      ! flag for Kalbach complex particle emission model
!   flagsurface     ! flag for surface effects in exciton model
!   mpreeqmode      ! designator for multiple pre - equilibrium model
!   pairmodel       ! model for preequilibrium pairing energy
!   pespinmodel     ! model for pre - equilibrium or compound spin distribution
!   phmodel         ! particle - hole state density model
!   preeqmode       ! designator for pre - equilibrium model
! Variables for fission
!   Rfiseps         ! ratio for limit for fission cross section per nucleus
!   fismodel        ! fission model alternative fission model for default barriers
!   fismodelalt     ! alternative fission model for default barriers
!   flagclass2      ! flag for class2 states in fission
!   flagffevap      ! flag for calculation of particle evaporation from fissi
!   flagffspin      ! flag to use spin distribution in initial FF population
!   flagfisout      ! flag for output of fission information
!   flagfispartdamp   ! flag for fission partial damping
!   flagfisfeed     ! flag for output of fission per excitation bin
!   flagfission     ! flag for fission
!   flaghbstate     ! flag for head band states in fission
!   flagoutfy       ! flag for output detailed fission yield calculation
!   fymodel         ! fission yield model, 1: Brosa 2: GEF
!   ffmodel         ! fission fragment model, 1: GEF 2: HF3D (Okumura) 3: SPY 4: Langevin-4D
!   pfnsmodel       ! PFNS  model, 1: Iwamoto 2: from FF decay
!   gefran          ! number of random events for GEF calculation
! Variables for astrophysics
!   flagastroex     ! flag for calculation of astrophysics reaction rate to f
!   flagastrogs     ! flag for calculation of astrophysics reaction rate with
!   nonthermlev     ! non - thermalized level in the calculation of astrophysic
! Variables for medical isotope production
!   Area            ! target area in cm^2
!   Eback           ! lower end of energy range in MeV for isotope
!   Ebeam           ! incident energy in MeV for isotope production
!   flagprod        ! flag for isotope production
!   Ibeam           ! beam current in mA for isotope production
!   radiounit       ! unit for radioactivity: Bq, kBq, MBq, Gbq, mCi, Ci or kCi
!   rhotarget       ! target material density
!   Tcool           ! cooling time per unit cooling time unit (y, d, h, m, s)
!   Tirrad          ! irradiation time per unit irradiation time unit (y, d, h, m, s)
!   unitTcool       ! cooling time unit (y, d, h, m, s)
!   unitTirrad      ! irradiation time unit (y, d, h, m, s)
!   yieldunit       ! unit for isotope yield: num (number), mug, mg, g, or kg
! Variables for gamma rays
!   flaggamma       ! flag for output of gamma - ray information
!   flagracap       ! flag for radiative capture model
!   flagupbend      ! flag for low-energy upbend of photon strength function
!   flagpsfglobal   ! flag for global photon strength functions only
!   flaggnorm       ! flag to normalize PSF to average radiative width
!   gammax          ! number of l - values for gamma multipolarity
!   ldmodelracap    ! level density model for direct radiative capture
!   strength        ! E1 strength function model
!   strengthM1      ! model for M1 gamma - ray strength function
! Variables for discrete levels
!   disctable       ! table with discrete levels
!   flagbestbr      ! flag to use only best set of branching ratios
!   flagpseudores   ! flag for using light nuclide discrete levels for resonances
!   flagelectron    ! flag for application of electron conversion coefficient
!   flaglevels      ! flag for output of discrete level information
!   nlevbin         ! number of excited levels for binary nucleus
!   nlevmax         ! maximum number of included discrete levels for target
!   nlevmaxres      ! maximum number of included discrete levels for residual nuclides
! Variables for level density
!   flagasys        ! flag for all level density parameters a from systematic
!   flagcolall      ! flag for collective enhancement of level density
!   flagcolldamp    ! flag for damping of coll. effects in eff. level density (without explicit coll. enh.)
!   flagctmglob     ! flag for global CTM model (no discrete level info)
!   flagdensity     ! flag for output of level densities
!   flagparity      ! flag for non - equal parity distribution
!   kvibmodel       ! model for vibrational enhancement
!   ldmodelCN       ! level density model for compound nucleus
!   ldmodelall      ! level density model for all nuclides
!   shellmodel      ! model for shell correction energies
!   spincutmodel    ! model for spin cutoff factor for ground state
! Variables for masses
!   flagexpmass     ! flag for using experimental nuclear mass if available
!   massmodel       ! model for theoretical nuclear mass
! Variables for OMP
!   flagriplomp     ! flag for RIPL OMP
! Variables for OMP
!   alphaomp        ! alpha optical model
!   deuteronomp     ! deuteron optical model
!   flagdisp        ! flag for dispersive optical model
!   flagincadj      ! flag for OMP adjustment on incident channel also
!   flagjlm         ! flag for using semi - microscopic JLM OMP
!   flagoutkd       ! flag for output of KD03 OMP parameters
!   flaglocalomp    ! flag for local (y) or global (n) optical model
!   flagompall      ! flag for new optical model calculation for all residual
!   flagomponly     ! flag to execute ONLY an optical model calculation
!   flagoutomp      ! flag for output of optical model parameters
!   flagriplrisk    ! flag for RIPL OMP outside mass validity range
!   flagsoukho      ! flag for Soukhovitskii OMP for actinides
!   jlmmode         ! option for JLM imaginary potential normalization
!   radialmodel     ! model for radial matter densities (JLM OMP only)
! Variables for energies
!   flagadd         ! flag for addition of discrete states to spectra flag
!   flagaddel       ! flag for addition of elastic peak to spectra
!   flagmulpre      ! flag for multiple pre - equilibrium calculation
!   flagpreeq       ! flag for pre - equilibrium calculation
!   flagwidth       ! flag for width fluctuation calculation
! Constants
!   parname         ! name of particle
!   parsym          ! symbol of particle
! Variables for reading input lines
!   inline    ! input line
!   nlines    ! number of input lines
!
! *** Declaration of local data
!
  implicit none
  character(len=1)  :: yesno        ! function to assign y or n to logical value
  character(len=12) :: rotstring    ! help variable
  character(len=12) :: sysstring    ! help variable
  integer           :: i            ! counter
  integer           :: type         ! particle type
!
! ************************** User input file ***************************
!
  write(*, '(/" ########## USER INPUT ##########")')
  write(*, '(/" USER INPUT FILE"/)')
  do i = 1, nlines
    write(*, '(1x, a)') trim(inline(i))
  enddo
!
! ********* All possible input parameters including defaults ***********
!
  write(*, '(/" USER INPUT FILE + DEFAULTS"/)')
  write(*, '(" Keyword           Value   Variable     Explanation"/)')
!
! 1. Four main keywords
!
  write(*, '(" #"/" # Four main keywords"/" #")')
  write(*, '(" projectile          ", a1, "     ptype0       type of incident particle")') ptype0
  write(*, '(" element            ", a2, "     Starget      symbol of target nucleus")') Starget
  write(*, '(" mass              ", i3, "     mass         mass number of target nucleus")') Atarget
  if (Ninc == 1 .and. .not. flaginitpop) then
    write(*, '(" energy           ", f8.3, " eninc        incident energy in MeV")') eninc(1)
  else
    write(*, '(" energy ", a14, "     energyfile   file with incident energies")') energyfile
  endif
!
! 2. Basic physical and numerical parameters
!
  write(*, '(" #"/" # Basic physical and numerical parameters")')
  write(*, '(" #")')
  write(*, '(" ejectiles", 7(1x, a1), "   outtype      outgoing particles")') (outtype(type), type = 0, 6)
  write(*, '(" maxz              ", i3, "     maxZ         maximal number of protons from the initial compound nucleus")') maxZ
  write(*, '(" maxn              ", i3, "     maxN         maximal number of neutrons from the initial compound nucleus")') &
 &  maxN
  write(*, '(" bins              ", i3, "     nbins        number of continuum excitation energy bins")') nbins0
  write(*, '(" equidistant         ", a1, "     flagequi     flag to use equidistant excitation bins instead of logarithmic ", &
 &  "bins")') yesno(flagequi)
  write(*, '(" equispec            ", a1, "     flagequispec flag to use equidistant bins for emission spectra")')  &
 &  yesno(flagequispec)
  write(*, '(" popmev              ", a1, "     flagpopmev   flag to use initial population per MeV instead of histograms")') &
 &  yesno(flagpopmev)
  write(*, '(" segment           ", i3, "     segment      number of segments to divide emission energy grid")') segment
  write(*, '(" maxlevelstar      ", i3, "     nlevmax      maximum number of included discrete levels for target")')  nlevmax
  write(*, '(" maxlevelsres      ", i3, "     nlevmaxres   maximum number of included discrete levels for residual nucleus")') &
 &  nlevmaxres
  do type = 0, 6
    write(*, '(" maxlevelsbin ", a1, "    ", i3, "     nlevbin      maximum number of included discrete levels for ", &
 &    a8, " channel")') parsym(type), nlevbin(type), parname(type)
  enddo
  write(*, '(" ltarget           ", i3, "     ltarget      excited level of target")') Ltarget
  write(*, '(" isomer          ", es9.2, " isomer       definition of isomer in seconds")') isomer
  write(*, '(" transpower        ", i3, "     transpower   power for transmission coefficient limit")') transpower
  write(*, '(" transeps        ", es9.2, " transeps     limit for transmission coefficient")') transeps
  write(*, '(" xseps           ", es9.2, " xseps        limit for cross sections")') xseps
  write(*, '(" popeps          ", es9.2, " popeps       limit for population cross section per nucleus")') popeps
  write(*, '(" Rfiseps         ", es9.2, " Rfiseps      ratio for limit for fission cross section per nucleus")') Rfiseps
  write(*, '(" elow            ", es9.2, " elow         minimal incident energy for nuclear model calculations")') eninclow
  write(*, '(" angles            ", i3, "     nangle       number of angles")') nangle
  write(*, '(" anglescont        ", i3, "     nanglecont   number of angles for continuum")') nanglecont
  write(*, '(" anglesrec         ", i3, "     nanglerec    number of recoil angles")') nanglerec
  write(*, '(" maxenrec          ", i3, "     maxenrec     number of recoil energies")') maxenrec
  write(*, '(" channels            ", a1, "     flagchannels flag for exclusive channels calculation")') yesno(flagchannels)
  write(*, '(" maxchannel         ", i2, "     maxchannel   maximal number of outgoing particles in", &
 &  " individual channel description")') maxchannel
  write(*, '(" micro               ", a1, "     flagmicro    flag for completely microscopic Talys calculation")') yesno(flagmicro)
  write(*, '(" best                ", a1, "     flagbest     flag to use best set of adjusted parameters")') yesno(flagbest)
  write(*, '(" bestbranch          ", a1, "     flagbestbr   flag to use flag to use only best set of branching ratios")') &
 &  yesno(flagbestbr)
  write(*, '(" bestend             ", a1, "     flagbestend  flag to put best set of parameters at end of input file")') &
 &  yesno(flagbestend)
  write(*, '(" relativistic        ", a1, "     flagrel      flag for relativistic kinematics")') yesno(flagrel)
  write(*, '(" recoil              ", a1, "     flagrecoil   flag for calculation of recoils")') yesno(flagrecoil)
  write(*, '(" labddx              ", a1, "     flaglabddx   flag for calculation of DDX in LAB system")') yesno(flaglabddx)
  write(*, '(" recoilaverage       ", a1, "     flagrecoilav flag for average velocity in recoil calculation")') yesno(flagrecoilav)
  write(*, '(" channelenergy       ", a1, "     flagEchannel flag for channel energy for emission spectrum")') yesno(flagEchannel)
  write(*, '(" reaction            ", a1, "     flagreaction flag for calculation of nuclear reactions")') yesno(flagreaction)
  write(*, '(" fit                 ", a1, "     flagfit      flag to use automatically fitted parameters")') &
 &  yesno(flagfit)
  write(*, '(" ngfit               ", a1, "     flagngfit    flag for using fitted (n,g) nuclear model parameters")') &
 &  yesno(flagngfit)
  write(*, '(" nffit               ", a1, "     flagnffit    flag for using fitted (n,f) nuclear model parameters")') &
 &  yesno(flagnffit)
  write(*, '(" nnfit               ", a1, "     flagnnfit    flag for using fitted (n,n), (n,2n) and (n,p) nuclear model", &
 &  " parameters")') yesno(flagnnfit)
  write(*, '(" nafit               ", a1, "     flagnafit    flag for using fitted (n,a) nuclear model parameters")') &
 &  yesno(flagnafit)
  write(*, '(" ndfit               ", a1, "     flagndfit    flag for using fitted (n,d) nuclear model parameters")') &
 &  yesno(flagndfit)
  write(*, '(" pnfit               ", a1, "     flagpnfit    flag for using fitted (p,n) nuclear model parameters")') &
 &  yesno(flagpnfit)
  write(*, '(" dnfit               ", a1, "     flagdnfit    flag for using fitted (d,n) nuclear model parameters")') &
 &  yesno(flagdnfit)
  write(*, '(" gnfit               ", a1, "     flaggnfit    flag for using fitted (g,n) nuclear model parameters")') &
 &  yesno(flaggnfit)
  write(*, '(" anfit               ", a1, "     flaganfit    flag for using fitted (a,n) nuclear model parameters")') &
 &  yesno(flaganfit)
  write(*, '(" gamgamfit           ", a1, "    flaggamgamfit flag for using fitted Gamma_gamma nuclear model parameters")') &
 &  yesno(flaggamgamfit)
  write(*, '(" macsfit             ", a1, "     flagmacsfit  flag for using fitted MACS nuclear model parameters")') &
 &  yesno(flagmacsfit)
  write(*, '(" astro               ", a1, "     flagastro    flag for calculation of astrophysics reaction rate")') yesno(flagastro)
  write(*, '(" astrogs             ", a1, "     flagastrogs  flag for", &
 &  " calculation of astrophysics reaction rate with target in ground state only")') yesno(flagastrogs)
  write(*, '(" astroex             ", a1, "     flagastroex  flag for", &
 &  " calculation of astrophysics reaction rate to long-lived excited states")') yesno(flagastroex)
  write(*, '(" nonthermlev       ", i3, "     nonthermlev  excited level non-thermalized in the calculation", &
 &  " of astrophysics rate")') nonthermlev
  write(*, '(" massmodel          ", i2, "     massmodel    model for theoretical nuclear mass")') massmodel
  write(*, '(" expmass             ", a1, "     flagexpmass  flag for using experimental nuclear mass if available")') &
 &  yesno(flagexpmass)
  write(*, '(" disctable          ", i2, "     disctable    table with discrete levels")') disctable
  write(*, '(" production          ", a1, "     flagprod     flag for isotope production")') yesno(flagprod)
  write(*, '(" outfy               ", a1, "     flagoutfy    flag for output detailed fission yield calculation")') &
 &  yesno(flagoutfy)
  write(*, '(" gefran        ", i7, "     gefran       number of random events for GEF calculation")') gefran
  write(*, '(" Estop            ", f8.3, " Estop        incident energy above which TALYS stops")') Estop
  write(*, '(" rpevap              ", a1, "     flagrpevap   flag for evaporation of residual products at high", &
 &  " incident energies")') yesno(flagrpevap)
  write(*, '(" maxZrp            ", i3, "     maxZrp       maximal number of protons from the initial", &
 &  " compound nucleus before residual evaporation")') maxZrp
  write(*, '(" maxNrp            ", i3, "     maxNrp       maximal number of neutrons from the initial", &
 &  " compound nucleus before residual evaporation")') maxNrp
  write(*, '(" user      ", a, t28, "user         user for this calculation")') trim(user)
  write(*, '(" source    ", a, t28, "source       source for this calculation")') trim(source)
  write(*, '(" format    ", a, t28, "oformat      format for output")') trim(oformat)
!
! Isotope production
!
  if (flagprod) then
    write(*, '(" #"/" # Isotope production "/" #")')
    write(*, '(" Ebeam            ", f8.3, " beam         incident energy in MeV for isotope production")') Ebeam
    write(*, '(" Eback            ", f8.3, " Eback        lower end of energy range in MeV for isotope production")') Eback
    write(*, '(" radiounit             ", a3, " radiounit    unit for radioactivity")') radiounit
    write(*, '(" yieldunit             ", a3, " yieldunit    unit for isotope yield")') yieldunit
    write(*, '(" Ibeam            ", f8.3, " Ibeam        beam current in mA")') Ibeam
    do i = 1, 5
      if (Tirrad(i) > 0) write(*, '(" Tirrad      ", i9, "     Tirrad       ", a1, " of irradiation time")') &
 &      Tirrad(i), unitTirrad(i)
    enddo
    write(*, '(" Area             ", f8.3, " Area         target area in cm^2")') Area
    do i = 1, 5
      if (Tcool(i) > 0) write(*, '(" Tcool       ", i9, "     Tcool        ", a1, " of cooling time")') Tcool(i), unitTcool(i)
    enddo
    write(*, '(" rho               ", f7.3, " rhotarget    target density [g/cm^3] ")') rhotarget
  endif
!
! 3. Optical model
!
  write(*, '(" #"/" # Optical model"/" #")')
  write(*, '(" localomp            ", a1, "     flaglocalomp flag for local (y) or global (n) optical model")') yesno(flaglocalomp)
  write(*, '(" dispersion          ", a1, "     flagdisp     flag for dispersive optical model")') yesno(flagdisp)
  write(*, '(" jlmomp              ", a1, "     flagjlm      flag for using semi-microscopic JLM OMP")') yesno(flagjlm)
  write(*, '(" pruitt              ", a1, "     flagpruitt   identifier for Pruitt parameters for KD03")') pruitt
  write(*, '(" pruittset         ", i3, "     pruittset    random set for Pruitt et al OMP")') pruittset
  write(*, '(" riplomp             ", a1, "     flagriplomp  flag for RIPL OMP")') yesno(flagriplomp)
  write(*, '(" riplrisk            ", a1, "     flagriplrisk flag for going outside RIPL mass validity range")') yesno(flagriplrisk)
  write(*, '(" optmodall           ", a1, "     flagompall   flag for new optical model calculation for all residual nuclei")') &
 &  yesno(flagompall)
  write(*, '(" incadjust           ", a1, "     flagincadj   flag for OMP adjustment on incident channel also")') yesno(flagincadj)
  write(*, '(" omponly             ", a1, "     flagomponly  flag to execute ONLY an optical model calculation")') &
 &  yesno(flagomponly)
  write(*, '(" autorot             ", a1, "     flagautorot  flag for automatic rotational coupled channels ", &
 &  "calculations for A > 150")') yesno(flagautorot)
  write(*, '(" spherical           ", a1, "     flagspher    flag to force spherical optical model")') yesno(flagspher)
  write(*, '(" soukho              ", a1, "     flagsoukho   flag for Soukhovitskii OMP for actinides")') yesno(flagsoukho)
  write(*, '(" coulomb             ", a1, "     flagcoulomb  flag for Coulomb excitation calculation with ECIS")') &
 &  yesno(flagcoulomb)
  write(*, '(" statepot            ", a1, "     flagstate    flag for optical model potential for each excited state")') &
 &  yesno(flagstate)
  write(*, '(" maxband           ", i3, "     maxband      highest vibrational band added to rotational model")') maxband
  write(*, '(" maxrot            ", i3, "     maxrot       number of included excited rotational levels")') maxrot
  sysstring = '            '
  i = - 1
  do type = 1, 6
    if (flagsys(type)) then
      i = i + 2
      write(sysstring(i:i), '(a1)') parsym(type)
    endif
  enddo
  write(*, '(" sysreaction  ", a12, " sysreaction  particles with reaction cross section from systematics")') sysstring
  rotstring = '            '
  i = - 1
  do type = 1, 6
    if (flagrot(type)) then
      i = i + 2
      write(rotstring(i:i), '(a1)') parsym(type)
    endif
  enddo
  write(*, '(" rotational   ", a12, " rotational   particles with possible rotational optical model")') rotstring
  write(*, '(" core              ", i3, "     core         even-even core for weakcoupling (-1 or 1)")') core
  write(*, '(" ecissave            ", a1, "     flagecissave flag for saving ECIS input and output files")') yesno(flagecissave)
  write(*, '(" eciscalc            ", a1, "     flageciscalc flag for new ECIS calculation for outgoing particles and", &
 &  " energy grid")') yesno(flageciscalc)
  write(*, '(" inccalc             ", a1, "     flaginccalc  flag for new ECIS calculation for incident channel")') &
 &  yesno(flaginccalc)
  write(*, '(" endfecis            ", a1, "     flagendfecis flag for new ECIS calculation for ENDF-6 files")') &
 &  yesno(flagendfecis)
  write(*, '(" radialmodel        ", i2, "     radialmodel  model for radial matter densities (JLM OMP only)")') radialmodel
  write(*, '(" jlmmode            ", i2, "     jlmmode      option for JLM imaginary potential normalization")') jlmmode
  write(*, '(" alphaomp           ", i2, "     alphaomp     alpha OMP (1=normal, 2= McFadden-Satchler,", &
 &  " 3-5= folding potential, 6,8= Avrigeanu, 7=Nolte)")') alphaomp
  write(*, '(" deuteronomp        ", i2, "     deuteronomp  deuteron OMP (1=normal, 2=Daehnick,", &
 &  " 3=Bojowald, 4=Han-Shi-Shen, 5=An-Cai)")') deuteronomp
!
! 4. Compound nucleus
!
  write(*, '(" #"/" # Compound nucleus"/" #")')
  if (Ninc > 1 .and. enincmin < ewfc .and. enincmax >= ewfc) then
    write(*, '(" widthfluc        ", f8.3, " ewfc         off-set incident energy for width fluctuation calculation")') ewfc
  else
    write(*, '(" widthfluc           ", a1, "     flagwidth    flag for width fluctuation calculation")') yesno(flagwidth)
  endif
  write(*, '(" widthmode          ", i2, "     wmode        designator for width fluctuation model")') wmode
  write(*, '(" WFCfactor          ", i2, "     WFCfactor    enhancement factor for WFC: 1: Original, 2: Ernebjerg and Herman", &
 &  " 3: Kawano")') WFCfactor
  write(*, '(" compound            ", a1, "     flagcomp     flag for compound nucleus model")') yesno(flagcomp)
  write(*, '(" fullhf              ", a1, "     flagfullhf   ", &
 &  "flag for full spin dependence of transmission coefficients")') yesno(flagfullhf)
  write(*, '(" eciscompound        ", a1, "     flageciscomp flag for compound nucleus calculation by ECIS")') yesno(flageciscomp)
  write(*, '(" cpang               ", a1, "     flagcpang    flag for", &
 &  " compound angular distribution calculation for incident charged particles")') yesno(flagcpang)
  if (Ninc > 1 .and. enincmin < eurr .and. enincmax >= eurr) then
    write(*, '(" urr              ", f8.3, " eurr         off-set incident energy for URR calculation")') eurr
  else
    write(*, '(" urr                 ", a1, "     flagurr      flag for URR calculation")') yesno(flagurr)
  endif
  write(*, '(" urrnjoy             ", a1, "     flagurrnjoy  flag for normalization of URR parameters with NJOY method")') &
 &  yesno(flagurrnjoy)
  write(*, '(" lurr              ", i3, "     lurr         maximal orbital angular momentum for URR")') lurr
!
! 5. Gamma emission
!
  write(*, '(" #"/" # Gamma emission"/" #")')
  write(*, '(" gammax             ", i2, "     gammax       number of l-values for gamma multipolarity")') gammax
  write(*, '(" strength           ", i2, "     strength     model for E1 gamma-ray strength function")') strength
  write(*, '(" strengthM1         ", i2, "     strengthM1   model for M1 gamma-ray strength function")') strengthM1
  write(*, '(" pseudoresonances    ", a1, "    flagpseudores flag for using light nuclide discrete levels for resonances")') &
 &  yesno(flagpseudores)
  write(*, '(" electronconv        ", a1, "     flagelectron flag for application of electron conversion coefficient")') &
 &  yesno(flagelectron)
  write(*, '(" racap               ", a1, "     flagracap    flag for radiative capture model")') yesno(flagracap)
  write(*, '(" ldmodelracap       ", i2, "     ldmodelracap level density model for direct radiative capture")') ldmodelracap
  write(*, '(" upbend              ", a1, "     flagupbend   flag for low-energy upbend of photon strength function")') &
 &  yesno(flagupbend)
  write(*, '(" psfglobal           ", a1, "    flagpsfglobal flag for global photon strength functions only")') yesno(flagpsfglobal)
  write(*, '(" gnorm               ", a1, "     flaggnorm    flag to normalize PSF to average radiative width")') yesno(flaggnorm)
!
! 6. Pre-equilibrium
!
  write(*, '(" #"/" # Pre-equilibrium"/" #")')
  if (Ninc > 1 .and. enincmin < epreeq .and. enincmax >= epreeq) then
    write(*, '(" preequilibrium   ", f8.3, " epreeq       on-set incident energy for preequilibrium calculation")') epreeq
  else
    write(*, '(" preequilibrium      ", a1, "     flagpreeq    flag for pre-equilibrium calculation")') yesno(flagpreeq)
  endif
  write(*, '(" preeqmode          ", i2, "     preeqmode    designator for pre-equilibrium model")') preeqmode
  if (Ninc > 1 .and. enincmin < emulpre .and. enincmax >= emulpre) then
    write(*, '(" multipreeq       ", f8.3, " emulpre      on-set incident energy for multiple preequilibrium")') emulpre
  else
    write(*, '(" multipreeq          ", a1, "     flagmulpre   flag for multiple pre-equilibrium calculation")') yesno(flagmulpre)
  endif
  write(*, '(" mpreeqmode         ", i2, "     mpreeqmode   designator for multiple pre-equilibrium model")') mpreeqmode
  write(*, '(" breakupmodel       ", i2, "     breakupmodel model for break-up reaction: 1. Kalbach 2. Avrigeanu")') breakupmodel
  write(*, '(" phmodel            ", i2, "     phmodel      particle-hole state density model")') phmodel
  write(*, '(" pairmodel          ", i2, "     pairmodel    designator for pre-equilibrium pairing model")') pairmodel
  write(*, '(" preeqspin           ", i1, "     pespinmodel  model for pre-equilibrium spin distribution")') pespinmodel
  write(*, '(" giantresonance      ", a1, "     flaggiant    flag for collective contribution from giant resonances")') &
 &  yesno(flaggiant0)
  write(*, '(" preeqsurface        ", a1, "     flagsurface  flag for surface effects in exciton model")') yesno(flagsurface)
  write(*, '(" preeqcomplex        ", a1, "     flagpecomp   flag for Kalbach complex particle emission model")') yesno(flagpecomp)
  write(*, '(" twocomponent        ", a1, "     flag2comp    flag for two-component pre-equilibrium model")') yesno(flag2comp)
  write(*, '(" ecisdwba            ", a1, "     flagecisdwba flag for new ECIS calculation for DWBA for MSD")') yesno(flagecisdwba)
  write(*, '(" onestep             ", a1, "     flagonestep  flag for continuum one-step direct only")') yesno(flagonestep)
!
! 7. Level densities
!
  write(*, '(" #"/" # Level densities"/" #")')
  write(*, '(" ldmodel            ", i2, "     ldmodelall   level density model")') ldmodelall
  write(*, '(" ldmodelCN          ", i2, "     ldmodelCN    level density model for compound nucleus")') ldmodelCN
  write(*, '(" shellmodel         ", i2, "     shellmodel   model for shell correction energies")') shellmodel
  write(*, '(" kvibmodel          ", i2, "     kvibmodel    model for vibrational enhancement")') kvibmodel
  write(*, '(" spincutmodel       ", i2, "     spincutmodel model for spin cutoff factor for ground state")') spincutmodel
  write(*, '(" asys                ", a1, "     flagasys     flag for all level density parameters a from systematics")') &
 &  yesno(flagasys)
  write(*, '(" parity              ", a1, "     flagparity   flag for non-equal parity distribution")') yesno(flagparity)
  write(*, '(" colenhance          ", a1, "     flagcolall   flag for collective enhancement of level density", &
 &  " for all nuclides")') yesno(flagcolall)
  write(*, '(" ctmglobal           ", a1, "     flagctmglob  flag for global CTM model (no discrete level info)")') &
 &  yesno(flagctmglob)
  write(*, '(" gshell              ", a1, "     flaggshell   flag for energy dependence of single particle", &
 &  " level density parameter g")') yesno(flaggshell)
  write(*, '(" colldamp            ", a1, "     flagcolldamp flag for damping of collective effects in effective", &
 &  " level density")') yesno(flagcolldamp)
!
! 8. Fission
!
  write(*, '(" #"/" # Fission"/" #")')
  write(*, '(" fission             ", a1, "     flagfission  flag for fission")') yesno(flagfission)
  write(*, '(" fismodel           ", i2, "     fismodel     fission model")') fismodel
  write(*, '(" fismodelalt        ", i2, "     fismodelalt  alternative fission model for default barriers")') fismodelalt
  write(*, '(" hbstate             ", a1, "     flaghbstate  flag for head band states in fission")') yesno(flaghbstate)
  write(*, '(" class2              ", a1, "     flagclass2   flag for class2 states in fission")') yesno(flagclass2)
  write(*, '(" fispartdamp         ", a1, "  flagfispartdamp flag for fission partial damping")') yesno(flagfispartdamp)
  write(*, '(" massdis             ", a1, "     flagmassdis  flag for calculation of fission fragment mass yields")') &
 &  yesno(flagmassdis)
  write(*, '(" ffevaporation       ", a1, "     flagffevap   flag for calculation of particle evaporation", &
 &  " from fission fragment mass yields")') yesno(flagffevap)
  write(*, '(" fisfeed             ", a1, "     flagfisfeed  flag for output of fission per excitation bin")') yesno(flagfisfeed)
  write(*, '(" fymodel             ", i1, "     fymodel      fission yield model, 1: Brosa 2: GEF")') fymodel
  write(*, '(" ffmodel             ", i1, "     ffmodel      fission fragment model, 1: GEF 2: HF3D (Okumura) 3: SPY", &
 &  " 4: Langevin-4D 0: own files")') ffmodel
  write(*, '(" pfnsmodel           ", i1, "     pfnsmodel    PFNS model, 1: Iwamoto 2: from FF decay")') pfnsmodel
  write(*, '(" ffspin              ", a1, "     flagffspin   flag to use spin distribution in initial FF population")') &
 &  yesno(flagffspin)
!
! 9. Output
!
  write(*, '(" #"/" # Output"/" #")')
  write(*, '(" outmain             ", a1, "     flagmain     flag for main output")') yesno(flagmain)
  write(*, '(" outbasic            ", a1, "     flagbasic    flag for output of basic information and results")') yesno(flagbasic)
  write(*, '(" outall              ", a1, "     flagoutall   flag for output of all data in main output file")') yesno(flagoutall)
  write(*, '(" outpopulation       ", a1, "     flagpop      flag for output of population")') yesno(flagpop)
  write(*, '(" outcheck            ", a1, "     flagcheck    flag for output of numerical checks")') yesno(flagcheck)
  write(*, '(" outlevels           ", a1, "     flaglevels   flag for output of discrete level information")') yesno(flaglevels)
  write(*, '(" outdensity          ", a1, "     flagdensity  flag for output of level densities")') yesno(flagdensity)
  write(*, '(" outomp              ", a1, "     flagoutomp   flag for output of optical model parameters")') yesno(flagoutomp)
  write(*, '(" outkd               ", a1, "     flagoutkd    flag for output of KD03 OMP parameters")') yesno(flagoutkd)
  write(*, '(" outdirect           ", a1, "     flagdirect   flag for output of direct reaction results")') yesno(flagdirect)
  write(*, '(" outinverse          ", a1, "     flaginverse  flag for output of transmission coefficients", &
 &  " and inverse reaction cross sections")') yesno(flaginverse)
  write(*, '(" outdecay            ", a1, "     flagdecay    flag for output of decay of each population bin")') yesno(flagdecay)
  write(*, '(" outtransenergy      ", a1, "     flagtransen  flag for output of transmission coefficients per energy")') &
 &  yesno(flagtransen)
  write(*, '(" outecis             ", a1, "     flagoutecis  flag for output of ECIS results")') yesno(flagoutecis)
  write(*, '(" outgamma            ", a1, "     flaggamma    flag for output of gamma-ray information")') yesno(flaggamma)
  write(*, '(" outpreequilibrium   ", a1, "     flagpeout    flag for output of pre-equilibrium results ")') yesno(flagpeout)
  write(*, '(" outfission          ", a1, "     flagfisout   flag for output of fission information")') yesno(flagfisout)
  write(*, '(" outdiscrete         ", a1, "     flagdisc     flag for output of discrete state cross sections")') yesno(flagdisc)
  write(*, '(" outspectra          ", a1, "     flagspec     flag for output of double-differential cross sections")') &
 &  yesno(flagspec)
  write(*, '(" outbinspectra       ", a1, "     flagbinspec  flag for output of emission spectrum per excitation bin")') &
 &  yesno(flagbinspec)
  write(*, '(" resonance           ", a1, "     flagres      flag for output of low energy resonance cross sections")') &
 &  yesno(flagres)
  write(*, '(" group               ", a1, "     flaggroup    flag for output of low energy groupwise cross sections")') &
 &  yesno(flaggroup)
  if (Ninc > 1 .and. enincmin < eadd .and. enincmax >= eadd) then
    write(*, '(" adddiscrete      ", f8.3, " eadd         on-set incident energy for addition of discrete peaks to spectra")') eadd
  else
    write(*, '(" adddiscrete         ", a1, "     flagadd      flag for addition of discrete states to spectra")') yesno(flagadd)
  endif
  if (Ninc > 1 .and. enincmin < eaddel .and. enincmax >= eaddel) then
    write(*, '(" addelastic       ", f8.3, " eaddel       on-set incident energy addition of elastic peak to spectra")')  eaddel
  else
    write(*, '(" addelastic          ", a1, "     flagaddel    flag for addition of elastic peak to spectra")') yesno(flagaddel)
  endif
  write(*, '(" outangle            ", a1, "     flagang      flag for output of angular distributions")') yesno(flagang)
  write(*, '(" outlegendre         ", a1, "     flaglegendre flag for output of Legendre coefficients")') yesno(flaglegendre)
  write(*, '(" ddxmode             ", i1, "     ddxmode      mode for double-differential cross sections")') ddxmode
  write(*, '(" outdwba             ", a1, "     flagoutdwba  flag for output of DWBA cross sections for MSD")') yesno(flagoutdwba)
  write(*, '(" outgamdis           ", a1, "     flaggamdis   flag for output of discrete gamma-ray intensities")') yesno(flaggamdis)
  write(*, '(" outexcitation       ", a1, "     flagexc      flag for output of excitation functions")') yesno(flagexc)
  write(*, '(" components          ", a1, "     flagcompo    flag for output of cross section components")') yesno(flagcompo)
  write(*, '(" endf                ", a1, "     flagendf     flag for information for ENDF-6 file")') yesno(flagendf)
  write(*, '(" endfdetail          ", a1, "     flagendfdet  flag for detailed ENDF-6 information per channel")') yesno(flagendfdet)
  write(*, '(" sacs                ", a1, "     flagsacs     flag for statistical analysis of cross sections")') yesno(flagsacs)
  write(*, '(" partable            ", a1, "     flagpartable flag for output of model parameters on separate file")') &
 &  yesno(flagpartable)
  write(*, '(" block               ", a1, "     flagblock    flag to block spectra, angle and gamma files")') yesno(flagblock)
  return
end subroutine inputout
! Copyright A.J. Koning 2021
