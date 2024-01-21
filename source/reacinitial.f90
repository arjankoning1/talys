subroutine reacinitial
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialization of arrays for various cross sections
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
! Variables for numerics
!   maxN              ! maximal number of neutrons away from initial compound nucleus
!   maxZ              ! maximal number of protons away from initial compound nucleus
! Variables for energies
!   Etotal            ! total energy of compound system (target + projectile)
! Variables for compound nucleus from target
!   Exinc             ! excitation energy of entrance bin
! Variables to normalize compound nucleus cross section
!   CNterm    ! compound nucleus formation cross section per spin
! Variables for incident channel
!   lmaxinc           ! maximal l - value for transm. coeff. for incident channel
!   maxA              ! maximal number of nucleons away from initial compound nucleus
! Variables to normalize compound nucleus cross section
!   J2beg             ! begin of J summation
!   J2end             ! end of J summation
! Variables for nuclides
!   primary           ! flag to designate primary (binary) reaction
! Variables for preequilibrium
!   ENHratio          !  breakup nucleons enhancing reaction cross
!   Spre              ! time-integrated strength of two-component exciton state
!   xsBF              ! nucleon inelastic breakup cross section
!   xsBFnuc           ! inelastic breakup enhancement brought by breakup neutrons
!   xsBUnuc           ! nucleon breakup cross section
!   xsEB              ! elastic breakup cross section
!   xspopph           ! population cross section per particle-hole configuration
!   xspopph2          ! population cross section per two-component particle-hole configuration
!   xspopnucT         ! total population cross section per nucleus including inel BU
!   xspreeqad         ! preequilibrium angular distribution per particle type and outgoing energy
!   xspreeqbu         ! preequilibrium cross section per particle type and outgoing energy for break-up
!   xspreeqdisc       ! preequilibrium cross section for discrete state
!   xspreeqdisctot    ! preequilibrium cross section summed over discrete states
!   xspreeqdiscsum    ! total preequilibrium cross section for discrete states
!   xspreeqJP         ! preeq. cross section per particle type, outgoing energy, J, P
!   xspreeqki         ! preequilibrium cross section per particle type and outgoing energy for knock-out
!   xspreeqps         ! preequilibrium cross section per particle type and outgoing energy for pickup and stripping
!   xspreeqtotbu      ! preequilibrium cross section per particle type for breakup
!   xspreeqtotki      ! preequilibrium cross section per particle type for knockout and inelastic
!   xspreeqtotps      ! preequilibrium cross section per particle type for pickup and stripping
!   xsstep            ! preeq. cross section per particle type, stage and outgoing E
!   xsstep2           ! two-component preequilibrium cross section
!   xssteptot         ! preequilibrium cross section per particle type and stage
!   wemission         ! emission rate per particle, exciton number and energy
!   wemission2        ! two-component emission rate
! Variables for MSD
!   Emsd            ! minimal outgoing energy for MSD calculation
!   msdstep         ! continuum n - step direct cross section
!   msdstep0        ! n - step cross section for MSD
!   msdstep1        ! continuum one - step direct cross section (unnormalized)
!   msdstepad       ! continuum n - step direct angular distribution
!   msdstepad0      ! n - step angular distribution for MSD
!   msdstepad1      ! continuum one - step direct angular distribution (unnormalized)
!   msdstepint      ! n - step direct cross section integrated over energy
!   msdstepintad    ! n - step direct angular distribution integrated over energy
!   msdsum          ! multi - step direct cross section summed over steps and integrated over energy
!   msdtot          ! multi - step direct cross section summed over steps
!   msdtotad        ! multi - step direct angular distribution summed over steps
!   msdtotintad     ! multi - step direct angular distribution summed over steps and integrated over energy
!   nangleint       ! number of possibilities to link intermediate angle to final angle
!   numJmsd         ! maximum spin for MSD
!   xscont          ! continuum one - step direct cross section
!   xscont1         ! continuum one - step direct cross section (unnormalized)
!   xscontad        ! continuum one - step direct angular distribution for MSD
!   xscontad1       ! continuum one - step direct angular distribution for MSD (unnormalized)
!   xsdw            ! DWBA angular distribution as a function of incident energy, outgoing energy, ang. mom. and angle
!   xsdwin          ! DWBA cross section as a function of incident energy, outgoing energy and angular momentum
! Variables for preequilibrium
!   Esurf             ! well depth for surface interaction
!   PP2           ! total strength
! Variables for giant resonances
!   collcontad    ! collective angular distribution in the continuum
!   eoutgr        ! emission energy
!   grcollad      ! giant resonance angular distribution
!   xscollcont    ! total collective cross section in the continuum
!   xsgrad        ! smoothed giant resonance angular distribution
!   xsgrcoll      ! giant resonance cross section
!   xsgrstate     ! smoothed giant resonance cross section per state
! Variables for ENDF data
!   e6           ! energies of ENDF-6 energy grid in MeV
!   xscompel6    ! compound elastic cross section
!   xsnonel6     ! non-elastic cross section
!   xselas6      ! total elastic cross section (neutrons only) for ENDF-6 file
!   xselassh6    ! shape elastic cross section (neutrons only) for ENDF-6 file
!   xsnon6       ! non-elastic cross section for ENDF-6 file
!   xsopt6       ! optical model reaction cross section for ENDF-6 file
!   xsreac6      ! reaction cross section for ENDF-6 file
!   xstot6       ! total cross section for ENDF-6 file
! Variables for incident channel
!   channelsum      ! sum over exclusive channel cross sections
!   cleg            ! compound nucleus Legendre coefficient
!   contrib         ! contribution to emission spectrum
!   directad        ! direct angular distribution
!   dleg            ! direct reaction Legendre coefficient
!   dorigin         ! origin of direct cross section (Direct or Preeq)
!   xsgr            ! total smoothed giant resonance cross section
!   multiplicity    ! particle multiplicity
!   partdecay       ! total decay per particle
!   popdecay        ! decay from population
!   preeqpop        ! pre-equilibrium population
!   preeqpopex      ! pre-equilibrium population
!   ruth            ! elastic/Rutherford ratio
!   elasni          ! nuclear+interference term
!   Tjlinc          ! transm. coeff. as a function of spin and l for inc. channel
!   Tlinc           ! transm. coeff. as a function of l for incident channel
!   xsabs           ! absorption cross section
!   xsbinary        ! cross section from initial compound to residual nucleus
!   xsbranch        ! branching ratio for isomeric cross section
!   xscollconttot   ! total collective cross section in the continuum
!   xscompcont      ! compound cross section for continuum
!   xsdirdisc       ! direct cross section for discrete state direct cross section
!   xsdirdisctot    ! direct cross section summed over discrete states
!   xsdirdiscsum    ! total direct cross section
!   xselasinc       ! total elastic cross section (neutrons only) for inc. channel
!   xsgrtot         ! total smoothed giant resonance cross section
!   xsgrsum         ! sum over giant resonance cross sections
!   xsmassprod      ! residual production cross section per mass unit
!   xsngnsum        ! sum over total (projectile, gamma - ejectile) cross section
!   xsoptinc        ! optical model reaction cross section for incident channel
!   xsparticle      ! total particle production cross section
!   xspop           ! population cross section
!   xspopex         ! population cross section summed over spin and parity
!   xspopexP        ! population cross section per parity
!   xspopnuc        ! population cross section per nucleus
!   xspopnucP       ! population cross section per nucleus per parity
!   xspreeq         ! preeq. cross section per particle type and outgoing energy
!   xspreeqsum      ! total preequilibrium cross section summed over particles
!   xspreeqtot      ! preequilibrium cross section per particle type
!   xsreacinc       ! reaction cross section for incident channel
!   xsresprod       ! total residual production ( = reaction)  cross section
!   xstotinc        ! total cross section (neutrons only) for incident channel
! Variables for spectra
!   buratio          ! break-up ratio
!   compspect        ! compound part of spectrum
!   Eaverage         ! average outgoing energy
!   eendout          ! last energy point of energy grid
!   espec            ! outgoing energy grid
!   preeqratio       ! pre-equilibrium ratio
!   preeqspect       ! multiple pre-equilibrium part of spectrum
!   xscompout        ! compound emission angular distribution
!   xscompoutad      ! compound emission angular distribution
!   xsdiscout        ! smoothed angular distribution for discrete state
!   xsdiscoutad      ! smoothed angular distribution for discrete state
!   xsmpreeqout      ! multiple preequilibrium angular distribution
!   xsmpreeqoutad    ! multiple preequilibrium angular distribution
!   xspreeqbuout     ! preequilibrium cross section for breakup
!   xspreeqkiout     ! preequilibrium cross section for knockout and inelastic
!   xspreeqout       ! preequilibrium angular distribution per particle type an
!   xspreeqoutad     ! preequilibrium angular distribution per particle type
!   xspreeqpsout     ! preequilibrium cross section for pickup and stripping
!   xssumout         ! cross section summed over mechanisms
!   xssumoutad       ! angular distribution summed over mechanisms
! Variables for excitation energy grid
!   deltaEx      ! excitation energy bin for population arrays
!   Ex           ! excitation energy
!   Exmax        ! maximum excitation energy for residual nucleus
!   Exmax0       ! maximum excitation energy for res. nucleus (incl. neg. en.)
!   fisfeedJP    ! fission contribution from excitation energy bin per J,P
!   maxex        ! maximum excitation energy bin for residual nucleus
!   maxJ         ! maximal J-value
!   nexmax       ! maximum excitation energy bin for residual nucleus
!   rhogrid      ! integrated level density
! Variables for direct capture
!   xsracappop    ! population cross section for radiative capture
!   xsracappopex  ! population cross section for radiative capture
! Variables for ENDF data
!   e6           ! energies of ENDF-6 energy grid in MeV
!   xscompel6    ! compound elastic cross section
!   xsnonel6     ! non-elastic cross section
!   xselas6      ! total elastic cross section (neutrons only) for ENDF-6 file
!   xselassh6    ! shape elastic cross section (neutrons only) for ENDF-6 file
!   xsopt6       ! optical model reaction cross section for ENDF-6 file
!   xsreac6      ! reaction cross section for ENDF-6 file
!   xstot6       ! total cross section for ENDF-6 file
! Variables for thermal cross sections
!   fnubar           ! nubar
!   fexclbranch      ! exclusive channel yield per isomer
!   fxsbinary        ! cross section from initial compound to residual n
!   fxsbranch        ! branching ratio for isomeric cross section
!   fxschaniso       ! channel cross section per isomer
!   fxschannel       ! channel cross section
!   fxscompdisc      ! compound cross section for discrete state
!   fxscompel        ! compound elastic cross section
!   fxscompnonel     ! total compound non-elastic cross section
!   fxsdirdisc       ! direct cross section for discrete state
!   fxsdirdiscsum    ! total direct cross section
!   fxsdisc          ! total cross section for discrete state
!   fxsdisctot       ! total cross section summed over discrete states
!   fxselasinc       ! total elastic cross section (neutrons only) for i
!   fxselastot       ! total elastic cross section (neutrons only) for i
!   fxsexclcont      ! exclusive single channel cross section for contin
!   fxsexclusive     ! exclusive single channel cross section
!   fxsgamchannel    ! gamma channel cross section
!   fxsgamdischan    ! discrete gamma channel cross section
!   fxsngn           ! total (n,gn) cross section
!   fxsnonel         ! non-elastic cross section for incident channel
!   fxspopex         ! population cross section summed over spin and par
!   fxspreeqsum      ! total preequilibrium cross section summed over pa
!   fxsracape        ! direct capture cross section
!   fxspopnuc        ! population cross section per nucleus
!   fxsratio         ! ratio of exclusive cross section over residual pr
!   fxsreacinc       ! reaction cross section for incident channel
!   fxstotinc        ! total cross section (neutrons only) for incident
! Variables for binary reactions
!   feedbinary       ! feeding from first compound nucleus
!   binemissum       ! integrated binary emission spectrum
!   Eaveragebin      ! average outgoing energy
!   xscompdisc       ! compound cross section for discrete state
!   xscompdisctot    ! compound cross section summed over discrete states
!   xscompel         ! compound elastic cross section
!   xscompnonel      ! total compound non-elastic cross section
!   xscompound       ! total compound cross section
!   xsconttot        ! total cross section for continuum
!   xsdircont        ! direct cross section for continuum
!   xsdirect         ! total direct cross section
!   xsdisc           ! total cross section for discrete state
!   xsdisctot        ! total cross section summed over discrete states
!   xselastot        ! total elastic cross section (shape + compound)
!   xsnonel          ! non-elastic cross section
!   xspopdir         ! direct population cross section per nucleus
!   xspopex0         ! binary population cross section
! Variables to prepare information for initial compound nucleus
!   enumhf     ! enumerator for compound nucleus formula
!   transjl    ! array for width fluctuation calculation
! Variables for binary emission spectra
!   binemis        ! emission spectra from initial compound nucleus
!   binnorm        ! normalization factor for binary spectra
!   xsbinemis      ! cross section for emission from first compound nucleus
!   xsbinemisad    ! angular distribution for emission from first compound nucleus
!   xscomp         ! compound elastic cross section
!   xscompad       ! compound emission angular distribution
!   xsemis         ! cross section for emission from compound nucleus
! Variables for angular distributions
!   cleg0      ! Legendre coefficient normalized to the first one
!   compad     ! compound angular distribution
!   discad     ! discrete state angular distribution
!   tleg       ! total Legendre coefficient
!   tlegnor    ! total Legendre coefficient normalized to 1
! Variables for astro initialization
!   Tastroinc         ! transmission coefficient for incident channel (Astrophysic
!   Tastroout         ! transmission coefficient for outgoing channel (Astrophysic
! Variables for energy grid, level density and transmission coefficients
!   lmaxhf      ! maximal l-value for transmission coefficients
!   nbintfis    ! number of bins
!   rho0        ! integrated level density
!   Tgam        ! gamma transmission coefficients
!   Tjlnex      ! transmission coefficients as a function of particle type, energy,
!   Tlnex       ! transmission coefficients as a function of particle type, energy
! Variables for fission transmission coefficients
!   denfis      ! fission level density
!   gamfis      ! fission width
!   rhofisA     ! integrated level density corresponding to tfisA
!   taufis      ! fission lifetime
!   tfis        ! fission transmission coefficient for Hill-Wheeler magnitude
!   tfisA       ! transmission coefficient for Hill-Wheeler magnitude
!   tfisdown    ! fission transmission coefficients
!   tfisup      ! fission transmission coefficients
! Variables for URR
!   sigurrc       ! (l,j) capture cross section for URR
!   sigurrf       ! (l,j) fission cross section for URR
!   sigurrs       ! (l,j) scattering cross section for URR
!   spot          ! potential scattering contribution
!   strengthl     ! l neutron strength function
!   strengthlj    ! (l,j) neutron strength function
!   xsurrN        ! URR cross section
!   xsurrT        ! URR cross section
!   urrwidth      ! channel width in URR
! Variables for multiple emission
!   Dmulti         ! depletion factor for multiple preequilibrium
!   Eaveragemul    ! average outgoing energy
!   Fcomp          ! compound population fraction per nucleus
!   Fdir           ! direct population fraction per nucleus
!   feedexcl       ! feeding terms from compound excitation ene
!   fisfeedex      ! fission contribution from excitation energy bin
!   Fpreeq         ! preequilibrium population fraction per nucleus
!   mcontrib       ! contribution to emission spectrum
!   mpecontrib     ! contribution to multiple pre-equilibrium emission spectr
!   popexcl        ! population cross section of bin just before decay
!   xsbinspec      ! emission spectrum from compound nucleus per bin
!   xsfeed         ! cross section from compound to residual nucleus
!   xsgamdis       ! discrete gamma-ray cross section
!   xsgamdistot    ! total discrete gamma-ray cross section
!   xsinitpop      ! initial population cross section
!   xsmpe          ! multiple-preequilibrium cross section per energy bin
!   xsmpeemis      ! multiple-preequilibrium emission spectrum from compound n
!   xsmpetot       ! total multiple-preequilibrium cross section
!   xsmpreeq       ! multiple pre-equilibrium emission spectrum
!   xsmpreeqad     ! multiple preequilibrium angular distribution
!   xsngn          ! total (projectile,gamma-ejectile) cross section
!   xsngnspec      ! total (projectile,gamma-ejectile) spectrum
!   xspartial      ! emitted cross section flux per energy bin
!   xspopcomp      ! compound population cross section per nucleus
!   xspopnuc0      ! population cross section per nucleus
!   xspoppreeq     ! preequilibrium population cross section per nucleus
! Variables for compound nucleus from target
!   Fnorm         ! multiplication factor
!   JmaxU         ! maximal total angular momentum
!   JminU         ! minimal total angular momentum
!   lmaxU         ! maximal orbital angular momentum
!   lminU         ! minimal orbital angular momentum
!   nulj          ! (l,j) number of degrees of freedom for URR calculation
!   Purrlj        ! (l,j) parity for URR calculation
!   Turrlj        !  transmission coefficient for URR calculation
!   Turrljinc     ! incident channel (l,j) transmission coefficient for URR ca
!   xsbinarylj    ! cross section from initial compound to residual nucleus
! Variables for isotope production
!   Erp             ! incident energy
!   Nenrp           ! number of incident energies for residual production cross
!   Tgrid           ! time
!   Tp              ! irradiation time with maximal yield per time unit
!   xsrp            ! residual production cross section in mb
! Variables for total cross sections
!   xsexclcont      ! exclusive single channel cross section for continuum
!   xsexclusive     ! exclusive single channel cross section
! Variables for exclusive channels
!   exclbranch        ! exclusive channel yield per isomer
!   Eavchannel        ! channel average energy
!   Especsum          ! total emission energy
!   fisstring         ! string for exclusive fission reaction channel
!   gamexcl           ! exclusive gamma cross section per excitation ener
!   gmult             ! continuum gamma multiplicity
!   specemis          ! exclusive emission contribution
!   xschancheck       ! integrated channel spectra
!   xschaniso         ! channel cross section per isomer
!   xschannel         ! channel cross section
!   xschannelsp       ! channel cross section spectra
!   xsexcl            ! exclusive cross section per excitation energy
!   xsfischancheck    ! integrated fission channel spectra
!   xsfischannel      ! fission channel cross section
!   xsfischannelsp    ! fission channel cross section spectra
!   xsgamchannel      ! gamma channel cross section
!   xsgamdischan      ! discrete gamma channel cross section
!   xsparcheck        ! total particle production cross section
!   xsratio           !  ratio of exclusive cross section over residual p
!   xsspeccheck       ! total particle production spectra
!   yieldchannel      ! relative yield
! Variables for energies
!   Ethrexcl       ! threshold incident energy for exclusive channel
!   idchannel      ! identifier for exclusive channel
!   Qexcl          ! Q-value for exclusive channel
!   reacstring     ! string for exclusive reaction channel
!   mulpreZN       ! logical for multiple pre-equilibrium per nucleus
!
! *** Declaration of local data
!
  implicit none
!
! *************** Initialization of primary reaction *******************
!
  primary = .true.
!
! *************** Initialize incident reaction arrays ******************
!
  lmaxinc = 0
!
! ********************** Initialization of energies ********************
!
  Exinc = Etotal
!
! *************** Initialize pre-equilibrium arrays ********************
!
  ENHratio = 0.
  Spre = 0.
  xsBF = 0.
  xsBFnuc = 0.
  xsBUnuc = 0.
  xsEB = 0.
  xspopph = 0.
  xspopph2 = 0.
  xspopnucT = 0.
  xspreeqad = 0.
  xspreeqbu = 0.
  xspreeqdisc = 0.
  xspreeqdisctot = 0.
  xspreeqJP = 0.
  xspreeqki = 0.
  xspreeqps = 0.
  xspreeqtotbu = 0.
  xspreeqtotki = 0.
  xspreeqtotps = 0.
  xsstep = 0.
  xsstep2 = 0.
  xssteptot = 0.
  PP2 = 0.
  wemission = 0.
  wemission2 = 0.
  xspreeqdiscsum = 0.
  Esurf = 0.
  collcontad = 0.
  eoutgr = 0.
  grcollad = 0.
  xscollcont = 0.
  xsgrad = 0.
  xsgrcoll = 0.
  xsgrstate = 0.
  Emsd = 0.
  msdstep = 0.
  msdstep0 = 0.
  msdstep1 = 0.
  msdstepad = 0.
  msdstepad0 = 0.
  msdstepad1 = 0.
  msdstepint = 0.
  msdstepintad = 0.
  msdsum = 0.
  msdtot = 0.
  msdtotad = 0.
  msdtotintad = 0.
  nangleint = 0
  xscont = 0.
  xscont1 = 0.
  xscontad = 0.
  xscontad1 = 0.
  xsdw = 0.
  xsdwin = 0.
  buratio = 0.
  compspect = 0.
  Eaverage = 0
  eendout = 0
  espec = 0.
  preeqratio = 0.
  preeqspect = 0.
  xscompout = 0.
  xscompoutad = 0.
  xsdiscout = 0.
  xsdiscoutad = 0.
  xsmpreeqout = 0.
  xsmpreeqoutad = 0.
  xspreeqbuout = 0.
  xspreeqkiout = 0.
  xspreeqout = 0.
  xspreeqoutad = 0.
  xspreeqpsout = 0.
  xssumout = 0.
  xssumoutad = 0.
!
! ********* Initialization of primary compound and binary arrays *******
!
  CNterm = 0.
  J2beg = 0
  J2end = 0
!
! ************* Initialization of total cross section arrays ***********
!
  maxA = maxZ + maxN
  channelsum = 0.
  cleg = 0.
  contrib = 0.
  directad = 0.
  dleg = 0.
  dorigin = '      '
  multiplicity = 0.
  partdecay = 0
  popdecay = 0
  preeqpop = 0.
  preeqpopex = 0.
  ruth = 0.
  elasni = 0.
  Tjlinc = 0.
  Tlinc = 0.
  xsabs = 0.
  xsbranch = 0.
  xscollcontJP = 0.
  xscollconttot = 0.
  xscompcont = 0.
  xsbinary = 0.
  xsdirdisc = 0.
  xsdirdisctot = 0.
  xsdirdiscsum = 0.
  xsgr = 0.
  xsgrsum = 0.
  xsgrtot = 0.
  xsmassprod = 0.
  xstotinc = 0.
  xsreacinc = 0.
  xsoptinc = 0.
  xsparticle = 0.
  xselasinc = 0.
  xsngnsum = 0.
  xspop = 0.
  xspopex = 0.
  xspopexP = 0.
  xspopnuc = 0.
  xspopnucP = 0.
  xspreeq = 0.
  xspreeqtot = 0.
  xspreeqsum = 0.
  xsresprod = 0.
  deltaEx = 0.
  Ex = 0.
  Exmax = 0.
  Exmax0 = 0.
  Exmax(0, 0) = Etotal
  Exmax0(0, 0) = Etotal
  fisfeedJP = 0.
  maxex = 0
  maxJ = numJ
  nexmax = -1
  rhogrid = 0.
  enumhf = 0.
  transjl = 0.
  xsracape=0.
  xsracapedisc=0.
  xsracapecont=0.
  xsracappop = 0.
  xsracappopex = 0.
  if (nin == 1) then
    fnubar = 0.
    fexclbranch = 0.
    fxsbinary = 0.
    fxsbranch = 0.
    fxschaniso = 0.
    fxschannel = 0.
    fxscompdisc = 0.
    fxscompel = 0.
    fxscompnonel = 0.
    fxsdirdisc = 0.
    fxsdirdiscsum = 0.
    fxsdisc = 0.
    fxsdisctot = 0.
    fxselasinc = 0.
    fxselastot = 0.
    fxsexclcont = 0.
    fxsexclusive = 0.
    fxsgamchannel = 0.
    fxsgamdischan = 0.
    fxsngn = 0.
    fxsnonel = 0.
    fxspopex = 0.
    fxspreeqsum = 0.
    fxsracape = 0.
    fxspopnuc = 0.
    fxsratio = 0.
    fxsreacinc = 0.
    fxstotinc = 0.
    e6 = 0.
    xscompel6 = 0.
    xsnonel6 = 0.
    xselas6 = 0.
    xselassh6 = 0.
    xsnon6 = 0.
    xsopt6 = 0.
    xsreac6 = 0.
    xstot6 = 0.
  endif
  feedbinary = 0.
  binemissum = 0.
  Eaveragebin = 0.
  xscompdisc = 0.
  xscompdisctot = 0.
  xscompound = 0.
  xsconttot = 0.
  xsdircont = 0.
  xsdirect = 0.
  xsdisc = 0.
  xsdisctot = 0.
  xspopdir = 0.
  xspopex0 = 0.
  xscompel = 0.
  xselastot = 0.
  xsnonel = 0.
  xscompnonel = 0.
  binemis = 0.
  binnorm = 0.
  xsbinemis = 0.
  xsbinemisad = 0.
  xscomp = 0.
  xscompad = 0.
  xsemis = 0.
  cleg0 = 0.
  compad = 0.
  discad = 0.
  tleg = 0.
  tlegnor = 0.
  Tastroinc = 0.
  Tastroout = 0.
  lmaxhf = 0
  nbintfis = 0
  rho0 = 0.
  Tgam = 0.
  Tjlnex = 0.
  Tlnex = 0.
  denfis = 0.
  gamfis = 0.
  rhofisA = 1.
  taufis = 0.
  tfis = 0.
  tfisA = 0.
  tfisdown = 0.
  tfisup = 0.
  sigurrc = 0.
  sigurrs = 0.
  sigurrf = 0.
  spot = 0.
  strengthl = 0.
  strengthlj = 0.
  urrwidth = 0.
  Dmulti = 0.
  Fcomp = 0.
  Fdir = 0.
  feedexcl = 0.
  fisfeedex = 0.
  Fpreeq = 0.
  mcontrib = 0.
  mpecontrib = 0.
  popexcl = 0.
  xsbinspec = 0.
  xsfeed = 0.
  xsgamdis = 0.
  xsgamdistot = 0.
  xsmpe = 0.
  xsmpeemis = 0.
  xsmpetot = 0.
  xsmpreeq = 0.
  xsmpreeqad = 0.
  xsngn = 0.
  xsngnspec = 0.
  xspartial = 0.
  xspopcomp = 0.
  if (.not.flagffruns) xsfistot0=0.
  if (.not. flagrpruns) xspopnuc0 = 0.
  xspoppreeq = 0.
  xsinitpop = 0.
  Fnorm = 1.
  nulj = 0
  Purrlj = 1
  Turrlj = 0.
  Turrljinc = 0.
  xsbinarylj = 0.
  Erp = 0.
  Nenrp = 0
  Tgrid = 0.
  Tp = 0
  xsrp = 0.
  xsexclcont = 0.
  xsexclusive = 0.
  exclbranch = 0.
  Eavchannel = 0.
  Especsum = 0.
  fisstring = '                  '
  gamexcl = 0.
  gmult = 0.
  specemis = 0.
  xschancheck = 0.
  xschaniso = 0.
  xschannel = 0.
  xschannelsp = 0.
  xsexcl = 0.
  xsfischancheck = 0.
  xsfischannel = 0.
  xsfischannelsp = 0.
  xsgamchannel = 0.
  xsgamdischan = 0.
  xsparcheck = 0.
  xsratio = 0.
  xsspeccheck = 0.
  yieldchannel = 0.
  mulpreZN = .false.
  Ethrexcl = 0.
  idchannel = -1
  Qexcl = 0.
  reacstring = '                  '
  return
end subroutine reacinitial
! Copyright A.J. Koning 2021
