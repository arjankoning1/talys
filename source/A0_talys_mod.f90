module A0_talys_mod
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : General module with all global variables
!
! Author    : Arjan Koning
!
! 2023-12-30: Original code
! 2025-02-21: Current version
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Definition of single and double precision variables
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, parameter :: sgl = selected_real_kind(6,37)   ! single precision kind
  integer, parameter :: dbl = selected_real_kind(15,307) ! double precision kind
!
!-----------------------------------------------------------------------------------------------------------------------------------
! All global dimensions
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, parameter :: memorypar=6                          ! memory parameter
  integer, parameter :: numpar=6                             ! number of particles
  integer, parameter :: numiso=20                            ! maximum number of isotopes per element
  integer, parameter :: numelem=124                          ! number of elements
  integer, parameter :: numelem0=124                         ! number of elements
  integer, parameter :: numl=60                              ! number of l values
  integer, parameter :: numlines=5000                        ! number of input lines
  integer, parameter :: numpop=1000                          ! number of population bins
  integer, parameter :: numenin=600                          ! number of incident energies
  integer, parameter :: numZ=2+2*memorypar                   ! maximum number of protons from initial compound nucleus
  integer, parameter :: numN=10+4*memorypar                  ! maximum number of neutrons from initial compound nucleus
  integer, parameter :: numZph=4                             ! maximum number of protons away from the initial compound nucleus
  integer, parameter :: numNph=8                             ! maximum number of neutrons away from the initial compound nucleus
  integer, parameter :: numbar=3                             ! number of fission barriers
  integer, parameter :: nummt=200                            ! number of MT numbers
  integer, parameter :: nummtall=1000                        ! number of MT numbers
  integer, parameter :: numgam=6                             ! maximum number of l-values for gamma multipolarity
  integer, parameter :: numrange=10                          ! number of energy ranges
  integer, parameter :: numadj=500                           ! maximum number of adjustable parameters
  integer, parameter :: numenadj=1000                        ! maximum number of energies for adjustable parameters
  integer, parameter :: numlev=40                            ! maximum number of discrete levels
  integer, parameter :: numisom=10                           ! number of isomers
  integer, parameter :: numresgrid=1000                      ! number of energies on resonance grid
  integer, parameter :: numflux=100                          ! number of integral experiments
  integer, parameter :: numfile=100                          ! maximum number of separate output files
  integer, parameter :: numlev2=200                          ! maximum number of levels
  integer, parameter :: numrotcc=4                           ! number of rotational deformation parameters
  integer, parameter :: numgamqrpa=300                       ! number of energies for QRPA strength function
  integer, parameter :: numTqrpa=11                          ! number of temperatures for QRPA strength functions
  integer, parameter :: numomp=500                           ! maximum number of lines in optical model file
  integer, parameter :: numompadj=13                         ! number of adjustable ranges for OMP
  integer, parameter :: numjlm=200                           ! maximum number of radial points
  integer, parameter :: numrot=700                           ! number of rotational states
  integer, parameter :: nummatchT=4000                       ! maximum number of energy points for T matching
  integer, parameter :: numdens=60                           ! number of energy points for tabulated level densities
  integer, parameter :: numdensracap=200                     ! number of energy points for tabulated level densities for direct cap.
  integer, parameter :: numen=260                            ! maximum number of outgoing energies
  integer, parameter :: numang=90                            ! maximum number of angles
  integer, parameter :: numangcont=36                        ! maximum number of angles for continuum
  integer, parameter :: numconf=72                           ! number of particle-hole combinations
  integer, parameter :: numexc=12                            ! number of excitons
  integer, parameter :: numJph=30                            ! maximum spin for particle-hole states
  integer, parameter :: numfact=6*numl                       ! number of terms for factorial logarithm
  integer, parameter :: numbins=20*(memorypar-1)             ! maximum number of continuum excitation energy bins
  integer, parameter :: numex=numlev+numbins                 ! maximum number of excitation energies
  integer, parameter :: numJ=40                              ! maximum J-value
  integer, parameter :: numenrec=4*(memorypar-1)             ! maximum number of recoil energies
  integer, parameter :: numangrec=9                          ! maximum number of recoil angles
  integer, parameter :: numendisc=400                        ! number of discrete outgoing energies
  integer, parameter :: numen2=numen+numendisc               ! maximum number of outgoing energies
  integer, parameter :: numenmsd=18                          ! maximum number of energy points for DWBA calculation for MSD
  integer, parameter :: nummsd=5                             ! maximum  number of MSD steps
  integer, parameter :: numbinfis=1000                       ! maximum number of bins for fission calculation
  integer, parameter :: numbeta=200                          ! number of points for WKB fission calculation
  integer, parameter :: numhill=20                           ! maximum number of Hill-Wheeler points
  integer, parameter :: numtrans=numl*12*(numex+1)+numhill+1 ! number of transmission coefficients
  integer, parameter :: nummold=32                           ! number of integration points for Moldauer calculation
  integer, parameter :: numgoe=50                            ! number of integration points for GOE calculation
  integer, parameter :: numT=30                              ! number of temperatures
  integer, parameter :: numZastro=4                          ! maximal number of protons away from initial CN for astroph. calcs
  integer, parameter :: numNastro=4                          ! maximal number of neutrons away from initial CN for astroph. calcs
  integer, parameter :: numZchan=6                           ! maximum number of outgoing protons in individual channel description
  integer, parameter :: numNchan=10                          ! maximum number of outgoing neutron in individual channel description
  integer, parameter :: numin=8                              ! maximum number of neutrons in channel description
  integer, parameter :: numip=4                              ! maximum number of protons in channel description
  integer, parameter :: numid=2                              ! maximum number of deuterons in channel description
  integer, parameter :: numit=1                              ! maximum number of tritons in channel description
  integer, parameter :: numih=1                              ! maximum number of helions in channel description
  integer, parameter :: numia=3                              ! maximum number of alphas in channel description
  integer, parameter :: numchantot=35*(memorypar-1)          ! maximum number of exclusive channels
  integer, parameter :: nummass=414                          ! number of masses
  integer, parameter :: nummass0=414                         ! number of masses
  integer, parameter :: numneu=nummass-numelem               ! number of neutrons
  integer, parameter :: numneu0=numneu                       ! number of neutrons
  integer, parameter :: numA=numZ+numN                       ! maximum number of nucleons from initial compound nucleus
  integer, parameter :: numnu=50                             ! number of neutrons from fission
  integer, parameter :: numenlow=20                          ! number of energies for low energy regime
  integer, parameter :: numtime=100                          ! number of time points
  integer, parameter :: numP=1000000                         ! number of points in reconstructed data file
  integer, parameter :: numpfns=300                          ! number of energies for PFNS grid
  integer, parameter :: numenout=1000                        ! number of outgoing energies
  integer, parameter :: numenrp=200                          ! number of incident energies for residual products
  integer, parameter :: numen6=memorypar*1700                ! number of energies for ENDF6 energy grid
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=1), dimension(-1:1)       :: cparity     ! parity (character)
  character(len=1), dimension(-1:numisom) :: isochar     ! symbol of isomer
  character(len=1), dimension(-1:numpar)  :: parsym      ! symbol of particle
  character(len=2), dimension(numelem)    :: nuc         ! symbol of nucleus
  character(len=8), dimension(-1:numpar)  :: parname     ! name of particle
  character(len=4), dimension(numiso)     :: natstring   ! string extension for file names
  integer                                 :: iso         ! counter for isotope
  integer, dimension(8)                   :: magic       ! magic numbers
  integer, dimension(0:numpar)            :: parA        ! mass number of particle
  integer, dimension(0:numpar)            :: parN        ! neutron number of particle
  integer, dimension(0:numpar)            :: parZ        ! charge number of particle
  integer, dimension(0:numpar)            :: spin2       ! 2 * spin of particle
  real(dbl)                               :: amu         ! atomic mass unit in MeV
  real(sgl)                               :: amu4pi2h2c2 ! amu/(4*pi*pi*clight*clight*hbar**2) in mb**-1.MeV**-1
  real(sgl)                               :: amupi2h3c2  ! amu/(pi*pi*clight*clight*hbar**3) in mb**-1.MeV**-2.s**-1
  real(sgl)                               :: avogadro    ! Avogadro's number
  real(sgl)                               :: clight      ! speed of light in vacuum in m/s
  real(sgl)                               :: deg2rad     ! conversion factor for degrees to radians
  real(sgl)                               :: e2          ! square of elementary charge in MeV.fm
  real(sgl)                               :: emass       ! electron mass in MeV/c^2
  real(dbl), dimension(0:numpar)          :: excmass     ! mass excess of particle in a.m.u.
  real(sgl)                               :: fislim      ! mass above which nuclide fissions
  real(sgl)                               :: fourpi      ! 4.*pi
  real(sgl)                               :: Emaxtalys   ! maximum acceptable energy for TALYS
  real(sgl)                               :: hbar        ! Planck's constant / 2.pi in MeV.s
  real(sgl)                               :: hbarc       ! hbar.c in MeV.fm
  real(sgl)                               :: kT          ! energy kT expressed in MeV corresponding to a temperature T9=1
  real(sgl)                               :: onethird    ! 1/3
  real(sgl)                               :: pardis      ! parity distribution
  real(dbl), dimension(0:numpar)          :: parmass     ! mass of particle in a.m.u.
  real(sgl), dimension(0:numpar)          :: parspin     ! spin of particle
  real(sgl)                               :: pi          ! pi
  real(sgl)                               :: pi2         ! pi**2
  real(sgl)                               :: pi2h2c2     ! 1/(pi*pi*clight*clight*hbar**2) in mb**-1.MeV**-2
  real(sgl)                               :: pi2h3c2     ! 1/(pi*pi*clight*clight*hbar**3) in mb**-1.MeV**-3.s**-1
  real(sgl)                               :: qelem       ! elementary charge in C
  real(sgl)                               :: rad2deg     ! conversion factor for radians to degrees
  real(sgl), dimension(0:2*numl)          :: sgn         ! sign
  real(sgl)                               :: sqrttwopi   ! sqrt(2.*pi)
  real(sgl)                               :: twopi       ! 2*pi
  real(sgl)                               :: twopihbar   ! 2*pi/hbar
  real(sgl)                               :: twothird    ! 2/3
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for reading input lines
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=132), dimension(numlines) :: inline  ! input line
  integer                                 :: nlines  ! number of input lines
  integer                                 :: nlines0 ! number of input lines
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for main input
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical          :: flagnatural ! flag for calculation of natural element
  character(len=1) :: ptype0      ! type of incident particle
  character(len=2) :: Starget     ! symbol of target nucleus
  integer          :: Arp         ! A of residual product
  integer          :: Atarget     ! mass number of target nucleus
  integer          :: Atarget0    ! mass number of target nucleus
  integer          :: Ainit       ! mass number of initial compound nucleus
  integer          :: k0          ! index of incident particle
  integer          :: Ltarget     ! excited level of target
  integer          :: Ninit       ! neutron number of initial compound nucleus
  integer          :: Ntarget     ! neutron number of target nucleus
  integer          :: Zinit       ! charge number of initial compound nucleus
  integer          :: Zrp         ! Z of residual product
  integer          :: Ztarget     ! charge number of target nucleus
  integer          :: Ztarget0    ! charge number of target nucleus
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for files
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=10)  :: date     ! date
  character(len=132) :: bestpath ! alternative directory for best values
  character(len=132) :: nulldev  ! null device
  character(len=132) :: path     ! directory containing files to be read
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for best files
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                                         :: flagbest      ! flag to use best set of adjusted parameters
  logical                                         :: flagbestend   ! flag to put best set of parameters at end of input file
  logical                                         :: flagrescue    ! flag for final rescue: normalization to data
  character(len=132), dimension(nummt,-1:numisom) :: rescuefile    ! file with incident energy dependent adjustment factors
  real(sgl), dimension(nummt,-1:numisom)          :: grescue       ! global multiplication factor for incident energy dependence
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for basic reactions
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                                         :: flagffruns    ! flag to designate subsequent evap. of fission products
  logical                                         :: flagrpruns    ! flag to designate that run is for residual product
  logical                                         :: flagang       ! flag for output of angular distributions
  logical                                         :: flagastro     ! flag for calculation of astrophysics reaction rate
  logical                                         :: flagbasic     ! flag for output of basic information and results
  logical                                         :: flagoutall    ! flag for output of all data in main output file
  logical                                         :: flagchannels  ! flag for exclusive channel calculation
  logical                                         :: flagEchannel  ! flag for channel energy for emission spectrum
  logical                                         :: flagendf      ! flag for information for ENDF-6 file
  logical                                         :: flagendfdet   ! flag for detailed ENDF-6 information per channel
  logical                                         :: flagendfecis  ! flag for new ECIS calculation for ENDF-6 files
  logical                                         :: flaglabddx    ! flag for calculation of DDX in LAB system
  logical                                         :: flagmassdis   ! flag for calculation of fission fragment mass distribution
  logical                                         :: flagmicro     ! flag for completely microscopic Talys calculation
  logical                                         :: flagpartable  ! flag for output of model parameters on separate file
  logical                                         :: flagreaction  ! flag for calculation of nuclear reactions
  logical                                         :: flagfit       ! flag to use automatically fitted parameters
  logical                                         :: flagngfit     ! flag for using fitted (n,g) nuclear model parameters
  logical                                         :: flagnffit     ! flag for using fitted (n,f) nuclear model parameters
  logical                                         :: flagnnfit     ! flag for using fitted (n,n'), etc. nuclear model parameters
  logical                                         :: flagnafit     ! flag for using fitted (n,a) nuclear model parameters
  logical                                         :: flagndfit     ! flag for using fitted (n,d) nuclear model parameters
  logical                                         :: flagpnfit     ! flag for using fitted (p,n) nuclear model parameters
  logical                                         :: flaggnfit     ! flag for using fitted (g,n) nuclear model parameters
  logical                                         :: flagdnfit     ! flag for using fitted (d,n) nuclear model parameters
  logical                                         :: flaganfit     ! flag for using fitted (a,n) nuclear model parameters
  logical                                         :: flaggamgamfit ! flag for using fitted Gamma_gamma nuclear model parameters
  logical                                         :: flagmacsfit   ! flag for using fitted MACS nuclear model parameters
  logical                                         :: flagrecoil    ! flag for calculation of recoils
  logical                                         :: flagrecoilav  ! flag for average velocity in recoil calculation
  logical                                         :: flagrel       ! flag for relativistic kinematics
  logical                                         :: flagrpevap    ! flag for eva. of residual products at high incident energies
  character(len=132)                              :: ompenergyfile ! file with energies for OMP calculation (ENDF files only)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for basic parameters
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                               :: flagequi     ! flag to use equidistant excitation instead of logarithmic bins
  logical                               :: flagequispec ! flagequispec: flag to use equidistant bins for emission spectra
  character(len=1), dimension(0:numpar) :: outtype      ! type of outgoing particles
  character(len=132)                    :: source       ! source of data
  character(len=132)                    :: oformat      ! format of data
  character(len=132)                    :: user         ! user of data
  integer                               :: Lisoinp      ! user assignment of target isomer number
  real(sgl)                             :: eninclow     ! minimal incident energy for nuclear model calculation
  real(sgl)                             :: isomer       ! definition of isomer in seconds
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for input energies
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                                      :: flaginitpop ! flag for initial population distribution
  logical                                      :: flagpopMeV  ! flag to use initial population per MeV instead of histogram
  character(len=132)                           :: energyfile  ! file with energies for OMP calculation
  integer                                      :: nin         ! counter for incident energy
  integer                                      :: nin0        ! counter for incident energy
  integer                                      :: npopE       ! number of energies for population distribution
  integer                                      :: npopJ       ! number of spins for population distribution
  integer                                      :: npopP       ! number of parities for population distribution
  integer                                      :: Ninc        ! number of incident energies
  real(sgl), dimension(0:numpop)               :: EdistE      ! excitation energy of population distribution
  real(sgl), dimension(0:numen6+2)             :: eninc       ! incident energy in MeV
  real(sgl)                                    :: deninc      ! incident energy increment
  real(sgl)                                    :: enincF      ! final incident energy
  real(sgl)                                    :: enincmax    ! maximum incident energy
  real(sgl)                                    :: enincmin    ! minimum incident energy
  real(sgl)                                    :: Estop       ! incident energy above which TALYS stops
  real(sgl), dimension(0:numpop)               :: PdistE      ! population distribution, spin-independent
  real(sgl), dimension(0:numpop, 0:numJ, -1:1) :: PdistJP     ! population distribution per spin and parity
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for numerics
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer   :: maxenrec   ! number of recoil energies
  integer   :: maxN       ! maximal number of neutrons away from initial compound nucleus
  integer   :: maxNrp     ! maximal number of neutrons away from the initial compound nucleus
  integer   :: maxZ       ! maximal number of protons away from initial compound nucleus
  integer   :: maxZrp     ! maximal number of protons away from the initial compound nucleus
  integer   :: nangle     ! number of angles
  integer   :: nanglecont ! number of angles for continuum
  integer   :: nbins0     ! number of continuum excitation energy bins
  integer   :: segment    ! help array for storing segment intersection points
  integer   :: maxchannel ! maximal number of outgoing particles in individual channel description
  integer   :: nanglerec  ! number of recoil angles
  integer   :: transpower ! power for transmission coefficient limit
  real(sgl) :: popeps     ! limit for population cross sections
  real(dbl) :: transeps   ! absolute limit for transmission coefficient
  real(sgl) :: xseps      ! limit for cross sections
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for masses
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                                          :: flagexpmass ! flag for using experimental nuclear mass if available
  character(len=132)                               :: massdir     ! directory with mass tables
  integer                                          :: massmodel   ! model for theoretical nuclear mass
  real(sgl), dimension(0:numZ+4,0:numN+4,0:numbar) :: beta2       ! deformation parameter
  real(sgl), dimension(0:numZ+4,0:numN+4)          :: massexcess  ! mass excess in MeV as read from user input file
  real(sgl), dimension(0:numZ+4,0:numN+4)          :: massnucleus ! mass of nucleus in amu as read from user input file
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for discrete levels
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                               :: flagbestbr   ! flag to use only best set of branching ratios
  logical                               :: flagelectron ! flag for application of electron conversion coefficient
  logical                               :: flagpseudores! flag for using light nuclide discrete levels for resonances
  logical                               :: flaglevels   ! flag for output of discrete level information
  character(len=132), dimension(0:numZ) :: deformfile   ! deformation parameter file
  character(len=132), dimension(0:numZ) :: levelfile    ! discrete level file
  integer                               :: disctable    ! table with discrete levels
  integer, dimension(0:numZ, 0:numN)    :: nlev         ! number of levels for nucleus
  integer, dimension(0:numpar)          :: nlevbin      ! number of excited levels for binary nucleus
  integer                               :: nlevmax      ! maximum number of included discrete levels for target nucleus
  integer                               :: nlevmaxres   ! maximum number of included discrete levels for residual nucleus
  real(sgl), dimension(0:numresgrid)    :: Eresgrid     ! energies on resonance grid
  real(sgl), dimension(0:numresgrid)    :: resgrid      ! resonance grid
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for level densities
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=132), dimension(0:numZ,0:numN)   :: densfile          ! tabulated level density file
  logical                                        :: filedensity       ! flag for level densities on separate files
  logical                                        :: flagasys          ! flag for all level density parameters a from systematics
  logical, dimension(0:numZ, 0:numN)             :: flagcol           ! flag for collective enhancement of level density
  logical                                        :: flagcolall        ! flag for collective enhancement of level density
  logical                                        :: flagcolldamp      ! flag for damp of coll eff in effec lev dens (no expl col en)
  logical                                        :: flagctmglob       ! flag for global CTM model (no discrete level info)
  logical                                        :: flagdensity       ! flag for output of level densities
  logical                                        :: flagparity        ! flag for non-equal parity distribution
  logical, dimension(0:numZ, 0:numN)             :: ldadjust          ! logical for energy-dependent level density adjustment
  integer                                        :: kvibmodel         ! model for vibrational enhancement
  integer, dimension(0:numZ, 0:numN)             :: ldmodel           ! level density model
  integer                                        :: ldmodelCN         ! level density model for compound nucleus
  integer                                        :: ldmodelall        ! level density model for all nuclides
  integer, dimension(0:numZ, 0:numN, 0:numbar)   :: Nlow              ! lowest discrete level for temperature matching
  integer, dimension(0:numZ, 0:numN, 0:numbar)   :: Ntop              ! highest discrete level for temperature matching
  integer                                        :: shellmodel        ! model for shell correction energies
  integer                                        :: spincutmodel      ! model for spin cutoff factor for ground state
  real(sgl), dimension(0:numZ,0:numN)            :: aadjust           ! adjustable factor for level density parameter
  real(sgl), dimension(0:numZ,0:numN)            :: alev              ! level density parameter
  real(sgl), dimension(0:numZ,0:numN)            :: alimit            ! asymptotic level density parameter
  real(sgl)                                      :: alphaldall        ! variable for level density
  real(sgl)                                      :: betaldall         ! variable for level density
  real(sgl)                                      :: gammashell1all    ! variable for level density
  real(sgl), dimension(0:numZ,0:numN)            :: alphald           ! alpha-constant for asymptotic level density parameter
  real(sgl), dimension(0:numZ,0:numN)            :: betald            ! beta-constant for asymptotic level density parameter
  real(sgl), dimension(0:numZ,0:numN,0:numbar)   :: cfermi            ! width of Fermi distribution for damping
  real(sgl)                                      :: cglobal           ! global constant to adjust tabulated level densities
  real(sgl), dimension(0:numZ,0:numN,0:numbar)   :: ctable            ! constant to adjust tabulated level densities
  real(sgl), dimension(0:numZ,0:numN,0:numbar)   :: ctableadjust      ! adjustable correction to adjust tabulated level densities
  real(sgl), dimension(0:numZ,0:numN)            :: D0                ! s-wave resonance spacing in eV
  real(sgl), dimension(0:numZ,0:numN,0:numbar)   :: deltaW            ! shell correction in nuclear mass
  real(sgl), dimension(0:numZ,0:numN,0:numbar)   :: E0                ! particle constant of temperature formula
  real(sgl), dimension(0:numZ,0:numN,0:numbar)   :: E0adjust          ! adjustable factor for E0
  real(sgl), dimension(0:numZ,0:numN,0:numbar)   :: Exmatch           ! matching point for Ex
  real(sgl), dimension(0:numZ,0:numN,0:numbar)   :: Exmatchadjust     ! adjustable factor for matching energy
  real(sgl), dimension(0:numZ,0:numN)            :: gammald           ! gamma-constant for asymptotic level density parameter
  real(sgl), dimension(0:numZ,0:numN)            :: gammashell1       ! gamma-constant for asymptotic level density parameter
  real(sgl)                                      :: gammashell2       ! gamma-constant for asymptotic level density parameter
  real(sgl), dimension(0:numZ,0:numN,0:numbar)   :: Krotconstant      ! normalization constant for rotational enhancement
  real(sgl), dimension(0:numZ,0:numN)            :: pair              ! pairing energy
  real(sgl)                                      :: pairconstant      ! constant for pairing energy systematics
  real(sgl)                                      :: pglobal           ! global constant to adjust tabulated level densities
  real(sgl), dimension(0:numZ, 0:numN, 0:numbar) :: Pshift            ! adjustable pairing shift
  real(sgl), dimension(0:numZ, 0:numN, 0:numbar) :: Pshiftadjust      ! adjustable correction to pairing shift
  real(sgl), dimension(0:numZ,0:numN)            :: Pshiftconstant    ! global constant for pairing shift
  real(sgl)                                      :: Pshiftconstantall ! variable for level density
  real(sgl), dimension(0:numZ,0:numN,0:numbar)   :: ptable            ! constant to adjust tabulated level densities
  real(sgl), dimension(0:numZ,0:numN,0:numbar)   :: ptableadjust      ! adjustable correction to adjust tabulated level densities
  real(sgl), dimension(0:numZ,0:numN,0:numbar)   :: Rclass2mom        ! norm. constant for moment of inertia for class 2 states
  real(sgl), dimension(0:numZ,0:numN)            :: Risomer           ! adjustable correction to level branching ratios
  real(sgl)                                      :: Rspincut          ! adjustable constant (global) for spin cutoff factor
  real(sgl), dimension(0:numZ,0:numN,0:numbar)   :: Rtransmom         ! norm. constant for moment of inertia for transition states
  real(sgl), dimension(0:numZ,0:numN,0:numbar)   :: s2adjust          ! adjustable constant (Z,A,barrier-dependent) for spin cutoff
  real(sgl), dimension(0:numZ,0:numN,0:numbar)   :: T                 ! temperature
  real(sgl), dimension(0:numZ,0:numN,0:numbar)   :: Tadjust           ! adjustable factor for temperature
  real(sgl), dimension(0:numZ,0:numN,0:numbar)   :: Ufermi            ! energy of Fermi distribution for damping
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for gamma rays
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical, dimension(0:numZ, 0:numN)                      :: gamadjust      ! logical for energy-dependent gamma adjustment
  logical                                                 :: flagracap      ! flag for radiative capture model
  logical                                                 :: filepsf        ! flag for photon strength functions on separate files
  logical                                                 :: flaggamma      ! flag for output of gamma-ray information
  logical                                                 :: flaggnorm      ! flag to normalize PSF to average radiative width
  logical                                                 :: flagupbend     ! flag for low-energy upbend of photon strength function
  logical                                                 :: flagpsfglobal  ! flag for global photon strength functions only
  integer                                                 :: gammax         ! number of l-values for gamma multipolarity
  integer                                                 :: ldmodelracap   ! level density model for direct radiative capture
  integer                                                 :: strength       ! E1 strength function model
  integer                                                 :: strengthM1     ! model for M1 gamma-ray strength function
  character(len=132), dimension(0:numZ,0:numN,0:1,numgam) :: Exlfile        ! tabulated gamma ray strength function
  real(sgl), dimension(-1:numpar)                         :: fisominit      ! correction factor for isospin forbidden transition
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam,2)        :: egr            ! energy of GR
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam,2)        :: egradjust      ! adjustable factor for energy of GR
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam,2)        :: epr            ! energy of PR
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam,2)        :: epradjust      ! adjustable factor for energy of PR
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam)          :: etable         ! constant to adjust tabulated strength functions
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam)          :: etableadjust   ! adjustable correction to adjust tabulated strength fun
  real(sgl), dimension(-1:numpar)                         :: fiso           ! correction factor for isospin forbidden transition
  real(sgl), dimension(-1:numpar)                         :: fisom          ! correction factor for isospin forbidden transition
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam)          :: ftable         ! constant to adjust tabulated strength functions
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam)          :: ftableadjust   ! adjustable correction to adjust tabulated strength fun
  real(sgl), dimension(0:numZ,0:numN)                     :: gamgam         ! total radiative width in eV
  real(sgl), dimension(0:numZ,0:numN)                     :: gamgamadjust   ! adjustable factor for radiative parameter
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam,2)        :: ggr            ! width of GR
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam,2)        :: ggradjust      ! adjustable factor for width of GR
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam,2)        :: gpr            ! width of PR
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam,2)        :: gpradjust      ! adjustable factor for width of PR
  real(sgl)                                               :: sfexpall       ! variable for spectrocopic factor
  real(sgl)                                               :: sfthall        ! variable for spectrocopic factor
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam,2)        :: sgr            ! strength of GR
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam,2)        :: sgradjust      ! adjustable factor for strength of GR
  real(sgl), dimension(0:numZ, 0:numN, 0:numlev)          :: spectfacexp    ! experimental spectroscopic factor
  real(sgl), dimension(0:numZ, 0:numN)                    :: spectfacth     ! theoretical spectroscopic factor
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam,2)        :: tpr            ! strength of PR
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam,2)        :: tpradjust      ! adjustable factor for strength of PR
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam,3)        :: upbend         ! properties of the low-energy upbend of given multipola
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam,3)        :: upbendadjust   ! properties of the low-energy upbend of given multipola
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam)          :: wtable         ! constant to adjust tabulated strength functions
  real(sgl), dimension(0:numZ,0:numN,0:1,numgam)          :: wtableadjust   ! adjustable correction to adjust tabulated strength fun
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables to set OMP parameters
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical, dimension(numpar)                                :: altomp        ! flag for alternative optical model
  logical                                                   :: flagdisp      ! flag for dispersive optical model
  logical                                                   :: flagincadj    ! flag for OMP adjustment on incident channel
  logical                                                   :: flagjlm       ! flag for using semi-microscopic JLM OMP
  logical                                                   :: flaglocalomp  ! flag for local (y) or global (n) optical model
  logical                                                   :: flagompall    ! flag for new optical model calculation for all nuclei
  logical                                                   :: flagomponly   ! flag to execute ONLY an optical model calculation
  logical                                                   :: flagriplomp   ! flag for RIPL OMP
  logical, dimension(0:numpar)                              :: ompadjustF    ! logical for local OMP adjustment
  logical, dimension(0:numpar)                              :: ompadjustp    ! flag for local optical model parameter adjustment
  logical                                                   :: flagoutomp    ! flag for output of optical model parameters
  logical                                                   :: flagoutkd     ! flag for output of KD03 OMP parameters
  logical                                                   :: flagsoukhoinp ! flag for Soukhovitskii OMP for actinides
  logical                                                   :: flagsoukho    ! flag for Soukhovitskii OMP for actinides
  logical                                                   :: flagriplrisk  ! flag for going outside RIPL mass validity range
  character(len=1)                                          :: pruitt        ! identifier for using Pruitt parameters for KD03
  character(len=132), dimension(0:numZph, 0:numNph, numpar) :: optmod        ! file with optical model parameters
  character(len=132), dimension(0:numZ)                     :: optmodfileN   ! optical model parameter file for neutrons
  character(len=132), dimension(0:numZ)                     :: optmodfileP   ! optical model parameter file for protons
  character(len=132), dimension(0:numZ)                     :: radialfile    ! radial matter density file
  integer                                                   :: alphaomp      ! alpha optical model
  integer                                                   :: deuteronomp   ! deuteron optical model
  integer                                                   :: jlmmode       ! option for JLM imaginary potential normalization
  integer, dimension(numpar)                                :: riplomp       ! RIPL OMP
  integer, dimension(0:numpar, numompadj)                   :: ompadjustN    ! number of energy ranges for local OMP adjustment
  integer                                                   :: pruittset     ! random set for Pruitt et al OMP
  integer                                                   :: radialmodel   ! model for radial matter densities (JLM OMP only)
  real(sgl)                                                 :: adepthcor     ! adjustable parameter for depth of DF alpha potential
  real(sgl)                                                 :: aradialcor    ! adjustable parameter for shape of DF alpha potential
  real(sgl), dimension(0:numpar)                            :: avadjust      ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: avdadjust     ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: avsoadjust    ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: awadjust      ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: awdadjust     ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: awsoadjust    ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: d1adjust      ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: d2adjust      ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: d3adjust      ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: Ejoin         ! joining energy for high energy OMP
  real(sgl)                                                 :: lvadjust      ! adjustable parameter for JLM OMP
  real(sgl)                                                 :: lv1adjust     ! adjustable parameter for JLM OMP
  real(sgl)                                                 :: lvsoadjust    ! adjustable parameter for JLM OMP
  real(sgl)                                                 :: lwadjust      ! adjustable parameter for JLM OMP
  real(sgl)                                                 :: lw1adjust     ! adjustable parameter for JLM OMP
  real(sgl)                                                 :: lwsoadjust    ! adjustable parameter for JLM OMP
  real(sgl), dimension(0:numpar,numompadj,numrange)         :: ompadjustD    ! depth of local OMP adjustment
  real(sgl), dimension(0:numpar,numompadj,numrange)         :: ompadjustE1   ! start energy of local OMP adjustment
  real(sgl), dimension(0:numpar,numompadj,numrange)         :: ompadjustE2   ! end energy of local OMP adjustment
  real(sgl), dimension(0:numpar,numompadj,numrange)         :: ompadjusts    ! variance of local OMP adjustment
  real(sgl)                                                 :: RprimeU       ! potential scattering radius
  real(sgl), dimension(0:numpar)                            :: rcadjust      ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: rvadjust      ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: rvdadjust     ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: rvsoadjust    ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: rwadjust      ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: rwdadjust     ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: rwsoadjust    ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: v1adjust      ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: v2adjust      ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: v3adjust      ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: v4adjust      ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: Vinfadjust    ! adj. factor for high E limit of real central pot.
  real(sgl), dimension(0:numpar)                            :: vso1adjust    ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: vso2adjust    ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: w1adjust      ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: w2adjust      ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: w3adjust      ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: w4adjust      ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: wso1adjust    ! adjustable factor for OMP (default 1.)
  real(sgl), dimension(0:numpar)                            :: wso2adjust    ! adjustable factor for OMP (default 1.)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for compound reactions
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=132)                              :: reslib       ! library with resonance parameters
  logical, dimension(0:numZ, 0:numN, -1:numpar)   :: adjustTJ     ! logical for energy-dependent TJ adjustment
  logical                                         :: flagcomp     ! flag for compound nucleus calculation
  logical                                         :: flageciscomp ! flag for compound nucleus calculation by ECIS
  logical                                         :: flagfullhf   ! flag for full spin dependence of transmission coefficients
  logical                                         :: flaggroup    ! flag for output of low energy groupwise cross sections
  logical                                         :: flagres      ! flag for output of low energy resonance cross sections
  logical                                         :: flagurr      ! flag for output of unresolved resonance parameters
  logical                                         :: flagurrnjoy  ! normalization of URR parameters with NJOY method
  integer                                         :: lenreslib    ! length of library name with resonance parameters
  integer                                         :: lurr         ! maximal orbital angular momentum for URR
  integer                                         :: wmode        ! designator for width fluctuation model
  integer                                         :: WFCfactor    ! enhancement factor for WFC: 1: Original, 2: Ernebjerg and Herman
  real(sgl)                                       :: eurr         ! off-set incident energy for URR calculation
  real(sgl)                                       :: ewfc         ! off-set incident energy for width fluctuation correction
  real(sgl), dimension(0:numZ, 0:numN)            :: skipCN       ! flag to skip compound nucleus in evaporation chain
  real(sgl), dimension(0:numZ, 0:numN, -1:numpar) :: TJadjust     ! adjustable factor for TJ (default 1.)
  real(sgl)                                       :: Tres         ! temperature for broadening low energy cross sections
  real(sgl), dimension(-1:numisom)                :: xsalphatherm ! thermal (n,a) cross section
  real(sgl), dimension(-1:numisom)                :: xscaptherm   ! thermal capture cross section
  real(sgl), dimension(-1:numisom)                :: xsptherm     ! thermal (n,p) cross section
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for direct reactions
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                      :: flagautorot  ! flag for automatic rotational coupled channels calculation
  logical                      :: flagcoulomb  ! flag for Coulomb excitation calculation with ECIS
  logical                      :: flagcpang    ! flag for compound angular distribution calculation
  logical                      :: flagdirect   ! flag for output of direct reaction results
  logical                      :: flagdisc     ! flag for output of discrete state cross sections
  logical                      :: flageciscalc ! flag for new ECIS calculation for outgoing channels
  logical                      :: flagecissave ! flag for saving ECIS input and output files
  logical                      :: flaggiant0   ! flag for collective contribution from giant resonances
  logical                      :: flaginccalc  ! flag for new ECIS calculation for incident channel
  logical                      :: flaglegendre ! flag for output of Legendre coefficients
  logical                      :: flagoutecis  ! flag for output of ECIS results
  logical, dimension(0:numpar) :: flagrot      ! flag for use of rotational optical model per outgoing particle
  logical                      :: flagspher    ! flag to force spherical optical model
  logical                      :: flagstate    ! flag for optical model potential for each excited state
  logical, dimension(0:numpar) :: flagsys      ! flag for reaction cross section from systematics
  logical                      :: flagtransen  ! flag for output of transmission coefficients per energy
  integer                      :: maxband      ! highest vibrational band added to rotational model
  integer                      :: maxrot       ! number of included excited rotational levels
  integer                      :: core         ! even-even core for weakcoupling (-1 or 1)
  real(sgl)                    :: eadd         ! on-set incident energy for addition of discrete states
  real(sgl)                    :: eaddel       ! on-set incident energy for addition of elastic peak
  real(sgl)                    :: elwidth      ! width of elastic peak in MeV
  real(sgl)                    :: soswitch     ! switch for deformed spin-orbit calculation
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for fission
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical, dimension(0:numZ, 0:numN)           :: fisadjust       ! logical for energy-dependent fission adjustment
  logical                                      :: flagclass2      ! flag for class2 states in fission
  logical                                      :: flagffevap      ! flag for calculation of particle evaporation per FF
  logical                                      :: flagffspin      ! flag to use spin distribution in initial FF population
  logical                                      :: flagfisfeed     ! flag for output of fission per excitation energy bin
  logical                                      :: flagfisout      ! flag for output of fission information
  logical                                      :: flagfispartdamp ! flag for fission partial damping
  logical                                      :: flagfission     ! flag for fission
  logical                                      :: flaghbstate     ! flag for head band states in fission
  logical                                      :: flagoutfy       ! flag for output detailed fission yield calculation
  character(len=132), dimension(0:numZ,0:numN) :: clas2file       ! file with class 2 transition states
  character(len=132), dimension(0:numZ,0:numN) :: hbtransfile     ! file with head band transition states
  character(len=132)                           :: yieldfile       ! fission yield file
  integer, dimension(0:numZ, 0:numN, numbar)   :: axtype          ! type of axiality of barrier
  integer                                      :: fismodel        ! fission model alternative fission model for default barriers
  integer                                      :: fismodelalt     ! alternative fission model for default barriers
  integer, dimension(0:numZ,0:numN)            :: fismodelx       ! fission model
  integer                                      :: fymodel         ! fission yield model, 1: Brosa 2: GEF 3: GEF + TALYS
  integer                                      :: ffmodel         ! fission fragment model
  integer                                      :: pfnsmodel       ! PFNS  model, 1: Iwamoto 2: from FF decay
  integer                                      :: gefran          ! number of random events for GEF calculation
  real(sgl), dimension(0:numZ,0:numN)          :: betafiscor      ! adjustable factor for fission path width
  real(sgl), dimension(0:numZ,0:numN)          :: betafiscoradjust! adjustable factor for fission path width
  real(sgl), dimension(0:numZ,0:numN,numbar)   :: fbaradjust      ! adjustable factor for fission parameters
  real(sgl), dimension(0:numZ,0:numN,numbar)   :: fbarrier        ! height of fission barrier
  real(sgl), dimension(0:numZ,0:numN,numbar)   :: fwidth          ! width of fission barrier
  real(sgl), dimension(0:numZ,0:numN,numbar)   :: fwidthadjust    ! adjustable factor for fission parameters
  real(sgl), dimension(0:numZ,0:numN,numbar)   :: bdamp           ! fission partial damping parameter
  real(sgl), dimension(0:numZ,0:numN,numbar)   :: bdampadjust     ! correction for fission partial damping parameter
  real(sgl)                                    :: Rfiseps         ! ratio for limit for fission cross section per nucleus
  real(sgl)                                    :: Cnubar1         ! adjustable parameter for nubar constant value
  real(sgl)                                    :: Cnubar2         ! adjustable parameter for nubar energy slope
  real(sgl)                                    :: Tmadjust        ! adjustable parameter for PFNS temperature
  real(sgl)                                    :: Fsadjust        ! adjustable parameter for PFNS scission fraction
  real(sgl)                                    :: Cbarrier        ! global multiplication factor for fission barrier of Sierk model
  real(sgl)                                    :: Rspincutff      ! adjustable parameter (global) for FF spin cutoff factor
  real(sgl), dimension(0:numZ,0:numN)          :: vfiscor         ! adjustable factor for fission path height
  real(sgl), dimension(0:numZ,0:numN)          :: vfiscoradjust   ! adjustable factor for fission path height
  real(sgl), dimension(0:numZ,0:numN,numbar)   :: widthc2         ! width of class2 states
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for preequilibrium
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                             :: flag2comp    ! flag for two-component pre-equilibrium model
  logical                             :: flagecisdwba ! flag for new ECIS calculation for DWBA for MSD
  logical                             :: flagoutdwba  ! flag for output of DWBA cross sections for MSD
  logical                             :: flaggshell   ! flag for energy dependence of p.h. level density parameter
  logical                             :: flagonestep  ! flag for continuum one-step direct only
  logical                             :: flagpecomp   ! flag for Kalbach complex particle emission model
  logical                             :: flagpeout    ! flag for output of pre-equilibrium results
  logical                             :: flagsurface  ! flag for surface effects in exciton model
  logical                             :: preeqadjust  ! logical for energy-dependent pre-eq adjustment
  logical                             :: preeqfirst   ! logical for first preequilibrium output
  integer                             :: breakupmodel ! model for break-up reaction: 1. Kalbach 2. Avrigeanu
  integer                             :: mpreeqmode   ! designator for multiple pre-equilibrium model
  integer                             :: msdbins      ! number of energy points for DWBA calculation for MSD
  integer                             :: pairmodel    ! model for preequilibrium pairing energy
  integer                             :: pespinmodel  ! model for pre-equilibrium or compound spin
  integer                             :: phmodel      ! particle-hole state density model
  integer                             :: preeqmode    ! designator for pre-equilibrium model
  real(sgl), dimension(0:numpar)      :: Cbreak       ! adjustable parameter for break-up reactions
  real(sgl), dimension(0:numpar)      :: Cknock       ! adjustable parameter for knockout reactions
  real(sgl), dimension(0:numpar)      :: Cstrip       ! adjustable parameter for stripping/pick-up reactions
  real(sgl)                           :: Emsdmin      ! minimal outgoing energy for MSD calculation
  real(sgl)                           :: emulpre      ! on-set incident energy for multiple preequilibrium
  real(sgl)                           :: epreeq       ! on-set incident energy for preequilibrium
  real(sgl)                           :: Esurf0       ! well depth for surface interaction
  real(sgl), dimension(0:numZ,0:numN) :: g            ! single-particle level density parameter
  real(sgl), dimension(0:numZ,0:numN) :: gadjust      ! adjustable factor for single - particle particle-hole states
  real(sgl), dimension(0:numZ,0:numN) :: gn           ! single-particle neutron level density parameter
  real(sgl), dimension(0:numZ,0:numN) :: gnadjust     ! adjustable factor for single-particle proton parameter
  real(sgl), dimension(0:numZ,0:numN) :: gp           ! single-particle proton level density parameter
  real(sgl), dimension(0:numZ,0:numN) :: gpadjust     ! adjustable factor for single-particle neutron parameter
  real(sgl)                           :: Kph          ! constant for single-particle level density parameter
  real(sgl)                           :: M2constant   ! constant for matrix element in exciton model
  real(sgl)                           :: M2limit      ! constant for asymptotic value for matrix element
  real(sgl)                           :: M2shift      ! constant for energy shift for matrix element
  real(sgl)                           :: Rgamma       ! adjustable parameter for pre-equilibrium gamma decay
  real(sgl)                           :: Rnunu        ! ratio for two-component matrix element
  real(sgl)                           :: Rnupi        ! ratio for two-component matrix element
  real(sgl)                           :: Rpinu        ! ratio for two-component matrix element
  real(sgl)                           :: Rpipi        ! ratio for two-component matrix element
  real(sgl)                           :: Rspincutpreeq ! adjustable constant (global) for preequilibrium spin cutoff factor
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for astrophysics
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical   :: flagastroex ! flag for calculation of astrophysics reaction with g.s. target
  logical   :: flagastrogs ! flag for calculation of astrophysics reaction to final states
  integer   :: nonthermlev ! non-thermalized level in the calculation of astrophysical rate
  integer   :: nTmax       ! effective number of temperatures for Maxwellian average
  real(sgl) :: astroE      ! energy, in MeV, for Maxwellian average
  real(sgl) :: astroT9     ! temperature, in 10^9 K, for Maxwellian average
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for medical isotope production
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                        :: flagprod   ! flag for isotope production
  character(len=1), dimension(5) :: unitTcool  ! cooling time unit (y,d,h,m,s)
  character(len=1), dimension(5) :: unitTirrad ! irradiation time unit (y,d,h,m,s)
  character(len=3)               :: radiounit  ! unit for radioactivity: Bq, kBq, MBq, Gbq, mCi,
  character(len=3)               :: yieldunit  ! unit for isotope yield: num (number), mug, mg, g, or kg
  integer, dimension(5)          :: Tcool      ! cooling time per unit cooling time unit (y,d,h,m,s)
  integer, dimension(5)          :: Tirrad     ! irradiation time per unit irradiation time unit
  real(sgl)                      :: Area       ! target area in cm^2
  real(sgl)                      :: Eback      ! lower end of energy range in MeV for isotope
  real(sgl)                      :: Ebeam      ! incident energy in MeV for isotope production
  real(sgl)                      :: Ibeam      ! beam current in mA for isotope production
  real(sgl)                      :: rhotarget  ! target material density
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for output
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical, dimension(0:numlev)            :: fileangle    ! designator for angular distributions on separate file
  logical                                 :: filechannels ! flag for exclusive channel cross sections on separate file
  logical, dimension(0:numlev)            :: filediscrete ! flag for discrete level cross sections on separate file
  logical                                 :: fileelastic  ! flag for elastic angular distribution on separate file
  logical                                 :: filefission  ! flag for fission cross sections on separate file
  logical                                 :: filegamdis   ! flag for gamma-ray intensities on separate file
  logical                                 :: filerecoil   ! flag for recoil spectra on separate file
  logical                                 :: fileresidual ! flag for residual production cross sections on separate file
  logical, dimension(0:numpar)            :: filespectrum ! designator for spectrum on separate file
  logical                                 :: filetotal    ! flag for total cross sections on separate file
  logical                                 :: flagbinspec  ! flag for output of emission spectrum per excitation bin
  logical                                 :: flagblock    ! flag to block spectra, angle and gamma files
  logical                                 :: flagcheck    ! flag for output of numerical checks
  logical                                 :: flagcompo    ! flag for output of cross section components
  logical                                 :: flagddx      ! flag for output of double-differential cross sections
  logical                                 :: flagdecay    ! flag for output of decay of each population bin
  logical                                 :: flagexc      ! flag for output of excitation functions
  logical                                 :: flaggamdis   ! flag for output of discrete gamma-ray intensities
  logical                                 :: flagintegral ! flag for calc. of effective cross section using integral data
  logical                                 :: flaginverse  ! flag for output of transmission coeff. and inverse cross sections
  logical                                 :: flagmain     ! flag for main output
  logical                                 :: flagpop      ! flag for output of population
  logical                                 :: flagsacs     ! flag for statistical analysis of cross sections
  logical                                 :: flagspec     ! flag for output of spectra
  character(len=132), dimension(numflux)  :: fluxname     ! name of integral spectrum
  character(len=132), dimension(numflux)  :: xsfluxfile   ! TALYS cross section file for integral data
  integer, dimension(0:numpar)            :: ddxacount    ! counter for double-differential cross section files
  integer, dimension(0:numpar)            :: ddxecount    ! counter for double-differential cross section files
  integer                                 :: ddxmode      ! mode for DDX: 0: None, 1: Ang. distr., 2: Spect. per angle, 3 both
  integer                                 :: Nflux        ! number of reactions with integral data
  real(sgl), dimension(0:numpar, numfile) :: fileddxa     ! designator for double-differential cross sections on separate fil
  real(sgl), dimension(0:numpar, numfile) :: fileddxe     ! designator for double-differential cross sections on separate fil
  real(sgl), dimension(numflux)           :: integralexp  ! experimental effective cross section
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for adjustment
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=132), dimension(numadj) :: adjustfile ! file for local adjustment
  character(len=132), dimension(numadj) :: adjustkey  ! keyword for local adjustment
  integer, dimension(numadj, 4)         :: adjustix   ! local adjustment index
  integer                               :: Nadjust    ! number of adjustable parameters
  integer, dimension(numadj)            :: nenadjust  ! number of tabulated energies of local adjustment
  real(sgl), dimension(numadj,4)        :: adjustpar  ! local adjustment parameters
  real(sgl), dimension(numadj,numenadj) :: Dadjust    ! tabulated depth of local adjustment
  real(sgl), dimension(numadj,numenadj) :: Eadjust    ! tabulated energy of local adjustment
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for existence arrays
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical, dimension(0:numpar,0:numpar,0:numlev)                      :: angexist     ! logical to determine existence of ang dis
  logical, dimension(0:numin,0:numip,0:numid,0:numit,0:numih,0:numia) :: chanexist    ! flag for existence of exclusive cross sect
  logical, dimension(0:numin,0:numip,0:numid,0:numit,0:numih,0:numia) :: chanfisexist ! flag for existence of exclusive fission cr
  logical, dimension(0:numin,0:numip,0:numid,0:numit,0:numih,0:numia,0:numlev) :: chanisoexist ! flag for existence of exclusive iso
  logical, dimension(0:numin,0:numip,0:numid,0:numit,0:numih,0:numia) :: chanopen     ! flag to open channel with first non-zero c
  logical, dimension(0:numin,0:numip,0:numid,0:numit,0:numih,0:numia) :: gamchanexist ! flag for existence of exclusive discrete g
  logical, dimension(0:numin,0:numip,0:numid,0:numit,0:numih,0:numia) :: recchanexist ! flag for existence of exclusive recoil spe
  logical, dimension(0:numin,0:numip,0:numid,0:numit,0:numih,0:numia) :: spfischanexist ! flag for existence of fission spectra
  logical, dimension(0:numin,0:numip,0:numid,0:numit,0:numih,0:numia) :: spchanexist  ! flag for existence of exclusive spectra
  logical, dimension(0:numpar)                        :: ddxexist1    ! flag for existence of DDX
  logical, dimension(0:numpar)                        :: ddxexist2    ! flag for existence of DDX
  logical, dimension(0:numpar)                        :: ddxexist3    ! flag for existence of DDX
  logical, dimension(0:numpar)                        :: ddxexist4    ! flag for existence of DDX
  logical                                             :: idnumfull    ! flag to designate maximum number of exclusive ch.
  logical, dimension(-1:numZ,-1:numN,-1:numisom)      :: prodexist    ! logical to determine existence of residual pr
  logical, dimension(0:numZ, 0:numN, -1:numisom)      :: Yexist       ! flag for existence of yield
  logical, dimension(nummass)                         :: fpaexist     ! flag for existence of fission product per mass unit
  logical, dimension(2, numelem, numneu, -1:1)        :: fpexist      ! flag for existence of fission product
  logical, dimension(0:numpar)                        :: nubarexist   ! flag for existence of nubar file
  logical, dimension(0:numZ, 0:numN)                  :: fisexist     ! flag for existence of fission cross section
  logical, dimension(0:numZ, 0:numN)                  :: tfisexist    ! flag for existence of fission transmission coefficients
  logical, dimension(0:numZ,0:numN,0:numlev,0:numlev) :: gamexist     ! flag for existence of gamma production cross section
  logical, dimension(0:numpar,0:numpar,0:numlev)      :: legexist     ! logical to determine existence of Legendre coefficients
  logical, dimension(0:numZ, 0:numN)                  :: recexist     ! flag for existence of recoils
  logical, dimension(0:numZ, 0:numN)                  :: rpexist      ! flag for existence of residual production cross section
  logical, dimension(0:numZ, 0:numN, 0:numlev)        :: rpisoexist   ! flag for existence of isomeric residual production cross se
  logical, dimension(0:numpar)                        :: spexist1     ! flag for existence of spectra
  logical, dimension(0:numpar)                        :: spexist2     ! flag for existence of spectra
  logical, dimension(-1:11,0:numl)                    :: urrexist     ! flag for existence of URR
  integer                                             :: opennum      ! total number of open channels
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for abundance
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer                      :: isonum  ! number of isotopes in element
  integer, dimension(numiso)   :: isotope ! isotope of natural element
  real(sgl), dimension(numiso) :: abun    ! natural abundance
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for nuclides
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical, dimension(-1:numpar)                :: parinclude ! logical to include outgoing particle
  logical, dimension(-1:numpar)                :: parskip    ! logical to skip outgoing particle
  logical, dimension(0:numZ, 0:numN)           :: invexist    ! logical to state necessity of new inverse cross section
  logical                                      :: primary     ! flag to designate primary (binary) reaction
  logical, dimension(0:numZ, 0:numN)           :: strucexist  ! flag to state whether structure info for this nucleus exists
  logical, dimension(0:numZ, 0:numN)           :: strucwrite  ! flag for output of nuclear structure info
  character(len=6)                             :: targetnuclide ! target nuclide
  character(len=6)                             :: targetnuclide0 ! target nuclide
  integer, dimension(0:numZ, 0:numN, 0:numpar) :: AA          ! mass number of residual nucleus
  integer, dimension(0:numZ, 0:numN, 0:numpar) :: Nindex      ! neutron number index for residual nucleus
  integer, dimension(0:numZ, 0:numN, 0:numpar) :: NN          ! neutron number of residual nucleus
  integer, dimension(0:numZ, 0:numN, 0:numpar) :: Zindex      ! charge number index for residual nucleus
  integer, dimension(0:numZ, 0:numN, 0:numpar) :: ZZ          ! charge number of residual nucleus
  integer                                      :: targetP     ! parity of target
  integer                                      :: targetspin2 ! 2 * spin of target
  real(sgl), dimension(0:numpar)               :: coulbar     ! Coulomb barrier
  real(sgl), dimension(0:numpar)               :: Q           ! Q-value
  real(sgl), dimension(numpar)                 :: Smyers      ! Myers-Swiatecki separation energy
  real(sgl), dimension(numT)                   :: T9          ! Temperature grid in 10**9 K
  real(sgl)                                    :: targetE     ! excitation energy of target
  real(dbl)                                    :: tarmass     ! mass of target nucleus
  real(sgl)                                    :: targetspin  ! spin of target
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for masses
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, dimension(0:numZ+4,0:numN+4)        :: gsparity ! ground state parity
  real(sgl), dimension(0:numZ+4,0:numN+4)      :: beta4    ! deformation parameters
  real(dbl), dimension(0:numZ+4,0:numN+4)      :: dumexc   ! theoretical mass excess from Duflo-Zuker formula
  real(dbl), dimension(0:numZ+4,0:numN+4)      :: expmass  ! flag for using experimental nuclear mass if available
  real(dbl), dimension(0:numZ+4,0:numN+4)      :: expmexc  ! experimental mass excess
  real(sgl), dimension(0:numZ+4,0:numN+4)      :: gsspin   ! ground state spin
  real(dbl), dimension(0:numZ+4,0:numN+4)      :: nucmass  ! mass of nucleus
  real(dbl), dimension(0:numZ,0:numN,0:numpar) :: redumass ! reduced mass
  real(dbl), dimension(0:numZ,0:numN,0:numpar) :: specmass ! specific mass for residual nucleus
  real(dbl), dimension(0:numZ+4,0:numN+4)      :: thmass   ! theoretical mass
  real(dbl), dimension(0:numZ+4,0:numN+4)      :: thmexc   ! theoretical mass excess
  real(dbl), dimension(0:numZ,0:numN,0:numpar) :: S ! separation energy
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for levels
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=1), dimension(0:numZ,0:numN,0:numlev,0:numlev) :: bassign     ! flag for assignment of branching ratio
  character(len=18), dimension(0:numZ,0:numN,0:numlev)         :: ENSDF       ! string from original ENSDF discrete level file
  character(len=1), dimension(0:numZ,0:numN,0:numlev2)         :: eassign     ! flag for assignment of energy
  character(len=1), dimension(0:numZ,0:numN,0:numlev2)         :: jassign     ! flag for assignment of spin
  character(len=1), dimension(0:numZ,0:numN,0:numlev2)         :: passign     ! flag for assignment of parity
  integer, dimension(0:numZ,0:numN,0:numlev,0:numlev)          :: branchlevel ! level to which branching takes place
  integer                                                      :: Liso        ! isomeric number of target
  integer, dimension(-1:numZ,-1:numN,0:numisom)                :: Lisomer     ! level number of isomer
  integer                                                      :: Ltarget0    ! excited level of target
  integer, dimension(0:numZ,0:numN,0:numlev)                   :: nbranch     ! number of branching levels
  integer, dimension(-1:numZ,-1:numN)                          :: Nisomer     ! number of isomers for this nuclide
  integer, dimension(0:numZ,0:numN)                            :: nlevmax2    ! maximum number of levels
  integer, dimension(0:numZ,0:numN)                            :: branchdone  ! flag for applying branching ratio normalization
  integer, dimension(0:numZ,0:numN,0:numlev2)                  :: levnum      ! number of level
  integer, dimension(0:numZ,0:numN,0:numlev2)                  :: parlev      ! parity of level
  real(sgl), dimension(0:numZ,0:numN,0:numlev,0:numlev)        :: branchratio ! gamma-ray branching ratio to level
  real(sgl), dimension(0:numZ,0:numN,0:numlev,0:numlev)        :: conv        ! conversion coefficient
  real(sgl), dimension(0:numZ,0:numN,0:numlev2)                :: edis        ! energy of level
  real(sgl), dimension(0:numZ,0:numN,0:numlev2)                :: jdis        ! spin of level
  real(sgl), dimension(0:numZ,0:numN,0:numlev2)                :: tau         ! lifetime of state in seconds
  real(sgl), dimension(0:numZ,0:numN,0:numlev2)                :: tauripl     ! lifetime of state in seconds from RIPL
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for deformation parameters
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=1), dimension(0:numZ,0:numN)           :: colltype   ! type of collectivity (D, V or R)
  character(len=1), dimension(0:numZ,0:numN)           :: deftype    ! deformation length (D) or parameter (B)
  character(len=1), dimension(0:numZ,0:numN,0:numlev2) :: leveltype  ! type of level (rotational (R) or vibrational (V)
  integer, dimension(0:numZ, 0:numN, 0:numlev2)        :: indexlevel ! level index
  integer, dimension(0:numZ, 0:numN, 0:numlev2)        :: indexcc    ! level index for coupled channel
  integer, dimension(0:numZ, 0:numN, 0:numlev2)        :: iphonon    ! phonon (1 or 2)
  integer, dimension(0:numZ, 0:numN, 0:numlev2)        :: Kband      ! magnetic quantum number
  integer, dimension(0:numZ, 0:numN, 0:numlev2)        :: lband      ! angular momentum
  integer, dimension(0:numZ, 0:numN)                   :: ndef       ! number of collective levels
  integer, dimension(0:numZ, 0:numN)                   :: nrot       ! number of deformation parameters for rotational nucleus
  integer, dimension(0:numZ, 0:numN, 0:numlev2)        :: pcore      ! parity of level of core nucleus
  integer, dimension(0:numZ, 0:numN, 0:numlev2)        :: vibband    ! band number of level
  real(sgl), dimension(0:numZ, 0:numN, 0:numlev2)      :: jcore      ! spin of level of core nucleus
  real(sgl), dimension(0:numZ, 0:numN, 0:numlev2)      :: deform     ! deformation parameter
  real(sgl), dimension(0:numZ, 0:numN, 0:numlev2)      :: defpar     ! deformation parameter
  real(sgl), dimension(0:numZ, 0:numN, 0:numbar)       :: Irigid     ! rigid body value of moment of inertia
  real(sgl), dimension(0:numZ, 0:numN)                 :: Irigid0    ! undeformed rigid body value of moment of inertia
  real(sgl), dimension(0:numZ, 0:numN, 0:numrotcc)     :: rotpar     ! deformation parameters for rotational nucleus
  real(sgl), dimension(0:3,2)                          :: betagr     ! deformation parameter for giant resonance
  real(sgl), dimension(0:3,2)                          :: Egrcoll    ! energy of giant resonance
  real(sgl), dimension(0:3,2)                          :: Ggrcoll    ! width of giant resonance
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for resonance parameters
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, dimension(0:numZ,0:numN)            :: Nrr      ! number of resonances
  real(sgl), dimension(0:numZ,0:numN)          :: dD0      ! uncertainty in D0
  real(sgl), dimension(0:numZ,0:numN)          :: dgamgam  ! uncertainty in gamgam
  real(sgl)                                    :: Eavres   ! average resonance energy
  real(sgl), dimension(0:numZ, 0:numN)         :: D0theo   ! mean s-wave resonance spacing
  real(sgl), dimension(0:numZ, 0:numN)         :: D0global ! global mean s-wave resonance spacing
  real(sgl), dimension(0:numZ, 0:numN)         :: dD0global! uncertainty in global mean s-wave resonance spacing
  real(sgl), dimension(0:numZ, 0:numN)         :: D1theo   ! mean p-wave resonance spacing
  real(sgl), dimension(0:numl)                 :: Dl       ! mean resonance spacing per l value
  real(sgl), dimension(0:numl, 0:numJ)         :: Dlj      ! mean resonance spacing per J,l value
  real(sgl), dimension(0:numZ, 0:numN, 0:numl) :: gamgamth ! theoretical total radiative width
  real(sgl), dimension(0:numZ, 0:numN)         :: swaveth  ! theoretical strength function for s-wave
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for gamma-ray strength functions
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical, dimension(0:numZ,0:numN,0:1,numgam)                         :: qrpaexist  ! flag for existence of tabulated QRPA PSF
  integer, dimension(0:numpar, 0:numen)                                :: lmax       ! maximal l-value for transmission coefficients
  integer, dimension(0:numZ,0:numN,0:1,numgam)                         :: ngr        ! number of GR
  integer                                                              :: nTqrpa     ! number of temperatures for QRPA
  real(sgl), dimension(0:numZ,0:numN,0:numgamqrpa,0:1,numgam)          :: eqrpa      ! energy grid for QRPA strength function
  real(sgl), dimension(0:numZ,0:numN,0:numgamqrpa,numTqrpa,0:1,numgam) :: fqrpa      ! tabulated QRPA strength functi
  real(sgl), dimension(250)                                            :: gamkopecky ! rad. width in eV by spline fit of Kopecky
  real(sgl), dimension(numgam)                                         :: kgr        ! constant for gamma-ray strength function
  real(sgl), dimension(numTqrpa)                                       :: Tqrpa      ! temperature for QRPA
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for OMP
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl) :: Fav   ! adjustable factor for OMP (default 1.)
  real(sgl) :: Favd  ! adjustable factor for OMP (default 1.)
  real(sgl) :: Favso ! adjustable factor for OMP (default 1.)
  real(sgl) :: Faw   ! adjustable factor for OMP (default 1.)
  real(sgl) :: Fawd  ! adjustable factor for OMP (default 1.)
  real(sgl) :: Fawso ! adjustable factor for OMP (default 1.)
  real(sgl) :: Fd1   ! adjustable factor for OMP (default 1.)
  real(sgl) :: Fd2   ! adjustable factor for OMP (default 1.)
  real(sgl) :: Fd3   ! adjustable factor for OMP (default 1.)
  real(sgl) :: Frc   ! adjustable factor for OMP (default 1.)
  real(sgl) :: Frv   ! adjustable factor for OMP (default 1.)
  real(sgl) :: Frvd  ! adjustable factor for OMP (default 1.)
  real(sgl) :: Frvso ! adjustable factor for OMP (default 1.)
  real(sgl) :: Frw   ! adjustable factor for OMP (default 1.)
  real(sgl) :: Frwd  ! adjustable factor for OMP (default 1.)
  real(sgl) :: Frwso ! adjustable factor for OMP (default 1.)
  real(sgl) :: Fv1   ! adjustable factor for OMP (default 1.)
  real(sgl) :: Fv2   ! adjustable factor for OMP (default 1.)
  real(sgl) :: Fv3   ! adjustable factor for OMP (default 1.)
  real(sgl) :: Fvso1 ! adjustable factor for OMP (default 1.)
  real(sgl) :: Fvso2 ! adjustable factor for OMP (default 1.)
  real(sgl) :: Fv4   ! adjustable factor for OMP (default 1.)
  real(sgl) :: Fw1   ! adjustable factor for OMP (default 1.)
  real(sgl) :: Fw2   ! adjustable factor for OMP (default 1.)
  real(sgl) :: Fw3   ! adjustable factor for OMP (default 1.)
  real(sgl) :: Fw4   ! adjustable factor for OMP (default 1.)
  real(sgl) :: Fwso1 ! adjustable factor for OMP (default 1.)
  real(sgl) :: Fwso2 ! adjustable factor for OMP (default 1.)
!
! opticaln
!
  real(sgl) :: av   ! real volume diffuseness
  real(sgl) :: avd  ! real surface diffuseness
  real(sgl) :: avso ! real spin-orbit diffuseness
  real(sgl) :: aw   ! imaginary volume diffuseness
  real(sgl) :: awd  ! imaginary surface diffuseness
  real(sgl) :: awso ! imaginary spin-orbit diffuseness
  real(sgl) :: rc   ! Coulomb radius
  real(sgl) :: rv   ! real volume radius
  real(sgl) :: rvd  ! real surface radius
  real(sgl) :: rw   ! imaginary volume radius
  real(sgl) :: rvso ! real spin-orbit radius
  real(sgl) :: rwd  ! imaginary surface radius
  real(sgl) :: rwso ! imaginary spin-orbit radius
  real(sgl) :: v    ! real volume depth
  real(sgl) :: vd   ! real surface depth
  real(sgl) :: w    ! imaginary volume depth
  real(sgl) :: vso  ! real spin-orbit depth
  real(sgl) :: wd   ! imaginary surface depth
  real(sgl) :: wso  ! imaginary spin-orbit depth
!
! omppar
!
  logical                                                    :: flagompejec! flag for OMP for ejectile equal to projectile
  logical, dimension(0:numZ,0:numN,numpar)                   :: disp      ! flag for dispersive optical model
  logical, dimension(0:numZ,0:numN,numpar)                   :: ompglobal ! flag for use of global optical model
  integer, dimension(0:numZ,0:numN,numpar)                   :: omplines  ! number of lines in optical model file
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: av0       ! diffuseness for real volume OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: avd0      ! diffuseness for surface OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: avso0     ! diffuseness for real spin-orbit OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: d1        ! parameter for imaginary surface OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: d2        ! parameter for imaginary surface OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: d3        ! parameter for imaginary surface OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: ef        ! Fermi energy
  real(sgl), dimension(0:numZph,0:numNph,numpar,0:numomp)    :: eomp      ! energies on optical model file
  real(sgl), dimension(numpar,10)                            :: Eompbeg0  ! upper energy of KD03 OMP
  real(sgl), dimension(numpar,10)                            :: Eompbeg1  ! lower energy of alternative OMP
  real(sgl), dimension(numpar,10)                            :: Eompend0  ! lower energy of KD03 OMP
  real(sgl), dimension(numpar,10)                            :: Eompend1  ! upper energy of alternative
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: rc0       ! Coulomb radius
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: rv0       ! radius for real volume OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: rvd0      ! radius for surface OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: rvso0     ! radius for real spin-orbit OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: v1        ! parameter for real volume OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: v2        ! parameter for real volume OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: v3        ! parameter for real volume OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: v4        ! parameter for real volume OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: vso1      ! parameter for real spin-orbit OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: vso2      ! parameter for real spin-orbit OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: w1        ! parameter for imaginary volume OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: w2        ! parameter for imaginary volume OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: w3        ! parameter for imaginary volume OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: w4        ! parameter for imaginary volume OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: wso1      ! parameter for imaginary spin-orbit OMP
  real(sgl), dimension(0:numZ,0:numN,numpar)                 :: wso2      ! parameter for imaginary spin-orbit OMP
  real(sgl), dimension(2)                                    :: V0        ! V at zero MeV
  real(sgl), dimension(2)                                    :: Vjoin     ! V at joining energy
  real(sgl), dimension(0:numZph,0:numNph,numpar,0:numomp,19) :: vomp      ! optical model parameters from file
  real(sgl), dimension(2)                                    :: Wjoin     ! W at joining energy
!
! KD03 global parameters
!
  real(sgl)   :: rv_0      !
  real(sgl)   :: rv_A      !
  real(sgl)   :: av_0      !
  real(sgl)   :: av_A      !
  real(sgl)   :: v1_0      !
  real(sgl)   :: v1_asymm  !
  real(sgl)   :: v1_A      !
  real(sgl)   :: v2_0      !
  real(sgl)   :: v2_A      !
  real(sgl)   :: v3_0      !
  real(sgl)   :: v3_A      !
  real(sgl)   :: v4_0      !
  real(sgl)   :: w1_0      !
  real(sgl)   :: w1_A      !
  real(sgl)   :: w2_0      !
  real(sgl)   :: w2_A      !
  real(sgl)   :: rd_0      !
  real(sgl)   :: rd_A      !
  real(sgl)   :: ad_0      !
  real(sgl)   :: ad_A      !
  real(sgl)   :: d1_0      !
  real(sgl)   :: d1_asymm  !
  real(sgl)   :: d2_0      !
  real(sgl)   :: d2_A      !
  real(sgl)   :: d2_A2     !
  real(sgl)   :: d2_A3     !
  real(sgl)   :: d3_0      !
  real(sgl)   :: rso_0     !
  real(sgl)   :: rso_A     !
  real(sgl)   :: aso_0     !
  real(sgl)   :: vso1_0    !
  real(sgl)   :: vso1_A    !
  real(sgl)   :: vso2_0    !
  real(sgl)   :: wso1_0    !
  real(sgl)   :: wso2_0    !
  real(sgl)   :: rc_0      !
  real(sgl)   :: rc_A      !
  real(sgl)   :: rc_A2     !
!
! spr
!
  real(sgl)                    :: Rprime    ! potential scattering radius
  real(sgl), dimension(0:numl) :: Sstrength ! s,p,d,etc-wave strength function
!
! JLM
!
  logical, dimension(0:numZ, 0:numN, numpar) :: jlmexist ! flag for existence of tabulated radial matter density
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for JLM
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(0:numZ, 0:numN, 6)         :: normjlm  ! JLM potential normalization factors
  real(sgl), dimension(0:numZ, 0:numN, numjlm, 6) :: potjlm   ! JLM potential depth values
  real(sgl), dimension(0:numZ, 0:numN, numjlm)    :: radjlm   ! radial points for JLM potential
  real(sgl), dimension(0:numZ, 0:numN, numjlm, 6) :: rhojlmn  ! density for neutrons
  real(sgl), dimension(0:numZ, 0:numN, numjlm, 6) :: rhojlmp  ! density for protons
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for ECIS
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                              :: flagecisinp ! flag for existence of ecis input file
  logical                              :: flaginvecis ! logical for calculating inverse channel OMP
  logical                              :: legendre    ! logical for output of Legendre coefficients
  character(len=1)                     :: tarparity   ! parity of target nucleus
  character(len=1), dimension(numlev2) :: Plevel      ! parity of level
  character(len=50)                    :: ecis1       ! 50 input flags ('T' or 'F') for ECIS
  character(len=50)                    :: ecis2       ! 50 input flags ('T' or 'F') for ECIS
  character(len=72)                    :: title       ! title of ECIS input file
  integer, dimension(numlev2)          :: iband       ! band number of level
  integer, dimension(numlev2)          :: idvib       ! identifier for existence of vibrational state inside rotational model
  integer, dimension(numlev2)          :: iph         ! help variable
  integer                              :: iqm         ! largest order of deformation
  integer                              :: iqmax       ! maximum l-value of multipole expansion
  integer                              :: iterm       ! number of iterations
  integer, dimension(numlev2)          :: Jband       ! angular momentum
  integer, dimension(numlev2)          :: Kmag        ! magnetic quantum number
  integer                              :: Nband       ! number of vibrational bands
  integer                              :: ncoll       ! number of nuclear states
  integer                              :: njmax       ! maximal number of j-values in ECIS
  integer                              :: npp         ! number of optical potentials
  integer                              :: nrad        ! number of radial points
  integer                              :: Nrotbeta    ! number of deformation parameters for rotational nucleus
  real(sgl)                            :: angbeg      ! first angle
  real(sgl)                            :: angend      ! last angle
  real(sgl)                            :: anginc      ! angle increment
  real(sgl)                            :: d2disp      ! constant for imaginary potential
  real(sgl)                            :: d3disp      ! constant for imaginary potential
  real(sgl)                            :: efer        ! Fermi energy
  real(sgl), dimension(numlev2)        :: Elevel      ! energy of level
  real(sgl)                            :: hint        ! integration step size h
  real(sgl), dimension(numlev2)        :: Jlevel      ! spin of level
  real(dbl)                            :: projmass    ! mass of projectile
  real(sgl)                            :: prodZ       ! product of charges of projectile and target nucleus
  real(dbl)                            :: resmass     ! mass of residual nucleus
  real(sgl)                            :: rmatch      ! matching radius
  real(sgl), dimension(numrot)         :: rotbeta     ! deformation parameters for rotational nucleus
  real(sgl)                            :: spin        ! spin of incident particle
  real(sgl)                            :: tarspin     ! spin of target nucleus
  real(sgl), dimension(numlev2)        :: vibbeta     ! vibrational deformation parameter
  real(sgl)                            :: w2disp      ! constant for imaginary potential
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for fission parameters
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(0:numZ,0:numN,numbar,0:numlev) :: efisc2hb ! energy of class2 states
  real(sgl), dimension(0:numZ,0:numN,numbar,0:numlev) :: efistrhb ! energy of head band transition states
  real(sgl), dimension(0:numZ, 0:numN, numbar)        :: fecont   ! start of continuum energy
  real(sgl), dimension(0:numZ,0:numN,numbar,0:numlev) :: jfisc2hb ! spin of class2 states
  real(sgl), dimension(0:numZ,0:numN,numbar,0:numlev) :: jfistrhb ! spin of head band transition states
  real(sgl), dimension(0:numZ, 0:numN, numbar)        :: minertia ! moment of inertia of fission barrier deformation
  real(sgl), dimension(0:numZ, 0:numN, numbar)        :: minertc2 ! moment of inertia for class2 states
  integer, dimension(0:numZ, 0:numN)                  :: nclass2  ! number of sets of class2 states
  integer, dimension(0:numZ, 0:numN, numbar)          :: nfisc2hb ! number of class2 states for barrier
  integer, dimension(0:numZ, 0:numN)                  :: nfisbar  ! number of fission barrier parameters
  integer, dimension(0:numZ, 0:numN, numbar)          :: nfistrhb ! number of head band transition states for barrier
  integer, dimension(0:numZ,0:numN,numbar,0:numlev)   :: pfisc2hb ! parity of class2 states
  integer, dimension(0:numZ,0:numN,numbar,0:numlev)   :: pfistrhb ! parity of head band transition states
!
! rotband
!
  real(sgl), dimension(0:numZ, 0:numN, numbar, 0:numrot) :: efistrrot ! energy of rotational transition states for barrier
  real(sgl), dimension(0:numZ, 0:numN, numbar, 0:numrot) :: jfistrrot ! spin of rotational transition states for barrier
  integer, dimension(0:numZ, 0:numN, numbar)             :: nfistrrot ! number of rotational transition states for barrier
  integer, dimension(0:numZ, 0:numN, numbar, 0:numrot)   :: pfistrrot ! parity of rotational transition states for barrier
!
! rotclass2
!
  real(sgl), dimension(0:numZ, 0:numN, numbar, 0:numrot) :: efisc2rot  ! energy of class2 rotational transition states
  real(sgl), dimension(0:numZ, 0:numN, numbar)           :: Emaxclass2 ! maximum energy for class2 states
  real(sgl), dimension(0:numZ, 0:numN, numbar, 0:numrot) :: jfisc2rot  ! spin of class2 rotational transition states
  integer, dimension(0:numZ, 0:numN, numbar)             :: nfisc2rot  ! number of class2 rotational transition states
  integer, dimension(0:numZ, 0:numN, numbar, 0:numrot)   :: pfisc2rot  ! parity of class2 rotational transition states
!
! fisdata
!
  real(sgl), dimension(7, 7)    :: barcof ! parameter values for barrier heights at l=0
  real(sgl), dimension(5, 4)    :: l20cof ! parameter values for l20 belonging to 20% barrier height
  real(sgl), dimension(5, 4)    :: l80cof ! parameter values for l80 belonging to 80% barrier height
  real(sgl), dimension(6, 4)    :: lmxcof ! parameter values for lmax, l-value where barrier disappears
  real(sgl), dimension(5, 6, 4) :: egscof ! parameter values for rotating ground state energy
  real(sgl), dimension(6, 11)   :: x1b    ! parameter for RLDM
  real(sgl), dimension(6, 11)   :: x2b    ! parameter for RLDM
  real(sgl), dimension(10, 20)  :: x3b    ! parameter for RLDM
  real(sgl), dimension(6, 11)   :: x1h    ! parameter for RLDM
  real(sgl), dimension(6, 11)   :: x2h    ! parameter for RLDM
  real(sgl), dimension(10, 20)  :: x3h    ! parameter for RLDM
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for WKB
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, dimension(0:2*numbar)                          :: iiextr    ! WKB counter
  integer                                                 :: nbeta     ! number of beta values
  integer                                                 :: nbinswkb  ! integration step for WKB calculation
  integer                                                 :: nextr     ! WKB counter
  real(sgl), dimension(numbeta)                           :: betafis   ! fission path width
  real(sgl), dimension(0:numZ, 0:numN, 0:numbins, numbar) :: Twkb      ! transmission coefficient of WKB potential
  real(sgl), dimension(0:numZ, 0:numN, 0:numbins, numbar) :: Twkbdir   ! transmission coefficient of WKB potential
  real(sgl), dimension(0:numZ, 0:numN, 0:numbins, numbar) :: Twkbphase ! transmission coefficient of WKB potential
  real(sgl), dimension(0:numZ, 0:numN, 0:numbins, numbar) :: Twkbtrans ! transmission coefficient of WKB potential
  real(sgl), dimension(0:numZ, 0:numN, 0:numbins)         :: Uwkb      ! energy of WKB potential
  real(sgl), dimension(numbeta)                           :: vfis      ! adjustable factor for fission path height
  real(sgl), dimension(2*numbar)                          :: Vheight   ! height of WKB potential
  real(sgl), dimension(2*numbar)                          :: Vpos      ! position of WKB potential
  real(sgl), dimension(2*numbar)                          :: Vwidth    ! width of WKB potential
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for level density
!-----------------------------------------------------------------------------------------------------------------------------------
!
! densitypar
!
  logical, dimension(0:numZ, 0:numN)             :: ldparexist ! flag for existence of tabulated level density
  integer, dimension(0:numZ, 0:numN, 0:numbar)   :: Nlast      ! last discrete level
  real(sgl), dimension(0:numZ, 0:numN, 0:numbar) :: aldcrit    ! critical level density parameter
  real(sgl), dimension(0:numZ, 0:numN, 0:numbar) :: delta      ! energy shift
  real(sgl), dimension(0:numZ, 0:numN)           :: delta0     ! systematical pairing energy
  real(sgl), dimension(0:numZ, 0:numN, 0:numbar) :: Dcrit      ! critical determinant
  real(sgl), dimension(0:numZ, 0:numN, 0:numbar) :: Econd      ! condensation energy
  real(sgl), dimension(0:numZ, 0:numN, 0:numbar) :: Scrit      ! critical entropy
  real(sgl), dimension(0:numZ, 0:numN)           :: Tcrit      ! critical temperature
  real(sgl), dimension(0:numZ, 0:numN, 0:numbar) :: Ucrit      ! critical U
!
! densitymatch
!
  integer                                        :: NLo         ! lowest discrete level for temperature matching
  integer                                        :: NP          ! highest discrete level for temperature matching
  real(sgl)                                      :: E0save      ! E0 value saved for matching routine
  real(sgl), dimension(0:numZ, 0:numN, 0:numbar) :: Ediscrete   ! energy of middle of discrete level region
  real(sgl)                                      :: EL          ! lower matching level energy
  real(sgl)                                      :: EP          ! higher matching level energy
  real(sgl)                                      :: Exmemp      ! empirical estimate for matching point for Ex
  real(sgl), dimension(nummatchT)                :: logrho      ! logarithm of level density
  real(sgl), dimension(0:numZ, 0:numN, 0:numbar) :: scutoffdisc ! spin cutoff factor for discrete level region
  real(sgl), dimension(nummatchT)                :: temprho     ! temperature
  real(sgl)                                      :: Tmemp       ! empirical estimate for temperature
!
! densitytable
!
  logical, dimension(0:numZ, 0:numN, 0:numbar)                       :: ldexist     ! flag for existence of level density
  integer, dimension(0:numZ, 0:numN)                                 :: nendens     ! number of energies for level density grid
  real(sgl), dimension(0:numdens)                                    :: edens       ! energy grid for tabulated level density
  real(sgl), dimension(0:numZ, 0:numN)                               :: Edensmax    ! maximum energy on level density table
  real(dbl), dimension(0:numZ,0:numN,0:numdens,0:numJ,-1:1,0:numbar) :: ldtable     ! level density from table
  real(dbl), dimension(0:numZ,0:numN,0:numdens,0:numbar)             :: ldtottable  ! total level density from table
  real(dbl), dimension(0:numZ,0:numN,0:numdens,-1:1,0:numbar)        :: ldtottableP ! total level density per parity from table
!
! densitycum
!
  real(sgl), dimension(0:numZ,0:numN)            :: CED0       ! C/E of D0
  real(sgl), dimension(0:numZ,0:numN)            :: CGD0       ! C/G of D0
  real(sgl), dimension(0:numZ,0:numN)            :: chi2D0     ! chi2 of D0
  real(sgl), dimension(0:numZ,0:numN)            :: FrmsD0     ! Frms of D0
  real(sgl), dimension(0:numZ,0:numN)            :: ErmsD0     ! Erms of D0
  real(dbl), dimension(0:numZ,0:numN)            :: chi2lev    ! chi2 of discrete levels
  real(dbl), dimension(0:numZ,0:numN)            :: Frmslev    ! Frms of discrete levels
  real(dbl), dimension(0:numZ,0:numN)            :: Ermslev    ! Erms of discrete levels
  real(dbl), dimension(0:numZ,0:numN)            :: avdevlev   ! average deviation from  discrete levels
  real(dbl), dimension(0:numZ,0:numN,0:numlev2)  :: Ncum       ! number of cumulative levels (integral of level density)
  real(dbl), dimension(0:numZ,0:numN,0:numlev2)  :: rhoexp     ! level density of experimental discrete levels
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for particle-hole density tables
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical, dimension(0:numZ,0:numN,0:numexc,0:numexc)                   :: phexist1   ! flag for existence of p-h state den. table
  logical, dimension(0:numZ,0:numN,0:numexc,0:numexc,0:numexc,0:numexc) :: phexist2   ! flag for existence of p-h state density
  integer, dimension(numconf)                                           :: hhtable    ! hole number from table
  integer, dimension(numconf)                                           :: hnutable   ! neutron hole number from table
  integer, dimension(numconf)                                           :: hpitable   ! proton hole number from table
  integer                                                               :: nenphdens  ! number of energies for p-h state den. grid
  integer                                                               :: Nphconf1   ! number of 1-component p-h configurations
  integer                                                               :: Nphconf2   ! number of 2-component p-h configurations
  integer, dimension(numconf)                                           :: pnutable   ! neutron particle number from table
  integer, dimension(numconf)                                           :: ppitable   ! proton particle number from table
  integer, dimension(numconf)                                           :: pptable    ! particle number from table
  real(sgl)                                                             :: Ephdensmax ! maximum energy on p-h state denity table
  real(sgl), dimension(0:1,0:1,0:numexc,0:numexc,0:numdens)                   :: phtable1   ! p-h state density
  real(sgl), dimension(0:1,0:1,0:numexc,0:numexc,0:numexc,0:numexc,0:numdens) :: phtable2   ! p-h state density
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for decay data
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, dimension(-1:numZ, 0:numN, -1:numisom)    :: rtyp      ! type of beta decay, beta-: 1 , beta+: 2 (from ENDF format)
  integer, dimension(-1:numZ, 0:numN, -1:numisom, 5) :: Td        ! half life per time unit
  real(sgl)                                          :: daysec    ! number of seconds in a day
  real(sgl)                                          :: hoursec   ! number of seconds in an hour
  real(sgl), dimension(-1:numZ, 0:numN, -1:numisom)  :: lambda    ! decay rate per isotope
  real(sgl)                                          :: minutesec ! number of seconds in a minute
  real(sgl), dimension(-1:numZ, 0:numN, -1:numisom)  :: Thalf     ! half life of nuclide in sec.
  real(sgl)                                          :: yearsec   ! number of seconds in a year
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for energy grid
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=7)                              :: ecisstatus ! status of ECIS file
  integer, dimension(0:numpar)                  :: ebegin     ! first energy point of energy grid
  integer, dimension(0:numpar)                  :: eendmax    ! last energy point of energy grid
  integer                                       :: maxen      ! total number of energies
  integer, dimension(nummt, -1:numisom)         :: Nrescue    ! number of energies for adjustment factors
  real(sgl), dimension(0:numang)                :: angle      ! angle in degrees
  real(sgl), dimension(0:numangcont)            :: anglecont  ! angle in degrees for continuum
  real(sgl), dimension(0:numpar)                :: coullimit  ! energy limit for charged particle OMP calculation
  real(sgl), dimension(nummt, -1:numisom)       :: Crescue    ! adjustment factor for this incident energy
  real(sgl), dimension(0:numen)                 :: deltaE     ! energy bin around outgoing energies
  real(sgl)                                     :: E1v        ! energy at end of 1/v region
  real(sgl), dimension(0:numen)                 :: Ebottom    ! bottom of outgoing energy bin
  real(sgl), dimension(0:numen)                 :: egrid      ! outgoing energy grid
  real(sgl)                                     :: Einc       ! incident energy in MeV
  real(sgl), dimension(nummt,-1:numisom,numen6) :: Erescue    ! energy grid for adjustment factors
  real(sgl), dimension(0:numen)                 :: Etop       ! top of outgoing energy bin
  real(sgl), dimension(nummt,-1:numisom,numen6) :: frescue    ! adjustment factor
  real(sgl)                                     :: translimit ! limit for transmission coefficient
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for energies
!-----------------------------------------------------------------------------------------------------------------------------------
!
! energies
!
  logical                                   :: flagadd     ! flag for addition of discrete states to spectra flag
  logical                                   :: flagaddel   ! flag for addition of elastic peak to spectra
  logical                                   :: flagcompang ! flag for compound angular distribution calculation
  logical                                   :: flaggiant   ! flag for collective contribution from giant resonances
  logical                                   :: flagmulpre  ! flag for multiple pre-equilibrium calculation
  logical                                   :: flagpreeq   ! flag for pre-equilibrium calculation
  logical                                   :: flagwidth   ! flag for width fluctuation calculation
  logical, dimension(0:numZ, 0:numN)        :: mulpreZN    ! logical for multiple pre-equilibrium per nucleus
  integer, dimension(0:numpar)              :: eend        ! last energy point of energy grid
  integer                                   :: eendhigh    ! last energy point for energy grid for any particle
  integer                                   :: nbins       ! number of continuum excitation energy bins
  integer, dimension(0:numpar)              :: nendisc     ! last discrete bin
  real(sgl)                                 :: Einc0       ! incident energy in MeV
  real(sgl), dimension(0:numpar, 0:numlev2) :: eoutdis     ! outgoing energy of discrete state reaction
  real(sgl)                                 :: Etotal      ! total energy of compound system (target + projectile)
  real(sgl)                                 :: eninccm     ! center-of-mass incident energy in MeV
  real(sgl)                                 :: speceps     ! limit for cross section spectra
  real(sgl)                                 :: wavenum     ! wave number
!
! exgrid
!
  real(dbl), dimension(0:numZ,0:numN,0:numlev) :: Ethresh ! threshold incident energy for residual nucleus
  real(dbl), dimension(0:numZ,0:numN,0:numlev) :: Qres    ! Q-value for residual nucleus
!
! grid
!
  integer                                :: Ninclow  ! number of incident energies below Elow
!
! channels
!
  character(len=18), dimension(0:numchantot)   :: reacstring     ! string for exclusive reaction channel
  integer, dimension(0:numchantot)             :: idchannel      ! identifier for exclusive channel
  real(dbl), dimension(0:numchantot, 0:numlev) :: Ethrexcl       ! threshold incident energy for exclusive channel
  real(dbl), dimension(0:numchantot, 0:numlev) :: Qexcl          ! Q-value for exclusive channel
!
! reaction codes
!
  integer, dimension(nummtall)           :: MTchannel
  character(len=16), dimension(nummtall) :: reactionstring
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for inverse channel data
!-----------------------------------------------------------------------------------------------------------------------------------
!
! inverse
!
  character(len=13) :: csfile    ! file with inverse reaction cross sections
  character(len=13) :: transfile ! file with transmission coefficients
!
! inverseread
!
  real(sgl), dimension(0:numpar, 0:numen, -1:1, 0:numl) :: Tjl    ! transmission coefficient per particle, energy, spin and l-value
  real(sgl), dimension(0:numpar, 0:numen, 0:numl)       :: Tl     ! transmission coefficients per particle, energy and l-value
  real(sgl), dimension(0:numpar, 0:numen)               :: xselas ! total elastic cross section (shape + compound)
  real(sgl), dimension(0:numpar, 0:numen)               :: xsopt  ! optical model reaction cross section
  real(sgl), dimension(0:numpar, 0:numen)               :: xsreac ! reaction cross section
  real(sgl), dimension(0:numpar, 0:numen)               :: xstot  ! total cross section (neutrons only)
!
! inversenorm
!
  real(sgl), dimension(0:numpar)        :: threshnorm ! normalization factor at threshold
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for preequilibrium initialization
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, parameter                               :: numparx=numexc/2 ! maximum number of particles
  integer                                          :: maxexc           ! maximal exciton number
  integer                                          :: maxJph           ! maximal spin for particle-hole states
  integer                                          :: maxpar           ! maximal particle number
  real(sgl), dimension(-1:numparx+1, -1:numparx+1) :: Apauli           ! two-component Pauli blocking correction factor
  real(sgl), dimension(-1:numparx+1,- 1:numparx+1, -1:numparx+1, -1:numparx+1) :: Apauli2 ! two-component Pauli blocking co
  real(sgl)                                        :: Efermi           ! depth of Fermi well
  real(sgl), dimension(0:numexc, 0:numexc)         :: ncomb            ! n!/(n!(n-k)!)
  real(sgl), dimension(0:numexc)                   :: nfac             ! n!
  real(sgl), dimension(2, 2, numparx)              :: Rblann           ! Blann's factor
  real(sgl), dimension(0:numexc, 0:numJ)           :: RnJ              ! spin distribution for particle-hole states
  real(sgl), dimension(0:numexc)                   :: RnJsum           ! (2J+1)*sum over spin distributions
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for exciton model initialization
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(0:numpar, 0:numparx) :: Qfactor ! Q-factor for neutron/proton distinction
  real(sgl), dimension(0:numpar)            :: wfac    ! factor for emission rate
  real(sgl), dimension(2, -200:10*numen)    :: wvol    ! absorption part of the optical potential averaged over the volume
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for initial compound nucleus
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(numfact) :: logfact ! factorial logarithm
  integer                       :: nmold   ! number of points for Gauss-Laguerre integration
  integer                       :: ngoep   ! number of points for Gauss-Legendre integration
  integer                       :: ngoes   ! number of points for Gauss-Legendre integration
  integer                       :: ngoet   ! number of points for Gauss-Legendre integration
  integer                       :: wpower  ! power used for rho*(t**wpower)
  real(sgl), dimension(nummold) :: xmold   ! variables for Gauss-Laguerre integration
  real(sgl), dimension(numgoe)  :: xgoep   ! variables for Gauss-Legendre integration
  real(sgl), dimension(numgoe)  :: xgoes   ! variables for Gauss-Legendre integration
  real(sgl), dimension(numgoe)  :: xgoet   ! variables for Gauss-Legendre integration
  real(sgl), dimension(nummold) :: wmold   ! variables for Gauss-Laguerre integration
  real(sgl), dimension(numgoe)  :: wgoep   ! variables for Gauss-Laguerre integration
  real(sgl), dimension(numgoe)  :: wgoes   ! variables for Gauss-Laguerre integration
  real(sgl), dimension(numgoe)  :: wgoet   ! variables for Gauss-Laguerre integration
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for incident channel
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl)                                               :: channelsum   ! sum over exclusive channel cross sections
  real(sgl), dimension(0:numpar,0:numlev,0:3*numl)        :: cleg         ! compound nucleus Legendre coefficient
  real(sgl), dimension(0:numpar, 0:numex)                 :: contrib      ! contribution to emission spectrum
  real(sgl), dimension(0:numpar, 0:numlev2, 0:numang)     :: directad     ! direct angular distribution
  real(sgl), dimension(0:numpar, 0:numlev2, 0:3*numl)     :: dleg         ! direct reaction Legendre coefficient
  character(len=6), dimension(0:numpar,0:numlev2)         :: dorigin      ! origin of direct cross section (Direct or Preeq)
  integer                                                 :: lmaxinc      ! maximal l-value for transm. coeff. for incident channel
  integer                                                 :: maxA         ! maximal number of nucleons away from initial CN
  real(sgl), dimension(0:numpar)                          :: multiplicity ! particle multiplicity
  real(dbl), dimension(-1:numpar,-1:1)                    :: partdecay    ! total decay per particle and parity
  real(dbl), dimension(-1:numpar)                         :: partdecaytot ! total decay per particle
  real(dbl), dimension(-1:numpar,0:numex,0:numJ,-1:1)     :: popdecay     ! decay from population
  real(sgl), dimension(0:numZ, 0:numN, 0:numex, 0:numJ, -1:1) :: preeqpop ! pre-equilibrium population cross section
  real(sgl), dimension(0:numZ, 0:numN, 0:numex)           :: preeqpopex   ! pre-equilibrium population c.s. summed over J and P
  real(sgl), dimension(0:numang)                          :: ruth         ! elastic/Rutherford ratio
  real(sgl), dimension(0:numang)                          :: elasni       ! nuclear+interference term
  real(sgl), dimension(-1:1,0:numl)                       :: Tjlinc       ! transm. coeff. of spin and l for incident channel
  real(sgl), dimension(0:numl)                            :: Tlinc        ! transm. coeff. as a function of l for incident channel
  real(sgl)                                               :: xsabs        ! absorption cross section
  real(sgl), dimension(-1:numpar)                         :: xsbinary     ! cross section from initial compound to residual nucleus
  real(sgl), dimension(0:numZ, 0:numN, 0:numlev)          :: xsbranch     ! branching ratio for isomeric cross section
  real(sgl), dimension(0:numpar, 0:numJph, -1:1, 0:numen) :: xscollcontJP ! total collective cross section per spin and parity
  real(sgl), dimension(0:numpar)                          :: xscollconttot! total collective cross section in the continuum
  real(sgl), dimension(0:numpar)                          :: xscompcont   ! compound cross section for continuum
  real(sgl)                                               :: xscoupled    ! inelastic cross section from coupled channels
  real(sgl), dimension(0:numpar, 0:numlev2)               :: xsdirdisc    ! direct cross section for discrete state
  real(sgl)                                               :: xsdirdiscsum ! total direct cross section
  real(sgl), dimension(0:numpar)                          :: xsdirdisctot ! direct cross section summed over discrete state
  real(sgl)                                               :: xselasinc    ! total elastic c.s. (neutrons only) for incident channel
  real(sgl), dimension(0:numpar, 0:numen)                 :: xsgr         ! total smoothed giant resonance cross section
  real(sgl)                                               :: xsgrsum      ! sum over giant resonance cross sections
  real(sgl), dimension(0:numpar)                          :: xsgrtot      ! total smoothed giant resonance cross section
  real(sgl), dimension(0:numA)                            :: xsmassprod   ! residual production cross section per mass unit
  real(sgl)                                               :: xsngnsum     ! sum over total (projectile,gamma-ejectile) cross section
  real(sgl)                                               :: xsoptinc     ! optical model reaction c.s. for incident channel
  real(sgl), dimension(0:numpar)                          :: xsparticle   ! total particle production cross section
  real(dbl), dimension(0:numZ,0:numN,0:numex,0:numJ,-1:1) :: xspop        ! population cross section
  real(dbl), dimension(0:numZ, 0:numN, 0:numex)           :: xspopex      ! population cross section summed over spin and parity
  real(dbl), dimension(0:numZ, 0:numN, 0:numex, -1:1)     :: xspopexP     ! population cross section per parity
  real(dbl), dimension(0:numZ, 0:numN)                    :: xspopnuc     ! population cross section per nucleus
  real(dbl), dimension(0:numZ, 0:numN, -1:1)              :: xspopnucP    ! population cross section per nucleus per parity
  real(sgl), dimension(0:numpar, 0:numen)                 :: xspreeq      ! preeq. c.s. per particle typ and outgoing energy
  real(sgl)                                               :: xspreeqsum   ! total preequilibrium cross section summed over particles
  real(sgl), dimension(0:numpar)                          :: xspreeqtot   ! preequilibrium cross section per particle type
  real(sgl)                                               :: xsreacinc    ! reaction cross section for incident channel
  real(sgl)                                               :: xsresprod    ! total residual production (= reaction) cross section
  real(sgl)                                               :: xstotinc     ! total cross section (neutrons only) for incident channel
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for excitation energy grid
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, dimension(0:numZ, 0:numN)                        :: maxex     ! maximum excitation energy bin for excited nucleus
  integer, dimension(0:numZ, 0:numN, 0:numex)               :: maxJ      ! maximal J-value
  integer, dimension(0:numpar)                              :: nexmax    ! maximum excitation energy bin for excited nucleus
  real(sgl), dimension(0:numZ, 0:numN, 0:numex)             :: deltaEx   ! excitation energy bin for population arrays
  real(sgl), dimension(0:numZ, 0:numN, 0:numex+1)           :: Ex        ! excitation energy
  real(sgl), dimension(0:numZ, 0:numN)                      :: Exmax     ! maximum excitation energy for excited nucleus
  real(sgl), dimension(0:numZ, 0:numN)                      :: Exmax0    ! maximum excitation energy (inc. negative energies)
  real(sgl), dimension(0:numZ,0:numN,0:numex+1,0:numJ,-1:1) :: fisfeedJP ! fission contribution from excitation energy bin per J, P
  real(dbl), dimension(0:numZ,0:numN,0:numex,0:numJ,-1:1)   :: rhogrid   ! integrated level density
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for recoil initialization
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(0:numangcont)     :: angcontmax    ! maximum of angular bin
  real(sgl), dimension(0:numangcont)     :: angcontmin    ! minimum of angular bin
  real(sgl), dimension(0:numangrec)      :: angrecmin     ! array for lower bin angle values
  real(sgl), dimension(0:numangrec)      :: angrecmax     ! array for upper bin angle values
  real(sgl), dimension(0:2*numangcont+1) :: cosangcontmax ! cosine of maximum of angular bin
  real(sgl), dimension(0:2*numangcont+1) :: cosangcontmin ! cosine of minimum of angular bin
  real(sgl), dimension(0:numang)         :: cosangmax     ! cosine of maximum of angular bin
  real(sgl), dimension(0:numang)         :: cosangmin     ! cosine of minimum of angular bin
  real(sgl), dimension(0:2*numangrec+1)  :: cosrecmin     ! array for cosine of lower bin angle values
  real(sgl), dimension(0:2*numangrec+1)  :: cosrecmax     ! array for cosine of upper bin angle values
  real(sgl), dimension(0:numang)         :: dcosang       ! width of cosine bin width of cosine bin
  real(sgl), dimension(0:2*numangcont+1) :: dcosangcont   ! width of cosine bin
  real(sgl), dimension(0:2*numangrec+1)  :: dcosangrec    ! array for deltacos
  real(sgl), dimension(0:2*numangcont+1) :: sinangcontmax ! sine of maximum of angular bin
  real(sgl), dimension(0:2*numangcont+1) :: sinangcontmin ! sine of minimum of angular bin
  real(sgl), dimension(0:numang)         :: sinangmax     ! sine of maximum of angular bin
  real(sgl), dimension(0:numang)         :: sinangmin     ! sine of minimum of angular bin
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for giant resonances
!-----------------------------------------------------------------------------------------------------------------------------------
!
! directread
!
  real(sgl), dimension(0:numpar,0:3,2,0:numangcont) :: grcollad ! giant resonance angular distribution
  real(sgl), dimension(0:numpar, 0:3, 2)            :: xsgrcoll ! giant resonance cross section
!
! giant
!
  real(sgl), dimension(0:numpar, 0:numen, 0:numangcont) :: collcontad ! collective angular distribution in the continuum
  real(sgl), dimension(0:numpar, 0:3, 2)                :: eoutgr     ! emission energy
  real(sgl), dimension(0:numpar, 0:numen)               :: xscollcont ! total collective cross section in the continuum
  real(sgl), dimension(0:numpar, 0:numen, 0:numangcont) :: xsgrad     ! smoothed giant resonance angular distribution
  real(sgl), dimension(0:numpar,0:3,2,0:numen)          :: xsgrstate  ! smoothed giant resonance cross section per state
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for direct capture initialization
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer                                         :: ispect       ! model for discrete levels
  integer                                         :: nlevexpracap ! number of experimental levels in the final nucleus
  integer, dimension(0:numZ, 0:numN)              :: nlevracap    ! number of levels in the final nucleus
  integer                                         :: racopt       ! OMP for radiative capture
  real(sgl)                                       :: avncap2      ! real volume diffuseness for JLM
  real(dbl), dimension(numdensracap,0:numJ)       :: chglnegj     ! help variable
  real(dbl), dimension(numdensracap,0:numJ)       :: chglposj     ! help variable
  real(sgl), dimension(0:numZ, 0:numN, 0:numdens) :: edensphjp    ! energy grid of ph spin- and parity-dependent level density
  real(sgl), dimension(numjlm)                    :: jlmracap2    ! JLM potential for direct capture
  real(dbl), dimension(0:numZ,0:numN,0:numdens,0:numJph,-1:1) :: phdensjp ! ph spin- and parity-dependent level density from table
  real(dbl), dimension(0:numZ,0:numN,0:numdens)   :: phdenstot    ! total ph level density from table
  real(sgl)                                       :: rvncap2      ! real volume radius for JLM
  real(sgl), dimension(0:numZ, 0:numN, 0:numex)   :: spectfac     ! spectroscopic factor
  real(sgl)                                       :: vncap2       ! real volume depth for JLM
  real(sgl), dimension(numenin)                   :: xsracap      ! direct radiative capture cross section
  real(sgl), dimension(numenin, 0:1, numgam)      :: xsracapEM    ! direct-semidirect radiative capture cross section
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for direct capture
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl)                                 :: xsracape     ! direct radiative capture cross section
  real(sgl)                                 :: xsracapecont ! direct radiative capture continuum cross section
  real(sgl)                                 :: xsracapedisc ! direct radiative capture discrete cross section
  real(sgl), dimension(0:numex,0:numJ,-1:1) :: xsracappop   ! population cross section for radiative capture
  real(sgl), dimension(0:numex)             :: xsracappopex ! population cross section for radiative capture
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for preequilibrium
!-----------------------------------------------------------------------------------------------------------------------------------
!
! preeq
!
  real(sgl) :: Esurf  ! well depth for surface interaction
  integer   :: h0     ! initial hole number
  integer   :: hnu0   ! initial neutron hole number
  integer   :: hpi0   ! initial proton hole number
  integer   :: p0     ! initial particle number
  integer   :: pnu0   ! initial neutron number
  integer   :: ppi0   ! initial proton number
  real(sgl) :: xsflux ! cross section flux
!
! exciton
!
  real(sgl)                                               :: Ecomp     ! total energy of composite system
  real(sgl), dimension(0:numpar, 0:numen, 0:numJph, -1:1) :: xspreeqJP ! preeq. cross section per particle type, outgoing E, J
  real(sgl), dimension(0:numpar, 0:numparx, 0:numen)      :: xsstep    ! preeq. cross section per particle type, stage and outgoi
!
! exciton2
!
  real(sgl),  dimension(0:numpar, 0:numparx, 0:numparx, 0:numen) :: xsstep2 ! two-component preequilibrium cross section
!
! emissionrate
!
  real(sgl), dimension(0:numpar, 0:numparx, 0:numparx, 0:numen)  :: wemission ! emission rate per particle, exciton number
!
! emissionrate2
!
  real(sgl), dimension(0:numpar, 0:numparx, 0:numparx, 0:numparx, 0:numparx, 0:numen) :: wemission2 ! two-component emission
!
! lifetime2
!
  real(sgl), dimension(0:numparx, 0:numparx, 0:numparx, 0:numparx) :: Spre ! time-integrated strength of two-component exciton
!
! msdplusmsc
!
  real(sgl), dimension(0:numpar, 0:numen2, 0:numangcont) :: xspreeqad ! preequilibrium angular distribution per particle
!
! preeqcomplex
!
  real(sgl), dimension(0:numpar, 0:numen) :: xspreeqbu ! breakup cross section per particle type and outgoing energy
  real(sgl), dimension(0:numpar, 0:numen) :: xspreeqki ! knockout cross section per particle type and outgoing energy
  real(sgl), dimension(0:numpar, 0:numen) :: xspreeqps ! pickup/stripping cross section per particle type and outgoing E
!
! breakup
!
  real(sgl) :: Ca        ! effective Coulomb barrier
  real(sgl) :: Deff      ! effective target-projectile separation
  real(sgl) :: Ecent     ! centroid energy for emission spectrum
  real(sgl) :: Sab       ! separation energy for projectile
!
! breakupAVR
!
  logical                                                    :: breakupexist ! logical for break up file
  real(sgl)                                                  :: ebubin       ! outgoing breakup nucleon energy bin for integration
  real(sgl), dimension(0:numpar, 0:numZ, 0:numN, 0:numenout) :: ENHratio     ! breakup nucleons enhancing reaction cross section rat
  real(sgl), dimension(0:numpar)                             :: xsEB         ! elastic breakup cross section
  real(sgl), dimension(0:numpar)                             :: xsBF         ! nucleon inelastic breakup cross section
  real(sgl), dimension(0:numZ, 0:numN)                       :: xsBFnuc      ! inelastic breakup enhancement brought by break-up
  real(sgl), dimension(0:numpar)                             :: xsBUnuc      ! nucleon breakup cross section
  real(sgl), dimension(0:numZ, 0:numN)                       :: xspopnucT    ! total population c.s. per nucleus including inel BU
!
! preeqcorrect
!
  real(sgl), dimension(0:numpar, 0:numlev2)  :: xspreeqdisc    ! preequilibrium cross section for discrete state
  real(sgl)                                  :: xspreeqdiscsum !  total preequilibrium cross section for discrete states
  real(sgl), dimension(0:numpar)             :: xspreeqdisctot ! preequilibrium cross section summed over discrete state
!
! preeqtotal
!
  real(sgl), dimension(0:numpar)             :: xspreeqtotbu ! preequilibrium cross section per particle type for breakup
  real(sgl), dimension(0:numpar)             :: xspreeqtotki ! preequilibrium cross section per particle type for knockout
  real(sgl), dimension(0:numpar)             :: xspreeqtotps ! preequilibrium cross section per particle type for pickup
  real(sgl), dimension(0:numpar, 0:numparx)  :: xssteptot    ! preequilibrium cross section per particle type and stage
!
! population
!
  real(sgl)                                                               :: preeqnorm  ! preequilibriu
  real(sgl), dimension(0:numZph, 0:numNph, 0:numex, 0:numparx, 0:numparx) :: xspopph    ! population cross section per particle-
  real(sgl), dimension(0:numZph, 0:numNph, 0:numex, 0:numparx, 0:numparx, 0:numparx, 0:numparx) :: xspopph2   ! population cross sec
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for exciton model
!-----------------------------------------------------------------------------------------------------------------------------------
!
! emissionrate
!
  real(sgl), dimension(0:numpar, 0:numparx, 0:numparx) :: wemispart ! emission rate per particle and exciton number
  real(sgl), dimension(0:numparx, 0:numparx)           :: wemistot  ! total emission rate per exciton number
!
! lifetime
!
  real(sgl), dimension(0:numparx, 0:numparx) :: depletion ! depletion factor at each stage
  real(sgl), dimension(0:numparx, 0:numparx) :: tauexc    ! lifetime of exciton state
!
! matrix
!
  real(sgl)                 :: M2      ! square of matrix element
  real(sgl)                 :: M2nunu  ! square of neutron-neutron matrix element
  real(sgl)                 :: M2nupi  ! square of neutron-proton matrix element
  real(sgl)                 :: M2pipi  ! square of proton-proton matrix element
  real(sgl)                 :: M2pinu  ! square of proton-neutron matrix element
  real(sgl), dimension(0:2) :: Wompfac ! adjustable constant for OMP based transition rates
!
! exchange2
!
  real(sgl), dimension(0:numparx, 0:numparx, 0:numparx, 0:numparx) :: Gnupi   ! two-component branching ratio
  real(sgl), dimension(0:numparx, 0:numparx, 0:numparx, 0:numparx) :: Gnuplus ! two-component branching ratio
  real(sgl), dimension(0:numparx, 0:numparx, 0:numparx, 0:numparx) :: Gpinu   ! two-component branching ratio
  real(sgl), dimension(0:numparx, 0:numparx, 0:numparx, 0:numparx) :: Gpiplus ! two-component branching ratio
  real(sgl), dimension(0:numparx, 0:numparx, 0:numparx, 0:numparx) :: Lexc    ! exchange term
  real(sgl), dimension(0:numparx, 0:numparx, 0:numparx, 0:numparx) :: tauexc2 ! lifetime of two-component exciton state
!
! emissionrate2
!
  real(sgl), dimension(0:numpar, 0:numparx, 0:numparx, 0:numparx, 0:numparx) :: wemispart2 ! two-component emission rate per exciton
  real(sgl), dimension(0:numparx, 0:numparx, 0:numparx, 0:numparx)           :: wemistot2  ! total two-comp. emis. rate per exciton
!
! lifetime2
!
  real(sgl), dimension(0:numparx, 0:numparx, 0:numparx, 0:numparx) :: PP2  ! total strength
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for MSD
!-----------------------------------------------------------------------------------------------------------------------------------
!
! msdinit
!
  integer, parameter               :: numJmsd=8 ! maximum spin for MSD
  real(sgl)                        :: dEmsd     ! energy bin for MSD
  real(sgl), dimension(0:numenmsd) :: Emsd      ! minimal outgoing energy for MSD calculation
  integer                          :: maxJmsd   ! maximal spin for MSD calculation
  integer                          :: maxmsd    ! number of MSD steps
  integer                          :: msdbins2  ! number of energy points for MSD calculation
!
! interangle
!
  integer, dimension(0:numangcont, 0:numangcont, 0:numangcont) :: nangleint ! number of possibilities to link intermedi
!
! dwbaecis
!
  real(sgl) :: betamsd ! deformation parameter
  real(sgl) :: Emsdin  ! incident MSD energy
  real(sgl) :: Emsdout ! outgoing MSD energy
  real(sgl) :: Exmsd   ! excitation energy for MSD energy grid
!
! dwbaread
!
  real(sgl), dimension(0:numenmsd, 0:numenmsd, 0:numJmsd, 0:numangcont, 0:2) :: xsdw   ! DWBA angular distribution per angle, incide
  real(sgl), dimension(0:numenmsd, 0:numenmsd, 0:numJmsd, 0:2)               :: xsdwin ! DWBA c.s. per incident E, outgoing E, and J
!
! onecontinuumA
!
  real(sgl), dimension(0:numpar, 0:numpar, 0:numenmsd, 0:numenmsd)               :: xscont1   ! continuum one-step direct c.s.
  real(sgl), dimension(0:numpar, 0:numpar, 0:numenmsd, 0:numenmsd, 0:numangcont) :: xscontad1 ! continuum one-step direct angular di
!
! onestepA
!
  real(sgl), dimension (0:numpar, 0:numen)              :: msdstep1   ! continuum one-step direct cross section (unnormalized)
  real(sgl), dimension(0:numpar, 0:numen, 0:numangcont) :: msdstepad1 ! continuum one-step direct angular distribution
!
! onestepB
!
  real(sgl), dimension(0:numpar, nummsd, 0:numen)               :: msdstep   ! continuum n-step direct cross section
  real(sgl), dimension(0:numpar, nummsd, 0:numen, 0:numangcont) :: msdstepad ! continuum n-step direct angular distribution
!
! onecontinuumB
!
  real(sgl), dimension(0:numpar, 0:numpar, 0:numenmsd, 0:numenmsd)               :: xscont     ! cont. one-step direct c.s. for MSD
  real(sgl), dimension(0:numpar, 0:numpar, 0:numenmsd, 0:numenmsd, 0:numangcont) :: xscontad   ! continuum one-step direct angular d
  real(sgl), dimension(0:numpar, nummsd, 0:numenmsd)                             :: msdstep0   ! n-step cross section for MSD
  real(sgl), dimension(0:numpar, nummsd, 0:numenmsd, 0:numangcont)               :: msdstepad0 ! n-step angular distribution for MSD
!
! msdtotal
!
  real(sgl)                                             :: msdall       ! total multi-step direct cross section
  real(sgl), dimension(0:numpar, nummsd)                :: msdstepint   ! n-step direct cross section integrated over energy
  real(sgl), dimension(0:numpar, nummsd, 0:numangcont)  :: msdstepintad ! n-step direct angular distribution integrated over ener
  real(sgl), dimension(0:numpar)                        :: msdsum       ! multi-step direct c.s. summed over steps, E-integrated
  real(sgl), dimension(0:numpar, 0:numen)               :: msdtot       ! multi-step direct cross section summed over steps
  real(sgl), dimension(0:numpar, 0:numen, 0:numangcont) :: msdtotad     ! multi-step direct angular dist. summed over steps
  real(sgl), dimension(0:numpar, 0:numangcont)          :: msdtotintad  ! multi-step direct ang. dist. summed over steps, E-integrat
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables to prepare information for initial compound nucleus
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer                                            :: tnum     ! counter for width fluctuation calculation
  integer                                            :: tNinc    ! counter for width fluctuation calculation
  real(dbl)                                          :: denomhf  ! denominator for compound nucleus formula
  real(dbl)                                          :: feed     ! feeding term for compound nucleus
  real(dbl), dimension(0:5, numtrans)                :: transjl  ! array for width fluctuation calculation
  real(dbl)                                          :: fiswidth ! fission width
  real(dbl), dimension(0:numpar,0:numex,0:numJ,-1:1) :: enumhf   ! enumerator for compound nucleus formula
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for ECIS calculation of compound cross sections (reference only)
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, parameter                   :: numcomp=100 ! number of nuclides for ECIS compound calculation
  character(len=1), dimension(numcomp) :: pcomp       ! parity of level
  integer                              :: ncont       ! number of continua
  integer                              :: nsp1        ! number of uncoupled states and continua
  integer                              :: nsp2        ! number of uncoupled states with angular distribution
  integer, dimension(0:numcomp)        :: typecomp    ! particle type
  real(sgl), dimension(0:numpar)       :: aldcomp     ! level density parameter with indices (Z,N)
  real(sgl)                            :: bz1         ! elastic enhancement factor
  real(sgl), dimension(0:numpar)       :: E0comp      ! constant of temperature formula
  real(dbl), dimension(numcomp)        :: ejeccomp    ! mass of projectile
  real(sgl), dimension(0:numcomp)      :: elevelcomp  ! energy of level
  real(sgl), dimension(0:numpar)       :: Excomp      ! matching Ex
  real(sgl), dimension(numcomp)        :: jcomp       ! spin of level
  real(dbl), dimension(numcomp)        :: masscomp    ! mass of nucleus with indices (Z,N)
  real(sgl), dimension(numcomp)        :: prodZcomp   ! product of charges
  real(sgl), dimension(numcomp)        :: spincomp    ! spin of incident particle
  real(sgl), dimension(0:numpar)       :: tempcomp    ! nuclear temperature
  real(sgl)                            :: tgo         ! slow s-wave neutron gamma width/spacing
  real(sgl), dimension(0:numpar)       :: Umcomp      ! matching point for U (excitation energy - pairing)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for compound nucleus normalization
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer                            :: J2beg    ! begin of J summation
  integer                            :: J2end    ! end of J summation
  integer                            :: pardif   ! difference between target and compound nucleus parity
  real(sgl)                          :: CNfactor ! factor for compound nucleus cross section: pi/[ k**2 (2s+1)(2I+1) ]
  real(sgl), dimension(-1:1, 0:numJ) :: CNterm   ! compound nucleus formation cross section per spin
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for compound nucleus from target
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, dimension(0:numl)                      :: JmaxU      ! maximal total angular momentum
  integer, dimension(0:numl)                      :: JminU      ! minimal total angular momentum
  integer                                         :: lmaxU      ! maximal orbital angular momentum
  integer                                         :: lminU      ! minimal orbital angular momentum
  integer, dimension(-1:numpar, 0:numl, 0:numJ  ) :: nulj       ! (l,j) number of degrees of freedom for URR calculation
  integer, dimension(0:numl, 0:numJ)              :: Purrlj     ! (l,j) parity for URR calculation
  integer                                         :: tnumi      ! counter for width fluctuation calculation
  integer                                         :: tnumo      ! counter for width fluctuation calculation
  real(sgl)                                       :: dExinc     ! excitation energy bin for mother nucleus
  real(sgl)                                       :: Exinc      ! excitation energy of entrance bin
  real(sgl), dimension(-1:numpar)                 :: Fnorm      ! multiplication factor
  real(sgl), dimension(-1:numpar, 0:numl, 0:numJ) :: Turrlj     ! transmission coefficient for URR calculation
  real(sgl), dimension(0:numl, 0:numJ)            :: Turrljinc  ! incident channel (l,j) transm. coefficient for URR calculation
  real(sgl), dimension(-1:numpar,0:numl,0:numJ)   :: xsbinarylj ! (l,j) cross section from initial compound to residual nucleus
  real(sgl)                                       :: Wab        ! width fluctuation factor
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for width fluctuation
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(numgoe, numgoe, numgoe) :: agoe1     ! variable for GOE triple integral calculation
  real(sgl), dimension(numgoe, numgoe, numgoe) :: agoe2     ! variable for GOE triple integral calculation
  real(sgl), dimension(numgoe, numgoe, numgoe) :: agoe3     ! variable for GOE triple integral calculation
  real(sgl), dimension(numgoe, numgoe, numgoe) :: agoe4     ! variable for GOE triple integral calculation
  real(sgl), dimension(numgoe, numgoe, numgoe) :: agoe5     ! variable for GOE triple integral calculation
  real(sgl), dimension(numgoe, numgoe, numgoe) :: agoe6     ! variable for GOE triple integral calculation
  real(sgl), dimension(numgoe, numgoe, numgoe) :: agoe7     ! variable for GOE triple integral calculation
  real(sgl), dimension(numgoe, numgoe, numgoe) :: agoe8     ! variable for GOE triple integral calculation
  real(sgl), dimension(numtrans)               :: freedom   ! number of degrees of freedom
  real(sgl), dimension(nummold)                :: prodwidth ! product of widths
  real(sgl)                                    :: sgoe1     ! variable for GOE triple integral calculation
  real(sgl)                                    :: sgoe2     ! variable for GOE triple integral calculation
  real(sgl)                                    :: sgoe3     ! variable for GOE triple integral calculation
  real(sgl)                                    :: sgoe4     ! variable for GOE triple integral calculation
  real(sgl)                                    :: sgoe5     ! variable for GOE triple integral calculation
  real(sgl)                                    :: sumhrtw   ! variable for HRTW calculation
  real(dbl), dimension(numtrans)               :: tjlav     ! array for width fluctuation calculation
  real(sgl), dimension(numtrans)               :: vhrtw     ! variable for HRTW calculation
  real(sgl), dimension(numtrans)               :: whrtw     ! variable for HRTW calculation
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for energy grid, level density and transmission coefficients
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, dimension(0:numpar,0:numex)               :: lmaxhf   ! maximal l-value for transmission coefficients
  integer, dimension(numbar)                         :: nbintfis ! number of bins
  real(sgl), dimension(numbinfis,numbar)             :: eintfis  ! excitation energy for fission
  real(dbl), dimension(0:numpar,0:numex,0:numJ,-1:1) :: rho0     ! integrated level density
  real(dbl), dimension(numbinfis,0:numJ,-1:1,numbar) :: rhofis   ! integrated level density corresponding to tfisA
  real(sgl), dimension(0:numZ,0:numN)                :: discfactor! correction for discrete level weight for NL > NT
  real(sgl), dimension(0:numex,0:numgam,0:1)         :: Tgam     ! gamma transmission coefficients
  real(sgl), dimension(0:numpar,0:numex,-1:1,0:numl) :: Tjlnex   ! transmission coefficients for particle, energy, spin and l
  real(sgl), dimension(0:numpar,0:numex,0:numl)      :: Tlnex    ! transmission coefficients for particle, emergy and l
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for fission transmission coefficients
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(dbl), dimension(0:numJ,-1:1)           :: denfis   ! fission level density
  real(dbl), dimension(0:numJ,-1:1)           :: gamfis   ! fission width
  real(dbl), dimension(0:numJ,-1:1,0:numhill) :: rhofisA  ! integrated level density corresponding to tfisA
  real(dbl), dimension(0:numJ,-1:1)           :: taufis   ! fission lifetime
  real(dbl), dimension(0:numJ,-1:1)           :: tfis     ! fission transmission coefficient for Hill-Wheeler magnitude
  real(dbl), dimension(0:numJ,-1:1,0:numhill) :: tfisA    ! transmission coefficient for Hill-Wheeler magnitude
  real(dbl), dimension(0:numJ,-1:1)           :: tfisdown ! fission transmission coefficients
  real(dbl), dimension(0:numJ,-1:1)           :: tfisup   ! fission transmission coefficients
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for astro
!-----------------------------------------------------------------------------------------------------------------------------------
!
! astro
!
  integer                                                             :: maxNastro      ! max number of protons from CN for astro
  integer                                                             :: maxZastro      ! max number of protons from CN for astro
  real(dbl), dimension(0:numZastro, 0:numNastro, numT)                :: macsastro      ! thermonuclear cross section
  real(dbl), dimension(0:numZastro, 0:numNastro, numT, 0:numlev)      :: macsastroex    ! thermonuclear reaction cross section
  real(dbl), dimension(numT)                                          :: macsastrofis   ! thermonuclear reaction cross for fission
  real(dbl), dimension(numT)                                          :: macsastroracap ! thermonuclear cross section for dir cap
  real(dbl), dimension(numT)                                          :: partf          ! integrated partition function
  real(dbl), dimension(0:numZastro, 0:numNastro, numT)                :: rateastro      ! thermonuclear reaction rate factor
  real(dbl), dimension(0:numZastro, 0:numNastro, numT, 0:numlev)      :: rateastroex    ! thermonuclear reaction rate to a given Ex
  real(dbl), dimension(numT)                                          :: rateastrofis   ! thermonuclear reaction rate for fission
  real(dbl), dimension(numT)                                          :: rateastroracap ! thermonuclear reaction rate for dir cap
  real(dbl), dimension(0:numZastro, 0:numNastro, 0:numenin)           :: xsastro        ! astrophysical cross section
  real(dbl), dimension(0:numZastro, 0:numNastro, 0:numenin, 0:numlev) :: xsastroex      ! astrophysical cross section to a given Ex
  real(dbl), dimension(0:numenin)                                     :: xsastrofis     ! astrophysical fission cross section
!
! astroprepare
!
  real(dbl)                                     :: rhoastrotot ! total level density for astrophysical case
  real(dbl), dimension(0:1,0:numex,0:numJ,-1:1) :: Tastroinc   ! transmission coefficient for incident channel for astrophysical c
  real(dbl), dimension(0:1,0:numpar,0:numex,0:numJ,-1:1) :: Tastroout   ! transmission coefficient for outgoing channel for astrophy
  real(dbl)                                     :: Tastrotot   ! total transmission coefficient for astrophysical case
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for URR
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                                       :: flagurrendf ! flag for URR info to ENDF
  real(sgl)                                     :: Rprime0     ! potential scattering radius
  real(sgl), dimension(0:numl,0:numJ)           :: sigurrc     ! (l,j) capture cross section for URR
  real(sgl), dimension(0:numl,0:numJ)           :: sigurrf     ! (l,j) fission cross section for URR
  real(sgl), dimension(0:numl,0:numJ)           :: sigurrs     ! (l,j) scattering cross section for URR
  real(sgl), dimension(0:numl)                  :: spot        ! potential scattering contribution
  real(sgl), dimension(0:numl,0:numJ)           :: strengthlj  ! (l,j) neutron strength function
  real(sgl), dimension(0:numl)                  :: strengthl   ! l neutron strength function
  real(sgl), dimension(4)                       :: xsurrN      ! URR cross section
  real(sgl), dimension(4)                       :: xsurrT      ! URR cross section
  real(sgl), dimension(-1:numpar,0:numl,0:numJ) :: urrwidth    ! channel width in URR
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for binary reactions
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(0:numpar)                          :: binemissum    ! integrated binary emission spectrum
  real(sgl), dimension(0:numpar,0:numex)                  :: feedbinary    ! feeding from first compound nucleus
  real(sgl), dimension(0:numZ,0:numN,0:numex,0:numJ,-1:1) :: sfactor       ! spin factor
  real(sgl), dimension(0:numpar)                          :: Eaveragebin   ! average outgoing energy
  real(sgl), dimension(0:numpar, 0:numlev)                :: xscompdisc    ! compound cross section for discrete state
  real(sgl), dimension(0:numpar)                          :: xscompdisctot ! compound cross section summed over discrete states
  real(sgl)                                               :: xscompel      ! compound elastic cross section
  real(sgl)                                               :: xscompnonel   ! total compound non-elastic cross section
  real(sgl), dimension(0:numpar)                          :: xscompound    ! total compound cross section
  real(sgl), dimension(0:numpar)                          :: xsconttot     ! total cross section for continuum
  real(sgl), dimension(0:numpar)                          :: xsdircont     ! direct cross section for continuum
  real(sgl), dimension(0:numpar)                          :: xsdirect      ! total direct cross section
  real(sgl), dimension(0:numpar, 0:numlev)                :: xsdisc        ! total cross section per discrete state
  real(sgl), dimension(0:numpar)                          :: xsdisctot     ! total cross section summed over discrete states
  real(sgl)                                               :: xselastot     ! total elastic cross section (shape + compound)
  real(sgl)                                               :: xsnonel       ! non-elastic cross section
  real(sgl), dimension(0:numZ, 0:numN)                    :: xspopdir      ! direct population cross section per nucleu
  real(sgl), dimension(0:numpar, 0:numlev)                :: xspopex0      ! binary population cross section
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for binary emission spectra
!-----------------------------------------------------------------------------------------------------------------------------------
!
! binaryspectra
!
  real(sgl), dimension(0:numpar, 0:numen)               :: xsbinemis   ! cross section for emission from first compound nucleus
  real(sgl), dimension(0:numpar, 0:numen, 0:numangcont) :: xsbinemisad ! angular distribution for emission from first compound nucle
  real(sgl), dimension(0:numpar, 0:numen)               :: xscomp      ! compound elastic cross section
  real(sgl), dimension(0:numpar, 0:numen, 0:numangcont) :: xscompad    ! compound emission angular distribution
!
! binemission
!
  real(sgl), dimension(0:numpar, 0:numex, 0:numen)      :: binemis ! emission spectrum from initial compound nucleus
  real(sgl), dimension(0:numpar)                        :: binnorm ! normalization factor for binary spectrum
  real(sgl), dimension(0:numpar, 0:numen)               :: xsemis  ! cross section for emission from compound nucleus
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for angular distributions
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(0:numpar,0:numlev,0:3*numl)   :: cleg0   ! Legendre coefficient normalized to the first one
  real(sgl), dimension(0:numpar, 0:numlev, 0:numang) :: compad  ! compound angular distribution
  real(sgl), dimension(0:numpar, 0:numlev, 0:numang) :: discad  ! discrete state angular distribution
  real(sgl), dimension(0:numpar,0:numlev,0:3*numl)   :: tleg    ! total Legendre coefficient
  real(sgl), dimension(0:numpar,0:numlev,0:3*numl)   :: tlegnor ! total Legendre coefficient normalized to 1
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for multiple emission
!-----------------------------------------------------------------------------------------------------------------------------------
!
! multiple
!
  real(sgl), dimension(0:numex)                       :: Dmulti     ! depletion factor for multiple preequilibrium
  real(sgl), dimension(0:numZ,0:numN)                 :: Fcomp      ! compound population fraction per nucleus
  real(sgl), dimension(0:numZ,0:numN)                 :: Fdir       ! direct population fraction per nucleus
  real(sgl), dimension (0:numZchan,0:numNchan,0:numpar,0:numex+1,0:numex+1) :: feedexcl   ! feeding terms from compound emission
  real(sgl), dimension(0:numZ,0:numN,0:numex+1)       :: fisfeedex  ! fission contribution from excitation energy bin
  real(sgl), dimension(0:numZ,0:numN)                 :: Fpreeq     ! preequilibrium population fraction per nucleus
  real(sgl), dimension(0:numpar,0:numex+1,0:numex+1)  :: mcontrib   ! contribution to emission spectrum
  real(sgl), dimension(0:numpar,0:numex+1,0:numex+1)  :: mpecontrib ! contribution to multiple pre-equilibrium emission
  real(sgl), dimension(0:numZ,0:numN,0:numex+1)       :: popexcl    ! population cross section of bin just before decay
  real(sgl), dimension(0:numpar, 0:numex+1, 0:numen)  :: xsbinspec  ! emission spectrum from compound nucleus per bin
  real(sgl), dimension(0:numZ-2,0:numN-2,-1:numpar)   :: xsfeed     ! cross section from compound to residual nucleus
  real(sgl), dimension(0:numpar,0:numex+1)            :: xsmpe      ! multiple-preequilibrium cross section per energy bin
  real(sgl), dimension(0:numpar, 0:numen)             :: xsmpeemis  ! multiple-preequilibrium emission spectrum from comp. nucleus
  real(sgl), dimension(0:numpar)                      :: xsmpetot   ! total multiple-preequilibrium cross section
  real(sgl), dimension(0:numpar,0:numen)              :: xsmpreeq   ! multiple pre-equilibrium emission spectrum
  real(sgl), dimension(0:numpar,0:numen,0:numangcont) :: xsmpreeqad ! multiple preequilibrium angular distribution
  real(sgl), dimension(-1:numpar)                     :: xsngn      ! total (projectile,gamma-ejectile) cross section
  real(sgl), dimension(0:numpar,0:numen)              :: xsngnspec  ! total (projectile,gamma-ejectile) spectrum
  real(sgl), dimension(0:numpar, 0:numex+1)           :: xspartial  ! emitted cross section flux per energy bin
  real(sgl), dimension(0:numZ,0:numN)                 :: xspopcomp  ! compound population cross section per nucleus
  real(sgl), dimension(0:numZ,0:numN)                 :: xspoppreeq ! preequilibrium population cross section per nucleus
  real(dbl), dimension(numelem,nummass)               :: xspopnuc0  ! population cross section per nucleus
!
! cascade
!
  real(sgl), dimension(0:numZ,0:numN,0:numlev,0:numlev) :: xsgamdis    ! discrete gamma-ray cross section
  real(sgl), dimension(0:numZ, 0:numN)                  :: xsgamdistot ! total discrete gamma-ray cross section
!
! excitation
!
  real(sgl)              :: xsinitpop ! initial population cross section
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for recoil
!-----------------------------------------------------------------------------------------------------------------------------------
!
! recoilinit
!
  integer, dimension(0:numpar)                                       :: iejlab     ! number of ejectile lab bins
  integer                                                            :: irecinit   ! counter
  real(sgl)                                                          :: angcm      ! CM angle with respect to the LAB
  real(sgl), dimension(0:numpar,0:numen2,0:2*numangcont+1)           :: areaejlab  ! total surface of LAB ddx bins
  real(sgl), dimension(0:numZ,0:numN,0:numenrec,0:2*numangrec+1)     :: areareclab ! total surface of LAB ddx bins
  real(sgl), dimension(0:numpar,0:numen2,0:numangcont)               :: ddxejlab   ! array with the ddx spectrum of light particle
  real(sgl), dimension(0:numZ,0:numN,0:numex,0:numenrec,0:numangrec) :: ddxrec     ! array with the lab double diff. xs of the
  real(sgl), dimension(0:numZ,0:numN,0:numex)                        :: ddxrectot  ! array with the total recoil flux in exc. bin
  real(sgl), dimension(0:numpar,0:numen2)                            :: dEejlab    ! width of ejectile lab bin
  real(sgl), dimension(0:numpar,0:numen2)                            :: Eejlab     ! center of ejectile lab bin
  real(sgl), dimension(0:numpar,0:numen2)                            :: Eejlabmax  ! maximum energy of ejectile lab bin
  real(sgl), dimension(0:numpar,0:numen2)                            :: Eejlabmin  ! minimum energy of ejectile lab bin
  real(sgl), dimension(0:numZ,0:numN,0:numenrec)                     :: Erec       ! recoil energy
  real(sgl)                                                          :: Erecinit   ! first compound nucleus recoil energy
  real(sgl), dimension(0:numZ,0:numN,0:numenrec)                     :: Erecmax    ! minimal energy limit of recoil bin
  real(sgl), dimension(0:numZ,0:numN,0:numenrec)                     :: Erecmin    ! minimal energy limit of recoil bin
  real(sgl), dimension(0:numZ,0:numN)                                :: recoilint  ! total recoil integrated over spectrum
  real(sgl), dimension(0:numZ,0:numN,0:numex)                        :: specrecoil ! recoil spectrum
  real(sgl), dimension(0:numpar,0:numen2)                            :: xsejlab    ! LAB ejectile spectrum
  real(sgl), dimension(0:numpar)                                     :: xsejlabint ! LAB energy-integrated spectrum
!
! cm2lab
!
  real(sgl)              :: cosejcm1    ! CM recoil energy cosine corresponding to (E1,ang1)
  real(sgl)              :: cosejcm2    ! CM recoil energy cosine corresponding to (E1,ang2)
  real(sgl)              :: cosejlab11  ! LAB ejectile angle cosine corresponding to (E1,ang1)
  real(sgl)              :: cosejlab12  ! LAB ejectile angle cosine corresponding to (E1,ang2)
  real(sgl)              :: cosejlab21  ! LAB ejectile angle cosine corresponding to (E2,ang1)
  real(sgl)              :: cosejlab22  ! LAB ejectile angle cosine corresponding to (E2,ang2)
  real(sgl)              :: cosreclab11 ! LAB recoil angle cosine corresponding to (E1,ang1)
  real(sgl)              :: cosreclab12 ! LAB recoil angle cosine corresponding to (E1,ang2)
  real(sgl)              :: cosreclab21 ! LAB recoil angle cosine corresponding to (E2,ang1)
  real(sgl)              :: cosreclab22 ! LAB recoil angle cosine corresponding to (E2,ang2)
  real(sgl)              :: Eejcm1      ! lower limit of CM ejectile energy bin
  real(sgl)              :: Eejcm2      ! upper limit of CM ejectile energy bin
  real(sgl)              :: Eejlab11    ! LAB ejectile energy corresponding to (E1,ang1)
  real(sgl)              :: Eejlab12    ! LAB ejectile energy corresponding to (E1,ang2)
  real(sgl)              :: Eejlab21    ! LAB ejectile energy corresponding to (E2,ang1)
  real(sgl)              :: Eejlab22    ! LAB ejectile energy corresponding to (E2,ang2)
  real(sgl)              :: ejectmass   ! ejectile mass
  real(sgl)              :: Ereclab11   ! LAB recoil energy corresponding to (E1,ang1)
  real(sgl)              :: Ereclab12   ! LAB recoil energy corresponding to (E1,ang2)
  real(sgl)              :: Ereclab21   ! LAB recoil energy corresponding to (E2,ang1)
  real(sgl)              :: Ereclab22   ! LAB recoil energy corresponding to (E2,ang2)
  real(sgl)              :: recoilmass  ! recoil mass
  real(sgl)              :: sinejcm1    ! CM recoil energy cosine corresponding to (E2,ang1)
  real(sgl)              :: sinejcm2    ! CM recoil energy cosine corresponding to (E2,ang2)
  real(sgl)              :: sinejlab11  ! LAB ejectile angle sine corresponding to (E1,ang1)
  real(sgl)              :: sinejlab12  ! LAB ejectile angle sine corresponding to (E1,ang2)
  real(sgl)              :: sinejlab21  ! LAB ejectile angle sine corresponding to (E2,ang1)
  real(sgl)              :: sinejlab22  ! LAB ejectile angle sine corresponding to (E2,ang2)
  real(sgl)              :: sinreclab11 ! LAB recoil angle sine corresponding to (E1,ang1)
  real(sgl)              :: sinreclab12 ! LAB recoil angle sine corresponding to (E1,ang2)
  real(sgl)              :: sinreclab21 ! LAB recoil angle sine corresponding to (E2,ang1)
  real(sgl)              :: sinreclab22 ! LAB recoil angle sine corresponding to (E2,ang2)
  real(sgl)              :: vcm         ! compound nucleus velocity
  real(sgl)              :: vejcm1      ! velocity corresponding to Eejcm1
  real(sgl)              :: vejcm2      ! velocity corresponding to Eejcm2
  real(sgl)              :: vreccm1     ! recoil velocity corresponding to Eejcm1
  real(sgl)              :: vreccm2     ! recoil velocity corresponding to Eejcm2
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for exclusive channels
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=18), dimension(0:numchantot)             :: fisstring      ! string for exclusive fission reaction channel
  integer                                                :: idnum          ! counter for exclusive channel
  real(sgl), dimension(0:numchantot, 0:numpar)           :: Eavchannel     ! channel average energy
  real(sgl), dimension(0:numchantot)                     :: Especsum       ! total emission energy
  real(sgl), dimension(0:numchantot, 0:numlev)           :: exclbranch     ! exclusive channel yield per isomer
  real(sgl), dimension(0:numchantot,0:numex+1)           :: gamexcl        ! exclusive gamma cross section per excitation energy
  real(sgl), dimension(0:numchantot)                     :: gmult          ! continuum gamma multiplicity
  real(sgl), dimension(0:numen)                          :: specemis       ! exclusive emission contribution
  real(sgl), dimension(0:numchantot)                     :: xschancheck    ! integrated channel spectrum
  real(sgl), dimension(0:numchantot, 0:numlev)           :: xschaniso      ! channel cross section per isomer
  real(sgl), dimension(0:numchantot)                     :: xschannel      ! channel cross section
  real(sgl), dimension(0:numchantot)                     :: yieldchannel   ! relative yield
  real(sgl), dimension(0:numchantot, 0:numpar, 0:numen)  :: xschannelsp    ! channel cross section spectra
  real(sgl), dimension(0:numchantot,0:numex+1)           :: xsexcl         ! exclusive cross section per excitation energy
  real(sgl), dimension(0:numchantot)                     :: xsfischannel   ! fission channel cross section
  real(sgl), dimension(0:numchantot, 0:numpar, 0:numen)  :: xsfischannelsp ! fission channel spectrum
  real(sgl), dimension(0:numchantot)                     :: xsgamchannel   ! gamma channel cross section
  real(sgl), dimension(0:numchantot, 0:numlev, 0:numlev) :: xsgamdischan   ! discrete gamma channel cross section
  real(sgl), dimension(0:numchantot)                     :: xsfischancheck ! integrated fission channel spectrum
  real(sgl), dimension(0:numpar)                         :: xsparcheck     ! total particle production cross section
  real(sgl), dimension(0:numchantot)                     :: xsratio        ! ratio of exclusive c.s. over residual production c.s.
  real(sgl), dimension(0:numpar, 0:numen)                :: xsspeccheck    ! total particle production spectrum
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for total cross sections
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(0:numpar) :: xsexclcont   ! exclusive single channel cross section for continuum
  real(sgl), dimension(0:numpar) :: xsexclusive  ! exclusive single channel cross section
  real(sgl)                      :: xsfistot     ! total fission cross section
  real(sgl)                      :: xsfistot0    ! total fission cross section
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for spectra
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, dimension(0:numpar)                           :: eendout       ! last energy point of energy grid
  real(sgl), dimension(0:numpar, 0:numen2)               :: buratio       ! break-up ratio
  real(sgl), dimension(numen)                            :: compspect     ! compound part of spectrum
  real(sgl), dimension(0:numpar)                         :: Eaverage      ! average outgoing energy
  real(sgl), dimension(0:numZ,0:numN,0:numpar)           :: Eaveragemul   ! average outgoing energy for multiple emission
  real(sgl), dimension(0:numpar, 0:numen2)               :: espec         ! outgoing energy grid
  real(sgl), dimension(0:numpar, 0:numen2)               :: preeqratio    ! pre-equilibrium ratio
  real(sgl), dimension(numen)                            :: preeqspect    ! multiple pre-equilibrium part of spectrum
  real(sgl), dimension(0:numpar, 0:numen2)               :: xscompout     ! compound emission cross section
  real(sgl), dimension(0:numpar, 0:numen2, 0:numangcont) :: xscompoutad   ! compound emission angular distribution
  real(sgl), dimension(0:numpar, 0:numen2)               :: xsdiscout     ! smoothed cross section for discrete state
  real(sgl), dimension(0:numpar, 0:numen2, 0:numangcont) :: xsdiscoutad   ! smoothed angular distribution for discrete state
  real(sgl), dimension(0:numpar, 0:numen2)               :: xsmpreeqout   ! multiple preequilibrium angular cross section
  real(sgl), dimension(0:numpar, 0:numen2, 0:numangcont) :: xsmpreeqoutad ! multiple preequilibrium angular distribution
  real(sgl), dimension(0:numpar, 0:numen2)               :: xspreeqbuout  ! preequilibrium cross section for breakup
  real(sgl), dimension(0:numpar, 0:numen2)               :: xspreeqkiout  ! preequilibrium cross section for knockout and inelastic
  real(sgl), dimension(0:numpar, 0:numen2)               :: xspreeqout    ! preequilibrium cross section per particle and out E
  real(sgl), dimension(0:numpar, 0:numen2, 0:numangcont) :: xspreeqoutad  ! preequilibrium angular distribution per particle and out
  real(sgl), dimension(0:numpar, 0:numen2)               :: xspreeqpsout  ! preequilibrium cross section for pickup and stripping
  real(sgl), dimension(0:numpar, 0:numen2)               :: xssumout      ! cross section summed over mechanisms
  real(sgl), dimension(0:numpar, 0:numen2, 0:numangcont) :: xssumoutad    ! angular distribution summed over mechanisms
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for mass distribution
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer                                         :: Aff          ! mass number of fission fragment
  integer                                         :: Zff          ! charge number of fission fragment
  integer                                         :: NEpfns       ! number of energies for PFNS
  real(sgl), dimension(nummass)                   :: disa         ! normalised fission fragment mass yield per excitation energy bin
  real(sgl), dimension(nummass)                   :: disacor      ! normalised fission product mass yield per excitation energy bin
  real(sgl), dimension(nummass, numelem)          :: disaz        ! normalised fission fragment isotope yield per exc. energy bin
  real(sgl), dimension(nummass, numelem)          :: disazcor     ! normalised fission product isotope yield per exc. energy bin
  real(sgl), dimension(0:numpar)                  :: Eavpfns      ! average energy of prompt fission neutrons spectrum
  real(sgl), dimension(numelem, numneu)           :: dExcff       ! width of excitation energy of fission fragment
  real(sgl), dimension(numen2)                    :: Epfns        ! energy of PFNS
  real(sgl), dimension(numen2)                    :: dEpfns       ! delta energy of PFNS
  real(sgl), dimension(0:numpar)                  :: Epfnsaverage ! average energy of PFNS
  real(sgl), dimension(numelem, numneu)           :: Excff        ! excitation energy of fission fragment
  real(sgl)                                       :: excfis       ! excitation energy at fission
  real(sgl)                                       :: fpeps        ! ratio for limit for fission product cross section
  real(sgl), dimension(numelem, numneu, 0:1)      :: fpratio      ! fission product isomeric ratio
  real(sgl), dimension(0:numpar, 0:numen2)        :: maxpfns      ! maximum energy of prompt fission neutrons spectrum
  real(sgl), dimension(0:numpar, nummass)         :: nuA          ! nu per A
  real(sgl), dimension(0:numpar, numelem, numneu) :: nuZA         ! nu per Z,A
  real(sgl), dimension(0:numpar, nummass)         :: EaverageA    ! average emission energy per A
  real(sgl), dimension(0:numpar, numelem, numneu) :: EaverageZA   ! average emission energy per (Z,A)
  real(sgl), dimension(0:numpar)                  :: nubar        ! average nu
  real(sgl), dimension(0:numpar, 0:numnu)         :: Pdisnu       ! prompt fission neutrons distribution
  real(sgl), dimension(0:numpar)                  :: Pdisnuav     ! average prompt fission neutrons distribution
  real(sgl), dimension(0:numpar, 0:numen2)        :: pfns         ! prompt fission neutron spectrum
  real(sgl), dimension(0:numpar, 0:numen2)        :: pfnscm       ! prompt fission neutron spectrum in CM
  real(sgl), dimension(numelem, numneu)           :: TKE          ! total kinetic energy
  real(sgl), dimension(nummass)                   :: xsApost      ! post-neutron emission corrected cross section
  real(sgl), dimension(nummass)                   :: xsApre       ! pre-neutron emission cross section
  real(sgl), dimension(numelem, numneu, 0:numlev) :: xsfpex       ! excitation energy spectrum per fission fragment
  real(sgl)                                       :: xstotpost    ! post-neutron emission fission product cross section
  real(sgl)                                       :: xstotpre     ! pre-neutron emission fission product cross section
  real(sgl), dimension(numelem, numneu)           :: xsZApost     ! post-neutron emission corrected isotopic cross section
  real(sgl), dimension(numelem, numneu)           :: xsZApre      ! pre-neutron emission isotopic cross section
  real(sgl), dimension(nummass)                   :: yieldApost   ! post-neutron emission corrected fission yield
  real(sgl), dimension(nummass)                   :: yieldApre    ! pre-neutron emission fission yield
  real(sgl), dimension(numelem, numneu, 0:1)      :: yieldfpex    ! fission yield per isomer
  real(sgl)                                       :: yieldtotpost ! post-neutron emission fission product yield
  real(sgl)                                       :: yieldtotpre  ! pre-neutron emission fission product yield
  real(sgl), dimension(numelem, numneu)           :: yieldZApost  ! post-neutron emission corrected isotopic yield
  real(sgl), dimension(numelem, numneu)           :: yieldZApre   ! pre-neutron emission isotopic yield
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for Brosa model
!-----------------------------------------------------------------------------------------------------------------------------------
!
! brosafy
!
  integer                 :: numtemp ! running variable denoting the numtempth temperature
  real(sgl), dimension(9) :: bf      ! Brosa barrier height
  real(sgl), dimension(9) :: bfsplin ! Brosa barrier height spline fit parameters
  real(sgl), dimension(9) :: hw      ! Brosa barrier width
  real(sgl), dimension(9) :: hwsplin ! Brosa barrier width spline fit parameters
!
! neck
!
  real(sgl) :: cur  ! Brosa parameter for neck rupture
  real(sgl) :: c0   ! Brosa curvature of neck
  real(sgl) :: totl ! Brosa parameter for neck rupture
  real(sgl) :: di   ! Brosa nucleon number density
  real(sgl) :: rest ! Brosa parameter for neck rupture
  real(sgl) :: r1   ! Brosa parameter for neck rupture
  real(sgl) :: r2   ! Brosa parameter for neck rupture
  real(sgl) :: r3   ! Brosa parameter for neck rupture
  real(sgl) :: z1   ! Brosa parameter for neck rupture
  real(sgl) :: z2   ! Brosa parameter for neck rupture
  real(sgl) :: z3   ! Brosa parameter for neck rupture
  real(sgl) :: vtot ! Brosa total volume
  real(sgl) :: rt   ! Brosa parameter for neck rupture
  real(sgl) :: rp   ! Brosa parameter for neck rupture
  real(sgl) :: rpt  ! Brosa parameter for neck rupture
  real(sgl) :: amm  ! Brosa parameter for neck rupture
  real(sgl) :: zee  ! Brosa parameter for neck rupture
  real(sgl) :: ess  ! Brosa parameter for neck rupture
  real(sgl) :: aaa  ! Brosa parameter for neck rupture
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for normalization
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(numen6) :: xseladjust  ! elastic cross section adjustment
  real(sgl), dimension(numen6) :: xsnonadjust ! nonelastic cross section adjustment
  real(sgl), dimension(numen6) :: xstotadjust ! total cross section adjustment
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for thermal cross sections
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(numenlow,0:numchantot,0:numlev)          :: fexclbranch   ! exclusive channel yield per isomer
  real(sgl), dimension(numenlow,0:numpar)                       :: fnubar        ! nubar
  real(sgl), dimension(numenlow,0:numpar)                       :: fxsbinary     ! cross section from initial compound nucleus
  real(sgl), dimension(numenlow,0:numZ,0:numN,0:numlev)         :: fxsbranch     ! branching ratio for isomeric cross section
  real(sgl), dimension(numenlow,0:numchantot,0:numlev)          :: fxschaniso    ! channel cross section per isomer
  real(sgl), dimension(numenlow,0:numchantot)                   :: fxschannel    ! channel cross section
  real(sgl), dimension(numenlow,0:numpar,0:numlev)              :: fxscompdisc   ! compound cross section for discrete state
  real(sgl), dimension(numenlow)                                :: fxscompel     ! compound elastic cross section
  real(sgl), dimension(numenlow)                                :: fxscompnonel  ! total compound non-elastic cross sec
  real(sgl), dimension(numenlow,0:numpar,0:numlev)              :: fxsdirdisc    ! direct cross section for discrete state
  real(sgl), dimension(numenlow)                                :: fxsdirdiscsum ! total direct cross section
  real(sgl), dimension(numenlow,0:numpar,0:numlev)              :: fxsdisc       ! total cross section for discrete state
  real(sgl), dimension(numenlow,0:numpar)                       :: fxsdisctot    ! total cross section summed over disc
  real(sgl), dimension(numenlow)                                :: fxselasinc    ! total elastic cross section (neutrons only)
  real(sgl), dimension(numenlow)                                :: fxselastot    ! total elastic cross section (neutrons only)
  real(sgl), dimension(numenlow,0:numpar)                       :: fxsexclcont   ! excl. single channel cross section for continuum
  real(sgl), dimension(numenlow,0:numpar)                       :: fxsexclusive  ! exclusive single channel cross section
  real(sgl), dimension(numenlow,0:numchantot)                   :: fxsgamchannel ! gamma channel cross section
  real(sgl), dimension(numenlow,0:numchantot,0:numlev,0:numlev) :: fxsgamdischan ! discrete gamma channel cross section
  real(sgl), dimension(numenlow,-1:numpar)                      :: fxsngn        ! total (projectile,gamma-ejectile) cross section
  real(sgl), dimension(numenlow)                                :: fxsnonel      ! non-elastic cross section for incident channel
  real(sgl), dimension(numenlow,0:numZ,0:numN,0:numlev)         :: fxspopex      ! pop. cross section summed over spin and parity
  real(sgl), dimension(numenlow,0:numZ,0:numN)                  :: fxspopnuc     ! population cross section per nucleus
  real(sgl), dimension(numenlow)                                :: fxspreeqsum   ! total preequilibrium cross section section
  real(sgl), dimension(numenlow)                                :: fxsracape     ! direct capture cross section
  real(sgl), dimension(numenlow,0:numchantot)                   :: fxsratio      ! ratio of excl. c.s. over res. production c.s.
  real(sgl), dimension(numenlow)                                :: fxsreacinc    ! reaction cross section for incident channel
  real(sgl), dimension(numenlow)                                :: fxstotinc     ! total cross section (neutrons only)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for isotope production
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, dimension(-1:numZ,-1:numN,-1:numisom)           :: Nenrp        ! number of incident energies for residual prod
  real(sgl)                                                :: targetdx     ! effective thickness of target
  real(sgl)                                                :: Vtar         ! active target volume
  real(sgl)                                                :: Mtar         ! active target mass
  real(sgl)                                                :: projnum      ! number of incident particles [s^-1]
  real(sgl)                                                :: heat         ! produced heat
  real(sgl), dimension(-1:numZ, -1:numN, -1:numisom)       :: prate        ! production rate per isotope
  real(sgl), dimension(-1:numZ,-1:numN,-1:numisom,numenrp) :: Erp          ! incident energy
  real(sgl), dimension(-1:numZ,-1:numN,-1:numisom,numenrp) :: xsrp         ! residual production cross section in mb
  integer                                                  :: Ntime        ! number of time points
  integer, dimension(0:numZ,0:numN,-1:numisom)             :: Tmaxactivity ! time of maximum activity of produced isoto
  integer, dimension(0:numZ,0:numN,-1:numisom,5)           :: Tp           ! irradiation time with maximal yield per time unit
  real(sgl)                                                :: Ntar0        ! number of original target atoms
  real(sgl), dimension(0:numtime)                          :: Tgrid        ! time
  real(sgl)                                                :: Tir          ! irradiation time per unit
  real(sgl)                                                :: Tco          ! cooling time per unit
  real(sgl), dimension(0:numZ,0:numN,-1:numisom,0:numtime) :: Niso         ! number of isotopes produced after irradiation
  real(sgl), dimension(0:numZ,0:numN,-1:numisom,0:numtime) :: activity     ! activity of produced isotope in MBq
  real(sgl), dimension(0:numZ,0:numN,-1:numisom,0:numtime) :: yield        ! yield of produced isotope in MBq/(mA.h)
  real(sgl), dimension(0:numZ,0:numN,-1:numisom,0:numtime) :: Nisorel      ! fraction of number of produced isotopes per ele
  real(sgl), dimension(0:numZ, 0:numtime)                  :: Nisotot      ! number of elemental isotopes produced after irr
  real(sgl), dimension(0:numZ,0:numN,-1:numisom)           :: Tmax         ! irradiation time with maximal yield
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for ENDF data
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer                      :: nen6      ! total number of energies
  real(dbl), dimension(numen6) :: e6        ! energies of ENDF-6 energy grid in MeV
  real(sgl), dimension(numen6) :: xscompel6 ! compound elastic cross section
  real(sgl), dimension(numen6) :: xsnonel6  ! non-elastic cross section
  real(sgl), dimension(numen6) :: xselas6   ! total elastic cross section (neutrons only) for ENDF-6 file
  real(sgl), dimension(numen6) :: xselassh6 ! shape elastic cross section (neutrons only) for ENDF-6 file
  real(sgl), dimension(numen6) :: xsnon6    ! non-elastic cross section for ENDF-6 file
  real(sgl), dimension(numen6) :: xsopt6    ! optical model reaction cross section for ENDF-6 file
  real(sgl), dimension(numen6) :: xsreac6   ! reaction cross section for ENDF-6 file
  real(sgl), dimension(numen6) :: xstot6    ! total cross section (neutrons only) for ENDF-6 file
end module A0_talys_mod
! Copyright A.J. Koning 2025
