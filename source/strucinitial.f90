subroutine strucinitial
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialization of arrays for various structure parameters
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
! Variables for levels
!   bassign        ! flag for assignment of branching ratio
!   branchlevel    ! level to which branching takes place
!   branchratio    ! gamma-ray branching ratio to level
!   conv           ! conversion coefficient
!   edis           ! energy of level
!   ENSDF          ! string from original ENSDF discrete level file
!   eassign        ! flag for assignment of energy
!   jassign        ! flag for assignment of spin
!   jdis           ! spin of level
!   Lisomer        ! level number of isomer
!   nbranch        ! number of branching levels
!   Nisomer        ! number of isomers for this nuclide
!   nlevmax2       ! maximum number of levels
!   parlev         ! parity of level
!   passign        ! flag for assignment of parity
!   tau            ! lifetime of state in seconds
!   tauripl        ! lifetime of state in seconds from RIPL
! Variables for energy grid
!   angle        ! angle in degrees
!   anglecont    ! angle in degrees for continuum
!   Crescue      ! adjustment factor for this incident energy
!   coullimit    ! energy limit for charged particle OMP calculation
!   deltaE       ! energy bin around outgoing energies
!   ebegin       ! first energy point of energy grid
!   Ebottom      ! bottom of outgoing energy bin
!   egrid        ! outgoing energy grid
!   eendmax      ! last energy point of energy grid
!   Erescue      ! energy grid for adjustment factors
!   Etop         ! top of outgoing energy bin
!   frescue      ! adjustment factor
!   Nrescue      ! number of energies for adjustment factors
! Variables for URR
!   flagurrendf     ! flag for URR info to ENDF
! Variables for main input
!   Ltarget         ! excited level of target
! Variables for incident channel
!   lmaxinc         ! maximal l - value for transm. coeff. for incident channel
! Variables for ECIS
!   flaginvecis     ! logical for calculating inverse channel OMP
! Variables for levels
!   Ltarget0        ! excited level of target
! Variables for existence libraries
!   angexist        ! flag for existence of angular distributions
!   chanexist       ! flag for existence of exclusive cross section
!   chanfisexist    ! flag for existence of exclusive fission cros
!   chanisoexist    ! flag for existence of exclusive iso
!   chanopen        ! flag to open channel with first non-zero cross s
!   ddxexist1       ! flag for existence of DDX
!   ddxexist2       ! flag for existence of DDX
!   ddxexist3       ! flag for existence of DDX
!   ddxexist4       ! flag for existence of DDX
!   fisexist        ! flag for existence of fission cross section
!   fpaexist        ! flag for existence of fission product per mass unit
!   fpexist         ! flag for existence of fission product
!   gamchanexist    ! flag for existence of exclusive discrete gam
!   gamexist        ! flag for existence of gamma production cross section
!   idnumfull       ! flag to designate maximum number of exclusive ch.
!   legexist        ! flag for existence of Legendre coefficients
!   nubarexist      ! flag for existence of nubar file
!   prodexist       ! logical to determine existence of residual production
!   recchanexist    ! flag for existence of recoil spectra
!   recexist        ! flag for existence of recoils
!   rpexist         ! flag for existence of residual production cross section
!   rpisoexist      ! flag for existence of isomeric residual production cross section
!   spchanexist     ! flag for existence of exclusive spectra
!   spfischanexist   ! flag for existence of exclusive fission spectra
!   spexist1        ! flag for existence of spectra
!   spexist2        ! flag for existence of spectra
!   urrexist        ! flag for existence of URR
!   Yexist          ! flag for existence of yield
! Variables for energies
!   Ethresh        ! threshold incident energy for residual nucleus
!   Qres           ! Q-value for residual nucleus
!   flagcompang    ! flag for compound angular distribution calculation
!   flaggiant      ! flag for collective contribution from giant resonances
!   flagmulpre     ! flag for multiple pre - equilibrium calculation
!   flagpreeq      ! flag for pre - equilibrium calculation
!   flagwidth      ! flag for width fluctuation calculation
! Variables for masses
!   beta4       ! deformation parameters
!   dumexc      ! theoretical mass excess from Duflo-Zuker formula
!   expmexc     ! experimental mass excess
!   expmass     ! flag for using experimental nuclear mass if available
!   gsparity    ! ground state parity
!   gsspin      ! ground state spin
!   nucmass     ! mass of nucleus
!   redumass    ! reduced mass
!   S           ! separation energy
!   specmass    ! specific mass for residual nucleus
!   thmass      ! theoretical mass
!   thmexc      ! theoretical mass excess
! Variables for deformation parameters
!   betagr        ! deformation parameter for giant resonance
!   colltype      ! type of collectivity (D, V or R)
!   deform        ! deformation parameter
!   defpar        ! deformation parameter
!   deftype       ! deformation length (D) or parameter (B)
!   Egrcoll       ! energy of giant resonance
!   Ggrcoll       ! width of giant resonance
!   indexcc       ! level index for coupled channel
!   indexlevel    ! level index
!   iphonon       ! phonon (1 or 2)
!   Irigid        ! rigid body value of moment of inertia
!   Irigid0       ! undeformed rigid body value of moment of inertia
!   jcore         ! spin of level of core nucleus
!   Kband         ! magnetic quantum number
!   lband         ! angular momentum
!   leveltype     ! type of level (rotational (R) or vibrational (V))
!   ndef          ! number of collective levels
!   nrot          ! number of deformation parameters for rotational nucleus
!   pcore         ! parity of level of core nucleus
!   rotpar        ! deformation parameters for rotational nucleus
!   vibband       ! band number of level
! Variables for resonance parameters
!   D0theo      ! mean s-wave resonance spacing
!   D1theo      ! mean p-wave resonance spacing
!   dD0         ! uncertainty in D0
!   dgamgam     ! uncertainty in gamgam
!   Dl          ! mean resonance spacing per l value
!   Dlj         ! mean resonance spacing per J,l value
!   gamgamth    ! theoretical total radiative width
!   Nrr         ! number of resonances
!   swaveth     ! theoretical strength function for s-wave
!  Variables for gamma-ray strength functions
!   eqrpa         ! energy grid for QRPA strength function
!   fqrpa         ! tabulated QRPA strength function
!   gamkopecky    ! radiative width in eV by spline fit of Kopecky
!   kgr           ! constant for gamma-ray strength function
!   lmax          ! maximal l-value for transmission coefficients
!   ngr           ! number of GR
!   qrpaexist     ! flag for existence of tabulated QRPA strength func.
!   Tqrpa         ! temperature for QRPA
! Variables for level density
!   aldcrit         ! critical level density parameter
!   Dcrit           ! critical determinant
!   delta           ! energy shift
!   delta0          ! systematical pairing energy
!   Econd           ! condensation energy
!   edens           ! energy grid for tabulated level densities
!   Edensmax        ! maximum energy on level density table
!   Ediscrete       ! energy of middle of discrete level region
!   ldexist         ! flag for existence of level density table
!   ldparexist      ! flag for existence of tabulated level density
!   ldtable         ! level density from table
!   ldtottable      ! total level density per parity from table
!   ldtottableP     ! total level density per parity from table
!   logrho          ! logarithm of level density
!   nendens         ! number of energies for level density grid
!   Nlast           ! last discrete level
!   Scrit           ! critical entropy
!   scutoffdisc     ! spin cutoff factor for discrete level region
!   sfactor         ! spin factor
!   Tcrit           ! critical temperature
!   temprho         ! temperature
!   Ucrit           ! critical U
! Variables for OMP
!   av0          ! diffuseness for real volume OMP
!   avd0         ! diffuseness for surface OMP
!   avso0        ! diffuseness for real spin-orbit OMP
!   d1           ! parameter for imaginary surface OMP
!   d2           ! parameter for imaginary surface OMP
!   d3           ! parameter for imaginary surface OMP
!   disp         ! flag for dispersive optical model
!   ef           ! Fermi energy
!   eomp         ! energies on optical model file
!   Eompbeg0     ! upper energy of KD03 OMP
!   Eompbeg1     ! lower energy of alternative OMP
!   Eompend0     ! lower energy of KD03 OMP
!   Eompend1     ! upper energy of alternative
!   jlmexist     ! flag for existence of tabulated radial matter density
!   ompglobal    ! flag for use of global optical model
!   omplines     ! number of lines in optical model file
!   rc0          ! Coulomb radius
!   Rprime       ! potential scattering radius
!   rv0          ! radius for real volume OMP
!   rvd0         ! radius for surface OMP
!   rvso0        ! radius for real spin-orbit OMP
!   Sstrength    ! s,p,d,etc-wave strength function
!   V0           ! V at zero MeV
!   v1           ! parameter for real volume OMP
!   v2           ! parameter for real volume OMP
!   v3           ! parameter for real volume OMP
!   Vjoin        ! V at joining energy
!   vomp         ! optical model parameters from file
!   vso1         ! parameter for real spin-orbit OMP
!   vso2         ! parameter for real spin-orbit OMP
!   w1           ! parameter for imaginary volume OMP
!   w2           ! parameter for imaginary volume OMP
!   w3           ! parameter for imaginary volume OMP
!   w4           ! parameter for imaginary volume OMP
!   Wjoin        ! W at joining energy
!   wso1         ! parameter for imaginary spin-orbit OMP
!   wso2         ! parameter for imaginary spin-orbit OMP
! Variables for OMP
!   normjlm      ! JLM potential normalization factors
!   potjlm       ! JLM potential depth values
!   radjlm       ! radial points for JLM potential
!   rhojlmn      ! density for neutrons
!   rhojlmp      ! density for protons
! Variables for decay data
!   lambda    ! decay rate per isotope
!   rtyp      ! type of beta decay, beta-: 1 , beta+: 2 (from ENDF format)
!   Td        ! half life per time unit
!   Thalf     ! half life of nuclide in sec.
! Variables for fission parameters
!   barcof        ! parameter values for barrier heights at l=0
!   efisc2hb      ! energy of class2 states
!   efisc2rot     ! energy of class2 rotational transition states
!   efistrhb      ! energy of head band transition states
!   efistrrot     ! energy of rotational transition states
!   egscof        ! parameter values for rotating ground state energy
!   Emaxclass2    ! maximum energy for class2 states
!   fecont        ! start of continuum energy
!   jfisc2hb      ! spin of class2 states
!   jfisc2rot     ! spin of class2 rotational transition states
!   jfistrrot     ! spin of rotational transition states
!   jfistrhb      ! spin of head band transition states
!   l20cof        ! parameter values for l20 belonging to 20% barrier height
!   l80cof        ! parameter values for l80 belonging to 80% barrier height
!   lmxcof        ! parameter values for lmax, l-value where barrier disappears
!   minertia      ! moment of inertia of fission barrier deformation
!   minertc2      ! moment of inertia for class2 states
!   nclass2       ! number of sets of class2 states
!   nfisbar       ! number of fission barrier parameters
!   nfisc2hb      ! number of class2 states for barrier
!   nfisc2rot     ! number of class2 rotational transition states for barrier
!   nfistrhb      ! number of head band transition states for barrier
!   nfistrrot     ! number of rotational transition states for barrier
!   pfisc2hb      ! parity of class2 states
!   pfisc2rot     ! parity of class2 rotational transition states
!   pfistrhb      ! parity of head band transition states
!   pfistrrot     ! parity of rotational transition states
! Variables for WKB
!   betafis    ! fission path width
!   iiextr     ! WKB counter
!   vfis       ! adjustable factor for fission path height
!   Vheight    ! height of WKB potential
!   Vpos       ! position of WKB potential
!   Vwidth     ! width of WKB potential
! Variables for direct capture initialization
!   xsracap       ! direct radiative capture cross section
!   xsracapEM     ! direct-semidirect radiative capture cross section as
! Variables for particle-hole state densities
!   Ephdensmax    ! maximum energy on p - h state density table
!   hhtable       ! hole number from table
!   hnutable      ! neutron hole number from table
!   hpitable      ! proton hole number from table
!   nenphdens     ! number of energies for p - h state density gri
!   Nphconf1      ! number of 1 - component p - h configurations
!   Nphconf2      ! number of 2 - component p - h configurations
!   phexist1      ! flag for existence of p-h state density tabl
!   phexist2      ! flag for existence of p-h state density tabl
!   phtable1      ! p-h state density from table
!   phtable2      ! p-h state density from table
!   ppitable      ! proton particle number from table
!   pnutable      ! neutron particle number from table
!   pptable       ! particle number from table
! Variables for normalization
!   xseladjust     ! elastic cross section adjustment
!   xsnonadjust    ! nonelastic cross section adjustment
!   xstotadjust    ! total cross section adjustment
!
! *** Declaration of local data
!
  implicit none
  character(len=2)    :: phstring1(14)
  character(len=4)    :: phstring2(72)
  character(len=132)  :: denfile
  logical             :: lexist
  integer             :: A
  integer             :: i
  integer             :: k
  integer             :: l
  integer             :: N
  integer             :: nen
  integer             :: nex
  integer             :: Nix
  integer             :: Z
  integer             :: Zix
  real(sgl)           :: Eout
  real(sgl)           :: Eeps
  real(sgl)           :: degrid
!
! *********** Initialization of nuclear structure arrays ***************
!
  angle = 0.
  anglecont = 0.
  coullimit = 0.
  Crescue = 1.
  deltaE = 0.
  ebegin = 0
  Ebottom = 0.
  egrid = 0.
  eendmax = 0
  Erescue = 0.
  Etop = 0.
  frescue = 0.
  Nrescue = 0
  bassign = ' '
  branchlevel = 0
  branchratio = 0.
  conv = 0.
  eassign = ' '
  jassign = ' '
  passign = ' '
  levnum = 0
  parlev = 1
  edis = 0.
  jdis = 0.
  tau = 0.
  tauripl = 0.
  ENSDF = ' '
  Lisomer = 0
  nbranch = 0
  Nisomer = 0
  nlevmax2 = 0
  branchdone = 0
  flagwidth = .false.
  flagpreeq = .false.
  flagcompang = .false.
  flaggiant = .false.
  flagmulpre = .false.
  flaginvecis = .true.
  preeqfirst = .true.
  Ethresh = 0.
  Qres = 0.
  flagurrendf = .false.
  JmaxU = 0
  JminU = numJ
  lmaxU = 0
  lminU = numl
  lmaxinc = 0
  Ltarget0 = Ltarget
  chanexist = .false.
  chanfisexist = .false.
  chanisoexist = .false.
  chanopen = .false.
  ddxexist1 = .false.
  ddxexist2 = .false.
  ddxexist3 = .false.
  ddxexist4 = .false.
  legexist = .false.
  angexist = .false.
  fisexist = .false.
  tfisexist = .false.
  gamchanexist = .false.
  gamexist = .false.
  idnumfull = .false.
  opennum = -1
  prodexist = .false.
  recchanexist = .false.
  recexist = .false.
  rpexist = .false.
  rpisoexist = .false.
  spchanexist = .false.
  spfischanexist = .false.
  spexist1 = .false.
  spexist2 = .false.
  breakupexist=.false.
  urrexist = .false.
  Yexist = .false.
  if (nin0 == 0) then
    fpexist = .false.
    fpaexist = .false.
    nubarexist = .false.
  endif
  beta4 = 0.
  dumexc = 0.
  expmexc = 0.
  expmass = 0.
  gsparity = 1
  nucmass = 0.
  redumass = 0.
  S = 0.
  specmass = 0.
  thmass = 0.
  thmexc = 0.
  do Nix = 0, numN + 4
    do Zix = 0, numZ + 4
      Z = Zinit - Zix
      N = Ninit - Nix
      A = Z + N
      if (mod(A, 2) == 0) then
        gsspin(Zix, Nix) = 0.
      else
        gsspin(Zix, Nix) = 0.5
      endif
    enddo
  enddo
  betagr = 0.
  colltype = ' '
  deform = 0.
  defpar = 0.
  deftype = 'B'
  Egrcoll = 0.
  Ggrcoll = 0.
  indexcc = 0
  indexlevel = 0
  iphonon = 0
  Irigid = 0.
  Irigid0 = 0.
  jcore = 0.
  Kband = 0
  lband = 0
  leveltype = 'D'
  ndef = 0
  nrot = 0
  pcore = 1
  rotpar = 0.
  vibband = 0
  D0theo = 0.
  D1theo = 0.
  dD0 = 0.
  D0global = 0.
  dD0global = 0.
  Ncum = 0.
  rhoexp = 0.
  dgamgam = 0.
  Dl = 0.
  Dlj = 0.
  gamgamth = 0.
  Nrr = 0
  swaveth = 0.
  aldcrit = 0.
  Dcrit = 0.
  delta = 0.
  delta0 = 0.
  Econd = 0.
  edens = 0.
  do nex = 1, 20
    edens(nex) = 0.25 * nex
  enddo
  do nex = 21, 30
    edens(nex) = 5. + 0.5 * (nex - 20)
  enddo
  do nex = 31, 40
    edens(nex) = 10. + nex - 30
  enddo
  edens(41) = 22.5
  edens(42) = 25.
  edens(43) = 30.
  do nex = 44, 60
    edens(nex) = 30. + 10. * (nex - 43)
  enddo
  nendens = 60
  edensmax = 200.
  do Nix = 0, numN
    do Zix = 0, numZ
      if (ldmodel(Zix, Nix) == 4) then
        nendens(Zix, Nix) = 55
        Edensmax(Zix, Nix) = 150.
      endif
    enddo
  enddo
  Ediscrete = 0.
  ldexist = .false.
  ldparexist = .false.
  ldtable = 0.
  ldtottable = 0.
  ldtottableP = 0.
  Nlast = 0
  Scrit = 0.
  scutoffdisc = 1.
  Tcrit = 0.
  Ucrit = 0.
  sfactor = 0.
  eqrpa = 0.
  fqrpa = 0.
  qrpaexist = .false.
  Tqrpa = 0.
  ngr = 1
  do l = 1, numgam
    kgr(l) = pi2h2c2 / (2 * l + 1.)
  enddo
  rhojlmn = 0.
  normjlm = 1.
  potjlm = 0.
  radjlm = 0.
  threshnorm = 1.
  lambda = 0.
  rtyp = 0
  prate = 0.
  Td = 0
  Thalf = 1.e30
  nfistrhb = 0
  nfisc2hb = 0
  minertia = 0.
  fecont = 0.
  minertc2 = 0.
  nfisbar = 0
  nclass2 = 0
  nfistrrot = 0
  nfisc2rot = 0
  Emaxclass2 = 0.
  pfistrhb = 1
  pfisc2hb = 1
  efistrhb = 0.
  jfistrhb = 0.
  efisc2hb = 0.
  jfisc2hb = 0.
  pfistrrot = 1
  efistrrot = 0.
  jfistrrot = 0.
  pfisc2rot = 1
  efisc2rot = 0.
  jfisc2rot = 0.
  betafis = 0.
  vfis = 0.
  Vpos = 0.
  Vheight = 0.
  Vwidth = 0.
  phexist2 = .false.
  phexist1 = .false.
  phtable2 = 0.
  phtable1 = 0.
  hhtable = 0
  hnutable = 0
  hpitable = 0
  pptable = 0
  pnutable = 0
  ppitable = 0
  RnJ = 0.
  Ephdensmax = 200.
  nenphdens = 60
  if (phmodel == 2) then
    denfile = trim(path)//'density/ph/Fe.ph'
    inquire (file = denfile, exist = lexist)
    if (lexist) then
      Nphconf2 = 72
      Nphconf1 = 14
      open (unit = 2, file = denfile, status = 'old')
      read(2, '(/////, 9x, 72(a4, 5x), 1x, 14(a2, 7x))') (phstring2(i), i = 1, 72), (phstring1(k), k = 1, 14)
      do i = 1, Nphconf2
        read(phstring2(i), '(4i1)') ppitable(i), hpitable(i), pnutable(i), hnutable(i)
      enddo
      do i = 1, Nphconf1
        read(phstring1(i), '(2i1)') pptable(i), hhtable(i)
      enddo
      close (unit=2)
    else
      Nphconf2 = 0
      Nphconf1 = 0
    endif
  endif
  eintfis = 0.
  rhofis = 0.
  ompglobal = .false.
  ef = 0.
  rc0 = 0.
  rv0 = 0.
  av0 = 0.
  v1 = 0.
  v2 = 0.
  v3 = 0.
  w1 = 0.
  w2 = 0.
  w3 = 0.
  w4 = 0.
  rvd0 = 0.
  avd0 = 0.
  d1 = 0.
  d2 = 0.
  d3 = 0.
  rvso0 = 0.
  avso0 = 0.
  vso1 = 0.
  vso2 = 0.
  wso1 = 0.
  wso2 = 0.
  if (v1adjust(0) == 1.) v1adjust(0) = v1adjust(k0)
  if (v2adjust(0) == 1.) v2adjust(0) = v2adjust(k0)
  if (v3adjust(0) == 1.) v3adjust(0) = v3adjust(k0)
  if (v4adjust(0) == 1.) v4adjust(0) = v4adjust(k0)
  if (rvadjust(0) == 1.) rvadjust(0) = rvadjust(k0)
  if (avadjust(0) == 1.) avadjust(0) = avadjust(k0)
  if (w1adjust(0) == 1.) w1adjust(0) = w1adjust(k0)
  if (w2adjust(0) == 1.) w2adjust(0) = w2adjust(k0)
  if (w3adjust(0) == 1.) w3adjust(0) = w3adjust(k0)
  if (w4adjust(0) == 1.) w4adjust(0) = w4adjust(k0)
  if (rwadjust(0) == 1.) rwadjust(0) = rwadjust(k0)
  if (awadjust(0) == 1.) awadjust(0) = awadjust(k0)
  if (rvdadjust(0) == 1.) rvdadjust(0) = rvdadjust(k0)
  if (avdadjust(0) == 1.) avdadjust(0) = avdadjust(k0)
  if (d1adjust(0) == 1.) d1adjust(0) = d1adjust(k0)
  if (d2adjust(0) == 1.) d2adjust(0) = d2adjust(k0)
  if (d3adjust(0) == 1.) d3adjust(0) = d3adjust(k0)
  if (rwdadjust(0) == 1.) rwdadjust(0) = rwdadjust(k0)
  if (awdadjust(0) == 1.) awdadjust(0) = awdadjust(k0)
  if (vso1adjust(0) == 1.) vso1adjust(0) = vso1adjust(k0)
  if (vso2adjust(0) == 1.) vso2adjust(0) = vso2adjust(k0)
  if (rvsoadjust(0) == 1.) rvsoadjust(0) = rvsoadjust(k0)
  if (avsoadjust(0) == 1.) avsoadjust(0) = avsoadjust(k0)
  if (wso1adjust(0) == 1.) wso1adjust(0) = wso1adjust(k0)
  if (wso2adjust(0) == 1.) wso2adjust(0) = wso2adjust(k0)
  if (rwsoadjust(0) == 1.) rwsoadjust(0) = rwsoadjust(k0)
  if (awsoadjust(0) == 1.) awsoadjust(0) = awsoadjust(k0)
  if (rcadjust(0) == 1.) rcadjust(0) = rcadjust(k0)
  flagompejec = .false.
  disp = .false.
  jlmexist = .false.
  omplines = 0
  eomp = 0.
  vomp = 0.
  xseladjust = 0.
  xsnonadjust = 0.
  xstotadjust = 0.
  Epfns = 0.
  Eout = 0.
  degrid = 0.001
  NEpfns = 0
  nen = 0
  do
    Eout = Eout + degrid
    Eeps = Eout + 1.e-4
    if (Eeps > 50.) exit
    if (nen == numpfns) exit
    nen = nen + 1
    epfns(nen) = Eout
    if (Eeps > 0.01) degrid = 0.01
    if (Eeps > 0.2) degrid = 0.02
    if (Eeps > 0.5) degrid = 0.05
    if (Eeps > 3.) degrid = 0.1
    if (Eeps > 10.) degrid = 0.5
  enddo
  NEpfns = nen
  dEpfns = 0.
  dEpfns(1)=Epfns(1)
  do nen = 2, NEpfns-1
    dEpfns(nen) = 0.5 * (Epfns(nen+1) - Epfns(nen-1))
  enddo
  dEpfns(NEpfns) = 0.5 * (Epfns(NEpfns) - Epfns(NEpfns-1))
  call reactioncodes
  return
end subroutine strucinitial
! Copyright A.J. Koning 2021
