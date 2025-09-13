subroutine endfecis
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : ECIS calculation for incident particle on ENDF-6 energy grid
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
!   sgl             ! single precision kind
! All global variables
!   numl            ! number of l values
! Variables for direct reactions
!   flagcoulomb     ! flag for Coulomb excitation calculation with ECIS
!   flagoutecis     ! flag for output of ECIS results
!   flagspher       ! flag to force spherical optical model
!   flagstate       ! flag for optical model potential for each excited state
!   maxband         ! highest vibrational band added to rotational model
!   soswitch        ! switch for deformed spin - orbit calculation
! Variables for basic reaction
!   flagendfecis    ! flag for new ECIS calculation for ENDF - 6 files
!   flagrel         ! flag for relativistic kinematics
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
!   Ztarget         ! charge number of target nucleus
! Variables for energy grid
!   coullimit       ! energy limit for charged particle OMP calculation
!   ecisstatus      ! status of ECIS file
! Variables for nuclides
!   coulbar         ! Coulomb barrier
!   Nindex          ! neutron number index for residual nucleus
!   tarmass         ! mass of target nucleus
!   Zindex          ! charge number index for residual nucleus
! Variables for files
!   nulldev         ! null device
! Constants
!   cparity         ! parity (character)
!   onethird        ! 1 / 3
!   parmass         ! mass of particle in a.m.u.
!   parspin         ! spin of particle
!   parZ            ! charge number of particle
! Variables for deformation parameters
!   colltype        ! type of collectivity (D, V or R)
!   defpar          ! deformation parameter
!   deftype         ! deformation length (D) or parameter (B)
!   indexlevel      ! level index
!   iphonon         ! phonon (1 or 2)
!   Kband           ! magnetic quantum number
!   lband           ! angular momentum
!   leveltype       ! type of level (rotational (R) or vibrational (V))
!   ndef            ! number of collective levels
!   nrot            ! number of deformation parameters for rotational nucleus
!   rotpar          ! deformation parameters for rotational nucleus
!   vibband         ! band number of level
! Variables for levels
!   edis            ! energy of level
!   jdis            ! spin of level
!   parlev          ! parity of level
! Variables for ECIS
!   angbeg          ! first angle
!   angend          ! last angle
!   anginc          ! angle increment
!   d2disp          ! constant for imaginary potential
!   d3disp          ! constant for imaginary potential
!   ecis1           ! 50 input flags ('T' or 'F') for ECIS
!   ecis2           ! 50 input flags ('T' or 'F') for ECIS
!   efer            ! Fermi energy
!   Elevel          ! energy of level
!   flagecisinp     ! flag for existence of ecis input file
!   iband           ! band number of level
!   idvib           ! identifier for existence of vibrational state inside rotational model
!   iph             ! help variable
!   iqm             ! largest order of deformation
!   iqmax           ! maximum l - value of multipole expansion
!   iterm           ! number of iterations
!   Jband           ! angular momentum
!   Jlevel          ! spin of level
!   Kmag            ! magnetic quantum number
!   legendre        ! logical for output of Legendre coefficients
!   Nband           ! number of vibrational bands
!   ncoll           ! number of nuclear states
!   njmax           ! maximal number of j - values in ECIS
!   npp             ! number of optical potentials
!   nrad            ! number of radial points
!   Nrotbeta        ! number of deformation parameters for rotational nucleus
!   Plevel          ! parity of level
!   prodZ           ! product of charges of projectile and target nucleus
!   projmass        ! mass of projectile
!   resmass         ! mass of residual nucleus
!   rmatch          ! matching radius
!   rotbeta         ! deformation parameters for rotational nucleus
!   spin            ! spin of incident particle
!   tarparity       ! parity of target nucleus
!   tarspin         ! spin of target nucleus
!   title           ! title of ECIS input file
!   vibbeta         ! vibrational deformation parameter
!   w2disp          ! constant for imaginary potential
! Variables for optical model
!   jlmexist        ! flag for existence of tabulated radial matter density
! Variables for optical model
!   d2              ! parameter for imaginary surface OMP
!   d3              ! parameter for imaginary surface OMP
!   disp            ! flag for dispersive optical model
!   ef              ! Fermi energy
!   w2              ! parameter for imaginary volume OMP
! Variables for ENDF data
!   e6              ! energies of ENDF - 6 energy grid in MeV
!   nen6            ! total number of energies
!
! *** Declaration of local data
!
  implicit none
  logical            :: jlmloc         ! flag for JLM OMP
  logical            :: rotational     ! flag for rotational input
  logical            :: vibrational    ! flag for vibrational input
  character(len=132) :: outfile        ! output file
  integer            :: i              ! level
  integer            :: i1             ! value
  integer            :: ii             ! counter
  integer            :: nen            ! energy counter
  integer            :: Nix            ! neutron number index for residual nucleus
  integer            :: Zix            ! charge number index for residual nucleus
  real(sgl)          :: e              ! energy
!
! ********************** Set ECIS input parameters *********************
!
! Specific ECIS flags:
!
  legendre = .false.
  rmatch = 0.
  projmass = parmass(k0)
  spin = parspin(k0)
  resmass = tarmass
  prodZ = real(Ztarget * parZ(k0))
  Nband = 0
  angbeg = 0.
  anginc = 180.
  angend = 180.
!
! Loop over energies on ENDF-6 energy grid.
!
  if (flagendfecis) open (unit = 9, file = 'ecisendf.inp', status = 'replace')
  Zix = Zindex(0, 0, k0)
  Nix = Nindex(0, 0, k0)
!
! Standard ECIS inputs for phenomenological optical potentials
!
!   outgoing particle, if available
!   rotational model
!
! Some input flags for ECIS are energy dependent for the rotational model so ecis1 will be defined inside the energy loop.
!
  ecis1 = 'FFFFFTFFFFFFFFFFFFFFFFFFTFFTFFFFFFFFFFFFFFFFFFFFFF'
  ecis2 = 'FFFFFFFFTFFFFFFFTTTFTTTFTFFFFFFFFFFFFFFFFFFFTFFFFF'
!
! 1. Spherical nucleus
!
  jlmloc = .false.
  if (colltype(Zix, Nix) == 'S' .or. flagspher) then
    rotational = .false.
    vibrational = .false.
    title = 'Spherical optical model                           '
    ncoll = 1
    npp = 1
    iterm = 1
    idvib(1) = 1
    Elevel(1) = 0.
    tarspin = 0.
    tarparity = '+'
    if (jlmexist(Zix, Nix, k0)) then
      ecis1(7:7) = 'T'
      ecis1(15:15) = 'T'
      ecis1(29:29) = 'T'
      ecis1(41:41) = 'T'
      rmatch = 18.
      nrad = 182
      jlmloc = .true.
    endif
  else
!
! 2. Deformed nucleus
!
    iterm = 0
    tarspin = jdis(Zix, Nix, 0)
    tarparity = cparity(parlev(Zix, Nix, 0))
    i1 = 0
    do i = 1, ndef(Zix, Nix)
      ii = indexlevel(Zix, Nix, i)
      if (leveltype(Zix, Nix, ii) /= 'V' .and. leveltype(Zix, Nix, ii) /= 'R') cycle
      if (colltype(Zix, Nix) == 'R' .and. vibband(Zix, Nix, i) > maxband) cycle
      i1 = i1 + 1
      idvib(i1) = vibband(Zix, Nix, i)
      Elevel(i1) = edis(Zix, Nix, ii)
      Jlevel(i1) = jdis(Zix, Nix, ii)
      Plevel(i1) = cparity(parlev(Zix, Nix, ii))
      iph(i1) = iphonon(Zix, Nix, i)
      if (iph(i1) == 2) ecis1(2:2) = 'T'
      Jband(i1) = lband(Zix, Nix, i)
      iband(i1) = vibband(Zix, Nix, i)
      Kmag(i1) = Kband(Zix, Nix, i)
      vibbeta(i1) = defpar(Zix, Nix, i)
      Nband = max(Nband, iband(i1))
    enddo
    ncoll = i1
    if (flagstate) then
      npp = ncoll
    else
      npp = 1
    endif
    ecis1(12:12) = 'T'
!
! 2a. Vibrational model
!
    if (colltype(Zix, Nix) == 'V') then
      rotational = .false.
      vibrational = .true.
      title = 'Vibrational optical model                         '
      do i = 1, ncoll
        idvib(i) = 0
      enddo
    else
!
! 2b. Rotational model
!
      rotational = .true.
      vibrational = .false.
      ecis1(1:1) = 'T'
      Nrotbeta = nrot(Zix, Nix)
      do i = 1, Nrotbeta
        rotbeta(i) = rotpar(Zix, Nix, i)
      enddo
      if (colltype(Zix, Nix) == 'R') then
        title = 'Symmetric rotational optical model                '
        iqm = 2 * Nrotbeta
      else
        title = 'Asymmetric rotational optical model               '
        ecis1(3:3) = 'T'
        iqm = 2 * (Nrotbeta - 1)
      endif
      iqmax = 8
    endif
  endif
!
! **************** ECIS input files for several energies ***************
!
  if (deftype(Zix, Nix) == 'B') ecis1(6:6) = 'F'
  if (flagrel) ecis1(8:8) = 'T'
  if (disp(Zix, Nix, k0)) then
    ecis1(10:10) = 'T'
    efer = ef(Zix, Nix, k0)
    w2disp = w2(Zix, Nix, k0)
    d3disp = d3(Zix, Nix, k0)
    d2disp = d2(Zix, Nix, k0)
  endif
  flagecisinp = .false.
  do nen = 1, nen6
    e = real(e6(nen))
    if (k0 > 1 .and. e < coullimit(k0)) cycle
!
! We use a simple formula to estimate the required number of j-values:
!    njmax=2.4*k*R;  R=1.25*A**1/3 ; k=0.22*sqrt(m(in amu)E(in MeV))
! and we always take a minimum of njmax=20.
!
    njmax = int(2.4 * 1.25 * (real(Atarget) **onethird) * 0.22 * sqrt(projmass * e))
    njmax = max(njmax, 20)
    njmax = min(njmax, numl - 2)
    if (jlmloc) njmax = 1600
!
! *************** Calculate optical potential parameters ***************
!
! optical: subroutine for determination of optical potential
!
    call optical(Zix, Nix, k0, e)
!
! ******************* Write ECIS input file ****************************
!
!   iterations in ECIS
! ecisinput  : subroutine to create ECIS input file
!
! For rotational nuclei, the switch at soswitch MeV needs to be made according to Pascal Romain.
!
    if (colltype(Zix, Nix) == 'R' .and. .not.flagspher) then
      if (k0 > 1 .and. flagcoulomb) ecis1(11:11) = 'T'
      if ((k0 == 1 .and. e <= soswitch) .or. (k0 > 1 .and. e <= coulbar(k0))) then
        ecis1(13:13) = 'F'
        ecis1(21:21) = 'T'
        ecis1(42:42) = 'T'
      else
        ecis1(13:13) = 'T'
        ecis1(21:21) = 'F'
        ecis1(42:42) = 'F'
      endif
      if (k0 > 1 .and. e <= 0.05 * coulbar(k0) .and. e <= 2. * Elevel(ncoll)) e = 0.1 * Elevel(ncoll)
      if (flagrel) ecis1(8:8) = 'T'
    endif
    flagecisinp = .true.
    call ecisinput(Zix, Nix, k0, e, rotational, vibrational, jlmloc)
  enddo
  if ( .not. flagendfecis) return
  if ( .not. flagecisinp) then
    close (unit = 9, status = ecisstatus)
    return
  endif
  write(9, '("fin")')
  close (unit = 9)
!
! ************ ECIS calculation for outgoing energies ******************
!
! ecist      : subroutine ecis, adapted for TALYS
!
  if (flagoutecis) then
    outfile = 'ecisendf.out '
  else
    outfile = nulldev
  endif
  call ecist('ecisendf.inp ',outfile,'ecis.endfcs  ','ecis.endfin  ','null         ','null         ','null         ', &
 &  'null         ')
  open (unit = 9, file = 'ecisendf.inp', status = 'unknown')
  close (unit = 9, status = ecisstatus)
  return
end subroutine endfecis
! Copyright A.J. Koning 2021

