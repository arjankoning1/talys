subroutine incidentecis
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : ECIS calculation for incident energy
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
! All global variables
!   numjlm         ! maximum number of radial points
!   numl           ! number of l values
! Variables for OMP
!   flagoutomp     ! flag for output of optical model parameters
! Variables for direct reactions
!   flagcoulomb    ! flag for Coulomb excitation calculation with ECIS
!   flaginccalc    ! flag for new ECIS calculation for incident
!   flagoutecis    ! flag for output of ECIS results
!   flagspher      ! flag to force spherical optical model
!   flagstate      ! flag for optical model potential for each
!   maxband        ! highest vibrational band added to rotation
!   soswitch       ! switch for deformed spin - orbit calculation
! Variables for basic reaction
!   flagrel        ! flag for relativistic kinematics
! Variables for numerics
!   nangle         ! number of angles
! Variables for main input
!   Atarget        ! mass number of target nucleus
!   k0             ! index of incident particle
!   Ztarget        ! charge number of target nucleus
! Variables for energy grid
!   ecisstatus     ! status of ECIS file
!   Einc           ! incident energy in MeV
! Variables for nuclides
!   AA             ! mass number of residual nucleus
!   coulbar        ! Coulomb barrier
!   Nindex         ! neutron number index for residual nucleus
!   tarmass        ! mass of target nucleus
!   Zindex         ! charge number index for residual nucleus
!   ZZ             ! charge number of residual nucleus
! Constants
!   cparity        ! parity (character)
!   nuc            ! symbol of nucleus
!   onethird       ! 1 / 3
!   parmass        ! mass of particle in a.m.u.
!   parname        ! name of particle
!   parspin        ! spin of particle
!   parZ           ! charge number of particle
! Variables for files
!   nulldev        ! null device
! Variables for deformation parameters
!   colltype       ! type of collectivity (D, V or R)
!   defpar         ! deformation parameter
!   deftype        ! deformation length (D) or parameter (B)
!   indexlevel     ! level index
!   iphonon        ! phonon (1 or 2)
!   Kband          ! magnetic quantum number
!   lband          ! angular momentum
!   leveltype      ! type of level (rotational (R) or vibrational (V)
!   ndef           ! number of collective levels
!   nrot           ! number of deformation parameters for rotational
!   rotpar         ! deformation parameters for rotational nucleus
!   vibband        ! band number of level
! Variables for levels
!   edis           ! energy of level
!   jdis           ! spin of level
!   parlev         ! parity of level
! Variables for ECIS
!   angbeg         ! first angle
!   angend         ! last angle
!   anginc         ! angle increment
!   d2disp         ! constant for imaginary potential
!   d3disp         ! constant for imaginary potential
!   ecis1          ! 50 input flags ('T' or 'F') for ECIS
!   ecis2          ! 50 input flags ('T' or 'F') for ECIS
!   efer           ! Fermi energy
!   Elevel         ! energy of level
!   hint           ! integration step size h
!   iband          ! band number of level
!   idvib          ! identifier for existence of vibrational state inside rotational
!   iph            ! help variable
!   iqm            ! largest order of deformation
!   iqmax          ! maximum l - value of multipole expansion
!   iterm          ! number of iterations
!   Jband          ! angular momentum
!   Jlevel         ! spin of level
!   Kmag           ! magnetic quantum number
!   legendre       ! logical for output of Legendre coefficients
!   Nband          ! number of vibrational bands
!   ncoll          ! number of nuclear states
!   njmax          ! maximal number of j - values in ECIS
!   npp            ! number of optical potentials
!   nrad           ! number of radial points
!   Nrotbeta       ! number of deformation parameters for rotational nucleus
!   Plevel         ! parity of level
!   prodZ          ! product of charges of projectile and target nucleus
!   projmass       ! mass of projectile
!   resmass        ! mass of residual nucleus
!   rmatch         ! matching radius
!   rotbeta        ! deformation parameters for rotational nucleus
!   spin           ! spin of incident particle
!   tarparity      ! parity of target nucleus
!   tarspin        ! spin of target nucleus
!   title          ! title of ECIS input file
!   vibbeta        ! vibrational deformation parameter
!   w2disp         ! constant for imaginary potential
! Variables for optical model
!   ef             ! Fermi energy
!   jlmexist       ! flag for existence of tabulated radial matter density
! Variables for JLM
!   normjlm        ! JLM potential normalization factors
!   potjlm         ! JLM potential depth values
!   radjlm         ! radial points for JLM potential
! Variables for optical model
!   av             ! real volume diffuseness
!   avd            ! real surface diffuseness
!   avso           ! real spin - orbit diffuseness
!   aw             ! imaginary volume diffuseness
!   awd            ! imaginary surface diffuseness
!   awso           ! imaginary spin - orbit diffuseness
!   d2             ! parameter for imaginary surface OMP
!   d3             ! parameter for imaginary surface OMP
!   disp           ! flag for dispersive optical model
!   ef             ! Fermi energy
!   rc             ! Coulomb radius
!   rv             ! real volume radius
!   rvd            ! real surface radius
!   rvso           ! real spin - orbit radius
!   rw             ! imaginary volume radius
!   rwd            ! imaginary surface radius
!   rwso           ! imaginary spin - orbit radius
!   v              ! real volume depth
!   vd             ! real surface depth
!   vso            ! real spin - orbit depth
!   w              ! imaginary volume depth
!   w2             ! parameter for imaginary volume OMP
!   wd             ! imaginary surface depth
!   wso            ! imaginary spin - orbit depth
!
! *** Declaration of local data
!
  implicit none
  logical            :: jlmloc         ! flag for JLM OMP
  logical            :: rotational     ! flag for rotational input
  logical            :: vibrational    ! flag for vibrational input
  character(len=132) :: inelfile       ! file for inelastic scattering cross sections
  character(len=132) :: outfile        ! output file
  character(len=21) :: optfile    ! file with optical potential
  character(len=13) :: Estr
  character(len=132) :: topline    ! topline
  character(len=15) :: col(5)     ! header
  character(len=15) :: un(5)     ! header
  character(len=80) :: quantity   ! quantity
  integer            :: A              ! mass number of target nucleus
  integer            :: i              ! counter
  integer            :: i1             ! value
  integer            :: ii             ! counter
  integer            :: Ncol
  integer            :: Nix            ! neutron number index for residual nucleus
  integer            :: Z              ! charge number of target nucleus
  integer            :: Zix            ! charge number index for residual nucleus
  real(sgl)          :: Ein            ! incident energy
!
! ********************** Set ECIS input parameters *********************
!
! Specific ECIS flags:
!
  if (flaginccalc) open (unit = 9, file = 'ecisinc.inp', status = 'unknown')
  legendre = .true.
  Ein = Einc
  hint = 0.
  rmatch = 0.
!
! We use a simple formula to estimate the required number of j-values:
!    njmax=2.4*k*R;  R=1.25*A**1/3 ; k=0.22*sqrt(m(in amu)E(in MeV))
! and we always take a minimum of njmax=20.
!
  projmass = parmass(k0)
  njmax = int(2.4 * 1.25 * (real(Atarget) **onethird) * 0.22 * sqrt(projmass * Ein))
  njmax = max(njmax, 20)
  njmax = min(njmax, numl)
  spin = parspin(k0)
  resmass = tarmass
  prodZ = real(Ztarget * parZ(k0))
  Nband = 0
  if (k0 == 1) then
    angbeg = 0.
  else
    angbeg = 0.00001
  endif
  anginc = 180. / nangle
  angend = 180.
  Zix = Zindex(0, 0, k0)
  Nix = Nindex(0, 0, k0)
!
! Standard ECIS inputs for phenomenological optical potentials
!
  ecis1 = 'FFFFFTFFFFFFFFFFFFFFFFFFTFFTFFFFFFFFFFFFFFFFFFFFFF'
  ecis2 = 'FFFFFFFFTFFFTTTFTTTFTTTFTFFFFFFFFFFFFFFFFFFFTFFFFF'
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
    inelfile = 'null         '
    if (jlmexist(Zix, Nix, k0)) then
      ecis1(7:7) = 'T'
      ecis1(15:15) = 'T'
      ecis1(29:29) = 'T'
      ecis1(41:41) = 'T'
      hint = 0.1
      rmatch = 18.
      nrad = 182
      njmax = 1600
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
!
! Some input flags for ECIS are energy dependent for the rotational model.
! For rotational nuclei, the switch at soswitch MeV needs to be made according to Pascal Romain.
!
      if (k0 > 1 .and. flagcoulomb) ecis1(11:11) = 'T'
      if ((k0 == 1 .and. Ein <= soswitch) .or. (k0 > 1 .and. Ein <= coulbar(k0))) then
        ecis1(13:13) = 'F'
        ecis1(21:21) = 'T'
        ecis1(42:42) = 'T'
      else
        ecis1(13:13) = 'T'
        ecis1(21:21) = 'F'
        ecis1(42:42) = 'F'
      endif
      if (k0 > 1 .and. Ein <= 0.05 * coulbar(k0) .and. Ein <= 2. * Elevel(ncoll)) Ein = 0.1 * Elevel(ncoll)
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
    inelfile = 'ecis.incin '
  endif
!
! ************** Calculate optical potential parameters ****************
!
! optical   : subroutine for determination of optical potential
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
  call optical(Zix, Nix, k0, Ein)
!
! Output of optical model parameters, if requested.
!
  if (flagoutomp) then
    Z = ZZ(0, 0, k0)
    A = AA(0, 0, k0)
    if (jlmexist(0, 0, k0)) then
      if (k0 <= 2) then
        call mom(Zix, Nix, dble(prodZ), dble(Ein))
        write(*, '(/" +++++++++ JLM OPTICAL MODEL POTENTIAL FOR INCIDENT CHANNEL ++++++++++",/)')
      else
        call foldalpha(Zix, Nix, Ein)
        write(*, '(/" +++++++++ DOUBLE FOLDING OPTICAL MODEL POTENTIAL FOR INCIDENT CHANNEL ++++++++++",/)')
      endif
      Estr=''
      write(Estr,'(es13.6)') Einc
      un = 'MeV'
      col(1) = 'radius'
      un(1) = 'fm'
      col(2) = 'V'
      col(3) = 'W'
      col(4) = 'Vso'
      col(5) = 'Wso'
      Ncol = 5
      optfile='optE0000.000.'//parsym(k0)
      write(optfile(5:12), '(f8.3)') Einc
      write(optfile(5:8), '(i4.4)') int(Einc)
      quantity='optical potential'
      open (unit=1, file=optfile, status='unknown')
      topline=trim(targetnuclide)//' '//parname(k0)//' optical potential at '//Estr//' MeV'
      call write_header(topline,source,user,date,oformat)
      call write_target
      call write_datablock(quantity,Ncol,numjlm,col,un)
      do i = 1, numjlm
        write(1, '(5es15.6)') radjlm(Zix, Nix, i), normjlm(Zix, Nix, 1) * potjlm(Zix, Nix, i, 1), &
 &        normjlm(Zix, Nix, 2) * potjlm(Zix, Nix, i, 2), normjlm(Zix, Nix, 5) * potjlm(Zix, Nix, i, 5), &
 &        normjlm(Zix, Nix, 6) * potjlm(Zix, Nix, i, 6)
      enddo
      close(unit = 1)
      call write_outfile(optfile,flagoutall)
    else
      write(*, '(/" +++++++++ OPTICAL MODEL PARAMETERS FOR INCIDENT CHANNEL ++++++++++")')
      write(*, '(/11x, a8, " on ", i3, a2/)') parname(k0), A, nuc(Z)
      write(*, '("  Energy", 5x, "V", 5x, "rv", 4x, "av", 4x, "W", 5x, "rw", &
 &      4x, "aw", 4x, "Vd", 3x, "rvd", 3x, "avd", 4x, "Wd", 3x, "rwd", 3x, "awd", 3x, "Vso", 3x, "rvso", 2x, "avso", &
 &      2x, "Wso", 3x, "rwso", 2x, "awso", 2x, "rc", 5x, "Ef" /)')
      write(*, '(1x, f8.3, 1x, 6(f6.2, f6.3, f6.3), f6.3, f8.3)') &
 &      Ein, v, rv, av, w, rw, aw, vd, rvd, avd, wd, rwd, awd, vso, rvso, avso, wso, rwso, awso, rc, ef(Zix, Nix, k0)
    endif
  endif
  if ( .not. flaginccalc) return
!
! ******************* Write ECIS input file ****************************
!
! ecisinput: subroutine to create ECIS input file
!
  call ecisinput(Zix, Nix, k0, Ein, rotational, vibrational, jlmloc)
  write(9, '("fin")')
  close (unit = 9)
  legendre = .false.
!
! ************ ECIS calculation for incident energy ********************
!
! ecist      : subroutine ecis, adapted for TALYS
!
  if (flagoutecis) then
    outfile = 'ecisinc.out  '
  else
    outfile = nulldev
  endif
  call ecist('ecisinc.inp  ', outfile, 'ecis.inccs   ', inelfile, 'ecis.inctr   ', 'ecis.incang  ', 'ecis.incleg  ')
  open (unit = 9, file = 'ecisinc.inp', status = 'unknown')
  close (unit = 9, status = ecisstatus)
  return
end subroutine incidentecis
! Copyright A.J. Koning 2021
