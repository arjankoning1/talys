subroutine directecis
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : ECIS calculation of direct cross section
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
!   numl           ! number of l values
!   numlev2        ! maximum number of levels
! Variables for direct reactions
!   flagoutecis    ! flag for output of ECIS results
!   flagstate      ! flag for optical model potential for each excited state
! Variables for basic reaction
!   flagrel        ! flag for relativistic kinematics
! Variables for numerics
!   nangle         ! number of angles
!   nanglecont     ! number of angles for continuum
! Variables for main input
!   Atarget        ! mass number of target nucleus
!   k0             ! index of incident particle
!   Ztarget        ! charge number of target nucleus
! Variables for discrete levels
!   nlev           ! number of levels for nucleus
! Variables for energy grid
!   ecisstatus     ! status of ECIS file
!   Einc           ! incident energy in MeV
! Variables for energies
!   eninccm        ! center - of - mass incident energy in MeV
!   eoutdis        ! outgoing energy of discrete state reaction
! Variables for nuclides
!   AA             ! mass number of residual nucleus
!   Nindex         ! neutron number index for residual nucleus
!   parskip        ! logical to skip outgoing particle
!   Q              ! Q - value
!   Zindex         ! charge number index for residual nucleus
! Variables for files
!   nulldev        ! null device
! Constants
!   cparity        ! parity (character)
!   onethird       ! 1 / 3
!   parA           ! mass number of particle
!   parmass        ! mass of particle in a.m.u.
!   parspin        ! spin of particle
!   parZ           ! charge number of particle
!   sgn            ! sign
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
!   flagecisinp    ! flag for existence of ecis input file
!   hint           ! integration step size h
!   iband          ! band number of level
!   idvib          ! identifier for existence of vibrational state inside rotational model
!   iph            ! help variable
!   iterm          ! number of iterations
!   Jband          ! angular momentum
!   Jlevel         ! spin of level
!   Kmag           ! magnetic quantum number
!   legendre       ! logical for output of Legendre coefficients
!   Nband          ! number of vibrational bands
!   ncoll          ! number of nuclear states
!   njmax          ! maximal number of j - values in ECIS
!   npp            ! number of optical potentials
!   Plevel         ! parity of level
!   projmass       ! mass of projectile
!   prodZ          ! product of charges of projectile and target nucleus
!   resmass        ! mass of residual nucleus
!   rmatch         ! matching radius
!   spin           ! spin of incident particle
!   tarparity      ! parity of target nucleus
!   tarspin        ! spin of target nucleus
!   title          ! title of ECIS input file
!   vibbeta        ! vibrational deformation parameter
!   w2disp         ! constant for imaginary potential
! Variables for deformation parameters
!   betagr         ! deformation parameter for giant resonance
!   deform         ! deformation parameter
!   deftype        ! deformation length (D) or parameter (B)
!   Egrcoll        ! energy of giant resonance
!   jcore          ! spin of level of core nucleus
!   pcore          ! parity of level of core nucleus
! Variables for levels
!   edis           ! energy of level
!   jdis           ! spin of level
!   parlev         ! parity of level
! Variables for masses
!   nucmass        ! mass of nucleus
! Variables for optical model
!   d2             ! parameter for imaginary surface OMP
!   d3             ! parameter for imaginary surface OMP
!   disp           ! flag for dispersive optical model
!   ef             ! Fermi energy
!   rv             ! real volume radius
!   w2             ! parameter for imaginary volume OMP
!
! *** Declaration of local data
!
  implicit none
  logical            :: jlmloc         ! flag for JLM OMP
  logical            :: rotational     ! flag for rotational input
  logical            :: vibrational    ! flag for vibrational input
  character(len=132) :: outfile        ! output file
  integer            :: A              ! mass number of target nucleus
  integer            :: i              ! counter
  integer            :: l              ! multipolarity
  integer            :: Nix            ! neutron number index for residual nucleus
  integer            :: odd            ! odd (1) or even (0) nucleus
  integer            :: type           ! particle type
  integer            :: Zix            ! charge number index for residual nucleus
  real(sgl)          :: vibfactor      ! term to prevent DWBA divergence at very high energy
!
! ********************** Set ECIS input parameters *********************
!
! Specific ECIS flags:
!
  open (unit = 9, file = 'ecisdisc.inp', status = 'replace')
  rotational = .false.
  vibrational = .true.
  jlmloc = .false.
  legendre = .true.
  title = 'Direct discrete cross sections by DWBA            '
  ecis1 = 'FFFFFTFFFFFFFFFFFFFFFFFFFFFTFTFFFFFFFFFFFFFFFFFFFF'
  ecis2 = 'FFFFFFFFTFFFFTTFTTTFTTTFTFFFFFFFFFFFFFFFFFFFTFFFFF'
  if (flagrel) ecis1(8:8) = 'T'
  ncoll = 2
  iterm = 1
  if (flagstate) then
    npp = ncoll
  else
    npp = 1
  endif
  hint = 0.
  rmatch = 0.
!
! We use a simple formula to estimate the required number of j-values:
!    njmax=2.4*k*R;  R=1.25*A**1/3 ; k=0.22*sqrt(m(in amu)E(in MeV))
! and we always take a minimum of njmax=20.
!
  projmass = parmass(k0)
  njmax = int(2.4 * 1.25 * (real(Atarget) **onethird) * 0.22 * sqrt(projmass * Einc))
  njmax = max(njmax, 20)
  njmax = min(njmax, numl)
  spin = parspin(k0)
  prodZ = real(Ztarget * parZ(k0))
  Nband = 1
  if (k0 == 1) then
    angbeg = 0.
  else
    angbeg = 0.00001
  endif
  angend = 180.
!
! ******************* Write ECIS input files ***************************
!
  flagecisinp = .false.
  do type = k0, k0
    if (parskip(type)) cycle
    Zix = Zindex(0, 0, type)
    Nix = Nindex(0, 0, type)
    A = AA(0, 0, type)
    if (deftype(Zix, Nix) == 'B') ecis1(6:6) = 'F'
    if (disp(Zix, Nix, k0)) then
      ecis1(10:10) = 'T'
      efer = ef(Zix, Nix, k0)
      w2disp = w2(Zix, Nix, k0)
      d3disp = d3(Zix, Nix, k0)
      d2disp = d2(Zix, Nix, k0)
    endif
    odd = mod(A, 2)
    resmass = nucmass(Zix, Nix)
    idvib(1) = 0
    idvib(2) = 0
    tarspin = 0.
    tarparity = '+'
    anginc = 180. / nangle
!
! 1. Direct collective states
!
! ecisinput: subroutine to create ECIS input file
!
    vibfactor = 1. - max((Einc - 200.) / 1600., 0.)
    do i = 0, numlev2
      if (i == 0 .and. type == k0) cycle
      if (deform(Zix, Nix, i) == 0.) cycle
      if (eoutdis(type, i) <= 0.) cycle
      Elevel(2) = edis(Zix, Nix, i) - Q(type)
      if (eninccm <= Elevel(2) + 0.1 * parA(type)) cycle
      if (odd == 0) then
        Jlevel(2) = jdis(Zix, Nix, i)
        Plevel(2) = cparity(parlev(Zix, Nix, i))
      else
        if (nlev(Zix, Nix) == 1 .or. (i == 1 .and. jcore(Zix, Nix, i) == 0)) then
          Jlevel(2) = 2.
          Plevel(2) = '+'
        else
          Jlevel(2) = jcore(Zix, Nix, i)
          Plevel(2) = cparity(pcore(Zix, Nix, i))
        endif
      endif
      iband(2) = 1
      Jband(1) = int(Jlevel(2))
      Kmag(1) = 0
      iph(2) = 1
      vibbeta(1) = deform(Zix, Nix, i) * vibfactor
      flagecisinp = .true.
      call ecisinput(Zix, Nix, type, Einc, rotational, vibrational, jlmloc)
    enddo
!
! 2. Giant resonance states
!
    if (type == k0) then
      anginc = 180. / nanglecont
      do l = 0, 3
        do i = 1, 2
          if (betagr(l, i) == 0.) cycle
          Elevel(2) = Egrcoll(l, i)
          if (eninccm <= Elevel(2) + 0.1 * parA(type)) cycle
          Jlevel(2) = real(l)
          iband(2) = 1
          Jband(1) = int(Jlevel(2))
          Kmag(1) = 0
          iph(2) = 1
          Plevel(2) = cparity(int(sgn(l)))
          vibbeta(1) = betagr(l, i)
          if (deftype(Zix, Nix) == 'D') vibbeta(1) = vibbeta(1) * rv * real(Atarget) **onethird
          flagecisinp = .true.
          call ecisinput(Zix, Nix, type, Einc, rotational, vibrational, jlmloc)
        enddo
      enddo
    endif
  enddo
  if ( .not. flagecisinp) then
    close (unit = 9, status = ecisstatus)
    return
  endif
  write(9, '("fin")')
  close (unit = 9)
  legendre = .false.
!
! ************ ECIS calculation for discrete levels ********************
!
! ecist      : subroutine ecis, adapted for TALYS
!
  if (flagoutecis) then
    outfile = 'ecisdisc.out '
  else
    outfile = nulldev
  endif
  call ecist('ecisdisc.inp ', outfile, 'ecis.dircs   ', 'ecis.dirin   ', 'null         ', 'ecis.dirang  ', 'ecis.dirleg  ')
  open (unit = 9, file = 'ecisdisc.inp', status = 'unknown')
  close (unit = 9, status = ecisstatus)
  return
end subroutine directecis
! Copyright A.J. Koning 2021
