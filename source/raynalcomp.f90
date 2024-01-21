subroutine raynalcomp
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : ECIS calculation of compound cross sections (reference only)
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
! Variables for level density
!   alev           ! level density parameter
!   E0             ! particle constant of temperature formula
!   Exmatch        ! matching point for Ex
!   pair           ! pairing energy
!   T              ! temperature
! Variables for direct reactions
!   flagoutecis    ! flag for output of ECIS results
! Variables for basic reaction
!   flagrel        ! flag for relativistic kinematics
! Variables for numerics
!   nangle         ! number of angles
! Variables for main input
!   Atarget        ! mass number of target nucleus
!   k0             ! index of incident particle
!   Ztarget        ! charge number of target nucleus
! Variables for energy grid
!   Einc           ! incident energy in MeV
! Variables for energies
!   eninccm        ! center - of - mass incident energy in MeV
!   flagwidth      ! flag for width fluctuation calculation
! Variables for excitation energy grid
!   deltaEx        ! excitation energy bin for population arrays
! Variables for nuclides
!   Nindex         ! neutron number index for residual nucleus
!   parinclude     ! logical to include outgoing particle
!   parskip        ! logical to skip outgoing particle
!   Q              ! Q - value
!   tarmass        ! mass of target nucleus
!   Zindex         ! charge number index for residual nucleus
!   ZZ             ! charge number of residual nucleus
! Variables for files
!   nulldev        ! null device
! Constants
!   cparity        ! parity (character)
!   onethird       ! 1 / 3
!   parmass        ! mass of particle in a.m.u.
!   parN           ! neutron number of particle
!   parspin        ! spin of particle
!   parZ           ! charge number of particle
! Variables for resonance parameters
!   swaveth        ! theoretical strength function for s - wave
! Variables for ECIS
!   angbeg         ! first angle
!   angend         ! last angle
!   anginc         ! angle increment
!   ecis1          ! 50 input flags ('T' or 'F') for ECIS
!   ecis2          ! 50 input flags ('T' or 'F') for ECIS
!   hint           ! integration step size h
!   iterm          ! number of iterations
!   ncoll          ! number of nuclear states
!   njmax          ! maximal number of j - values in ECIS
!   npp            ! number of optical potentials
!   projmass       ! mass of projectile
!   prodZ          ! product of charges of projectile and target nucleus
!   resmass        ! mass of residual nucleus
!   rmatch         ! matching radius
!   spin           ! spin of incident particle
!   tarparity      ! parity of target nucleus
!   title          ! title of ECIS input file
! Variables for levels
!   edis           ! energy of level
!   jdis           ! spin of level
!   parlev         ! parity of level
! Variables for level density
!   Nlast          ! last discrete level
! Variables for masses
!   nucmass        ! mass of nucleus
! Variables for OMP
!   disp           ! flag for dispersive optical model
! Variables for ECIS calculation of compound cross sections (reference only)
!   aldcomp        ! level density parameter with indices (Z, N)
!   bz1            ! elastic enhancement factor
!   E0comp         ! constant of temperature formula
!   ejeccomp       ! mass of projectile
!   elevelcomp     ! energy of level
!   Excomp         ! Matching Ex
!   jcomp          ! spin of level
!   masscomp       ! mass of nucleus with indices (Z, N)
!   ncont          ! number of continua
!   nsp1           ! number of uncoupled states and continua
!   nsp2           ! number of uncoupled states with angular distribution
!   pcomp          ! parity of level
!   prodZcomp      ! product of charges
!   spincomp       ! spin of incident particle
!   tempcomp       ! nuclear temperature
!   tgo            ! slow s - wave neutron gamma width / spacing
!   typecomp       ! particle type
!   Umcomp         ! matching point for U (excitation energy - pairing)
!
! *** Declaration of local data
!
  implicit none
  integer   :: ildens ! counter for continua
  integer   :: ilevel ! level counter
  integer   :: nex    ! excitation energy bin of compound nucleus
  integer   :: Nix    ! neutron number index for residual nucleus
  integer   :: NLmax  ! last discrete level
  integer   :: type   ! particle type
  integer   :: Z      ! charge number of target nucleus
  integer   :: Zix    ! charge number index for residual nucleus
  real(sgl) :: econt  ! first continuum energy
  real(sgl) :: ethrcm ! threshold energy
!
! In addition to a calculation by TALYS, a compound nucleus run by ECIS is performed.
! The results will however not be used for TALYS but are just for comparison.
!
! ********************** Set ECIS input parameters *********************
!
! Specific ECIS flags:
!
  open (unit = 1, file = 'eciscomp.inp', status = 'replace')
  title = 'Compound cross sections by ECIS                   '
  ecis1 = 'FFTFFTFFFFFFFFTFFFFFFFFFFFFTFFFFFFFFFFFFFFFFFFFFFF'
  ecis2 = 'FFFFFFFFTFFFFTTFTTTFTTTFTFFFFFTFTTFFFFFFFFFFTFFFFF'
  if (flagrel) ecis1(8:8) = 'T'
  Zix = Zindex(0, 0, k0)
  Nix = Nindex(0, 0, k0)
  if (disp(Zix, Nix, k0)) ecis1(10:10) = 'T'
!
! Gamma emission
!
  if (parinclude(0)) ecis2(36:36) = 'T'
!
! Width fluctuations
!
  if (flagwidth) then
    bz1 = 0.
  else
    bz1 = 1.
    ecis2(37:37) = 'F'
  endif
!
! ******************** Enumerate discrete levels ***********************
!
  ilevel = 0
Loop1:  do type = 1, 6
    if (parskip(type)) cycle
    Zix = Zindex(0, 0, type)
    Nix = Nindex(0, 0, type)
    Z = ZZ(0, 0, type)
    ethrcm = eninccm + Q(type)
    NLmax = min(14, Nlast(Zix, Nix, 0))
    do nex = 0, NLmax
      if (type == k0 .and. nex == 0) cycle
      if (edis(Zix, Nix, nex) > ethrcm) cycle Loop1
      ilevel = ilevel + 1
      typecomp(ilevel) = type
      elevelcomp(ilevel) = edis(Zix, Nix, nex)
      jcomp(ilevel) = jdis(Zix, Nix, nex)
      if (parlev(Zix, Nix, nex) == 1) then
        pcomp(ilevel) = '+'
      else
        pcomp(ilevel) = '-'
      endif
      spincomp(ilevel) = parspin(type)
      ejeccomp(ilevel) = parmass(type)
      masscomp(ilevel) = nucmass(Zix, Nix)
      prodZcomp(ilevel) = parZ(type) * Z
    enddo
  enddo Loop1
!
! ************** Add continua for photons and particles ****************
!
! 1. Photons
!
! swaveth     : theoretical strength function for s-wave
!
  tgo = 0.
  if (parinclude(0)) then
    Zix = 0
    Nix = 0
    aldcomp(0) = alev(Zix, Nix)
    Umcomp(0) = max(Exmatch(Zix, Nix, 0) - pair(Zix, Nix), 0.)
    tempcomp(0) = T(Zix, Nix, 0)
    E0comp(0) = E0(Zix, Nix, 0)
    Excomp(0) = Umcomp(0) + pair(Zix, Nix)
    tgo = swaveth(Zix, Nix)
  endif
!
! 2. Particles
!
  ildens = 0
  do type = 1, 6
    if (parskip(type)) cycle
    Zix = Zindex(0, 0, type)
    Nix = Nindex(0, 0, type)
    ethrcm = eninccm + Q(type)
    econt = edis(Zix, Nix, NLmax) + 0.5 * deltaEx(Zix, Nix, NLmax + 1)
    if (econt < ethrcm) then
      ilevel = ilevel + 1
      typecomp(ilevel) = type
      elevelcomp(ilevel) = econt
      jcomp(ilevel) = 0.
      pcomp(ilevel) = '+'
      spincomp(ilevel) = parspin(type)
      ejeccomp(ilevel) = parmass(type)
      masscomp(ilevel) = nucmass(Zix, Nix)
      prodZcomp(ilevel) = parZ(type) * Z
      ildens = ildens + 1
      aldcomp(ildens) = alev(Zix, Nix)
      Umcomp(ildens) = max(Exmatch(Zix, Nix, 0) - pair(Zix, Nix), 0.)
      tempcomp(ildens) = T(Zix, Nix, 0)
      E0comp(ildens) = E0(Zix, Nix, 0)
      Excomp(ildens) = Umcomp(ildens) + pair(Zix, Nix)
    endif
  enddo
!
! ************************ ECIS parameters *****************************
!
  elevelcomp(0) = 0.
  typecomp(0) = k0
  ncoll = 1
  iterm = 1
  npp = ilevel + 1
  hint = 0.
  rmatch = 0.
  nsp1 = ilevel
  nsp2 = ilevel - ildens
  ncont = ildens
!
! We use a simple formula to estimate the required number of j-values:
!    njmax=2.4*k*R;  R=1.25*A**1/3 ; k=0.22*sqrt(m(in amu)E(in MeV))
! and we always take a minimum of njmax=20.
!
  projmass = parmass(k0)
  njmax = int(2.4 * 1.25 * (real(Atarget) **onethird) * 0.22 * sqrt(projmass * Einc))
  njmax = max(njmax, 20)
  njmax = min(njmax, numl)
  tarparity = cparity(parlev(parZ(k0), parN(k0), 0))
  spin = parspin(k0)
  resmass = tarmass
  prodZ = real(Ztarget * parZ(k0))
  if (k0 == 1) then
    angbeg = 0.
  else
    angbeg = 0.00001
  endif
  anginc = 180. / nangle
  angend = 180.
!
! ******************* Write ECIS input file ****************************
!
! eciscompound      : subroutine to create ECIS input file for compound cross section
!
  call eciscompound
  write(1, '("fin")')
  close (unit = 1)
!
! ************************** ECIS calculation **************************
!
! ecist      : subroutine ecis, adapted for TALYS
!
  if (flagoutecis) then
    call ecist('eciscomp.inp ', 'eciscomp.out ', 'ecis.comcs   ', 'ecis.comin   ', 'null         ', &
      'ecis.comang  ', 'ecis.comleg  ')
  else
    call ecist('eciscomp.inp ', nulldev, 'ecis.comcs   ', 'ecis.comin   ', 'null         ', &
      'ecis.comang  ', 'ecis.comleg  ')
  endif
  return
end subroutine raynalcomp
! Copyright A.J. Koning 2021
