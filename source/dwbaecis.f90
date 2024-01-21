subroutine dwbaecis
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : ECIS calculations of DWBA for MSD
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
!   flagoutecis     ! flag for output of ECIS results
! Variables for basic reaction
!   flagrel         ! flag for relativistic kinematics
! Variables for numerics
!   nanglecont      ! number of angles for continuum
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
! Variables for preequilibrium
!   flagecisdwba    ! flag for new ECIS calculation for DWBA for MSD
!   flagonestep     ! flag for continuum one - step direct only
!   flagoutdwba     ! flag for output of DWBA cross sections for MSD
! Variables for energy grid
!   ecisstatus      ! status of ECIS file
!   Einc            ! incident energy in MeV
! Variables for nuclides
!   Nindex          ! neutron number index for residual nucleus
!   parskip         ! logical to skip outgoing particle
!   Zindex          ! charge number index for residual nucleus
! Variables for files
!   nulldev         ! null device
! Constants
!   onethird        ! 1 / 3
! Variables for ECIS
!   angbeg          ! first angle
!   angend          ! last angle
!   anginc          ! angle increment
!   ecis1           ! 50 input flags ('T' or 'F') for ECIS
!   ecis2           ! 50 input flags ('T' or 'F') for ECIS
!   hint            ! integration step size h
!   iterm           ! number of iterations
!   ncoll           ! number of nuclear states
!   njmax           ! maximal number of j - values in ECIS
!   npp             ! number of optical potentials
!   projmass        ! mass of projectile
!   rmatch          ! matching radius
!   title           ! title of ECIS input file
! Variables for masses
!   S               ! separation energy
!   specmass        ! specific mass for residual nucleus
! Variables for MSD
!   betamsd         ! deformation parameter
!   Emsd            ! minimal outgoing energy for MSD calculation
!   Emsdin          ! incident MSD energy
!   Emsdout         ! outgoing MSD energy
!   Exmsd           ! excitation energy for MSD energy grid
!   maxJmsd         ! maximal spin for MSD calculation
!   msdbins2        ! number of energy points for MSD calculation
!
! *** Declaration of local data
!
  implicit none
  logical            :: lexist     ! logical to determine existence
  character(len=132) :: outfile    ! output file
  integer            :: ii         ! counter
  integer            :: itype      ! help variable
  integer            :: nen1       ! energy counter
  integer            :: nen1end    ! help variable
  integer            :: nen2       ! energy counter
  integer            :: Nix        ! neutron number index for residual nucleus
  integer            :: type       ! particle type
  integer            :: Zix        ! charge number index for residual nucleus
  real(sgl)          :: QQ         ! Q-value
!
! ********************** Set ECIS input parameters *********************
!
  open (unit = 9, file = 'ecisdwba.inp', status = 'replace')
  title = 'DWBA cross sections for MSD                       '
  ecis1 = 'FFFFFFFFFFFTFFFFFFFFFFFFFFFTFTFFFFFFFFFFFFFFFFFFFF'
  ecis2 = 'FFFFFFFFTFFFFTFFTTTFTTTFTFFFFFFFFFFFFFFFFFFFTFFFFF'
  if (flagrel) ecis1(8:8) = 'T'
  ncoll = maxJmsd + 2
  iterm = 1
  npp = 3
  hint = 0.
  rmatch = 15.
!
! We use a simple formula to estimate the required number of j-values:
!    njmax=2.4*k*R;  R=1.25*A**1/3 ; k=0.22*sqrt(m(in amu)E(in MeV))
! and we always take a minimum of njmax=20.
!
  njmax = int(2.4 * 1.25 * (real(Atarget) **onethird) * 0.22 * sqrt(projmass * Einc))
  njmax = max(njmax, 20)
  njmax = min(njmax, numl)
  betamsd = 0.02
  if (k0 == 1) then
    angbeg = 0.
  else
    angbeg = 0.00001
  endif
  anginc = 180. / nanglecont
  angend = 180.
!
! ***************** Loop over possible one-step reactions **************
!
! ecisdwbamac  : subroutine to create ECIS input file for macroscopic DWBA calculation for MSD
! dwbaread     : subroutine to read ECIS results for DWBA for MSD
! dwbaout      : subroutine for output of ECIS results for DWBA for MSD
! dwbaint      : subroutine to interpolate DWBA cross sections for MSD
! onecontinuumA: subroutine for unnormalized one-step direct cross sections for MSD
! onestepA     : subroutine for unnormalized one-step direct cross sections for outgoing energy grid
!
  if (flagoutdwba) write( * , '(/" ++++++++++ DWBA CROSS SECTIONS FOR MSD ++++++++++")')
!
! *************************** Macroscopic MSD **************************
!
! 2. First exchange one-step reaction for multi-step
!
  do ii = 1, 2
    if (ii == 2) then
      inquire (file = 'ecis.msdin', exist = lexist)
      if ( .not. lexist) then
        write(*, '(" TALYS-error: The first calculation of a run should always be done with ecissave y and ecisdwba y")')
        stop
      endif
      open (unit = 8, file = 'ecis.msdang', status = 'unknown')
      open (unit = 10, file = 'ecis.msdin', status = 'unknown')
    endif
    itype = k0
    do type = 1, 2
      if (parskip(type)) cycle
      Zix = Zindex(0, 0, type)
      Nix = Nindex(0, 0, type)
      QQ = S(0, 0, itype) - S(0, 0, type)
      if (flagonestep) then
        nen1end = 0
      else
        nen1end = msdbins2
      endif
      do nen1 = 0, nen1end, 2
        Emsdin = real(specmass(Zix, Nix, itype) * Emsd(nen1))
        do nen2 = nen1 + 2, msdbins2, 2
          Emsdout = Emsd(nen2)
          Exmsd = Emsdin - Emsdout + QQ
          if (Exmsd < 0..or.(Exmsd + 0.1) >= Emsdin) cycle
          if (ii == 1 .and. flagecisdwba) call ecisdwbamac(itype, type)
          if (ii == 2) then
            call dwbaread(nen1, nen2)
            if (flagoutdwba) call dwbaout(itype, type, nen1, nen2)
          endif
        enddo
      enddo
      if (ii == 2) then
        call dwbaint
        if ( .not. flagonestep) call onecontinuumA(itype, type)
        call onestepA(type)
      endif
    enddo
    if (.not. flagonestep) then
!
! 3. Inelastic one-step reaction for multi-step
!
      do type = 1, 2
        if (parskip(type)) cycle
        if (type == k0) cycle
        itype = type
        Zix = Zindex(0, 0, type)
        Nix = Nindex(0, 0, type)
        do nen1 = 0, msdbins2, 2
          Emsdin = real(specmass(Zix, Nix, itype) * Emsd(nen1))
          do nen2 = nen1 + 2, msdbins2, 2
            Emsdout = Emsd(nen2)
            Exmsd = Emsdin - Emsdout
            if (Exmsd < 0..or.(Exmsd + 0.1) >= Emsdin) cycle
            if (ii == 1 .and. flagecisdwba) call ecisdwbamac(itype, type)
            if (ii == 2) then
              call dwbaread(nen1, nen2)
              if (flagoutdwba) call dwbaout(itype, type, nen1, nen2)
            endif
          enddo
        enddo
        if (ii == 2) then
          call dwbaint
          call onecontinuumA(itype, type)
        endif
      enddo
!
! 4. Second exchange one-step reaction for multi-step
!
      do itype = 1, 2
        if (parskip(itype)) cycle
        if (itype == k0) cycle
        type = k0
        Zix = Zindex(0, 0, type)
        Nix = Nindex(0, 0, type)
        QQ = S(0, 0, itype) - S(0, 0, type)
        do nen1 = 0, msdbins2, 2
          Emsdin = real(specmass(Zix, Nix, itype) * Emsd(nen1))
          do nen2 = nen1 + 2, msdbins2, 2
            Emsdout = Emsd(nen2)
            Exmsd = Emsdin - Emsdout + QQ
            if (Exmsd < 0..or.(Exmsd + 0.1) >= Emsdin) cycle
            if (ii == 1 .and. flagecisdwba) call ecisdwbamac(itype, type)
            if (ii == 2) then
              call dwbaread(nen1, nen2)
              if (flagoutdwba) call dwbaout(itype, type, nen1, nen2)
            endif
          enddo
        enddo
        if (ii == 2) then
          call dwbaint
          call onecontinuumA(itype, type)
        endif
      enddo
    endif
    if (ii == 1 .and. flagecisdwba) then
      write(9, '("fin")')
      close (unit = 9)
!
! **************** ECIS calculation for DWBA for MSD *******************
!
! ecist      : subroutine ecis, adapted for TALYS
!
      if (flagoutecis) then
        outfile = 'ecisdwba.out '
      else
        outfile = nulldev
      endif
      call ecist('ecisdwba.inp ', outfile, 'ecis.msdcs   ', 'ecis.msdin   ', 'null         ', &
        'ecis.msdang  ', 'null         ')
      open (unit = 9, file = 'ecisdwba.inp', status = 'unknown')
      close (unit = 9, status = ecisstatus)
    endif
    if (ii == 2) then
      close (unit = 8, status = ecisstatus)
      close (unit = 10, status = ecisstatus)
    endif
  enddo
  return
end subroutine dwbaecis
! Copyright A.J. Koning 2021
