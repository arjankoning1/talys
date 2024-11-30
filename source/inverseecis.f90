subroutine inverseecis(Zcomp, Ncomp)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : ECIS calculation for outgoing particles and energy grid
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
!   numjlm          ! maximum number of radial points
!   numl            ! number of l values
! Variables for OMP
!   flagoutomp      ! flag for output of optical model parameters
! Variables for direct reactions
!   flageciscalc    ! flag for new ECIS calculation for outgoing
!   flagoutecis     ! flag for output of ECIS results
!   flagrot         ! flag for use of rotational optical model p
!   flagstate       ! flag for optical model potential for each
!   maxband         ! highest vibrational band added to rotation
!   soswitch        ! switch for deformed spin - orbit calculation
! Variables for basic reaction
!   flagrel         ! flag for relativistic kinematics
! Variables for energy grid
!   ebegin          ! first energy point of energy grid
!   ecisstatus      ! status of ECIS file
!   eendmax         ! last energy point of energy grid
!   egrid           ! outgoing energy grid
! Variables for inverse channel data
!   csfile          ! file with inverse reaction cross sections
!   transfile       ! file with transmission coefficients
! Variables for nuclides
!   AA              ! mass number of residual nucleus
!   coulbar         ! Coulomb barrier
!   invexist        ! logical to state necessity of new inverse cross section
!   Nindex          ! neutron number index for residual nucleus
!   parskip         ! logical to skip outgoing particle
!   Zindex          ! charge number index for residual nucleus
!   ZZ              ! charge number of residual nucleus
! Variables for files
!   nulldev         ! null device
! Constants
!   cparity         ! parity (character)
!   nuc             ! symbol of nucleus
!   onethird        ! 1 / 3
!   parmass         ! mass of particle in a.m.u.
!   parname         ! name of particle
!   parspin         ! spin of particle
!   parZ            ! charge number of particle
! Variables for levels
!   edis            ! energy of level
!   jdis            ! spin of level
!   parlev          ! parity of level
! Variables for deformation parameters
!   colltype        ! type of collectivity (D, V or R)
!   defpar          ! deformation parameter
!   deftype         ! deformation length (D) or parameter (B)
!   indexlevel      ! level index
!   iphonon         ! phonon (1 or 2)
!   Kband           ! magnetic quantum number
!   lband           ! angular momentum
!   leveltype       ! type of level (rotational (R) or vibrational (V)
!   ndef            ! number of collective levels
!   nrot            ! number of deformation parameters for rotational
!   rotpar          ! deformation parameters for rotational nucleus
!   vibband         ! band number of level
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
!   flaginvecis     ! logical for calculating inverse channel OMP
!   hint            ! integration step size h
!   iband           ! band number of level
!   idvib           ! identifier for existence of vibrational state inside rotational
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
! Variables for JLM
!   rhojlmn         ! density for neutrons
!   rhojlmp         ! density for protons
! Variables for masses
!   nucmass         ! mass of nucleus
!   specmass        ! specific mass for residual nucleus
! Variables for optical model
!   flaompejec      ! flag for OMP for ejectile equal to projectile
!   av              ! real volume diffuseness
!   avd             ! real surface diffuseness
!   avso            ! real spin - orbit diffuseness
!   aw              ! imaginary volume diffuseness
!   awd             ! imaginary surface diffuseness
!   awso            ! imaginary spin - orbit diffuseness
!   d2              ! parameter for imaginary surface OMP
!   d3              ! parameter for imaginary surface OMP
!   disp            ! flag for dispersive optical model
!   ef              ! Fermi energy
!   rc              ! Coulomb radius
!   rv              ! real volume radius
!   rvd             ! real surface radius
!   rvso            ! real spin - orbit radius
!   rw              ! imaginary volume radius
!   rwd             ! imaginary surface radius
!   rwso            ! imaginary spin - orbit radius
!   v               ! real volume depth
!   vd              ! real surface depth
!   vso             ! real spin - orbit depth
!   w               ! imaginary volume depth
!   w2              ! parameter for imaginary volume OMP
!   wd              ! imaginary surface depth
!   wso             ! imaginary spin - orbit depth
!
! *** Declaration of local data
!
  implicit none
  logical            :: jlmloc         ! flag for JLM OMP
  logical            :: rotational     ! flag for rotational input
  logical            :: vibrational    ! flag for vibrational input
  character(len=3)   :: massstring !
  character(len=6)   :: finalnuclide !
  character(len=8)   :: ompfile !
  character(len=15)  :: col(20)    ! header
  character(len=15)  :: un(20)    ! header
  character(len=80)  :: quantity   ! quantity
  character(len=132) :: topline   ! topline
  character(len=132) :: outfile        ! output file
  integer            :: A              ! mass number of target nucleus
  integer            :: i              ! counter
  integer            :: i1             ! value
  integer            :: ii             ! counter
  integer            :: Ncomp          ! neutron number index for compound nucleus
  integer            :: nen            ! energy counter
  integer            :: Nix            ! neutron number index for residual nucleus
  integer            :: Ncol
  integer            :: type           ! particle type
  integer            :: Z              ! charge number of target nucleus
  integer            :: Zcomp          ! proton number index for compound nucleus
  integer            :: Zix            ! charge number index for residual nucleus
  real(sgl)          :: e              ! energy
!
! ********************** Set ECIS input parameters *********************
!
! Specific ECIS flags:
!
  flaginvecis = .true.
  legendre = .false.
  hint = 0.
  rmatch = 0.
  anginc = 180.
  angend = 180.
!
! Loop over all particle types and energies on standard energy grid.
!
  if (flagoutomp) then
    write(*, '(/" ######### OPTICAL MODEL PARAMETERS ##########",/)')
    if (jlmexist(Zcomp, Ncomp, 1) .or. jlmexist(Zcomp, Ncomp, 6)) then
      un=''
      col(1)='radius'
      un(1)='fm'
      col(2)='protons'
      col(3)='neutrons'
      Ncol=3
      quantity='radial density'
      topline=trim(targetnuclide)//' '//trim(quantity)//' for OMP'
      open (unit=1, file='radial.out', status='unknown')
      call write_header(topline,source,user,date,oformat)
      call write_target
      call write_datablock(quantity,Ncol,numjlm,col,un)
      do i = 1, numjlm
        write(1, '(3es15.6)') 0.1*real(i), rhojlmp(Zcomp, Ncomp, i, 1), rhojlmn(Zcomp, Ncomp, i, 1)
      enddo
      close(unit =1)
      call write_outfile('radial.out',flagoutall)
    endif
  endif
  flagecisinp = .false.
  flagompejec = .true.
  if (flageciscalc) open (unit = 9, file = 'ecisinv.inp', status = 'unknown')
  do type = 1, 6
    if (parskip(type)) cycle
    Zix = Zindex(Zcomp, Ncomp, type)
    Nix = Nindex(Zcomp, Ncomp, type)
    Z = ZZ(Zcomp, Ncomp, type)
    A = AA(Zcomp, Ncomp, type)
    if (type == 1) then
      angbeg = 0.
    else
      angbeg = 0.00001
    endif
    if (jlmexist(Zix, Nix, type) .and. (colltype(Zix, Nix) == 'S' .or. .not. flagrot(type))) then
      jlmloc = .true.
    else
      jlmloc = .false.
    endif
!
! Output of optical model parameters, if requested.
!
    if (flagoutomp .and. .not. jlmloc) then
      un='fm'
      col(1)='E'
      un(1)='MeV'
      col(2)='V'
      un(2)='MeV'
      col(3)='rv'
      col(4)='av'
      col(5)='W'
      un(5)='MeV'
      col(6)='rw'
      col(7)='aw'
      col(8)='Vd'
      un(8)='MeV'
      col(9)='rvd'
      col(10)='avd'
      col(11)='Wd'
      un(11)='MeV'
      col(12)='rwd'
      col(13)='awd'
      col(14)='Vso'
      un(14)='MeV'
      col(15)='rvso'
      col(16)='avso'
      col(17)='Wso'
      un(17)='MeV'
      col(18)='rwso'
      col(19)='awso'
      col(20)='rc'
      Ncol=20
      quantity='optical model parameters'
      massstring='   '
      write(massstring,'(i3)') A
      finalnuclide=trim(nuc(Z))//adjustl(massstring)
      topline=trim(finalnuclide)//' '//parname(type)//' '//trim(quantity)
      ompfile='omppar.'//parsym(type)
      open (unit=1, file=ompfile, status='replace')
      call write_header(topline,source,user,date,oformat)
      call write_residual(Z,A,finalnuclide)
      write(1,'("# parameters:")') 
      call write_char(2,'particle',parname(type))
      Nen =  eendmax(type) - ebegin(type) + 1
      call write_datablock(quantity,Ncol,Nen,col,un)
    endif
!
! Standard ECIS inputs for phenomenological optical potentials
!
! Some input flags for ECIS are energy dependent for the rotational model so ecis1 will be defined inside the energy loop.
!
    ecis1 = 'FFFFFTFFFFFFFFFFFFFFFFFFTFFTFFFFFFFFFFFFFFFFFFFFFF'
    ecis2 = 'FFFFFFFFTFFFTFFFTTTFTTTFTFFFFFFFFFFFFFFFFFFFTFFFFF'
    Nband = 0
!
! 1. Spherical nucleus
!
    jlmloc = .false.
    if (colltype(Zix, Nix) == 'S' .or. .not. flagrot(type)) then
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
      if (jlmexist(Zix, Nix, type)) then
        ecis1(7:7) = 'T'
        ecis1(15:15) = 'T'
        ecis1(29:29) = 'T'
        ecis1(41:41) = 'T'
        hint = 0.1
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
    if (disp(Zix, Nix, type)) then
      ecis1(10:10) = 'T'
      efer = ef(Zix, Nix, type)
      w2disp = w2(Zix, Nix, type)
      d3disp = d3(Zix, Nix, type)
      d2disp = d2(Zix, Nix, type)
    endif
    projmass = parmass(type)
    spin = parspin(type)
    resmass = nucmass(Zix, Nix)
    prodZ = real(Z * parZ(type))
    do nen = ebegin(type), eendmax(type)
      e = real(egrid(nen) / specmass(Zix, Nix, type))
!
! We use a simple formula to estimate the required number of j-values:
!    njmax=2.4*k*R;  R=1.25*A**1/3 ; k=0.22*sqrt(m(in amu)E(in MeV))
! and we always take a minimum of njmax=20.
!
      njmax = int(2.4 * 1.25 * (real(A) **onethird) * 0.22 * sqrt(projmass * e))
      njmax = max(njmax, 20)
      njmax = min(njmax, numl - 2)
      if (jlmloc) njmax = 1600
!
! *************** Calculate optical potential parameters ***************
!
! optical: subroutine for determination of optical potential
!
      call optical(Zix, Nix, type, e)
      if (flagoutomp .and. .not. jlmloc) then
        write(1, '(20es15.6)') &
 &        e, v, rv, av, w, rw, aw, vd, rvd, avd, wd, rwd, awd, vso, rvso, avso, wso, rwso, awso, rc
      endif
      if ( .not. flageciscalc) cycle
!
! ******************* Write ECIS input file ****************************
!
! ecisinput: subroutine to create ECIS input file
!
! For rotational nuclei, the switch at soswitch MeV needs to be made according to Pascal Romain.
!
      if (colltype(Zix, Nix) == 'R' .and. flagrot(type)) then
        if ((type == 1 .and. e <= soswitch) .or. (type > 1 .and. e <= coulbar(type))) then
          ecis1(13:13) = 'F'
          ecis1(21:21) = 'T'
          ecis1(42:42) = 'T'
        else
          ecis1(13:13) = 'T'
          ecis1(21:21) = 'F'
          ecis1(42:42) = 'F'
        endif
        if (type > 1 .and. e <= 0.05 * coulbar(type) .and. e <= 2. * Elevel(ncoll)) e = 0.1 * Elevel(ncoll)
        if (flagrel) ecis1(8:8) = 'T'
        if (disp(Zix, Nix, type)) ecis1(10:10) = 'T'
      endif
      flagecisinp = .true.
      call ecisinput(Zix, Nix, type, e, rotational, vibrational, jlmloc)
    enddo
    if (flagoutomp .and. .not. jlmloc) then
      close(unit =1)
      call write_outfile(ompfile,flagoutall)
    endif
  enddo
  flaginvecis = .false.
  flagompejec = .false.
  if ( .not. flageciscalc) return
  if ( .not. flagecisinp) then
    close (unit = 9)
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
    outfile = 'ecisinv.out  '
  else
    outfile = nulldev
  endif
  call ecist('ecisinv.inp  ', outfile, csfile, 'ecis.invin   ', transfile, 'null         ', 'null         ')
  invexist(Zcomp, Ncomp) = .true.
  open (unit = 9, file = 'ecisinv.inp', status = 'unknown')
  close (unit = 9, status = ecisstatus)
  return
end subroutine inverseecis
! Copyright A.J. Koning 2021
