subroutine omppar(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Optical model parameters
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
  use A1_error_handling_mod
!
! Definition of single and double precision variables
!   sgl             ! single precision kind
! All global variables
!   numen           ! maximum number of outgoing energies
!   numNph          ! maximum number of neutrons away from the initial compound nucleus
!   numomp          ! maximum number of lines in optical model file
!   numZph          ! maximum number of protons away from the initial compound nucleus
! Variables for OMP
!   flagriplomp     ! flag for RIPL OMP
! Variables for OMP
!   Ejoin           ! joining energy for high energy OMP
!   flagdisp        ! flag for dispersive optical model
!   flagjlm         ! flag for using semi - microscopic JLM OMP
!   flaglocalomp    ! flag for local (y) or global (n) optical model
!   flagriplrisk    ! flag for going outside RIPL mass validity range
!   optmod          ! file with optical model parameters
!   optmodfileN     ! optical model parameter file for neutrons
!   optmodfileP     ! optical model parameter file for protons
!   riplomp         ! RIPL OMP
! Variables for input energies
!   enincmax        ! maximum incident energy
! Variables for nuclides
!   AA              ! mass number of residual nucleus
! Variables for nuclides
!   Nindex          ! neutron number index for residual nucleus
!   NN              ! neutron number of residual nucleus
!   Zindex          ! charge number index for residual nucleus
!   ZZ              ! charge number of residual nucleus
! Constants
!   nuc             ! symbol of nucleus
!   onethird        ! 1 / 3
!   parsym          ! symbol of particle
!   parN            ! neutron number of particle
!   parZ            ! charge number of particle
!   twothird        ! 2 / 3
! Variables for files
!   path            ! directory containing files to be read
! Variables for deformation parameters
!   colltype        ! type of collectivity (D, V or R)
! Variables for optical model
!   av0             ! diffuseness for real volume OMP
!   avd0            ! diffuseness for surface OMP
!   avso0           ! diffuseness for real spin - orbit OMP
!   d1              ! parameter for imaginary surface OMP
!   d2              ! parameter for imaginary surface OMP
!   d3              ! parameter for imaginary surface OMP
!   disp            ! flag for dispersive optical model
!   ef              ! Fermi energy
!   eomp            ! energies on optical model file
!   Eompbeg0        ! upper energy of KD03 OMP
!   Eompbeg1        ! lower energy of alternative OMP
!   Eompend0        ! lower energy of KD03 OMP
!   Eompend1        ! upper energy of alternative
!   ompglobal       ! flag for use of global optical model
!   omplines        ! number of lines in optical model file
!   rc0             ! Coulomb radius
!   rv0             ! radius for real volume OMP
!   rvd0            ! radius for surface OMP
!   rvso0           ! radius for real spin - orbit OMP
!   v               ! real volume depth
!   V0              ! V at zero MeV
!   v1              ! parameter for real volume OMP
!   v2              ! parameter for real volume OMP
!   v3              ! parameter for real volume OMP
!   v4              ! parameter for real volume OMP
!   Vjoin           ! V at joining energy
!   vomp            ! optical model parameters from file
!   vso1            ! parameter for real spin - orbit OMP
!   vso2            ! parameter for real spin - orbit OMP
!   w               ! imaginary volume depth
!   w1              ! parameter for imaginary volume OMP
!   w2              ! parameter for imaginary volume OMP
!   w3              ! parameter for imaginary volume OMP
!   w4              ! parameter for imaginary volume OMP
!   Wjoin           ! W at joining energy
!   wso1            ! parameter for imaginary spin - orbit OMP
!   wso2            ! parameter for imaginary spin - orbit OMP
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist        ! logical to determine existence
  character(len=1)  :: omptype       ! type of optical model (spherical or coupled)
  character(len=1)  :: Rpart         !
  character(len=8)  :: ompchar       ! help variable
  character(len=132):: optmodfile    ! file with optical model parameters
  character(len=132):: ompfile       ! optical model parameter file
  character(len=132):: string        ! string
  integer           :: A             ! mass number of target nucleus
  integer           :: Abeg          ! first A to be included
  integer           :: Aend          ! last A to be included
  integer           :: At            !
  integer           :: i             ! counter
  integer           :: j             ! counter
  integer           :: kk            ! counter
  integer           :: ia            ! mass number from abundance table
  integer           :: iz            ! charge number
  integer           :: ii            ! counter
  integer           :: istat         ! logical for file access
  integer           :: k             ! designator for particle
  integer           :: N             ! neutron number of residual nucleus
  integer           :: NE            ! number of incident energies
  integer           :: nen           ! energy counter
  integer           :: Nix           ! neutron number index for residual nucleus
  integer           :: nomp          ! number of particles for optical model parameters
  integer           :: Rnum          !
  integer           :: Z             ! charge number of target nucleus
  integer           :: Zbeg          ! first Z to be included
  integer           :: Zend          ! last Z to be included
  integer           :: Zix           ! charge number index for residual nucleus
  integer           :: Zt            !
  real(sgl)         :: dEripl        !
  real(sgl)         :: dum           ! dummy value
  real(sgl)         :: e             ! energy
  real(sgl)         :: eta(2)        ! sign of asymmetry term
  real(sgl)         :: Ebeg          ! first energy of experimental energy grid
  real(sgl)         :: Eeps          ! help variable
  real(sgl)         :: eferm         ! Fermi energy
  real(sgl)         :: Efin          !
  real(sgl)         :: Er            ! resonance energy in LAB system
  real(sgl)         :: Eripl(numen)  !
!
! ************** Read optical model parameters from database ***********
!
!
  Z = ZZ(Zix, Nix, 0)
  N = NN(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  omptype = ' '
  if (flaglocalomp) then
    do k = 1, 2
      ompchar = parsym(k)//'-'//trim(nuc(Z))//'.omp'
      if (k == 1) then
        if (optmodfileN(Zix)(1:1) /= ' ') then
          ompfile = optmodfileN(Zix)
        else
          ompfile = trim(path)//'optical/neutron/'//ompchar
        endif
      else
        if (optmodfileP(Zix)(1:1) /= ' ') then
          ompfile = optmodfileP(Zix)
        else
          ompfile = trim(path)//'optical/proton/'//ompchar
        endif
      endif
      omptype = ' '
      inquire (file = ompfile, exist = lexist)
      if ( .not. lexist) cycle
      open (unit = 2, file = ompfile, status = 'old', iostat = istat)
      if (istat /= 0) call read_error(ompfile, istat)
Loop1: do
        read(2, '(4x, 2i4, 3x, a1)', iostat = istat) ia, nomp, omptype
        if (istat == -1) exit
        if (A == ia) then
          omptype = ' '
          do i = 1, nomp
            read(2, '(4x, f7.2, f8.3)') ef(Zix, Nix, k), rc0(Zix, Nix, k)
            read(2, '(2f8.3, f6.1, f10.4, f9.6, f6.1, f7.1)') rv0(Zix, Nix, k), &
              av0(Zix, Nix, k), v1(Zix, Nix, k), v2(Zix, Nix, k), v3(Zix, Nix, k), w1(Zix, Nix, k), w2(Zix, Nix, k)
            read(2, '(2f8.3, f6.1, f10.4, f7.2)') rvd0(Zix, Nix, k), &
              avd0(Zix, Nix, k), d1(Zix, Nix, k), d2(Zix, Nix, k), d3(Zix, Nix, k)
            read(2, '(2f8.3, f6.1, f10.4, f6.1, f7.1)') rvso0(Zix, Nix, k), &
              avso0(Zix, Nix, k), vso1(Zix, Nix, k), vso2(Zix, Nix, k), wso1(Zix, Nix, k), wso2(Zix, Nix, k)
            v4(Zix, Nix, k) = 7.e-9
            if (flagdisp) then
              if (nomp == 1 .or. flagjlm) then
                disp(Zix, Nix, k) = .false.
              else
                disp(Zix, Nix, k) = .true.
              endif
            else
              exit Loop1
            endif
          enddo
        else
          do i = 1, nomp
            do ii = 1, 4
              read(2, '()')
            enddo
          enddo
        endif
      enddo Loop1
      close (unit = 2)
!
! Reduce d1 parameter (of Wd) in case of coupled-channels, unless already specified in OMP parameterization.
!
      if (colltype(Zix, Nix) /= 'S' .and. omptype /= 'C') d1(Zix, Nix, k) = 0.85 * d1(Zix, Nix, k)
    enddo
  endif
!
! ************************* Global optical model ***********************
!
! KD03-Delaroche OMP, 
! A.J. Koning and J.P. Delaroche, 
! Local and global nucleon optical models from 1 keV to 200 MeV, Nucl. Phys. A713 (2003) 231.
!
! Option to use global OMP parameters from Uncertainty-quantified phenomenological optical potentials for 
! single-nucleon scattering, C. D. Pruitt, J. E. Escher, and R. Rahman, Phys. Rev. C 107, 014602 
!
! Test if local OMP has been assigned.
!
  eta(1) = -1.
  eta(2) = 1.
  do k= 1, 2
    if (rv0(Zix, Nix, k) == 0.) then
      if (flagdisp) disp(Zix, Nix, k) = .false.
      ompglobal(Zix, Nix, k) = .true.
      if (k == 1) ef(Zix, Nix, 1) = - 11.2814 + 0.02646 * A
      if (k == 2) ef(Zix, Nix, 2) = - 8.4075 + 0.01378 * A
      call kd03(k,pruitt)
      rv0(Zix, Nix, k) = rv_0 - rv_A * A **( - onethird)
      av0(Zix, Nix, k) = av_0 - av_A * A
      v1(Zix, Nix, k) = v1_0 + eta(k) * v1_asymm * real(N - Z) / A - v1_A * A
      v2(Zix, Nix, k) = v2_0 + eta(k) * v2_A * A
      v3(Zix, Nix, k) = v3_0 + eta(k) * v3_A * A
      v4(Zix, Nix, k) = v4_0
      w1(Zix, Nix, k) = w1_0 + w1_A * A
      w2(Zix, Nix, k) = w2_0 + w2_A * A
      rvd0(Zix, Nix, k) = rd_0 - rd_A * A **onethird
      avd0(Zix, Nix, k) = ad_0 + eta(k) * ad_A * A
      d1(Zix, Nix, k) = d1_0 + eta(k) * d1_asymm * real(N - Z) / A
      d2(Zix, Nix, k) = d2_0 + d2_A / (1. + exp((A - d2_A3) / d2_A2))
      d3(Zix, Nix, k) = d3_0
      rvso0(Zix, Nix, k) = rso_0 - rso_A * A **( - onethird)
      avso0(Zix, Nix, k) = aso_0
      vso1(Zix, Nix, k) = vso1_0 + vso1_A * A
      vso2(Zix, Nix, k) = vso2_0
      wso1(Zix, Nix, k) = wso1_0
      wso2(Zix, Nix, k) = wso2_0
    endif
!
! Reduce d1 parameter (of Wd) in case of coupled-channels, unless already specified in OMP parameterization.
!
    if (colltype(Zix, Nix) /= 'S' .and. omptype /= 'C') d1(Zix, Nix, k) = 0.85 * d1(Zix, Nix, k)
  enddo
  rc0(Zix, Nix, 1) = 0.
  rc0(Zix, Nix, 2) = rc_0 + rc_A * A **( - twothird) + rc_A2 * A **( - 5. / 3.)
!
! ************** Optical model parameters from RIPL ********************
!
  if (Zix <= numZph .and. Nix <= numNph) then
    if (flagriplomp) then
!
! Copy RIPL parameter files to working directory
!
      ompfile = trim(path)//'optical/ripl/om-parameter-u.dat'
      open (unit = 41, file = ompfile, status = 'unknown')
      open (unit = 42, file = 'om-parameter-u.dat', status = 'unknown')
      do
        read(41, '(a)', iostat = istat) string
        if (istat ==  -1) exit
        write(42, '(a)') trim(string)
      enddo
      close (unit = 42)
      close (unit = 41)
      ompfile = trim(path)//'optical/ripl/gs-mass-sp.dat'
      open (unit = 41, file = ompfile, status = 'unknown')
      open (unit = 42, file = 'gs-mass-sp.dat', status = 'unknown')
      do
        read(41, '(a)', iostat = istat) string
        if (istat ==  -1) exit
        write(42, '(a)') trim(string)
      enddo
      close (unit = 42)
      close (unit = 41)
!
! Set energy grid for OMP table
!
      Er = 0.
      nen = 0
      dEripl = 0.001
      do
        Er = Er + dEripl
        Eeps = Er + 1.e-4
        if (Er > enincmax + 12.) exit
        if (nen == numen) exit
        nen = nen + 1
        Eripl(nen) = Er
        if (Eeps > 0.01) dEripl = 0.01
        if (Eeps > 0.1) dEripl = 0.1
        if (Eeps > 4.) dEripl = 0.2
        if (Eeps > 10.) dEripl = 0.5
        if (Eeps > 30.) dEripl = 1.
        if (Eeps > 100.) dEripl = 2.
      enddo
      NE = nen
      do k = 1, 6
        if (riplomp(k) > 0 .and. Zix == parZ(k) .and. Nix == parN(k) .and. optmod(Zix, Nix, k) == ' ') then
!
! Test if RIPL number is correct and inside mass and energy ranges
!
          Zt = ZZ(0, 0, k)
          At = AA(0, 0, k)
          ompfile = trim(path)//'optical/ripl/om-index.txt'
          open (unit = 31, file = ompfile, status = 'old')
          read(31, '(//)')
          do
            read(31, '(i4, 3x, a1, 23x, i2, 1x, i2, 2x, i3, 1x, i3, 1x, f5.1, 1x, f5.1)', iostat = istat) &
 &            Rnum, Rpart, Zbeg, Zend, Abeg, Aend, Ebeg, Efin
            if (istat ==  -1) exit
            if (riplomp(k) == Rnum) then
              if (parsym(k) == Rpart .and. Zt >= Zbeg .and. Zt <= Zend .and. At >= Abeg .and. At <= Aend) then
                Eompbeg1(k, 1) = Ebeg
                Eompend1(k, 1) = Efin
                Eompbeg0(k, 1) = Ebeg * 0.8
                Eompend0(k, 1) = Efin * 1.2
                close (unit = 31)
                goto 190
              endif
            endif
          enddo
          close (unit = 31)
          write(*, '(" TALYS-warning: RIPL OMP ",i4," for ",i3,a2," out of range")') riplomp(k), At, nuc(Zt)
          if ( .not. flagriplrisk) stop
!
! Create input file for RIPL OMP retrieval code
!
  190     open (unit = 42, file = 'ominput.inp', status = 'unknown')
          write(42, * ) NE
          write(42, * ) (Eripl(i), i = 1, NE)
          write(42, * ) Zt, At, riplomp(k), - 2
          close (unit = 42)
!
! Run RIPL OMP retrieval code
!
          write( * , * ) "RIPL OMP for ", Zt, At
          call riplomp_mod(NE, Zt, At, Eripl, riplomp(k), - 2, 0)
          close (unit = 35)
!
! Retrieve  RIPL OMP parameters on energy grid
!
          e = 0.
          open(42, file = 'omp-table.dat', status = 'unknown')
  192     read(42, '(a)') string
          kk = index(string, 'Ef =')
          if (kk > 0) read(string(kk + 4:132), * ) ef(Zix, Nix, k)
          read(string(1:7), * , err = 192, end = 192) e
          if (e > 0.) then
            backspace 42
            do i = 1, NE
              read(42, * ) eomp(Zix, Nix, k, i), (vomp(Zix, Nix, k, i, j), j = 1, 3), dum, (vomp(Zix, Nix, k, i, j), j = 4, 18)
            enddo
            close(42)
            omplines(Zix, Nix, k) = NE
          else
            goto 192
          endif
        endif
      enddo
    endif
  endif
!
! ************** Optical model file from user input file ***************
!
  if (Zix <= numZph .and. Nix <= numNph) then
    do k = 1, 6
      optmodfile = optmod(Zix, Nix, k)
      if (optmodfile(1:1) /= ' ') then
        open (unit = 2, file = optmodfile, status = 'old', iostat = istat)
        if (istat /= 0) call read_error(optmodfile, istat)
        eferm = 0.
        read(2, '(a)', iostat = istat) string
        if (istat > 0) call read_error(optmodfile, istat)
        read(string, *, iostat = istat) iz, ia, omplines(Zix, Nix, k), eferm
        if (istat > 0) call read_error(optmodfile, istat)
        if (eferm > 0.) ef(Zix, Nix, k) = eferm
        if (omplines(Zix, Nix, k) > numomp) then
          write(*, '(" TALYS-error: number of lines in OMP file > ", i4)') numomp
          stop
        endif
        call range_integer_error(optmodfile, omplines(Zix, Nix, k), 0, numomp, index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', &
 &        index3 = k, name3 = 'k')
        eomp(Zix, Nix, k, 0) = 0.
        do nen = 1, omplines(Zix, Nix, k)
          read(2, * , iostat = istat) eomp(Zix, Nix, k, nen), (vomp(Zix, Nix, k, nen, ii), ii = 1, 19)
          if (istat == -1) exit
          if (istat /= 0) call read_error(optmodfile, istat)
        enddo
        close (unit = 2)
        if (omplines(Zix, Nix, k) == 1) then
          omplines(Zix, Nix, k) = 2
          eomp(Zix, Nix, k, 2) = eomp(Zix, Nix, k, 1)
          do ii = 1, 19
            vomp(Zix, Nix, k, 2, ii) = vomp(Zix, Nix, k, 1, ii)
          enddo
        endif
      endif
    enddo
  endif
!
! Set OMP parameters for extension up to 1 GeV
!
! optical : subroutine for determination of optical potential
!
  do k = 1, 2
    w3(Zix, Nix, k) = 25. - 0.0417 * A
    w4(Zix, Nix, k) = 250.
    if (Zix == Zindex(0, 0, k) .and. Nix == Nindex(0, 0, k) .and. enincmax > Ejoin(k)) then
      e = 0.
      call optical(Zix, Nix, k, e)
      V0(k) = v
      e = Ejoin(k)
      call optical(Zix, Nix, k, e)
      Vjoin(k) = v
      Wjoin(k) = w
    endif
  enddo
!
! Set energy ranges of alternative optical models
!
! Eompbeg0 <=  Eompbeg1 <=  Eompend1 <=  Eompend0
!
! Deuteron OMPs
!
  do i = 2, 5
    Eompbeg0(3, i) = 0.
    Eompbeg1(3, i) = 0.
    Eompend1(3, i) = 200.
    Eompend0(3, i) = 300.
  enddo
  Eompend1(3, 2) = 90.
  Eompend0(3, 2) = 150.
  Eompend1(3, 3) = 100.
  Eompend0(3, 3) = 150.
!
! Alpha OMPs
!
  do i = 2, 8
    Eompbeg0(6, i) = 0.
    Eompbeg1(6, i) = 0.
    Eompend1(6, i) = 200.
    Eompend0(6, i) = 300.
  enddo
  Eompend1(6, 2) = 25.
  Eompend0(6, 2) = 50.
  return
end subroutine omppar
! Copyright A.J. Koning 2021
