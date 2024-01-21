subroutine read_energies
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input energy variables
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
!   sgl               ! single precision kind
! All global variables
!   numenin        ! number of incident energies
!   numJ           ! maximum J - value
!   numpop         ! number of population bins
! Variables for input energies
!   deninc      ! incident energy increment
!   flaginitpop ! flag for initial population distribution
!   energyfile  ! file with energies for OMP calculation
!   nin         ! counter for incident energy
!   npopE       ! number of energies for population distribution
!   npopJ       ! number of spins for population distribution
!   npopP       ! number of parities for population distribution
!   Ninc        ! number of incident energies
!   EdistE      ! excitation energy of population distribution
!   eninc       ! incident energy in MeV
!   enincmax    ! maximum incident energy
!   enincmin    ! minimum incident energy
!   Estop       ! incident energy above which TALYS stops
!   PdistE      ! population distribution, spin-independent
!   PdistJP     ! population distribution per spin and parity
! Variables for main input
!   ptype0          ! type of incident particle
! Constants
!   Emaxtalys    ! maximum acceptable energy for TALYS
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   range_real_error    ! Test if real variable is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical      :: fexist    ! flag for energy grid
  logical      :: lexist    ! logical to determine existence
  integer      :: i         ! counter
  integer      :: istat     ! logical for file access
  integer      :: J         ! spin of level
  integer      :: k         ! counter
  integer      :: negrid    ! number of grid points
  integer      :: nen       ! energy counter
  integer      :: parity    ! parity
  integer      :: pbeg      ! help variable
  real(sgl)    :: E         ! incident energy
  real(sgl)    :: Ein       ! incident energy
  real(sgl)    :: etmp      ! help variable
!
! It is possible to define a population distribution as the initial state, through the keyword projectile 0.
! In that case, we assign a photon projectile.
!
  if (ptype0 == '0') then
    flaginitpop = .true.
  else
    flaginitpop = .false.
  endif
  npopE = 0
  npopJ = 0
  npopP = 1
  EdistE = 0.
  PdistE = 0.
  PdistJP = 0.
  Ninc = 0
  enincmax = 0.
  enincmin = 0.
!
! ******** Set energy grids ************8
!
  Ein = 0.
  fexist = .false.
  if (flaginitpop) then
!
! 1. Population distribution as the initial state
!
    inquire (file = energyfile, exist = lexist)
    if ( .not. lexist) then
      write(*, '(" TALYS-error: if projectile 0, specify a range of excitation energies in a file ", a)') trim(energyfile)
      stop
    endif
    open (unit = 2, file = energyfile, status = 'old', iostat = istat)
    if (istat /= 0) call read_error(energyfile, istat)
    read(2, * , iostat = istat) npopE, npopJ, npopP
    if (istat /= 0) call read_error(energyfile, istat)
    call range_integer_error(energyfile, npopE, 1, numpop)
    call range_integer_error(energyfile, npopJ, 0, numJ + 1)
    if (npopJ > 0) call range_integer_error(energyfile, npopP, 1, 2)
!
! Only excitation energy distribution (no spins)
!
    if (npopJ == 0) then
      do nen = 1, npopE
        read(2, * , iostat = istat) EdistE(nen), PdistE(nen)
        if (istat /= 0) call read_error(energyfile, istat)
      enddo
    else
!
! Spin-dependent excitation energy distribution (no total)
!
      if (npopP == 1) then
        pbeg = 1
      else
        pbeg = - 1
      endif
      do parity = pbeg, 1, 2
        do nen = 1, npopE
          read(2, * , iostat = istat) EdistE(nen), (PdistJP(nen, J, parity), J = 0, npopJ - 1)
          if (istat /= 0) call read_error(energyfile, istat)
        enddo
      enddo
      do nen = 1, npopE
        if (EdistE(nen) <= EdistE(nen - 1)) then
          if (EdistE(1) == 0.) cycle
          write(*, '(" TALYS-error: excitation energies must be given in ascending order, or the number", &
 &          " of population bins is not correct")')
          stop
        endif
      enddo
      if (npopP == 1) then
        do nen = 1, npopE
          do J = 0, npopJ - 1
            PdistJP(nen, J, 1) = 0.5 * PdistJP(nen, J, 1)
            PdistJP(nen, J, -1) = PdistJP(nen, J, 1)
          enddo
        enddo
      endif
    endif
    close (unit = 2)
    Ninc = 1
    eninc(1) = EdistE(npopE)
    enincmin = eninc(1)
    enincmax = eninc(1)
  else
!
! 2. Normal case of a projectile with a target nucleus
!
! A. If no incident energy is given in the input file, incident energies should be read from a file.
!
    eninc(0) = 0.
    if (eninc(1) == 0.) then
      inquire (file = energyfile, exist = lexist)
      if ( .not. lexist) then
        call incidentgrid(energyfile(1:14), fexist)
        if (fexist) then
          goto 300
        else
          write(*, '(" TALYS-error: give a single incident energy in the input file using the energy keyword ")')
          write(*, '(14x, "or specify a range of incident energies in a file ")')
          write(*, '(14x, "or give a correct name for a pre-defined energy grid ", a)') trim(energyfile)
          stop
        endif
      endif
      nen = 0
      open (unit = 2, file = energyfile, status = 'old')
      do
        read(2, * , iostat = istat) Ein
        if (istat == -1) exit
        if (istat /= 0) call read_error(energyfile, istat)
        if (Ein /= 0.) then
          nen = nen + 1
!
! There is a maximum number of incident energies
!
          call range_integer_error(energyfile, nen, 0, numenin)
          eninc(nen) = Ein
        endif
      enddo
      close (unit = 2)
      if (nen == 0) then
        write(*, '(" TALYS-error: there are no incident energies in file ", a)') trim(energyfile)
        stop
      endif
!
! Sort incident energies in ascending order and remove double points
!
      do i = 1, nen
        do k = 1, i
          if (eninc(i) >= eninc(k)) cycle
          etmp = eninc(i)
          eninc(i) = eninc(k)
          eninc(k) = etmp
        enddo
      enddo
      Ninc = nen
      do i = 1, nen - 1
        if (eninc(i) == eninc(i + 1)) then
          do k = i + 1, nen
            eninc(k) = eninc(k + 1)
          enddo
          Ninc = Ninc - 1
        endif
      enddo
!
! The minimum and maximum value of all the incident energies is determined.
!
      enincmin = eninc(1)
      enincmax = eninc(1)
      do nen = 2, Ninc
        enincmin = min(enincmin, eninc(nen))
        enincmax = max(enincmax, eninc(nen))
      enddo
    else
      if (enincF == 0.) then
!
! B1. Single value given in the user input file
!
        Ninc = 1
        enincmin = eninc(1)
        enincmax = eninc(1)
      else
!
! B2. Energy grid based on input values E1, E2, dE
!
        call range_real_error('final incident energy', eninc(1), 0., enincF)
        call range_real_error('energy increment', deninc, 0., Emaxtalys)
        if (deninc == 0.) then
          negrid = 10
          deninc = (enincF - eninc(1)) / (negrid - 1)
        endif
        nen = 1
        do
          nen = nen + 1
          call range_integer_error('number of incident energies', nen, 0, numenin)
          E = eninc(nen - 1) + deninc
          if (E < enincF - 1.e-4) then
            eninc(nen) = E
          else
            eninc(nen) = enincF
            Ninc = nen
            enincmin = eninc(1)
            enincmax = eninc(nen)
            exit
          endif
        enddo
      endif
    endif
!
! Remove incident energies above energy given by Estop
!
    nen = Ninc
    do i = 1, nen
      if (eninc(i) > Estop) then
        Ninc = i - 1
        exit
      endif
    enddo
  endif
!
! In case of built-in energy range, write an explicit 'energies' file
!
300 if (enincF > 0. .or. fexist) then
    open (unit = 2, file = 'energies', status = 'replace')
    do nen = 1, Ninc
      write(2, '(1p, g12.4)') eninc(nen)
    enddo
    close (unit = 2)
  endif
  return
end subroutine read_energies
! Copyright A.J. Koning 2021
