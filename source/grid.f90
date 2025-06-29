subroutine grid
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Energy and angle grid
!
! Author    : Arjan Koning and Stephane Hilaire
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
!   sgl              ! single precision kind
! All global variables
!   numen            ! maximum number of outgoing energies
!   numen6           ! number of energies for ENDF6 energy grid
!   numenlow         ! number of energies for low energy regime
!   numisom          ! number of isomers
!   nummt            ! number of MT numbers
! Variables for basic reaction
!   flagendf         ! flag for information for ENDF - 6 file
!   flagequispec     ! flag to use equidistant bins for emission spectra
!   flagendfecis     ! flag for new ECIS calculation for ENDF - 6 files
! Variables for best files
!   flagrescue       ! flag for final rescue: normalization to data
!   grescue          ! global multipl. factor for incident energy dep. adj. factors
!   rescuefile       ! file with incident energy dependent adjustment factors
! Variables for level density
!   D0               ! s - wave resonance spacing in eV
! Variables for numerics
!   nangle           ! number of angles
!   nanglecont       ! number of angles for continuum
!   segment          ! number of segments to divide emission energy grid
!   transpower       ! power for transmission coefficient limit
! Variables for direct reactions
!   flagecissave     ! flag for saving ECIS input and output files
!   flaginccalc      ! flag for new ECIS calculation for incident channel
! Variables for basic parameters
!   eninclow         ! minimal incident energy for nuclear model calculations
! Variables for input energies
!   eninc            ! incident energy in MeV
!   enincmax         ! maximum incident energy
!   Ninc           ! number of incident energies
! Variables for main input
!   Atarget          ! mass number of target nucleus
! Variables for energy grid
!   angle            ! angle in degrees
!   anglecont        ! angle in degrees for continuum
!   coullimit        ! energy limit for charged particle OMP calculation
!   Crescue          ! adjustment factor for this incident energy
!   deltaE           ! energy bin around outgoing energies
!   E1v              ! energy at end of 1 / v region
!   ebegin           ! first energy point of energy grid
!   Ebottom          ! bottom of outgoing energy bin
!   ecisstatus       ! status of ECIS file
!   eendmax          ! last energy point of energy grid
!   egrid            ! outgoing energy grid
!   Einc             ! incident energy in MeV
!   Erescue          ! energy grid for adjustment factors
!   Etop             ! top of outgoing energy bin
!   frescue          ! adjustment factor
!   maxen            ! total number of energies
!   Nrescue          ! number of energies for adjustment factors
!   translimit       ! limit for transmission coefficient
! Variables for energies
!   eend             ! last energy point of energy grid
!   Ninclow        ! number of incident energies below Elow
! Variables for nuclides
!   coulbar          ! Coulomb barrier
!   parskip          ! logical to skip outgoing particle
! Variables for resonance parameters
!   D0theo           ! mean s - wave resonance spacing
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical   :: lexist      ! logical to determine existence
  integer   :: iang        ! running variable for angle
  integer   :: is          ! isotope counter: -1=total, 0=ground state 1=isomer
  integer   :: istat       ! logical for file access
  integer   :: mhalf       ! help variable
  integer   :: mt          ! MT number
  integer   :: nen         ! energy counter
  integer   :: nen0        ! energy counter
  integer   :: type        ! particle type
  real(sgl) :: coulfactor  ! constant for Coulomb barrier
  real(sgl) :: dang        ! delta angle
  real(sgl) :: degrid      ! energy increment
  real(sgl) :: Eeps        ! help variable
  real(sgl) :: Elimit
  real(sgl) :: Eout        ! outgoing energy
  real(sgl) :: val         ! real value
!
! ************************ Basic outgoing energy grid ******************
!
! The basic outgoing energy grid we use is:
!
!   0.001, 0.002, 0.005    MeV
!   0.01,  0.02,  0.05     MeV
!   0.1-   2 MeV : dE= 0.1 MeV
!   2  -   4 MeV : dE= 0.2 MeV
!   4  -  20 MeV : dE= 0.5 MeV
!   20  -  40 MeV : dE= 1.0 MeV
!   40  - 200 MeV : dE= 2.0 MeV
!   200 - 300 MeV : dE= 5.0 MeV
!   above 300 MeV : dE=10.0 MeV
!
! This grid ensures that the calculation for outgoing energies of a few MeV (around the evaporation peak) is sufficiently precise,
! whereas at higher energies a somewhat coarser energy grid can be used.
! For the same reason, this energy grid is used for the calculation of transmission coefficients.
! The grid can be further subdivided with the segment keyword.
!
  Eout = 0.
  degrid = 0.001
  nen = 0
  Elimit = enincmax + S(0,0,k0) + targetE + 1.
  do
    Eout = Eout + degrid
    Eeps = Eout + 1.e-4
    if (Eeps > Elimit) exit
    if (nen == numen) exit
    nen = nen + 1
    egrid(nen) = Eout
    if (Eeps > 0.002) degrid = 0.003
    if (Eeps > 0.005) degrid = 0.005
    if (Eeps > 0.01) degrid = 0.01
    if (Eeps > 0.02) degrid = 0.03
    if (Eeps > 0.05) degrid = 0.05
    if (Eeps > 0.1) degrid = 0.1 / segment
    if (Eeps > 2.) degrid = 0.2 / segment
    if (Eeps > 4.) degrid = 0.5 / segment
    if (Eeps > 20.) degrid = 1. / segment
    if (Eeps > 40.) degrid = 2. / segment
    if (Eeps > 200.) degrid = 5. / segment
    if (Eeps > 300.) degrid = 10. / segment
  enddo
  maxen = nen
  call range_integer_error('number of energies', maxen, 1, numen)
!
! Equidistant energy grid for two spectrum regions e.g. for PFNS
!
      if (flagequispec) then
        maxen = numen - 2
        mhalf = numen / 2
        do nen = 1, mhalf
          egrid(nen) = 0.1 * nen
        enddo
        dEgrid = (enincmax + 12. - egrid(mhalf)) / mhalf
        do nen = mhalf+1, numen
          egrid(nen) = egrid(mhalf) + dEgrid * (nen-mhalf)
        enddo
      endif

!
! The widths of the bins around the emission energies are set.
!
  do nen = 2, maxen - 1
    deltaE(nen) = 0.5 * (egrid(nen + 1) - egrid(nen - 1))
    Etop(nen) = 0.5 * (egrid(nen) + egrid(nen + 1))
    Ebottom(nen) = 0.5 * (egrid(nen) + egrid(nen - 1))
  enddo
  deltaE(1) = egrid(1) + 0.5 * (egrid(2) - egrid(1))
  Etop(1) = deltaE(1)
  Ebottom(1) = 0.
  deltaE(0) = 0.
  deltaE(maxen) = 0.5 * (egrid(maxen) - egrid(maxen - 1))
  Etop(maxen) = egrid(maxen)
  Ebottom(maxen) = 0.5 * (egrid(maxen) + egrid(maxen - 1))
!
! ******************** Set lower limit for energy grid *****************
!
! For charged particles it is not necessary, or even numerically possible, to calculate transmission coefficients and cross sections
! for very low energies.
! Therefore, we relate their energy grids to the corresponding Coulomb barrier.
! The first outgoing energy is a factor of coulfactor lower than the Coulomb barrier.
! The last outgoing energy of the grid depends on the incident energy and is initialized later in the energies subroutine.
!
  ebegin(0) = 1
  ebegin(1) = 1
Loop1:  do type = 2, 6
    ebegin(type) = 0
    if (parskip(type)) cycle
    coulfactor = 0.01 * (1. + Atarget / 200.)
    coullimit(type) = coulfactor * coulbar(type)
    do nen = 1, maxen
      if (egrid(nen) > coullimit(type)) then
        ebegin(type) = nen
        cycle Loop1
      endif
    enddo
  enddo Loop1
!
! ********* Set upper limit for energy grid and other energies *********
!
! energies: subroutine for energies
!
  Einc = enincmax
  call energies
  do type = 0, 6
    if (parskip(type)) cycle
    eendmax(type) = eend(type)
  enddo
!
! *************** Determine number of low incident energies ************
!
! TALYS performs full nuclear reaction calculations for incident energies above Elow only.
! We keep track of the number of lower energies, for which simple empirical cross section estimates are made.
!
! locate   : subroutine to find value in ordered table
!
! For excitation functions that extend to very low incident energies, the energy at the end of the 1/v region and
! the energy at the end of the resonance region are also inserted.
!
  if (eninclow == 0.) then
    if (D0(0, 0) == 0.) then
      eninclow = min(D0theo(0, 0) * 1.e-6, 1.)
    else
      eninclow = min(D0(0, 0) * 1.e-6, 1.)
    endif
  endif
  if (Ninc >= numenlow - 2) eninclow = min(eninclow, eninc(numenlow - 2))
  eninclow = max(eninclow, 1.e-11)
  E1v = 0.2 * eninclow
  if (eninc(1) < E1v) then
    call locate(eninc, 1, Ninc, E1v, nen0)
    if (eninc(nen0) / E1v < 0.99) then
      Ninc = Ninc + 1
      do nen = Ninc, nen0 + 2, - 1
        eninc(nen) = eninc(nen - 1)
      enddo
      eninc(nen0 + 1) = E1v
    endif
    call locate(eninc, 1, Ninc, eninclow, nen0)
    if (eninc(nen0) / eninclow < 0.99) then
      Ninc = Ninc + 1
      do nen = Ninc, nen0 + 2, - 1
        eninc(nen) = eninc(nen - 1)
      enddo
      eninc(nen0 + 1) = eninclow
    endif
  endif
  Ninclow = 0
  do nen = 1, Ninc
    if (eninc(nen) < eninclow) Ninclow = Ninclow + 1
  enddo
!
! ************** Set limit for transmission coefficients ***************
!
  translimit = 1. / (10 **transpower)
!
! **************************** Basic angle grid ************************
!
! 1. Discrete angular distributions
!
  dang = 180. / nangle
  do iang = 0, nangle
    angle(iang) = iang * dang
  enddo
!
! 2. Continuum angular distributions
!
  dang = 180. / nanglecont
  do iang = 0, nanglecont
    anglecont(iang) = iang * dang
  enddo
!
! *** Open files with basic reaction information for incident channel **
!
  if (flagecissave) then
    ecisstatus = 'keep'
  else
    ecisstatus = 'delete'
  endif
  if ( .not. flaginccalc) then
    inquire (file = 'incident.cs', exist = lexist)
    if ( .not. lexist) then
      write(*, '(" TALYS-error: The first calculation of a run should always be done with ecissave y and inccalc y")')
      stop
    endif
    open (unit = 13, file = 'incident.cs', status = 'old', iostat = istat)
    if (istat /= 0) call read_error('incident.cs', istat)
    open (unit = 17, file = 'incident.tr', status = 'old', iostat = istat)
    if (istat /= 0) call read_error('incident.tr', istat)
    open (unit = 18, file = 'incident.ang', status = 'old', iostat = istat)
    if (istat /= 0) call read_error('incident.ang', istat)
    open (unit = 19, file = 'incident.leg', status = 'old', iostat = istat)
    if (istat /= 0) call read_error('incident.leg', istat)
    open (unit = 20, file = 'incident.in', status = 'old', iostat = istat)
    if (istat /= 0) call read_error('incident.in', istat)
  endif
  if (flagendf .and. .not. flagendfecis) then
    inquire (file = 'endf.cs', exist = lexist)
    if ( .not. lexist) then
      write(*, '(" TALYS-error: The first calculation of a run should always be done with ecissave y and endfecis y")')
      stop
    endif
    open (unit = 23, file = 'endf.cs', status = 'old', iostat = istat)
    if (istat /= 0) call read_error('endf.cs', istat)
  endif
!
! To fit data very precisely, a "normal" TALYS calculation may not suffice.
! Therefore, as a final "rescue" an option is included to read incident energy dependent adjustment factors for the most important
! cross sections.
! The MT number from the ENDF-6 format is used as index.
!
  if (flagrescue) then
    do mt = 1, nummt
      do is = - 1, numisom
        Crescue(mt, is) = 1.
        do nen = 1, numen6
          Erescue(mt, is, nen) = 0.
          frescue(mt, is, nen) = grescue(mt, is)
        enddo
        if (rescuefile(mt, is)(1:1) /= ' ') then
          open (unit = 2, file = rescuefile(mt, is), status = 'old', iostat = istat)
          if (istat /= 0) call read_error(rescuefile(mt, is), istat)
          nen = 1
          do
            read(2, * , iostat = istat) Erescue(mt, is, nen), val
            if (istat == -1) exit
            if (istat /= 0) call read_error(rescuefile(mt, is), istat)
            frescue(mt, is, nen) = grescue(mt, is) * val
            if (nen > 1) then
              if (Erescue(mt, is, nen) <= Erescue(mt, is, nen - 1)) then
                write(*, '(" TALYS-error: energies in rescuefile must be given in ascending order")')
                stop
              endif
            endif
            nen = nen + 1
            call range_integer_error(rescuefile(mt, is), nen, 0, numen6 + 1)
          enddo
          close (unit = 2)
          Nrescue(mt, is) = nen - 1
        endif
      enddo
    enddo
  endif
  return
end subroutine grid
! Copyright A.J. Koning 2021
