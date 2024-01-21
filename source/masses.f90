subroutine masses
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Nuclear masses
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
!   sgl            ! single precision kind
!   dbl            ! double precision kind
! All global wariables
!   numNastro      ! maximal number of neutrons away from initial CN for astroph. calcs
!   numZastro      ! maximal number of protons away from initial CN for astroph. calcs
! Variables for masses
!   beta2          ! deformation parameter
!   flagexpmass    ! flag for using experimental nuclear mass if available
!   massdir        ! directory with mass tables
!   massexcess     ! mass excess in MeV as read from user input file
!   massmodel      ! model for theoretical nuclear mass
!   massnucleus    ! mass of nucleus in amu as read from user input file
! Variables for numerics
!   maxN           ! maximal number of neutrons away from initial compound nucleus
!   maxZ           ! maximal number of protons away from initial compound nucleus
! Variables for main input
!   Ainit          ! mass number of initial compound nucleus
!   k0             ! index of incident particle
!   Ninit          ! neutron number of initial compound nucleus
!   Zinit          ! charge number of initial compound nucleus
! Variables for nuclides
!   parskip        ! logical to skip outgoing particle
! Variables for files
!   path           ! directory containing files to be read
! Constants
!   amu            ! atomic mass unit in MeV
!   nuc            ! symbol of nucleus
!   parmass        ! mass of particle in a.m.u.
!   parN           ! neutron number of particle
!   parZ           ! charge number of particle
! Variables for masses
!   beta4          ! deformation parameters
!   dumexc         ! theoretical mass excess from Duflo - Zuker formula
!   expmass        ! flag for using experimental nuclear mass if available
!   expmexc        ! experimental mass excess
!   gsparity       ! ground state parity
!   gsspin         ! ground state spin
!   nucmass        ! mass of nucleus
!   redumass       ! reduced mass
!   specmass       ! specific mass for residual nucleus
!   thmass         ! theoretical mass
!   thmexc         ! theoretical mass excess
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: flagduflo    ! flag to check whether Duflo-Zuker calculation is  required
  logical           :: lexist       ! logical to determine existence
  character(len=7)  :: masschar     ! help variable
  character(len=132):: massfile     ! mass file
  integer           :: A            ! mass number of target nucleus
  integer           :: Abegin       ! first A to be included
  integer           :: Aend         ! last A to be included
  integer           :: i            ! counter
  integer           :: ia           ! mass number from abundance table
  integer           :: istat        ! logical for file access
  integer           :: L            ! counter
  integer           :: N            ! neutron number of residual nucleus
  integer           :: Nbegin       ! first N to be included
  integer           :: Nend         ! maximal neutron number
  integer           :: Nix          ! neutron number index for residual nucleus
  integer           :: p            ! parity
  integer           :: type         ! particle type
  integer           :: Z            ! charge number of target nucleus
  integer           :: Zix          ! charge number index for residual nucleus
  real(sgl)         :: b2           ! beta2
  real(sgl)         :: b4           ! beta4
  real(sgl)         :: exc          ! mass excess
  real(sgl)         :: gs           ! ground state spin
  real(dbl)         :: expmass1     ! experimental mass
  real(dbl)         :: expmexc1     ! experimental mass excess
  real(dbl)         :: thmass1      ! theoretical mass
  real(dbl)         :: thmexc1      ! theoretical mass excess
!
! ************************ Read nuclear masses *************************
!
! We read both the experimental masses, from AME2020, and the theoretical masses from the masstable.
! The default option is to adopt the experimental nuclear mass, when available.
! We also read the experimental and theoretical mass excess, to enable a more precise calculation of separation energies.
! If we need separation energies and specific masses for nuclides up to (maxZ,maxN), we need nuclear masses up to (maxZ+4,maxN+4).
!
! There are also input options for theoretical mass models:
! massmodel 0: Duflo-Zuker
! massmodel 1: Moeller
! massmodel 2: Goriely HFB-Skyrme model
! massmodel 3: HFB-Gogny D1M model
! where, if massmodel 1, 2 or 3, massmodel 0 is used when no tabulated values are available.
! Also, with the input option expmass n  the use of experimental masses can be disabled.
!
  do Zix = 0, maxZ + 4
    Z = Zinit - Zix
    if (Z <= 0) cycle
    Nbegin = Ninit - maxN - 4
    Nend = Ninit
    Abegin = Z + Nbegin
    Aend = Z + Nend
    masschar = trim(nuc(Z))//'.mass'
    if (flagexpmass) then
      massfile = trim(path)//'masses/ame2020/'//masschar
      inquire (file = massfile, exist = lexist)
      if (lexist) then
        open (unit = 1, file = massfile, status = 'old', iostat = istat)
        if (istat /= 0) call read_error(massfile, istat)
        do
          read(1, '(4x, i4, 2f12.6)', iostat = istat) ia, expmass1, expmexc1
          if (istat == -1) exit
          if (istat /= 0) call read_error(massfile, istat)
          if (ia < Abegin) cycle
          if (ia > Aend) exit
          N = ia - Z
          Nix = Ninit - N
          expmass(Zix, Nix) = expmass1
          expmexc(Zix, Nix) = expmexc1
        enddo
        close (unit = 1)
      endif
    endif
    if (massmodel == 1) massfile = trim(path)//'masses/frdm/'//masschar
    if (massmodel == 0 .or. massmodel == 2) massfile = trim(path)//'masses/hfb/'//masschar
    if (massmodel == 3) massfile = trim(path)//'masses/hfbd1m/'//masschar
    if (massdir(1:1) /= ' ') then
      L = 0
      do i = 1, 132
        if (massdir(i:i) == ' ') exit
        L = L + 1
      enddo
      massfile = massdir(1:L)//'/'//masschar
    endif
    inquire (file = massfile, exist = lexist)
    if (lexist) then
      open (unit = 1, file = massfile, status = 'old', iostat = istat)
      if (istat /= 0) call read_error(massfile, istat)
      do
        read(1, '(4x, i4, 2f12.6, 2f8.4, 20x, f4.1, i2)', iostat = istat) ia, thmass1, thmexc1, b2, b4, gs, p
        if (istat == -1) exit
        if (istat /= 0) call read_error(massfile, istat)
        if (ia < Abegin) cycle
        if (ia > Aend) exit
        N = ia - Z
        Nix = Ninit - N
        if (massmodel /= 0) then
          thmass(Zix, Nix) = thmass1
          thmexc(Zix, Nix) = thmexc1
        endif
        if (beta2(Zix, Nix, 0) == 0.)  beta2(Zix, Nix, 0) = b2
        beta4(Zix, Nix) = b4
        gsspin(Zix, Nix) = gs
        gsparity(Zix, Nix) = p
      enddo
      close (unit = 1)
    endif
  enddo
  if (massmodel == 0) then
    flagduflo = .true.
  else
    flagduflo = .false.
  endif
  do Zix = 0, maxZ + 4
    do Nix = 0, maxN + 4
      A = Ainit - Zix - Nix
      if (massnucleus(Zix, Nix) /= 0.) then
        nucmass(Zix, Nix) = massnucleus(Zix, Nix)
        expmexc(Zix, Nix) = (massnucleus(Zix, Nix) - A) * amu
        thmexc(Zix, Nix) = expmexc(Zix, Nix)
        cycle
      endif
      if (massexcess(Zix, Nix) /= 0.) then
        expmexc(Zix, Nix) = massexcess(Zix, Nix)
        thmexc(Zix, Nix) = massexcess(Zix, Nix)
        nucmass(Zix, Nix) = A + massexcess(Zix, Nix) / amu
        cycle
      endif
      if (flagexpmass .and. expmass(Zix, Nix) /= 0.) then
        nucmass(Zix, Nix) = expmass(Zix, Nix)
      else
        nucmass(Zix, Nix) = thmass(Zix, Nix)
      endif
      if (nucmass(Zix, Nix) == 0.) flagduflo = .true.
    enddo
  enddo
!
! The target nucleus MUST be present in the masstable. This is to avoid unbound nuclei.
!
  if (nucmass(parZ(k0), parN(k0)) == 0. .and. .not.flagduflo) then
    write(*, '(" TALYS-error: Target nucleus not in masstable")')
    stop
  endif
!
! ********* Use analytical mass formula for remaining nuclei ***********
!
! If a residual nucleus is not in the experimental/theoretical mass table, or if massmodel=0,
! we use the analytical formula of Duflo-Zuker.
!
  if (flagduflo) then
Loop1: do Zix = 0, maxZ + 4
      Z = Zinit - Zix
      if (Z <= 0) cycle
      do Nix = 0, maxN + 4
        N = Ninit - Nix
        if (N <= 0) cycle Loop1
        A = Z + N
        if (nucmass(Zix, Nix) == 0. .and. expmass(Zix, Nix) == 0. .and. Zix <= numZastro .and. Nix <= numNastro) &
 &        write(*, '(" TALYS-warning: Duflo-Zuker mass for ", a, i3)') trim(nuc(Z)), A
        call duflo(N, Z, exc)
        dumexc(Zix, Nix) = exc
        thmass1 = A + exc / amu
        if (nucmass(Zix, Nix) == 0) nucmass(Zix, Nix) = thmass1
        if (expmass(Zix, Nix) == 0..and.massmodel == 0) nucmass(Zix, Nix) = thmass1
        if (massmodel == 0) thmexc(Zix, Nix) = dumexc(Zix, Nix)
      enddo
    enddo Loop1
  endif
!
! ********************* Reduced and specific masses ********************
!
  do Zix = 0, maxZ + 2
    do Nix = 0, maxN + 2
      do type = 0, 6
        if (parskip(type)) cycle
        specmass(Zix, Nix, type) = nucmass(Zix, Nix) / (nucmass(Zix, Nix) + parmass(type))
        redumass(Zix, Nix, type) = specmass(Zix, Nix, type) * parmass(type)
      enddo
    enddo
  enddo
  return
end subroutine masses
! Copyright A.J. Koning 2021
