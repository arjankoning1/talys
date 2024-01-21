subroutine deformpar(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deformation parameters
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
!   numbar         ! number of fission barriers
!   numlev         ! maximum number of discrete levels
!   numlev2        ! maximum number of levels
!   numrotcc       ! number of rotational deformation parameters
! Variables for masses
!   beta2          ! deformation parameter
! Variables for direct reactions
!   flagautorot    ! flag for automatic rotational coupled channels
!   flagspher      ! flag to force spherical optical model
!   maxband        ! highest vibrational band added to rotational model
!   maxrot         ! number of included excited rotational levels
! Variables for discrete levels
!   deformfile     ! deformation parameter file
!   disctable      ! table with discrete levels
! Variables for main input
!   k0             ! index of incident particle
! Variables for nuclides
!   NN             ! neutron number of residual nucleus
! Variables for nuclides
!   AA             ! mass number of residual nucleus
!   ZZ             ! charge number of residual nucleus
! Constants
!   amu            ! atomic mass unit in MeV
!   hbarc          ! hbar.c in MeV.fm
!   magic          ! magic numbers
!   onethird       ! 1 / 3
!   nuc            ! symbol of nucleus
!   parmass        ! mass of particle in a.m.u.
!   sgn            ! sign
! Variables for files
!   path           ! directory containing files to be read
! Variables for levels
!   edis           ! energy of level
!   jdis           ! spin of level
!   parlev         ! parity of level
! Variables for ECIS
!   spin           ! spin of incident particle
! Variables for deformation parameters
!   colltype       ! type of collectivity (D, V or R)
!   deform         ! deformation parameter
!   defpar         ! deformation parameter
!   deftype        ! deformation length (D) or parameter (B)
!   indexcc        ! level index for coupled channel
!   indexlevel     ! level index
!   iphonon        ! phonon (1 or 2)
!   Irigid         ! rigid body value of moment of inertia
!   Irigid0        ! undeformed rigid body value of moment of inertia
!   Kband          ! magnetic quantum number
!   lband          ! angular momentum
!   leveltype      ! type of level (rotational (R) or vibrational (V))
!   ndef           ! number of collective levels
!   nrot           ! number of deformation parameters for rotational nucleus
!   rotpar         ! deformation parameters for rotational nucleus
!   vibband        ! band number of level
! Variables for masses
!   beta4          ! deformation parameters
!
! *** Declaration of local data
!
  implicit none
  logical           :: first2               ! flag to determine first state of specific spin
  logical           :: first3               ! flag to determine first state of specific spin
  logical           :: first4               ! flag to determine first state of specific spin
  logical           :: lexist               ! logical to determine existence
  character(len=1)  :: colltype1            ! help variable
  character(len=1)  :: deftype1             ! deformation length (D) or parameter (B)
  character(len=1)  :: leveltype1           ! type of level (rotational (R) or vibrational (V))
  character(len=6)  :: defchar              ! help variable
  character(len=132):: deffile              ! deformation parameter file
  integer           :: A                    ! mass number of target nucleus
  integer           :: distance             ! number of nucleons to closest magic number
  integer           :: i                    ! counter
  integer           :: ia                   ! mass number from abundance table
  integer           :: ibar                 ! fission barrier
  integer           :: idef                 ! counter for deformation
  integer           :: ii                   ! counter
  integer           :: iirot                ! help variable
  integer           :: iphonon1             ! phonon (1 or 2)
  integer           :: irot                 ! counter for rotational band
  integer           :: istat                ! logical for file access
  integer           :: k                    ! designator for particle
  integer           :: k2                   ! counter
  integer           :: Kmag1                ! magnetic quantum number
  integer           :: lband1               ! angular momentum
  integer           :: N                    ! neutron number of residual nucleus
  integer           :: natpar               ! natural parity
  integer           :: ndisc                ! number of lines on discrete level file for nucleus
  integer           :: nex                  ! excitation energy bin of compound nucleus
  integer           :: Nix                  ! neutron number index for residual nucleus
  integer           :: nrotlev              ! number of rotational levels
  integer           :: odd                  ! odd (1) or even (0) nucleus
  integer           :: type                 ! particle type
  integer           :: vibband1             ! band number of level
  integer           :: Z                    ! charge number of target nucleus
  integer           :: Zix                  ! charge number index for residual nucleus
  real(sgl)         :: deform1(numrotcc)    ! deformation parameter
  real(sgl)         :: dspin                ! angular momentum increase for rotational band
  real(sgl)         :: R                    ! radius
!
! ************************ Read deformation parameters *****************
!
  Z = ZZ(Zix, Nix, 0)
  N = NN(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  colltype(Zix, Nix) = 'S'
  defchar = trim(nuc(Z))//'.def'
!
! 1. Inquire whether file is present
!
  if (deformfile(Zix)(1:1) /= ' ') then
    deffile = deformfile(Zix)
  else
    deffile = trim(path)//'deformation/'//defchar
  endif
  inquire (file = deffile, exist = lexist)
  if (lexist .and. disctable /= 3) then
    open (unit = 2, file = deffile, status = 'old')
!
! 2. Search for the isotope under consideration
!
    do
      read(2, '(4x, 2i4, 2(3x, a1))', iostat = istat) ia, ndisc, colltype1, deftype1
      if (istat == -1) exit
      if (A == ia) then
!
! Initialization
!
        do i = 1, numrotcc
          deform1(i) = 0.
        enddo
!
! 3. Read deformation parameters
!
! The coupling scheme is read from the nuclear structure database.
! Different actions are performed if the requested calculation is for a spherical (S), vibrational (V) or a rotational (R) nucleus.
!
        if (colltype1 /= 'R' .and. colltype1 /= 'A' .and. colltype1 /= 'V') colltype1 = 'S'
        if (flagspher .and. colltype1 == 'A') then
          iirot = 1
        else
          iirot = numrotcc
        endif
        if (flagspher) colltype1 = 'S'
        colltype(Zix, Nix) = colltype1
        deftype(Zix, Nix) = deftype1
        idef = 0
        irot = 0
        ii = 0
        do i = 1, ndisc
          read(2, '(i4, 3x, a1, 4i4, 4f9.5)') nex, leveltype1, vibband1, lband1, Kmag1, iphonon1, (deform1(k), k = 1, iirot)
          if (nex <= numlev2) then
            idef = idef + 1
            indexlevel(Zix, Nix, i) = nex
            if (leveltype1 == 'R') then
              irot = irot + 1
              if (irot > maxrot+1) leveltype1 = 'D'
            endif
            leveltype(Zix, Nix, nex) = leveltype1
            vibband(Zix, Nix, i) = vibband1
            if (colltype(Zix, Nix) == 'R' .and. vibband1 > maxband) leveltype1 = 'D'
            iphonon(Zix, Nix, i) = max(iphonon1, 1)
!
! Default: read deformation parameters from mass table
!
            if (leveltype1 == 'R' .and. rotpar(Zix, Nix, 1) == 0.) then
              if (nex == 0) then
                if (deform1(1) == 0.) then
                  deform1(1) = beta2(Zix, Nix, 0)
                  deform1(2) = beta4(Zix, Nix)
                endif
              endif
              do k = 1, numrotcc
                rotpar(Zix, Nix, k) = deform1(k)
                if (deform1(k) == 0.) then
                  nrot(Zix, Nix) = k - 1
                  goto 130
                endif
              enddo
            endif
            if (leveltype1 == 'V' .and. defpar(Zix, Nix, vibband1) == 0.) then
              lband(Zix, Nix, vibband1) = lband1
              Kband(Zix, Nix, vibband1) = Kmag1
              defpar(Zix, Nix, vibband1) = deform1(1)
            endif
    130     if (colltype(Zix, Nix) == 'S' .or. leveltype1 == 'D') then
!
! For DWBA, we only include natural parity states.
!
              natpar = int(sgn(int(jdis(Zix, Nix, nex))))
              if (parlev(Zix, Nix, nex) == natpar) deform(Zix, Nix, nex) = deform1(1)
              if (i > 1 .and. i <= numrotcc+1 .and. leveltype1 == 'R') deform(Zix, Nix, nex) = rotpar(Zix, Nix, i - 1)
            else
              ii = ii + 1
              indexcc(Zix, Nix, ii) = nex
            endif
          endif
        enddo
        ndef(Zix, Nix) = idef
        exit
      else
        do i = 1, ndisc
          read(2, '()')
        enddo
      endif
    enddo
    close (unit = 2)
  endif
!
! ******************** Default deformation parameters ******************
!
! Automatic assignment of rotational deformation parameters.
! Calculate distance to closest magic number as a measure for deformation.
!
  distance = 1000
  do k2 = 1, 8
    distance = min(abs(N - magic(k2)), distance)
    distance = min(abs(Z - magic(k2)), distance)
  enddo
  odd = mod(A, 2)
!
! Read rotational deformation parameters
!
  if (colltype(Zix, Nix) == 'S' .and. ( .not. flagspher) .and. A > 150 .and. distance >= 8 .and. flagautorot) then
    indexlevel(Zix, Nix, 1) = 0
    indexcc(Zix, Nix, 1) = 0
    leveltype(Zix, Nix, 0) = 'R'
    if (odd == 0) then
      dspin = 2.
    else
      dspin = 1.
    endif
    ndef(Zix, Nix) = maxrot + 1
Loop1: do i = 2, ndef(Zix, Nix)
      ii = i - 1
      spin = jdis(Zix, Nix, 0) + dspin * ii
      do nex = 1, numlev2
        if (spin == jdis(Zix, Nix, nex) .and. parlev(Zix, Nix, nex) == parlev(Zix, Nix, 0)) then
          indexlevel(Zix, Nix, i) = nex
          indexcc(Zix, Nix, i) = nex
          leveltype(Zix, Nix, nex) = 'R'
          cycle Loop1
        endif
      enddo
    enddo Loop1
    if (indexcc(Zix, Nix, ndef(Zix, Nix)) /= 0) colltype(Zix, Nix) = 'R'
    nrot(Zix, Nix) = 2
    rotpar(Zix, Nix, 1) = beta2(Zix, Nix, 0)
    rotpar(Zix, Nix, 2) = beta4(Zix, Nix)
    deftype(Zix, Nix) = 'B'
    close (unit = 2)
    if (odd == 0) goto 400
  endif
!
! Number of rotational levels should remain below maxrot
!
  nrotlev = 0
  do nex = 1, numlev2
    if (leveltype(Zix, Nix, nex) == 'R') then
      nrotlev = nrotlev + 1
      if (nrotlev > maxrot) leveltype(Zix, Nix, nex) = 'D'
    endif
  enddo
!
! Assign vibrational deformation parameters
!
! Systematics for first 2+, 3- and 4+ vibrational states, derived from individual deformation parameters.
! Also, we assign small deformation parameter to all remaining discrete levels.
!
  if (odd /= 0) goto 400
  if (colltype(Zix, Nix) /= 'S') then
    first2 = .false.
    first3 = .false.
    first4 = .false.
  else
    first2 = .true.
    first3 = .true.
    first4 = .true.
  endif
  type = 2 * Zix + Nix
  do k = 0, numlev
    if (k == 0 .and. type == k0) cycle
    if (colltype(Zix, Nix) /= 'S' .and. leveltype(Zix, Nix, k) == 'V') cycle
    if (colltype(Zix, Nix) == 'A' .and. leveltype(Zix, Nix, k) == 'R') cycle
    if (colltype(Zix, Nix) == 'R' .and. leveltype(Zix, Nix, k) == 'R') cycle
    if (jdis(Zix, Nix, k) == 0.) cycle
    if (first2 .and. jdis(Zix, Nix, k) == 2..and. parlev(Zix, Nix, k) == 1) then
      if (leveltype(Zix, Nix, k) /= 'R' .and. deform(Zix, Nix, k) == 0.) then
        deform(Zix, Nix, k) = 0.40 * exp( - 0.012 * A) + 0.025 * min(distance, 5)
        if (edis(Zix, Nix, k) <= 0.1) deform(Zix, Nix, k) = 0.02
        if (deftype(Zix, Nix) == 'D') deform(Zix, Nix, k) = deform(Zix, Nix, k) * 1.24 * (A **onethird)
      endif
      first2 = .false.
      cycle
    endif
    if (first3 .and. jdis(Zix, Nix, k) == 3..and. parlev(Zix, Nix, k) ==  - 1) then
      if (leveltype(Zix, Nix, k) /= 'R' .and. deform(Zix, Nix, k) == 0.) then
        deform(Zix, Nix, k) = 0.35 * exp( - 0.008 * A)
        if (edis(Zix, Nix, k) <= 0.1) deform(Zix, Nix, k) = 0.02
        if (deftype(Zix, Nix) == 'D') deform(Zix, Nix, k) = deform(Zix, Nix, k) * 1.24 * (A **onethird)
      endif
      first3 = .false.
      cycle
    endif
    if (first4 .and. jdis(Zix, Nix, k) == 4..and. parlev(Zix, Nix, k) == 1) then
      if (leveltype(Zix, Nix, k) /= 'R' .and. deform(Zix, Nix, k) == 0.) then
        deform(Zix, Nix, k) = 0.20 * exp( - 0.006 * A)
        if (edis(Zix, Nix, k) <= 0.1) deform(Zix, Nix, k) = 0.02
        if (deftype(Zix, Nix) == 'D') deform(Zix, Nix, k) = deform(Zix, Nix, k) * 1.24 * (A **onethird)
      endif
      first4 = .false.
      cycle
    endif
    if (deform(Zix, Nix, k) /= 0.) cycle
    natpar = int(sgn(int(jdis(Zix, Nix, k))))
    if (parlev(Zix, Nix, k) == natpar) deform(Zix, Nix, k) = 0.02
    if (deftype(Zix, Nix) == 'D') deform(Zix, Nix, k) = deform(Zix, Nix, k) * 1.24 * (A **onethird)
  enddo
!
! ************** Rigid body value for moment of inertia ****************
!
  400 R = 1.2 * A **onethird
  Irigid0(Zix, Nix) = 0.4 * R * R * A * parmass(1) * amu / (hbarc **2)
  do ibar = 0, numbar
    Irigid(Zix, Nix, ibar) = (1 + abs(beta2(Zix, Nix, ibar)) / 3.) * Irigid0(Zix, Nix)
  enddo
  return
end subroutine deformpar
! Copyright A.J. Koning 2021
