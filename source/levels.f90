subroutine levels(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Discrete levels
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
! All global variables
!   numisom        ! number of isomers
!   numJ           ! maximum J - value
!   numlev         ! maximum number of discrete levels
!   numlev2        ! maximum number of levels
! Variables for direct reactions
!   flagautorot    ! flag for automatic rotational coupled channels
! Variables for basic parameters
!   isomer         ! definition of isomer in seconds
!   Lisoinp        ! user assignment of target isomer number
! Variables for discrete levels
!   disctable      ! table with discrete levels
!   levelfile      ! discrete level file
!   nlev           ! number of levels for nucleus
! Variables for level density
!   Ntop           ! highest discrete level for temperature matching
!   Risomer        ! adjustable correction to level branching ratios
! Variables for masses
!   massmodel      ! model for theoretical nuclear mass
! Variables for main input
!   k0             ! index of incident particle
!   Ltarget        ! excited level of target
! Variables for nuclides
!   AA             ! mass number of residual nucleus
!   ZZ             ! charge number of residual nucleus
! Variables for files
!   path           ! directory containing files to be read
! Constants
!   nuc            ! symbol of nucleus
!   parN           ! neutron number of particle
!   parZ           ! charge number of particle
! Variables for levels
!   bassign        ! flag for assignment of branching ratio
!   branchlevel    ! level to which branching takes place
!   branchratio    ! gamma - ray branching ratio to level
!   conv           ! conversion coefficient
!   edis           ! energy of level
!   ENSDF          ! string from original ENSDF discrete level file
!   jassign        ! flag for assignment of spin
!   jdis           ! spin of level
!   Liso           ! isomeric number of target
!   Lisomer        ! level number of isomer
!   Ltarget0       ! excited level of target
!   nbranch        ! number of branching levels
!   Nisomer        ! number of isomers for this nuclide
!   nlevmax2       ! maximum number of levels
!   parlev         ! parity of level
!   passign        ! flag for assignment of parity
!   tau            ! lifetime of state in seconds
! Variables for masses
!   gsparity       ! ground state parity
!   gsspin         ! ground state spin
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist            ! logical to determine existence
  logical           :: brexist(0:numlev) ! logical to determine existence
  character(len=1)  :: bas(numlev)       ! symbol for assignment of branching ratio
  character(len=6)  :: levelchar         ! help variable
  character(len=132):: levfile           ! level file
  integer           :: A                 ! mass number of target nucleus
  integer           :: i                 ! counter
  integer           :: ia                ! mass number from abundance table
  integer           :: ii                ! counter
  integer           :: istat             ! logical for file access
  integer           :: j                 ! counter
  integer           :: k                 ! counter
  integer           :: klev(numlev)      ! level number
  integer           :: Lis               ! isomer number
  integer           :: N                 ! neutron number of residual nucleus
  integer           :: nb                ! help variable
  integer           :: lbr               ! help variable
  integer           :: Nix               ! neutron number index for residual nucleus
  integer           :: nlev2             ! number of excited levels for nucleus
  integer           :: nlevlines         ! number of lines on discrete level file for nucleus
  integer           :: nnn               ! number of levels in discrete level file
  integer           :: Z                 ! charge number of target nucleus
  integer           :: Zix               ! charge number index for residual nucleus
  real(sgl)         :: br(numlev)        ! branching ratio multiplied by initial flux
  real(sgl)         :: con(numlev)       ! conversion factor
  real(sgl)         :: sum               ! help variable
!
! ******************** Default nuclear levels **************************
!
! For any nuclide, we first assign default ground state spins and parities.
! If there is information in the discrete level file, this will of course be overwritten.
! The index 0 of edis, etc. represents the ground state, the index 1 the first excited state, etc.
!
  Z = ZZ(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  levnum(Zix, Nix, 0) = 0
  edis(Zix, Nix, 0) = 0.
  jdis(Zix, Nix, 0) = gsspin(Zix, Nix)
  parlev(Zix, Nix, 0) = gsparity(Zix, Nix)
!
! We also assign a default value to the first excited state, if it does not exist in the discrete level file
!
  levnum(Zix, Nix, 1) = 1
  edis(Zix, Nix, 1) = min(26. / A, 10.)
  jdis(Zix, Nix, 1) = jdis(Zix, Nix, 0) + 2.
  parlev(Zix, Nix, 1) = parlev(Zix, Nix, 0)
  jassign(Zix, Nix, 1) = 'J'
  eassign(Zix, Nix, 1) = 'E'
  passign(Zix, Nix, 1) = 'P'
  branchratio(Zix, Nix, 1, 1) = 1.
  nbranch(Zix, Nix, 1) = 1
  bassign(Zix, Nix, 1, 0) = 'B'
  conv(Zix, Nix, 1, 0) = 0.
  tau(Zix, Nix, 1) = 0.
!
! ************************ Read nuclear levels *************************
!
! Note that nlev is the number of excited levels, i.e. excluding the ground state.
!
! 1. Inquire whether file is present
!
  if (levelfile(Zix)(1:1) /= ' ') then
    levfile = levelfile(Zix)
  else
    levelchar = trim(nuc(Z))//'.lev'
    if (disctable == 1) levfile = trim(path)//'levels/final/'//levelchar
    if (disctable == 2) levfile = trim(path)//'levels/exp/'//levelchar
    if (disctable == 3) levfile = trim(path)//'levels/hfb/'//levelchar
    inquire (file = levfile, exist = lexist)
    if ( .not. lexist) return
  endif
!
! 2. Search for the isotope under consideration
!
  con = 0.
  br = 0.
  bas = 'B'
  nlev2 = 0
  nnn = 0
  open (unit = 2, file = levfile, status = 'old', iostat = istat)
  if (istat /= 0) call read_error(levfile, istat)
  do
   read(2, '(4x, i4, 2i5)', iostat = istat) ia, nlevlines, nnn
   if (istat == -1) exit
    if (A == ia) then
      nlev2 = min(nnn, numlev)
      if (nnn < 2) flagautorot = .false.
!
! 3. Read discrete level information
!
      do i = 0, nlev2
        read(2, '(4x, f11.6, f6.1, 3x, i2, i3, 18x, e10.3, 3a1, a18)') edis(Zix, Nix, i), jdis(Zix, Nix, i), parlev(Zix, Nix, i), &
 &        nb, tau(Zix, Nix, i), eassign(Zix, Nix, i), jassign(Zix, Nix, i), passign(Zix, Nix, i), ENSDF(Zix, Nix, i)
        do j = 1, nb
          read(2, '(29x, i3, f10.6, e10.3, 5x, a1)') klev(j), br(j), con(j), bas(j)
        enddo
!
! Branching ratio from input
!
        if (nbranch(Zix, Nix, i) > 0) then
          do j = 1, nbranch(Zix, Nix, i)
            conv(Zix, Nix, i, j) = con(j)
            bassign(Zix, Nix, i, j) = bas(j)
          enddo
        else
!
! Branching ratio from discrete level file
!
          ii = 0
          do j = 1, nb
            if (br(j) /= 0.) then
              ii = ii + 1
              branchlevel(Zix, Nix, i, ii) = klev(j)
              branchratio(Zix, Nix, i, ii) = br(j)
              conv(Zix, Nix, i, ii) = con(j)
              bassign(Zix, Nix, i, ii) = bas(j)
            endif
          enddo
          nbranch(Zix, Nix, i) = ii
        endif
!
! Spins beyond numJ are set to numJ
!
        jdis(Zix, Nix, i) = min(jdis(Zix, Nix, i), real(numJ))
      enddo
!
! Lifetimes below the isomeric definition are set to zero.
! The isomeric number is determined.
!
      do i = 0, nlev2
        tauripl(Zix, Nix, i) = tau(Zix, Nix, i)
        if (tau(Zix, Nix, i) < isomer) tau(Zix, Nix, i) = 0.
        levnum(Zix, Nix, i) = i
      enddo
      if (massmodel <= 1) then
        if (jassign(Zix, Nix, 0) == 'J' .and. passign(Zix, Nix, 0) == 'P') then
          jdis(Zix, Nix, 0) = gsspin(Zix, Nix)
          parlev(Zix, Nix, 0) = gsparity(Zix, Nix)
        endif
        if (nlev2 == 0) then
          jdis(Zix, Nix, 1) = jdis(Zix, Nix, 0) + 2.
          parlev(Zix, Nix, 1) = parlev(Zix, Nix, 0)
        endif
      endif
!
! Read extra levels which are used only for the level density matching problem or for direct reactions (deformation parameters).
! The branching ratios are not read for these higher levels.
!
      nlevmax2(Zix, Nix) = min(nnn, numlev2)
      do i = nlev2 + 1, nlevmax2(Zix, Nix)
        read(2, '(4x, f11.6, f6.1, 3x, i2, i3, 18x, e10.3, 3a1)') edis(Zix, Nix, i), jdis(Zix, Nix, i), parlev(Zix, Nix, i), &
 &        nb, tau(Zix, Nix, i), eassign(Zix, Nix, i), jassign(Zix, Nix, i), passign(Zix, Nix, i)
        jdis(Zix, Nix, i) = min(jdis(Zix, Nix, i), real(numJ))
        do j = 1, nb
          read(2, * )
        enddo
      enddo
      exit
    else
      do i = 1, nlevlines
        read(2, '()')
      enddo
    endif
  enddo
  close (unit = 2)
  nlev(Zix, Nix) = min(nlev(Zix, Nix), nlev2)
  nlev(Zix, Nix) = max(nlev(Zix, Nix), 1)
!
! The maximal value of Ntop is always given by the last discrete level of the discrete level file.
!
  nlevmax2(Zix, Nix) = min(nnn, numlev2)
  nlevmax2(Zix, Nix) = max(nlevmax2(Zix, Nix), 1)
  Ntop(Zix, Nix, 0) = min(nlevmax2(Zix, Nix), Ntop(Zix, Nix, 0))
!
! Check existence of excited level of target
!
  if (Lisoinp ==  - 1 .and. Ltarget /= 0) then
    if (Zix == parZ(k0) .and. Nix == parN(k0) .and. Ltarget > nlev(Zix, Nix)) then
      write(*, '(" TALYS-error: excited level of target does not exist")')
      stop
    endif
  endif
!
! Determine isomeric level number
!
  Lis = 0
  do i = 1, nlevmax2(Zix, Nix)
    if (Lis < numisom .and. tau(Zix, Nix, i) >= isomer) then
      Lis = Lis + 1
      Lisomer(Zix, Nix, Lis) = i
    endif
  enddo
  Nisomer(Zix, Nix) = Lis
!
! Determine isomeric number of target
!
  if (Zix == parZ(k0) .and. Nix == parN(k0)) then
    if (Lisoinp ==  -1) then
      Liso = 0
      if (Ltarget /= 0) then
        do i = 1, Lis
          if (Ltarget == Lisomer(Zix, Nix, i)) Liso = i
        enddo
      endif
    else
      Liso = Lisoinp
    endif
    if (Liso > 0) targetnuclide = trim(targetnuclide0) // isochar(min(Liso,numisom))
  endif
!
! Special treatment for isomers in the continuum. There are about 10 known isomers whose level number is larger than 30.
! To avoid wasting too much memory we renumber the isomer in the continuum to the last discrete level taken into account
! in the calculation.
!
  Lis = Nisomer(Zix, Nix) + 1
  do i = nlevmax2(Zix, Nix), nlev(Zix, Nix) + 1, - 1
    if (tau(Zix, Nix, i) > isomer .and. isomer >= 0.1) then
      Lis = Lis - 1
      N = nlev(Zix, Nix) - Nisomer(Zix, Nix) + Lis
      if (Lis >= 0 .and. N >= 0) then
        levnum(Zix, Nix, N) = min(i, 99)
        edis(Zix, Nix, N) = edis(Zix, Nix, i)
        jdis(Zix, Nix, N) = jdis(Zix, Nix, i)
        parlev(Zix, Nix, N) = parlev(Zix, Nix, i)
        tau(Zix, Nix, N) = tau(Zix, Nix, i)
        if (N /= i) tau(Zix, Nix, i) = 0.
        eassign(Zix, Nix, N) = ' '
        jassign(Zix, Nix, N) = ' '
        passign(Zix, Nix, N) = ' '
        if (Ltarget0 == Lisomer(Zix, Nix, Lis) .and. Zix == parZ(k0) .and. Nix == parN(k0) .and. tau(Zix, Nix, N) > isomer ) & 
 &        Ltarget = N
      endif
    endif
  enddo
!
! Adjust branching ratios for isomeric cross sections
!
  if (Lis > 0 .and. Risomer(Zix, Nix) /= 1.) then
    do i = 1, nlev2
      if (tau(Zix,Nix,i) >= isomer) then
        brexist = .false.
        brexist(i) = .true.
        do j = i+1, nlev2
          do k = 1, nbranch(Zix,Nix,j)
            lbr = branchlevel(Zix,Nix,j,k)
            if (brexist(lbr)) brexist(j) = .true.
          enddo
        enddo   
        do j = i+1, nlev2
          do k = 1, nbranch(Zix,Nix,j)
            lbr = branchlevel(Zix,Nix,j,k)
            if (brexist(lbr)) branchratio(Zix,Nix,j,k) = Risomer(Zix,Nix) * branchratio(Zix,Nix,j,k)
          enddo
        enddo
      endif
    enddo
    do i = 1, nlev2
      sum = 0.
      do k = 1, nbranch(Zix,Nix,i)
        sum = sum + branchratio(Zix,Nix,i,k)
      enddo 
      if (sum >  0.) then
        do k = 1, nbranch(Zix,Nix,i)
          branchratio(Zix,Nix,i,k) = branchratio(Zix,Nix,i,k) / sum
        enddo
      endif
    enddo
  endif
!
! Extract level information for resonances of light nuclides
!
  if (flagpseudores .and. Zix == parZ(k0) .and. Nix == parN(k0)) call pseudo_resonance
  return
end subroutine levels
! Copyright A.J. Koning 2021
