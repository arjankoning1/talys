subroutine fissionpar(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Fission parameters
!
! Author    : Stephane Hilaire, Marieke Duijvestijn and Arjan Koning
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
!   numbar         ! number of fission barriers
!   numlev         ! maximum number of discrete levels
! Variables for fission
!   axtype           ! type of axiality of barrier
!   betafiscor       ! adjustable factor for fission path width
!   betafiscoradjust ! adjustable factor for fission path width
!   clas2file        ! file with class 2 transition states
!   fbarrier         ! height of fission barrier
!   fismodelalt      ! alternative fission model for default barriers
!   fismodelx        ! fission model
!   flagclass2       ! flag for class2 states in fission
!   flaghbstate      ! flag for head band states in fission
!   fwidth           ! width of fission barrier
!   hbtransfile      ! file with head band transition states
!   vfiscor          ! adjustable factor for fission path height
!   vfiscoradjust    ! adjustable factor for fission path height
! Variables for numerics
!   nbins0         ! number of continuum excitation energy bins
! Variables for level density
!   Rclass2mom     ! norm. constant for moment of inertia for class 2 states
!   Rtransmom      ! norm. constant for moment of inertia for transition states
! Variables for nuclides
!   AA             ! mass number of residual nucleus
!   NN             ! neutron number of residual nucleus
!   ZZ             ! charge number of residual nucleus
! Variables for files
!   path           ! directory containing files to be read
! Constants
!   fislim         ! mass above which nuclide fissions
!   nuc            ! symbol of nucleus
! Variables for deformation parameters
!   Irigid         ! rigid body value of moment of inertia
! Variables for fission parameters
!   fecont         ! start of continuum energy
!   efisc2hb       ! energy of class2 states
!   efistrhb       ! energy of head band transition states
!   jfisc2hb       ! spin of class2 states
!   jfistrhb       ! spin of head band transition states
!   minertc2       ! moment of inertia for class2 states
!   minertia       ! moment of inertia of fission barrier deformation
!   nclass2        ! number of sets of class2 states
!   nfisbar        ! number of fission barrier parameters
!   nfisc2hb       ! number of class2 states for barrier
!   nfistrhb       ! number of head band transition states for barrier
!   pfisc2hb       ! parity of class2 states
!   pfistrhb       ! parity of head band transition states
! Variables for WKB
!   betafis        ! fission path width
!   nbeta          ! number of beta values
!   nbinswkb       ! integration step for WKB calculation
!   vfis           ! adjustable factor for fission path height
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist      ! logical to determine existence
  character(len=6)  :: fischar     ! help variable
  character(len=132):: c2file      ! file with class 2 states
  character(len=132):: fisfile     ! fission file
  character(len=132):: hbsfile     ! file with head band transition states
  integer           :: A           ! mass number of target nucleus
  integer           :: fislocal    ! fission model
  integer           :: i           ! counter
  integer           :: ia          ! mass number from abundance table
  integer           :: il          ! angular momentum
  integer           :: istat       ! logical for file access
  integer           :: j           ! counter
  integer           :: ib          ! index
  integer           :: modn        ! help variable
  integer           :: modz        ! help variable
  integer           :: N           ! neutron number of residual nucleus
  integer           :: nbar        ! number of fission barriers
  integer           :: Nix         ! neutron number index for residual nucleus
  integer           :: Z           ! charge number of target nucleus
  integer           :: Zix         ! charge number index for residual nucleus
  real(sgl)         :: bar1        ! inner and outer barrier heights
  real(sgl)         :: bar2        ! inner and outer barrier heights
  real(sgl)         :: bb          ! fission path parameter
  real(sgl)         :: egs         ! rotating ground state energy
  real(sgl)         :: esp         ! saddle point energy
  real(sgl)         :: hw1         ! inner and outer barrier curvatures
  real(sgl)         :: hw2         ! inner and outer barrier curvatures
  real(sgl)         :: lbar0       ! l-value for which bfis becomes zero
  real(sgl)         :: vv          ! fission path parameter
  real(sgl)         :: b22         ! help variable
  real(sgl)         :: b30         ! help variable
  real(sgl)         :: rmiu        ! moment of inertia parameter
!
! ****************** Read fission barrier parameters *******************
!
! Note that next to the chosen fission model (fismodel), there is always an alternative fission model (fismodelalt) which comes
! into play if fission parameters for the first choice model are not available.

! Determine whether barrier parameters have been provided in input
!
  nfisbar(Zix, Nix) = 0
  do i = 1, numbar
    if (fbarrier(Zix, Nix, i) /= 0..or.fwidth(Zix, Nix, i) /= 0.) nfisbar(Zix, Nix) = nfisbar(Zix, Nix) + 1
  enddo
   20 fislocal = fismodelx(Zix, Nix)
!
! Fission parameters from database
!
  Z = ZZ(Zix, Nix, 0)
  N = NN(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  fischar = trim(nuc(Z))//'.bar'
!
! Fismodel 1: Experimental parameters
!
! Fission barriers from database may have been overruled by user input.
! As starting point, we take the RIPL values as compiled by Maslov.
!
  if (fislocal == 1) then
    fisfile = trim(path)//'fission/barrier/'//fischar
    inquire (file = fisfile, exist = lexist)
    if (lexist) then
      open (unit = 2, file = fisfile, status = 'old', iostat = istat)
      if (istat /= 0) call read_error(fisfile, istat)
      do
        read(2, '(4x, i4, 4x, 1x, 2(f8.2), 5x, 2(f8.2))', iostat = istat) ia, bar1, hw1, bar2, hw2
        if (istat == -1) exit
        if (A == ia) then
          if (fbarrier(Zix, Nix, 1) == 0.) fbarrier(Zix, Nix, 1) = bar1
          if (fwidth(Zix, Nix, 1) == 0.) fwidth(Zix, Nix, 1) = hw1
          if (fbarrier(Zix, Nix, 2) == 0.) fbarrier(Zix, Nix, 2) = bar2
          if (fwidth(Zix, Nix, 2) == 0.) fwidth(Zix, Nix, 2) = hw2
          if (nfisbar(Zix, Nix) /= 3) nfisbar(Zix, Nix) = 2
          exit
        endif
      enddo
      close (unit = 2)
    endif
    if (nfisbar(Zix, Nix) == 0) fislocal = 2
  endif
!
! Fismodel 2: Mamdouh parameters
!
  if (fislocal == 2) then
    fisfile = trim(path)//'fission/mamdouh/'//fischar
    inquire (file = fisfile, exist = lexist)
    if (lexist) then
      open (unit = 2, file = fisfile, status = 'old', iostat = istat)
      if (istat /= 0) call read_error(fisfile, istat)
      do
        read(2, '(4x, i4, 2(24x, f8.2))', iostat = istat) ia, bar1, bar2
        if (istat == -1) exit
        if (A == ia) then
          if (fbarrier(Zix, Nix, 1) == 0.) fbarrier(Zix, Nix, 1) = bar1
          if (fbarrier(Zix, Nix, 2) == 0.) fbarrier(Zix, Nix, 2) = bar2
          if (fbarrier(Zix, Nix, 1) == 0..or.fbarrier(Zix, Nix, 2) == 0.) then
            nfisbar(Zix, Nix) = 1
          else
            nfisbar(Zix, Nix) = 2
          endif
        endif
      enddo
      close (unit = 2)
    endif
  endif
  if (fismodelx(Zix, Nix) <= 2 .and. nfisbar(Zix, Nix) == 0) fislocal = fismodelalt
!
! Fismodel 3: Sierk
!
! barsierk: subroutine for fission barrier heights, rotating gs energy and lbar0
! fisdata : subroutine to fit parameter values for reconstruction of fission barriers
!
! Empirical adjustment of fission barrier to globally fit subactinide fission
!
  if (fislocal == 3) then
    call fisdata
    il = 0
    nfisbar(Zix, Nix) = 1
    call barsierk(Z, A, il, bar1, egs, lbar0)
    if (fbarrier(Zix, Nix, 1) == 0.) fbarrier(Zix, Nix, 1) = Cbarrier * bar1
    if (fwidth(Zix, Nix, 1) == 0.) fwidth(Zix, Nix, 1) = 0.24
  endif
!
! Fismodel 4: Rotating Liquid Drop Model
!
! rldm: subroutine for saddle point energies, rotating gs energy
!
  if (fislocal == 4) then
    call fisdata
    il = 0
    nfisbar(Zix, Nix) = 1
    call rldm(Z, A, il, egs, esp)
    if (fbarrier(Zix, Nix, 1) == 0.) fbarrier(Zix, Nix, 1) = esp - egs
    if (fwidth(Zix, Nix, 1) == 0.) fwidth(Zix, Nix, 1) = 0.24
  endif
!
! Fismodel 5: WKB approximation
!
! Read the potential energy curve and call the WKB subroutine
!
! wkb       : subroutine for WKB approximation for fission
!
  nextr = 0
  iiextr(0)=1
  if (fislocal >= 5) then
    nfisbar(Zix, Nix) = 0
    fischar = trim(nuc(Z))//'.fis'
    if (fislocal == 5) fisfile = trim(path)//'fission/hfbpath/'//fischar
    if (fislocal == 6) fisfile = trim(path)//'fission/hfbpath_bskg3/'//fischar
    inquire (file = fisfile, exist = lexist)
    if (lexist) then
      open (unit = 2, file = fisfile, status = 'old')
300   if (fislocal.eq.5) then
        read(2, '(/11x, i4, 12x, i4//)', end = 330) ia, nbeta
      else
        read(2, '(/11x, i5, 12x, i4//)', end = 330) ia, nbeta
      endif
      if (A /= ia) then
        do i = 1, nbeta
          read(2, '()')
        enddo
        goto 300
      else
        do i = 1, nbeta
          if (fislocal.eq.5) then
            read(2, '(f10.3, 20x, f10.3)') bb, vv
!
! 6/12/2024: new mean rmiu value adapted to U236, so that action integral S~S(exp)=52.21
!               rmiu=0.0047*A**(5./3.)  --> S~39.5
!
            rmiu=rmiufiscor(Zix,Nix)*0.0063*A**(5./3.)
            betafis(i) = betafiscoradjust(Zix, Nix) * betafiscor(Zix, Nix) * bb
          else
            read(2,*) bb,b22,b30,vv,rmiu,ib
!
! for BSkG3 the deformation corresponds to the index i (rmiu coherently associated)
!
            betafis(i) = betafiscoradjust(Zix, Nix) * betafiscor(Zix, Nix) * float(i)
!
! ib=-1 corresponds to the ground state
! ib= 1 corresponds to identified maxima
! ib= 2 corresponds to identified minima
! ib= 9 corresponds to non-considered extrema (if more than 3 maxima exist)
!
            if (ib.eq.1.or.ib.eq.2) then
              nextr=nextr+1
              iiextr(nextr)=i
            endif
          endif
          vfis(i) = vfiscoradjust(Zix, Nix) * vfiscor(Zix, Nix) * vv
          rmiufis(i) = rmiufiscoradjust(Zix, Nix) * rmiufiscor(Zix, Nix) * rmiu
        enddo
        if (nbins0 == 0) then
          nbinswkb = 30
        else
          nbinswkb = nbins0
        endif
        if (fismodel == 6) then
          iiextr(nextr+1) = nbeta
          if (nextr > 5) write(*, '(" TALYS-warning: number of barriers ",i3," > 3 for Z = ",i3," A = ",i3)') nextr, Z, A
        endif
        call wkb(Z, A, Zix, Nix, nbar)
        nfisbar(Zix, Nix) = nbar
      endif
      goto 340
330   fismodelx(Zix, Nix) = fismodelalt
      axtype(Zix, Nix, 2) = 2
      close (unit = 2)
      goto 20
340   close (unit = 2)
    else
      fismodelx(Zix, Nix) = fismodelalt
      axtype(Zix, Nix, 2) = 2
      goto 20
    endif
  endif
!
! Read fission states
!
  modz = mod(Z, 2)
  modn = mod(N, 2)
  if (modz == 0) then
    if (modn == 0) then
      hbsfile = trim(path)//'fission/states/hbstates.ee'
      c2file = trim(path)//'fission/states/class2states.ee'
    else
      hbsfile = trim(path)//'fission/states/hbstates.eo'
      c2file = trim(path)//'fission/states/class2states.eo'
    endif
  else
    if (modn == 0) then
      hbsfile = trim(path)//'fission/states/hbstates.oe'
      c2file = trim(path)//'fission/states/class2states.oe'
    else
      hbsfile = trim(path)//'fission/states/hbstates.oo'
      c2file = trim(path)//'fission/states/class2states.oo'
    endif
  endif
!
! Use user-defined files for head band and class 2 transition states
!
  if (hbtransfile(Zix, Nix)(1:1) /= ' ') hbsfile = hbtransfile(Zix, Nix)
  if (clas2file(Zix, Nix)(1:1) /= ' ')  c2file = clas2file(Zix, Nix)
!
! Read head band transition states
!
  if (flaghbstate) then
    open (unit = 2, file = hbsfile, status = 'old', iostat = istat)
    if (istat /= 0) call read_error(c2file, istat)
    do i = 1, nfisbar(Zix, Nix)
      read(2, '(4x, i4, f8.3)', iostat = istat) nfistrhb(Zix, Nix, i), fecont(Zix, Nix, i)
      if (istat /= 0) cycle
      call range_integer_error(c2file, nfistrhb(Zix, Nix, i), 0, numlev, index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', &
 &        index3 = i, name3 = 'i')
      do j = 1, nfistrhb(Zix, Nix, i)
        read(2, '(4x, f11.6, f6.1, i5)', iostat = istat) efistrhb(Zix, Nix, i, j), jfistrhb(Zix, Nix, i, j), &
          pfistrhb(Zix, Nix, i, j)
        if (istat /= 0) cycle
      enddo
    enddo
    close (unit = 2)
  endif
!
! Class2 states
!
  if (flagclass2) then
    open (unit = 2, file = c2file, status = 'old', iostat = istat)
    if (istat /= 0) call read_error(c2file, istat)
    nclass2(Zix, Nix) = nfisbar(Zix, Nix) - 1
    do i = 1, nclass2(Zix, Nix)
      read(2, '(4x, i4)', iostat = istat) nfisc2hb(Zix, Nix, i)
      if (istat /= 0) cycle
      call range_integer_error(c2file, nfisc2hb(Zix, Nix, i), 0, numlev, index1 = Z, name1 = 'Z', index2 = A, name2 = 'A', &
 &        index3 = i, name3 = 'i')
      do j = 1, nfisc2hb(Zix, Nix, i)
        read(2, '(4x, f11.6, f6.1, i5)', iostat = istat) efisc2hb(Zix, Nix, i, j), jfisc2hb(Zix, Nix, i, j), &
 &        pfisc2hb(Zix, Nix, i, j)
        if (istat /= 0) cycle
      enddo
    enddo
    close (unit = 2)
  endif
!
! ************************* Default parameters *************************
!
  if (fwidth(Zix, Nix, 1) == 0.) fwidth(Zix, Nix, 1) = 1.
  if (fwidth(Zix, Nix, 2) == 0.) fwidth(Zix, Nix, 2) = 0.6
  if (nfisbar(Zix, Nix) == 1 .and. fbarrier(Zix, Nix, 1) == 0.) then
    fbarrier(Zix, Nix, 1) = fbarrier(Zix, Nix, 2)
    fwidth(Zix, Nix, 1) = fwidth(Zix, Nix, 2)
  endif
  do i = 1, numbar
    minertia(Zix, Nix, i) = Rtransmom(Zix, Nix, i) * Irigid(Zix, Nix, i)
    minertc2(Zix, Nix, i) = Rclass2mom(Zix, Nix, i) * Irigid(Zix, Nix, i)
  enddo
!
! ********** Rotational bands on transition and class2 states **********
!
! rotband   : subroutine to build rotational bands on transition states
! rotclass2 : subroutine to build rotational bands on class2 states
!
  call rotband(Zix, Nix)
  if (flagclass2) call rotclass2(Zix, Nix)
  return
end subroutine fissionpar
! Copyright A.J. Koning 2021
