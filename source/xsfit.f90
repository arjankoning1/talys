subroutine xsfit(Z, A)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Adjusted parameters to fit cross sections
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     25-10-2024   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! use A0_kinds_mod, only: & ! Definition of single and double precision variables
!              sgl, &          ! single precision kind
!              dbl             ! double precision kind
! use A0_talys_mod, only: & ! All global variables
!              alphaomp, &     ! alpha optical model
!              Cstrip, &       ! adjustable parameter for stripping / pick - up react
!              ctable, &       ! constant to adjust tabulated level densities
!              fbaradjust, &   ! adjustable factor for fission parameters
!              fisadjust, &    ! logical for energy - dependent fission adjustment
!              fismodel, &     ! fission model alternative fission model for default barriers
!              flagcol, &      ! flag for collective enhancement of level density
!              fwidth, &       ! width of fission barrier
!              fwidthadjust, & ! adjustable factor for fission parameters
!              gamadjust, &    ! logical for energy - dependent gamma adjustment
!              iso, &          ! counter for isotope
!              ldadjust, &     ! logical for energy - dependent level density adjus
!              ldmodel, &      ! level density model for direct radiative capture
!              Ltarget, &      ! excited level of target
!              nin, &          ! counter for incident energy
!              nuc, &          ! symbol of nucleus
!              ompadjustp, &   ! flag for local optical model parameter adjustmen
!              path, &         ! directory containing files to be read
!              ptable, &       ! constant to adjust tabulated level densities
!              rvadjust, &     ! adjustable factor for OMP (default 1.)
!              s2adjust, &     ! adjustable constant (Z, A, barrier - dependent) for
!              strength, &     ! E1 strength function model
!              vfis, &         ! adjustable factor for fission path height
!              vfiscor         ! adjustable factor for fission path height
!
! *** Declaration of local data
!
  implicit none

  integer, parameter :: numfit=11            !
  logical            :: first                ! tracks whether a TALYS run is the first one
  logical            :: flagassign           ! flag to assign value or not
  logical            :: lexist               ! logical to determine existence
  character(len=1)   :: yesno                ! y or n function
  character(len=1)   :: ch                   ! character
  character(len=2)   :: ele
  character(len=132)   :: colf                 !
  character(len=132)   :: proj                 !
  character(len=132)   :: elf                  !
  character(len=132) :: cval                 ! character value
  character(len=132) :: ffile(numfit)        !
  character(len=132) :: key                  ! keyword
  character(len=132) :: keystring(20)        !
  character(len=132) :: line                 ! input line
  character(len=132) :: parfile              !
  character(len=132) :: value                ! value or string
  character(len=132) :: word(40)             ! words on input line
  integer            :: A                    ! mass number of target nucleus
  integer            :: aomp                 !
  integer            :: class                ! input class
  integer            :: fis                  !
  integer            :: fitpar(numfit)       !
  integer            :: i                    ! counter
  integer            :: ia                   ! mass number from abundance table
  integer            :: ibar                 ! fission barrier
  integer            :: igr                  ! giant resonance
  integer            :: irad                 ! variable to indicate M(=0) or E(=1) radiation
  integer            :: istat                ! logical for file access
  integer            :: ival                 ! integer value
  integer            :: k                    ! designator for particle
  integer            :: keyix
  integer            :: ldm                  !
  integer            :: Lt                   !
  integer            :: lval                 ! multipolarity
  integer            :: Nix                  ! neutron number index for residual nucleus
  integer            :: Npar                 ! number of parameters
  integer            :: psf                  !
  integer            :: type                 ! particle type
  integer            :: Z                    ! charge number of target nucleus
  integer            :: Zix                  ! charge number index for residual nucleus
  real(sgl)          :: val                  ! real value
!
! ************************ Initialization ******************************
!
  fitpar = 0
  ffile = 'x'
  if (flagnnfit) fitpar(1) = 1
  if (flagngfit) fitpar(2) = 1
  if (flagnafit) fitpar(3) = 1
  if (flagnffit) fitpar(4) = 1
  if (flaganfit) fitpar(5) = 1
  if (flagdnfit) fitpar(6) = 1
  if (flagpnfit) fitpar(7) = 1
  if (flaggnfit) fitpar(8) = 1
  if (flaggamgamfit) fitpar(9) = 1
  if (flagmacsfit) fitpar(10) = 1
  if (flagndfit) fitpar(11) = 1
  ffile(1) = 'nn.par'
  ffile(2) = 'ng.par'
  ffile(3) = 'na.par'
  ffile(4) = 'nf.par'
  ffile(5) = 'an.par'
  ffile(6) = 'dn.par'
  ffile(7) = 'pn.par'
  ffile(8) = 'gn.par'
  ffile(9) = 'gamgam.par'
  ffile(10) = 'macs.par'
  ffile(11) = 'nd.par'
!
! ************************* Read parameters ****************************
!
  first = .true.
  do k = 1, numfit
    if (fitpar(k) == 0) cycle
    parfile = trim(path)//'best/fits/'//trim(ffile(k))
    inquire (file = parfile, exist = lexist)
    if ( .not. lexist) cycle
    Lt = 0
    open (unit = 2, file = parfile, status = 'old')
Loop1: do
      do
        read(2, '(a)', iostat = istat) line
        if (istat == - 1) exit Loop1
        key = 'projectile'
        keyix=index(line,trim(key))
        if (keyix > 0) read(line(keyix+len_trim(key):132),'(a)', iostat = istat) proj
        key = 'element'
        keyix=index(line,trim(key))
        if (keyix > 0) read(line(keyix+len_trim(key):132),'(a)', iostat = istat) elf
        key = 'mass'
        keyix=index(line,trim(key))
        if (keyix > 0) read(line(keyix+len_trim(key):132),*, iostat = istat) ia
        key = 'ldmodel'
        keyix=index(line,trim(key))
        if (keyix > 0) read(line(keyix+len_trim(key):132),*, iostat = istat) ldm
        key = 'colenhance'
        keyix=index(line,trim(key))
        if (keyix > 0) read(line(keyix+len_trim(key):132),'(a)', iostat = istat) colf
        key = 'strength'
        keyix=index(line,trim(key))
        if (keyix > 0) read(line(keyix+len_trim(key):132),*, iostat = istat) psf
        key = 'fismodel'
        keyix=index(line,trim(key))
        if (keyix > 0) read(line(keyix+len_trim(key):132),*, iostat = istat) fis
        key = 'alphaomp'
        keyix=index(line,trim(key))
        if (keyix > 0) read(line(keyix+len_trim(key):132),*, iostat = istat) aomp
        key = '##parameters'
        keyix=index(line,trim(key))
        if (keyix > 0) then
          read(line(keyix+len_trim(key):132),*, iostat = istat) Npar
          do i = 1, Npar
            read(2, '(a)') keystring(i)
          enddo
        endif
        key = '#####'
        keyix=index(line,trim(key))
        if (keyix > 0) exit
      enddo
      ele = trim(adjustl(elf))
      ch = ele(1:1)
      if (ch  >= 'a'  .and. ch <= 'z') ele(1:1) = achar(iachar(ch) - 32)
      if (trim(ele) /= trim(nuc(Z))) cycle
      if (ia /= A) cycle
      if (Lt /= Ltarget) cycle
      if (ldm /= ldmodel(0, 0)) cycle
      if (k /= 4 .and. trim(adjustl(colf)) /= yesno(flagcol(0, 0))) cycle
      if ((k == 2 .or. k == 8 .or. k == 9 .or. k==10) .and. psf /= strength) cycle
      if ((k == 3 .or. k == 5) .and. aomp /= alphaomp) cycle
      if (k == 4 .and. fis /= fismodel) cycle
      if (first) then
        open (unit = 1, file = 'adjust.dat', status = 'replace')
        first = .false.
        write(1, '("# Adjusted parameters from fitting")')
      else
        open (unit = 1, file = 'adjust.dat', status = 'old', position = 'append')
      endif
      do i = 1, Npar
        write(1, '(a)') trim(keystring(i))
      enddo
      close(1)
      do i = 1, Npar
        line = keystring(i)
        call getkeywords(line, word)
        key = word(1)
        value = word(2)
        ch = word(2)(1:1)
        Zix = 0
        Nix = 0
        type = 0
        lval = 0
        ibar = 0
        igr = 1
        if (key == 'rvadjust') then
          class = 6
          call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
          if (flagassign) then
            rvadjust(type) = val
            ompadjustp(type) = .true.
          endif
          cycle
        endif
        if (key == 'rwdadjust') then
          class = 6
          call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
          if (flagassign) then
            rwdadjust(type) = val
            ompadjustp(type) = .true.
          endif
          cycle
        endif
        if (key == 'awdadjust') then
          class = 6
          call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
          if (flagassign) then
            awdadjust(type) = val
            ompadjustp(type) = .true.
          endif
          cycle
        endif
        if (key == 'gadjust') then
          class = 1
          call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
          if (flagassign) gadjust(Zix, Nix) = val
          cycle
        endif
        if (key == 'wtable') then
          class = 5
          call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
          if (flagassign) then
            wtable(Zix, Nix, irad, lval) = val
            gamadjust(Zix, Nix) = .true.
          endif
          cycle
        endif
        if (key == 'fisbaradjust') then
          ibar = 1
          class = 3
          call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
          if (flagassign) then
            fbaradjust(Zix, Nix, ibar) = val
            fisadjust(Zix, Nix) = .true.
          endif
          cycle
        endif
        if (key == 'fishwadjust') then
          ibar = 1
          class = 3
          call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
          if (flagassign) then
            fwidthadjust(Zix, Nix, ibar) = val
            fisadjust(Zix, Nix) = .true.
          endif
          cycle
        endif
        if (key == 'vfiscor') then
          class = 1
          call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
          if (flagassign) vfiscor(Zix, Nix) = val
          cycle
        endif
        if (key == 'ctable') then
          class = 3
          call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
          if (flagassign) then
            ctable(Zix, Nix, ibar) = val
            ldadjust(Zix, Nix) = .true.
          endif
          cycle
        endif
        if (key == 'ptable') then
          class = 3
          call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
          if (flagassign) then
            ptable(Zix, Nix, ibar) = val
            ldadjust(Zix, Nix) = .true.
          endif
          cycle
        endif
        if (key == 'ctableadjust') then
          class = 3
          call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
          if (flagassign) then
            ctableadjust(Zix, Nix, ibar) = val
            ldadjust(Zix, Nix) = .true.
          endif
          cycle
        endif
        if (key == 'ptableadjust') then
          class = 3
          call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
          if (flagassign) then
            ptableadjust(Zix, Nix, ibar) = val
            ldadjust(Zix, Nix) = .true.
          endif
          cycle
        endif
        if (key == 's2adjust') then
          class = 3
          call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
          if (flagassign) s2adjust(Zix, Nix, ibar) = val
          cycle
        endif
        if (key == 'risomer') then
          class = 1
          call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
          if (flagassign) risomer(Zix, Nix) = val
          cycle
        endif
        if (key == 'cstrip') then
          class = 6
          call getvalues(class, word, Zix, Nix, type, ibar, irad, lval, igr, val, ival, cval, flagassign)
          if (flagassign) Cstrip(type) = val
          cycle
        endif
      enddo
      exit
    enddo Loop1
  enddo 
  return
end subroutine xsfit
! Copyright A.J. Koning 2023
