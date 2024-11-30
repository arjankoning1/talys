    subroutine resonance
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Reconstruction and broadening of resonance information
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
!   sgl          ! single precision kind
!   dbl          ! double precision kind
! All global variables
!   numP         ! number of points in reconstructed data file
! Variables for basic reaction
!   flagastro    ! flag for calculation of astrophysics reaction rate
! Variables for astrophysics
!   nTmax        ! effective number of temperatures for Maxwellian
! Variables for compound reactions
!   flaggroup    ! flag for output of low energy groupwise cross sections
!   reslib       ! library with resonance parameters
!   Tres         ! temperature for broadening low energy cross sections
! Variables for main input
!   Atarget      ! mass number of target nucleus
!   k0           ! index of incident particle
!   Ztarget      ! charge number of target nucleus
! Variables for nuclides
!   T9           ! Temperature grid in 10 **9 K
! Constants
!   isochar      ! symbol of isomer
!   nuc          ! symbol of nucleus
!   parN         ! neutron number of particle
!   parsym       ! symbol of particle
!   parZ         ! charge number of particle
! Variables for files
!   path         ! directory containing files to be read
! Variables for levels
!   Liso         ! isomeric number of target
! Variables for masses
!   redumass     ! reduced mass
!   specmass     ! specific mass for residual nucleus
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist           ! logical to determine existence
  character(len=9)  :: afile            ! TALYS file with mass number
  character(len=12) :: xsfile           ! file with channel cross sections
  character(len=12) :: rpfile           ! file with residual production cross sections
  character(len=12) :: pfile            ! parameter file
  character(len=72) :: resfile          ! file with residual production cross sections
  character(len=132):: infile           ! input file
  character(len=72) :: outfile          ! output file
  character(len=132) :: headstring(100) ! string for first part of filename
  character(len=132) :: line
  character(len=132) :: key
  character(len=80) :: string           ! line with parameter value
  integer           :: i                ! counter
  integer           :: ifile            ! file
  integer           :: istat            ! error code
  integer           :: iT               ! counter
  integer           :: j                ! counter
  integer           :: k                ! designator for particle
  integer           :: keyix
  integer           :: MF               ! MF-number
  integer           :: MT               ! MT-number
  integer           :: nen              ! energy counter
  integer           :: Nh               ! help variable
  integer           :: Nheader
  integer           :: Nix              ! neutron number index for residual nucleus
  integer           :: nlin             ! number of lines
  integer           :: NP0              ! number of points
  integer           :: NPr              ! number of points
  integer           :: NR               ! number of interpolation ranges
  integer           :: NT               ! help variable
  integer           :: Ntemp            ! number of temperatures
  integer           :: Ntot             ! number of nucleon units in exit channel
  integer           :: Zix              ! charge number index for residual nucleus
  real(sgl)         :: am               ! reduced mass
  real(sgl)         :: convfac          ! conversion factor for the average velocity at  temperature T9
  real(sgl)         :: E(numP)          ! incident energy
  real(sgl)         :: ee               ! energy
  real(sgl)         :: ee0              ! incident energy
  real(sgl)         :: ee1              ! help variable
  real(sgl)         :: Et(numP)         ! energy
  real(sgl)         :: fact0            ! reaction rate factor for photons
  real(sgl)         :: fex              ! help variable
  real(sgl)         :: smacs            ! total conversion factor for the average velocity at  temperature T9
  real(sgl)         :: temp             ! nuclear temperature
  real(sgl)         :: x(3)             ! help variable
  real(sgl)         :: xs(numP)         ! help variable
  real(sgl)         :: xst(numP)        ! help variable
  real(sgl)         :: y(3)             ! coordinates of the point to test
  real(dbl)         :: dE               ! help variable
  real(dbl)         :: fact             ! reaction rate factor for particles
  real(dbl)         :: rate             ! Maxwellian reaction rate
  real(dbl)         :: sum              ! help variable
  real(dbl)         :: term             ! help variable
  real(dbl)         :: term0            ! help variable
  real(dbl)         :: xxs              ! cross section
!
!
!
  write(*,'()')
  headstring=''
  if (Liso == 0) then
    afile = '000.res. '
  else
    afile = '000i.res.'
    write(afile(4:4), '(a1)') isochar(min(Liso,numisom))
  endif
  write(afile(1:3), '(i3.3)') Atarget
  resfile = 'n-'//trim(nuc(Ztarget))//trim(afile)//trim(reslib)
  infile = trim(path)//'resfiles/'//trim(reslib)//'/' //trim(resfile)
  inquire (file = infile, exist = lexist)
  if ( .not. lexist) then
    write(*, '(" TALYS-warning: ", a, " does not exist")') trim(infile)
    return
  endif
  open(unit=41,file=infile,status='unknown')
  open(unit=42,file=resfile,status='unknown')
  do
    read(41,'(a)',iostat = istat) string
    if (istat == -1) exit
    write(42,'(a)') string
  enddo
  close (41)
  close (42)
  if (flagastro) then
    fact0 = 3.951776e+17
    Zix = parZ(k0)
    Nix = parN(k0)
    am = redumass(Zix, Nix, k0)
    convfac = 2.45484e+5 / sqrt(am)
    open(unit = 2, file = 'astrorateres.g', status = 'replace')
    write(2, '("# Reaction rate for ", a, "(", a1, ",g)")') trim(targetnuclide), parsym(k0)
    write(2, '("#    T     kT[keV]     Rate       MACS")')
    Ntemp = nTmax
  else
    Ntemp = 1
  endif
Loop1:  do iT = 1, Ntemp
    if (flagastro) then
      temp = T9(iT) * 1.e9
    else
      temp = Tres
    endif
    if (flaggroup) then
      outfile = resfile(1:11)//'group'
    else
      if (temp > 0.) then
        outfile = resfile(1:11)//'point'
      else
        outfile = resfile(1:11)//'point0'
      endif
    endif
    call prepro(resfile, outfile, temp, flaggroup)
    open (unit = 1, file = outfile, status = 'old')
Loop2: do
      read(1, '(a80)', iostat = istat) string
      if (istat /= 0) exit
      read(string(71:72), '(i2)') MF
      if (MF == 3) then
        read(string(73:75), '(i3)', iostat = istat) MT
        if (istat /= 0) exit
        read(1, '(a80)', iostat = istat) string
        if (istat /= 0) exit
        read(string(45:55), '(i11)', iostat = istat) NR
        if (istat > 0) exit
        read(string(56:66), '(i11)', iostat = istat) NPr
        if (istat > 0) exit
        nlin = 1 + (NR - 1) / 3
        do i = 1, nlin
          read(1, '(a80)', iostat = istat) string
          if (istat /= 0) exit Loop2
        enddo
        nlin = 1 + (NPr - 1) / 3
        k = 0
        do i = 1, nlin
          read(1, '(6e11.6)', iostat = istat) (x(j), y(j), j = 1, 3)
          if (istat /= 0) exit loop2
          do j = 1, 3
            k = k + 1
            E(k) = x(j) * 1.e-6
            xs(k) = y(j) * 1.e3
          enddo
        enddo
        NPr = NPr - 1
        NP0 = NPr
        do k = NP0, 1, - 1
          if (xs(k) > 0.) then
            NPr = k
            exit
          endif
        enddo
        read(1, '(a80)', iostat = istat) string
        if (istat /= 0) exit Loop2
        if (MT == 1) xsfile = 'total.tot'
        if (MT == 2) xsfile = 'elastic.tot'
        if (MT == 18) xsfile = 'fission.tot'
        rpfile = 'rp000000.tot'
        if (MT == 102) then
          xsfile = 'xs000000.tot'
          write(rpfile(3:5), '(i3.3)') Ztarget
          write(rpfile(6:8), '(i3.3)') Atarget+1
        endif
        if (MT == 103)  then
          xsfile = 'xs010000.tot'
          write(rpfile(3:5), '(i3.3)') Ztarget-1
          write(rpfile(6:8), '(i3.3)') Atarget
        endif
        if (MT == 107) then
          xsfile = 'xs000001.tot'
          write(rpfile(3:5), '(i3.3)') Ztarget-2
          write(rpfile(6:8), '(i3.3)') Atarget-3
        endif
!
! Read TALYS output files
!
        do ifile = 1, 2
          if (ifile == 1) then
            pfile = xsfile
          else
            if (MT < 102) cycle
            pfile = rpfile
          endif
          inquire (file = pfile, exist = lexist)
          if ( .not. lexist) cycle Loop1
          open (unit = 2, file = pfile, status = 'old')
          k=0
          Nheader=0
          do
            read(2,'(a)',iostat = istat) line
            if (istat == -1) exit
            k=k+1
            headstring(k)=line
            key='entries'
            keyix=index(line,trim(key))
            if (keyix > 0) then
              read(line(keyix+len_trim(key)+2:80),*, iostat = istat) NT
              if (istat == -1) exit
              do i= 1, 2
                read(2,'(a)',iostat = istat) line
                k=k+1
                headstring(k)=line
              enddo
              Nheader=k
              do k = 1, NT
                read(2, *, iostat = istat) Et(k), xst(k)
                if (istat /= 0) exit
              enddo
              exit
            endif
          enddo
          close (2)
          Nh = NT + 1
          do k = 1, NT
            if (Et(k) > E(NPr)) then
              Nh = k
              exit
            endif
          enddo
          Ntot = NPr + NT - Nh + 1
          write(headstring(Nheader-2)(14:20), '(i7)') Ntot
          inquire (file = pfile, exist = lexist)
!
! Write final output files
!
          open (unit = 2, file = pfile, status = 'replace')
            do i = 1, Nheader
            write(2, '(a)') trim(headstring(i))
          enddo
          do k = 1, NPr
            write(2, '(2es15.6)') E(k), xs(k)
          enddo
          do k = Nh, NT
            write(2, '(2es15.6)') Et(k), xst(k)
          enddo
          close(2)
        enddo
      endif
    enddo Loop2
    close(unit = 1)
    if (flagastro) then
      fact = 3.7335e+10 / (sqrt(am * (temp **3)))
      sum = 0.
      term0 = 0.
      do nen = 1, Ntot
        ee = real(Et(nen) * specmass(Zix, Nix, k0))
        ee1 = ee
        if (nen > 1) then
          ee0 = real(Et(nen - 1) * specmass(Zix, Nix, k0))
        else
          ee0 = ee
        endif
        xxs = xst(nen) / 1000.
        fex = 11.605 * ee / temp
        if (fex > 80.) cycle
        term = fact * ee * xxs * exp( - fex)
        dE = (ee1 - ee0)
        if (nen > 1) sum = sum + (term + term0) / 2. * dE
        term0 = term
      enddo
      rate = sum
      smacs = rate / convfac / sqrt(temp)
      write(2, '(f8.4, f9.3, 2es12.5)') temp, temp*0.086173e3, rate, smacs
    endif
  enddo Loop1
  if (flagastro) close(unit = 2)
  return
end subroutine resonance
! Copyright A.J. Koning 2021
