subroutine integral
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculate effective cross section for integral spectrum
!
! Author    : Arjan Koning
!
! 2023-03-19: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl            ! single precision kind
! All global variables
!   numenin        ! number of incident energies
!   numflux        ! number of integral experiments
!   numlev         ! maximum number of discrete levels
!   numP           ! number of points in reconstructed data file
! Variables for output
!   fluxname       ! name of integral spectrum
!   integralexp    ! experimental effective cross section
!   Nflux          ! number of reactions with integral data
!   xsfluxfile     ! TALYS cross section file for integral data
! Variables for levels
!   Liso         ! isomeric number of target
! Variables for main input
!   Atarget        ! mass number of target nucleus
!   k0             ! index of incident particle
!   Starget        ! symbol of target nucleus
! Constants
!   parsym         ! symbol of particle
! Variables for files
!   path           ! directory containing files to be read
!
! *** Declaration of local data
!
  implicit none
  logical            :: lexist               ! logical to determine existence
  character(len=1)   :: isostring(numflux)   ! isomer string
  character(len=3)   :: Astring              ! mass string
  character(len=3)   :: ext                  ! file extension
  character(len=8)   :: reac                 ! reaction string
  character(len=8)   :: reacstr(numflux)     ! reaction string
  character(len=132) :: word(40)             ! word
  character(len=132) :: fluxfile             ! file with experimental integral spectrum
  character(len=132) :: xsfile               ! cross section file
  character(len=132) :: key
  character(len=200) :: line                 ! line
  integer            :: i                    ! counter
  integer            :: igr                  ! 
  integer            :: is                   ! isotope counter: -1=total, 0=ground state 1=isomer
  integer            :: istat                ! logical for file access
  integer            :: k                    ! designator for particle
  integer            :: keyix
  integer            :: L                    ! counter
  integer            :: nen                  ! energy counter
  integer            :: nen0                 ! energy counter
  integer            :: Nspec                ! number of spectral energies
  integer            :: Nxs                  ! number of incident energies
  real(sgl)          :: dum                  ! help variable
  real(sgl)          :: Ea1                  ! help variable
  real(sgl)          :: Eb1                  ! help variable
  real(sgl)          :: Efl                  ! energy of bin for flux
  real(sgl)          :: Eflux(0:numenin)     ! energy of bin for flux
  real(sgl)          :: Elow                 ! lower energy of bin for flux
  real(sgl)          :: Eup                  ! upper energy of bin for flux
  real(sgl)          :: Exs(0:numP)          ! energy of cross section file
  real(sgl)          :: fluxsum              ! check for conservation of flux per P,J,j,l
  real(sgl)          :: fspec(numenin)       ! spectrum values
  real(sgl)          :: ratio                ! if 0<ratio<1 then x is between xl1 and xl2
  real(sgl)          :: sp                   ! help variable
  real(sgl)          :: xs(0:numP)           ! help variable
  real(sgl)          :: xsa                  ! help variable
  real(sgl)          :: xsb                  ! help variable
  real(sgl)          :: xseff                ! effective cross section
  real(sgl)          :: xsexp                ! experimental effective cross section
  real(sgl)          :: xsf                  ! filename
!
! ********* Read reaction channels with integral cross sections ********
!
! This is for the case where a simple 'integral y' is given in the input file, i.e. no explicit information per case
!
  if (Nflux == 0 .and. Liso ==.0) then
    Astring = '   '
    write(Astring(1:3), '(i3.3)') Atarget
    xsfile = trim(path)//'integral/sacs/'//trim(Starget)//Astring//'.sacs'
    open (unit=2, file=xsfile, status='old', iostat = istat)
    if (istat == 0) then
      i = 1
      do
        read(2, '(a)', iostat = istat) line
        if (istat /= 0) exit
        call getkeywords(line(1:132), word)
        reacstr(i) = trim(word(3))
        isostring(i) = trim(word(4))
        fluxname(i) = trim(word(5))
        read(line(41:49), *) integralexp(i)
!       read(line(51:59), *) dintegralexp(i)
        i = i + 1
      enddo
      close (unit = 2)
      Nflux = i - 1
    endif
  endif
!
! Determine corresponding TALYS output file
!
  do i = 1,Nflux
    reac = ''
    ext = 'tot'
    if (trim(reacstr(i)) == 'n,g') reac = 'xs000000'
    if (trim(reacstr(i)) == 'n,p') reac = 'xs010000'
    if (trim(reacstr(i)) == 'n,d') reac = 'xs001000'
    if (trim(reacstr(i)) == 'n,t') reac = 'xs000100'
    if (trim(reacstr(i)) == 'n,h') reac = 'xs000010'
    if (trim(reacstr(i)) == 'n,a') reac = 'xs000001'
    if (trim(reacstr(i)) == 'n,2n') reac = 'xs200000'
    if (trim(reacstr(i)) == 'n,np') reac = 'xs110000'
    if (trim(reacstr(i)) == 'n,na') reac = 'xs100001'
    if (trim(reacstr(i)) == 'n,2a') reac = 'xs000002'
    if (trim(reacstr(i)) == 'n,2na') reac = 'xs200001'
    if (trim(reacstr(i)) == 'n,3n') reac = 'xs300000'
    if (trim(reacstr(i)) == 'n,4n') reac = 'xs400000'
    if (trim(reacstr(i)) == 'n,pa') reac = 'xs010001'
    if (trim(reacstr(i)) == 'n,da') reac = 'xs001001'
    if (trim(reacstr(i)) == 'n,nt') reac = 'xs100100'
    if (trim(reacstr(i)) == 'n,n2a') reac = 'xs100002'
    if (trim(reacstr(i)) == 'n,3na') reac = 'xs300001'
    if (isostring(i) == 'g') ext = 'L00'
    if (isostring(i) == 'm') then
      ext = 'L  '
      write(ext(2:3), '(i2.2)') L
      do L = 1,numlev
        write(ext(2:3), '(i2.2)') L
        inquire(file = reac// '.' // ext,exist = lexist)
        if (lexist) exit
      enddo
    endif
    xsfluxfile(i) = reac// '.' //ext
  enddo
!
! ********* Read integral spectrum from experimental database **********
!
  open (unit = 1, file = 'integral.dat', status = 'replace')
  write(1, '("# ", a1, " + ", a, ": Effective cross sections from integral data")') parsym(k0), trim(targetnuclide)
  write(1, '("# Channel      Spectrum        Eff. xs (b) Exp. xs (b)         Ratio")')
  do i = 1, Nflux
    fluxfile=trim(path)//'integral/spectra/'//trim(fluxname(i))//'.txt'
    open (unit = 2, file = fluxfile, status = 'old', iostat = istat)
    if (istat == 0) then
      read(2,'(/)')
      nen=0
      fluxsum=0
      do
        read(2,*,iostat=istat) igr,Eup,elow,dum,sp
        if (istat.ne.0) exit
        nen=nen+1
        fspec(nen)=sp
        Eflux(nen)=0.5*(Eup+Elow)*1.e-6
        fluxsum=fluxsum+sp
      enddo
      close (unit = 2)
      Nspec=nen
    else
      write(*, '(" TALYS-warning: integral spectrum file ", a, "does not exist")') trim(fluxfile)
      cycle
    endif
!
! *************** Read cross sections from TALYS output files **********
!
    is = 1
  200   open (unit = 3, file = xsfluxfile(i), status = 'old', iostat = istat)
    if (istat == 0) then
      Exs(0) = 0.
      xs(0) = 0.
      do
        read(3,'(a)',iostat = istat) line
        if (istat == -1) exit
        key='entries'
        keyix=index(line,trim(key))
        if (keyix > 0) then
          read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Nxs
          read(3,'(/)')
          do nen = 1, Nxs
            read(3, *, iostat=istat ) Exs(nen), xs(nen)
            if (istat == -1) exit
          enddo
          exit
        endif
      enddo
      close (unit = 3)
    else
      if (is < numlev) then
        is = is + 1
        do k = 1, 20
          if (xsfluxfile(i)(k:k) == 'L') then
            write(xsfluxfile(i)(k+1:k+2), '(i2.2)') is
            goto 200
          endif
        enddo
      endif
      write(*, '(" TALYS-warning: cross section file ", a, " does not exist")') trim(xsfluxfile(i))
      cycle
    endif
!
! *************** Calculate effective cross section by folding *********
!
! Interpolate cross section grid on the grid of the integral spectrum
!
    xseff = 0.
    do nen = 1, Nspec
      Efl = Eflux(nen)
      if (Efl >= Exs(0) .and. Efl <= Exs(Nxs)) then
        call locate(Exs, 0, Nxs, Efl, nen0)
        Ea1 = Exs(nen0)
        Eb1 = Exs(nen0 + 1)
        xsa = xs(nen0)
        xsb = xs(nen0 + 1)
        call pol1(Ea1, Eb1, xsa, xsb, Efl, xsf)
        xseff = xseff + xsf * fspec(nen)
      endif
    enddo
    xseff = 0.001 * xseff / fluxsum
    xsexp = integralexp(i)
    if (xsexp > 0.) then
      ratio = xseff / xsexp
    else
      ratio = 1.
    endif
    write(1, '(2a15, 2es12.5, f15.5)') xsfluxfile(i), fluxname(i), xseff, xsexp, ratio
  enddo
  close (unit = 1)
  return
end subroutine integral
! Copyright A.J. Koning 2021
