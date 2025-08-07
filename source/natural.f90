  subroutine natural
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculation for natural elements
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
!   numen2          ! maximum number of outgoing energies
!   numia           ! maximum number of alphas in channel description
!   numid           ! maximum number of deuterons in channel description
!   numih           ! maximum number of helions in channel description
!   numin           ! maximum number of neutrons in channel description
!   numip           ! maximum number of protons in channel description
!   numiso          ! maximum number of isotopes per element
!   numit           ! maximum number of tritons in channel description
!   numlev          ! maximum number of discrete levels
!   numN            ! maximum number of neutrons from initial compound nucleus
!   numtime         ! number of time points
!   numZ            ! maximum number of protons from initial compound nucleus
! Variables for output
!   filechannels    ! flag for exclusive channel cross sections on separate file
!   fileelastic     ! flag for elastic angular distribution on separate file
!   filefission     ! flag for fission cross sections on separate file
!   fileresidual    ! flag for residual production cross sections on separate file
!   filespectrum    ! designator for spectrum on separate file
!   filetotal       ! flag for total cross sections on separate file
! Variables for medical isotope production
!   flagprod        ! flag for isotope production
! Variables for basic reaction
!   flagmassdis     ! flag for calculation of fission fragment mass yields
!   flagrecoil      ! flag for calculation of recoils
! Variables for numerics
!   maxchannel      ! maximal number of outgoing particles in individual channel description
!   nangle          ! number of angles
! Variables for input energies
!   eninc           ! incident energy in MeV
!   Ninc          ! number of incident energies
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
!   Starget         ! symbol of target nucleus
!   Ztarget         ! charge number of target nucleus
! Variables for fission
!   flagffevap      ! flag for calculation of particle evaporation from fissi
! Variables for abundance
!   abun            ! Natural abundance
!   isonum          ! number of isotopes
!   isotope         ! isotope number of residual nucleus
! Constants
!   iso             ! counter for isotope
!   natstring       ! string extension for file names
!   nuc             ! symbol of nucleus
!   parname         ! name of particle
!   parsym          ! symbol of particle
! Variables for isotope production
!   Tgrid           ! time
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical            :: elexist                       ! logical to determine existence of elastic scattering file
  logical            :: fissionexist                  ! logical to determine existence of fission file
  logical            :: fyexist                       ! flag for existence of fission yields
  logical            :: lexist                        ! logical to determine existence
  logical            :: resexist                      ! logical to determine existence of residual production file
  logical            :: specexist                     ! logical to determine existence of spectrum file
  logical            :: xsexist                       ! flag for existence of cross section file
  character(len=3)   :: massstring
  character(len=6)  :: finalnuclide !
  character(len=10)  :: totstring                     ! string for totfile
  character(len=12)  :: Estr
  character(len=13)   :: prodfile                      ! file with total particle production cross sections
  character(len=15)  :: fisfile                       ! fission file
  character(len=15)  :: Yfile                         ! file with production yields
  character(len=16)  :: resfile                       ! file with residual production cross sections
  character(len=16)  :: xsfile                        ! file with channel cross sections
  character(len=20)  :: totfile                       ! file with total cross sections
  character(len=20)  :: tot0file                      ! file with total cross sections
  character(len=21)  :: discfile                      ! file with elastic scattering angular distribution
  character(len=21)  :: specfile                      ! file with composite particle spectra
  character(len=21)  :: fyfile                        ! file with fission yields
  character(len=29)  :: recfile                       ! file with recoil spectra
  character(len=132) :: str(9)                        ! input line
  character(len=132) :: key
  character(len=200) :: line
  character(len=18)  :: reaction   ! reaction
  character(len=15)  :: col(11)    ! header
  character(len=15)  :: un(11)    ! header
  character(len=80)  :: quantity   ! quantity
  character(len=132) :: topline   ! topline
  character(len=200) :: resline(20) 
  integer            :: Ncol       ! counter
  integer            :: abeg                          ! start of A loop
  integer            :: aend                          ! end of A loop
  integer            :: i                             ! level
  integer            :: ia                            ! mass number from abundance table
  integer            :: iang                          ! running variable for angle
  integer            :: id                            ! counter for deuterons
  integer            :: ih                            ! hole number
  integer            :: in                            ! counter for neutrons
  integer            :: ip                            ! particle number
  integer            :: istat                         ! logical for file access
  integer            :: it                            ! counter for tritons
  integer            :: iz                            ! charge number of residual nucleus
  integer            :: j                             ! counter
  integer            :: jend                          ! end of jloop
  integer            :: k                             ! designator for particle
  integer            :: k2                            ! counter
  integer            :: keyix
  integer            :: n1                            ! number of coordinate grid points
  integer            :: n2                            ! counter
  integer            :: indent
  integer            :: id2
  integer            :: id4
  integer            :: nen                           ! energy counter
  integer            :: nenen                         ! energy counter
  integer            :: neniso(numiso)                ! number of emission energies per isotope
  integer            :: nex                           ! excitation energy bin of compound nucleus
  integer            :: Nk
  integer            :: npart                         ! number of particles in outgoing channel
  integer            :: type                          ! particle type
  integer            :: zbeg                          ! start of Z loop
  integer            :: zend                          ! end of Z loop
  integer            :: MF
  integer            :: MT
  integer            :: ident
  real(sgl)          :: act                           ! activity
  real(sgl)          :: actnat(0:numtime)             ! activity for natural element
  real(sgl)          :: ang(0:nangle)                 ! angle
  real(sgl)          :: E                             ! incident energy
  real(sgl)          :: Ea                            ! start energy of local adjustment
  real(sgl)          :: Eb                            ! end energy of local adjustment
  real(sgl)          :: Efac                          ! help variable
  real(sgl)          :: en(numen2)                    ! incident energy
  real(sgl)          :: enspec(numiso, 0:7*numen2)    ! emission energy
  real(sgl)          :: enspecnat(0:7*numen2)         ! emission energy for natural element
  real(sgl)          :: entmp                         ! help variable
  real(sgl)          :: frac                          ! help variable
  real(sgl)          :: fracnat(0:numtime)            ! fraction for natural element
  real(sgl)          :: Nis                           ! number of produced isotopes
  real(sgl)          :: Nisnat(0:numtime)             ! number of produced isotopes for natural element
  real(sgl)          :: absum                         ! help variable
  real(sgl)          :: natbranch(0:numen2)           ! branching ratio for natural targets
  real(sgl)          :: br                            ! help variable
  real(sgl)          :: xs                            ! help variable
  real(sgl)          :: xs1                           ! help variable
  real(sgl)          :: xs1nat(0:numen2)              ! help variable
  real(sgl)          :: xs2                           ! help variable
  real(sgl)          :: xs2nat(0:numen2)              ! help variable
  real(sgl)          :: xsnat(0:numen2)               ! cross section for natural element
  real(sgl)          :: xsprodnat(numen2)             ! production cross sections for natural element
  real(sgl)          :: xsspec(numiso, 0:7*numen2, 5) ! differential cross section
  real(sgl)          :: xsspecnat(0:7*numen2, 5)      ! differential cross section for natural element
  real(sgl)          :: xst(10)                       ! help variable
  real(sgl)          :: xstotnat(10, numen2)          ! total cross sections for natural element
  real(sgl)          :: xsyieldnat(numen2)            ! yields for natural element
  real(sgl)          :: y                             ! coordinates of the point to test
  real(sgl)          :: Ynat(0:numtime)               ! yield for natural element
!
! **************** Create runs and directories per isotope *************
!
! talysinput   : subroutine for user input and defaults
! talysinitial : subroutine for initialization of nuclear structure
! talysreaction: subroutine with reaction models
!
! For mono-isotopic nuclides we are done.
!
  indent = 0
  id2 = indent + 2
  id4 = indent + 4
  if (isonum == 1) return
!
! do a full TALYS calculation for each isotope
!
  do iso = 2, isonum
    call talysinput
    call talysinitial
    call talysreaction
  enddo
!
! ******** Merge output files in results for natural elements **********
!
  Atarget=0
  targetnuclide=Starget//'0'
!
! 1. Total cross sections
!
  MF = 3
  en = 0.
  xsprodnat = 0.
  xsyieldnat = 0.
  xst = 0.
  xstotnat = 0.
  if (filetotal) then
    do i = 1, isonum
      totfile = 'all.tot'//natstring(i)
      inquire (file = totfile, exist = lexist)
      if (lexist) then
        open (2, file = totfile, status = 'old', iostat = istat)
        if (istat /= 0) call read_error(totfile, istat)
        do 
          read(2,'(a)',iostat = istat) line
          if (istat == -1) exit
          key='entries'
          keyix=index(line,trim(key))
          if (keyix > 0) then
            read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Ninc
            if (istat /= 0) call read_error(totfile, istat, eor = 'continue', eof = 'continue')
            read(2,'(/)')
            do k = 1, Ninc
              read(2, '(11e15.6)', iostat = istat) en(k), (xst(k2), k2 = 1, 10)
              if (istat /= 0) call read_error(totfile, istat, eor = 'continue', eof = 'continue')
              do k2 = 1, 10
                xstotnat(k2, k) = xstotnat(k2, k) + abun(i) * xst(k2)
              enddo
            enddo
            exit
          endif
        enddo
        close (unit = 2)
      endif
    enddo
    un = 'mb'
    col(1)='E'
    un(1)='MeV'
    col(2)='xs'
    col(2)='Non-elastic'
    col(3)='Elastic'
    col(4)='Total'
    col(5)='Compound_elast.'
    col(6)='Shape_elastic'
    col(7)='Reaction'
    col(8)='Compound_nonel.'
    col(9)='Direct'
    col(10)='Preequilibrium'
    col(11)='Direct_capture'
    Ncol=11
    open (1, file = 'all.tot', status = 'replace')
    quantity='cross section'
    reaction='('//parsym(k0)//',all)'
    topline=trim(targetnuclide)//trim(reaction)//' general '//trim(quantity)
    call write_header(indent,topline,source,user,date,oformat)
    call write_target(indent)
    call write_reaction(indent,reaction,0.d0,0.d0,0,0)
    call write_quantity(id2,quantity)
    call write_datablock(id2,Ncol,Ninc,col,un)
    do k = 1, Ninc
      write(1, '(11es15.6)') en(k), (xstotnat(k2, k), k2 = 1, 10)
    enddo
    close (unit = 1)
    do k = 1, numen2
      en(k) = 0.
      do k2 = 1, 10
        xst(k2) = 0.
        xstotnat(k2, k) = 0.
      enddo
    enddo
    if (k0 == 1) then
      jend = 4
    else
      jend = 3
    endif
    quantity='cross section'
    col(1)='E'
    col(2)='xs'
    Ncol=2
    do j = 1, jend
      if (j == 1) totstring = 'total'
      if (j == 2) totstring = 'nonelastic'
      if (j == 3) totstring = 'reaction'
      if (j == 4) totstring = 'elastic'
      tot0file = trim(totstring)//'.tot'
      do i = 1, isonum
        totfile = trim(tot0file) //natstring(i)
        inquire (file = totfile, exist = lexist)
        if (lexist) then
          open (2, file = totfile, status = 'old', iostat = istat)
          if (istat /= 0) call read_error(totfile, istat)
          do 
            read(2,'(a)',iostat = istat) line
            if (istat == -1) exit
            key='entries'
            keyix=index(line,trim(key))
            if (keyix > 0) then
              read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Ninc
              if (istat /= 0) call read_error(totfile, istat, eor = 'continue', eof = 'continue')
              read(2,'(/)')
              do k = 1, Ninc
                read(2, '(2es15.6)', iostat = istat) en(k), xst(j)
                if (istat /= 0) call read_error(totfile, istat, eor = 'continue', eof = 'continue')
                xstotnat(j, k) = xstotnat(j, k) + abun(i) * xst(j)
              enddo
              exit
            endif
          enddo
          close (unit = 2)
        endif
      enddo
      if (j == 1) reaction='('//parsym(k0)//',tot)'
      if (j == 2) reaction='('//parsym(k0)//',non)' 
      if (j == 3) reaction='('//parsym(k0)//',reac)' 
      if (j == 4) reaction='('//parsym(k0)//',el)' 
      MT = 0
      if (j == 1) MT = 1
      if (j == 2) MT = 3
      if (j == 4) MT = 2
      open (1, file = tot0file, status = 'replace')
      topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
      call write_header(indent,topline,source,user,date,oformat)
      call write_target(indent)
      call write_reaction(indent,reaction,0.d0,0.d0,MF,MT)
      call write_quantity(id2,quantity)
      call write_datablock(id2,Ncol,Ninc,col,un)
      do k = 1, Ninc
        write(1, '(2es15.6)') en(k), xstotnat(j, k)
      enddo
      close (unit = 1)
    enddo
!
! 2. Particle production cross sections
!
    do type = 0, 6
      do k = 1, Ninc
        xsprodnat(k) = 0.
        xsyieldnat(k) = 0.
      enddo
      do i = 1, isonum
        prodfile = ' prod.tot'//natstring(i)
        write(prodfile(1:1), '(a1)') parsym(type)
        inquire (file = prodfile, exist = lexist)
        if (lexist) then
          open (2, file = prodfile, status = 'old', iostat = istat)
          if (istat /= 0) call read_error(prodfile, istat)
          do 
            read(2,'(a)',iostat = istat) line
            if (istat == -1) exit
            key='entries'
            keyix=index(line,trim(key))
            if (keyix > 0) then
              read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Ninc
              if (istat /= 0) call read_error(prodfile, istat, eor = 'continue', eof = 'continue')
              read(2,'(/)')
              do k = 1, Ninc
                read(2, '(3es15.6)', iostat = istat) en(k), xs, y
                if (istat /= 0) call read_error(prodfile, istat, eor = 'continue', eof = 'continue')
                xsprodnat(k) = xsprodnat(k) + abun(i) * xs
                xsyieldnat(k) = xsyieldnat(k) + abun(i) * y
              enddo
              exit
            endif
          enddo
          close (unit = 2)
        endif
      enddo
      quantity='cross section'
      col(1)='E'
      col(2)='xs'
      col(3)='Multiplicity'
      un(3)=''
      Ncol=3
      reaction='('//parsym(k0)//',x'//parsym(type)//')'
      topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
      if (type == 0) MT = 202
      if (type == 1) MT = 201
      if (type > 1) MT = 201 + type
      open (1, file = prodfile(1:9), status = 'replace')
      call write_header(indent,topline,source,user,date,oformat)
      call write_target(indent)
      call write_reaction(indent,reaction,0.d0,0.d0,MF,MT)
      call write_quantity(id2,quantity)
      call write_datablock(id2,Ncol,Ninc,col,un)
      do k = 1, Ninc
        write(1, '(3es15.6)') en(k), xsprodnat(k), xsyieldnat(k)
      enddo
      close (unit = 1)
    enddo
  endif
!
! 3. Elastic scattering angular distribution
!
  if (fileelastic) then
    do k = 1, Ninc
      Estr=''
      write(Estr,'(es12.6)') eninc(k)
      xsnat = 0.
      xs1nat = 0.
      xs2nat = 0.
      elexist = .false.
      do i = 1, isonum
        discfile = 'nn        ang.L00'//natstring(i)
        write(discfile(1:2), '(2a1)') parsym(k0), parsym(k0)
        write(discfile(3:10), '(f8.3)') eninc(k)
        write(discfile(3:6), '(i4.4)') int(eninc(k))
        inquire (file = discfile, exist = lexist)
        if (lexist) then
          elexist = .true.
          open (2, file = discfile, status = 'old', iostat = istat)
          do 
            read(2,'(a)',iostat = istat) line
            if (istat == -1) exit
            key='entries'
            keyix=index(line,trim(key))
            if (keyix > 0) then
              read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Ninc
              if (istat /= 0) call read_error(discfile, istat, eor = 'continue', eof = 'continue')
              read(2,'(/)')
              do iang = 0, nangle
                read(2, '(4e15.6)', iostat = istat) ang(iang), xs, xs1, xs2
                if (istat /= 0) call read_error(discfile, istat, eor = 'continue', eof = 'continue')
                xsnat(iang) = xsnat(iang) + abun(i) * xs
                xs1nat(iang) = xs1nat(iang) + abun(i) * xs1
                xs2nat(iang) = xs2nat(iang) + abun(i) * xs2
              enddo
              exit
            endif
          enddo
          close (unit = 2)
        endif
      enddo
      if (elexist) then
        quantity='angular distribution'
        reaction='('//parsym(k0)//',el)'
        topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
        un = 'mb/sr'
        col(1)='Angle'
        un(1)='deg'
        col(2)='xs'
        col(3)='Direct'
        col(4)='Compound'
        Ncol=4
        open (1, file = discfile(1:16), status = 'replace')
        call write_header(indent,topline,source,user,date,oformat)
        call write_target(indent)
        call write_reaction(indent,reaction,0.D0,0.D0,4,2)
        call write_real(id2,'E-incident [MeV]',eninc(k))
        call write_quantity(id2,quantity)
        call write_datablock(id2,Ncol,nangle+1,col,un)
        do iang = 0, nangle
          write(1, '(4es15.6)') ang(iang), xsnat(iang), xs1nat(iang), xs2nat(iang)
        enddo
        close (unit = 1)
      endif
    enddo
  endif
!
! 4. Composite particle spectra
!
  un = 'mb/MeV'
  col(1)='E-out'
  un(1)='MeV'
  col(2)='xs'
  col(3)='Direct'
  col(4)='Preequilibrium'
  col(5)='Multiple_preeq'
  col(6)='Compound'
  Ncol=6
  quantity='emission spectrum'
  do k = 1, Ninc
    Estr=''
    write(Estr,'(es12.6)') eninc(k)
    do type = 1, 6
      if ( .not. filespectrum(type)) cycle
      specexist = .false.
      do k2 = 1, 7*numen2
        enspecnat(k2) = 0.
        do j = 1, 5
          xsspecnat(k2, j) = 0.
        enddo
      enddo
      do i = 1, isonum
        do k2 = 1, numen2
          enspec(i, k2) = 0.
          do j = 1, 5
            xsspec(i, k2, j) = 0.
          enddo
        enddo
!
! In general, the secondary spectra for the various isotopes are all different.
! Therefore, we first read the secondary energy grids and cross sections into memory.
!
        neniso(i) = 0
        specfile = parsym(type)//'spec0000.000.tot'//natstring(i)
        write(specfile(6:13), '(f8.3)') eninc(k)
        write(specfile(6:9), '(i4.4)') int(eninc(k))
        inquire (file = specfile, exist = lexist)
        if (lexist) then
          specexist = .true.
          open (2, file = specfile, status = 'old', iostat = istat)
          if (istat /= 0) call read_error(specfile, istat)
          do 
            read(2,'(a)',iostat = istat) line
            if (istat == -1) exit
            key='entries'
            keyix=index(line,trim(key))
            if (keyix > 0) then
              read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Ninc
              if (istat /= 0) call read_error(specfile, istat, eor = 'continue', eof = 'continue')
              read(2,'(/)')
              do k2 = 1, numen2
                read(2, '(6es15.6)', iostat = istat) enspec(i, k2), (xsspec(i, k2, j), j = 1, 5)
                if (istat /= 0) call read_error(specfile, istat, eor = 'continue', eof = 'continue')
                neniso(i) = neniso(i) + 1
              enddo
              exit
            endif
          enddo
          close (unit = 2)
        endif
      enddo
!
! Make one unifying energy grid by removing double energies and sorting the remaining energies.
!
      if (specexist) then
        nenen = 0
        do i = 1, isonum
Loop1:    do k2 = 1, neniso(i)
            E = enspec(i, k2)
            do nen = 1, nenen
              if (E == enspecnat(nen)) cycle Loop1
            enddo
            nenen = nenen + 1
            enspecnat(nenen) = E
          enddo Loop1
        enddo
        do n1 = 1, nenen
          do n2 = n1, nenen
            if (enspecnat(n1) <= enspecnat(n2)) cycle
            entmp = enspecnat(n1)
            enspecnat(n1) = enspecnat(n2)
            enspecnat(n2) = entmp
          enddo
        enddo
!
! Interpolation and construction of natural spectra.
!
        do nen = 1, nenen
          do i = 1, isonum
            do k2 = 1, neniso(i) - 1
              Ea = enspec(i, k2)
              Eb = enspec(i, k2 + 1)
              if (enspecnat(nen) >= Ea .and. enspecnat(nen) < Eb) then
                Efac = (enspecnat(nen) - Ea) / (Eb - Ea)
                do j = 1, 5
                  xs = xsspec(i, k2, j) + Efac * (xsspec(i, k2 + 1, j) - xsspec(i, k2, j))
                  xsspecnat(nen, j) = xsspecnat(nen, j) + abun(i) * xs
                enddo
                cycle
              endif
            enddo
          enddo
        enddo
        reaction='('//parsym(k0)//',x'//parsym(type)//')'
        topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
        open (1, file = specfile(1:17), status = 'replace')
        call write_header(indent,topline,source,user,date,oformat)
        call write_target(indent)
        call write_reaction(indent,reaction,0.D0,0.D0,6,5)
        call write_real(id2,'E-incident [MeV]',eninc(k))
        call write_quantity(id2,quantity)
        call write_datablock(id2,Ncol,nenen,col,un)
        do nen = 1, nenen
          write(1, '(6es15.6)') enspecnat(nen), (xsspecnat(nen, j), j = 1, 5)
        enddo
        close (unit = 1)
      endif
    enddo
  enddo
!
! 5. Residual production cross sections
!
  if (fileresidual) then
    quantity='cross section'
    col(1)='E'
    col(2)='xs'
    un(2)='mb'
    col(3)='Isomeric_ratio'
    un(3)=''
    reaction='('//parsym(k0)//',x)'
    zbeg = max(Ztarget - numZ - 2, 1)
    zend = Ztarget + 2
    abeg = max(isotope(1) - numN - 2, 1)
    aend = isotope(isonum) + 4
    do iz = zbeg, zend
      do ia = abeg, aend
!
! Total
!
        xsnat = 0.
        resexist = .false.
        do i = 1, isonum
          resfile = 'rp000000.tot'//natstring(i)
          write(resfile(3:8), '(2i3.3)') iz, ia
          inquire (file = resfile, exist = lexist)
          if (lexist) then
            resexist = .true.
            open (2, file = resfile, status = 'old', iostat = istat)
            if (istat /= 0) call read_error(resfile, istat)
            do 
              read(2,'(a)',iostat = istat) line
              if (istat == -1) exit
              key='entries'
              keyix=index(line,trim(key))
              if (keyix > 0) then
                read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Ninc
                if (istat /= 0) call read_error(resfile, istat, eor = 'continue', eof = 'continue')
                read(2,'(/)')
                do k = 1, Ninc
                  read(2, '(2es15.6)', iostat = istat) en(k), xs
                  if (istat /= 0) call read_error(resfile, istat, eor = 'continue', eof = 'continue')
                  xsnat(k) = xsnat(k) + abun(i) * xs
                enddo
                exit
              endif
            enddo
            close (unit = 2)
          endif
        enddo
        if (resexist) then
          Ncol=2
          massstring='   '
          write(massstring,'(i3)') ia
          finalnuclide=trim(nuc(iz))//trim(adjustl(massstring))
          topline=trim(targetnuclide)//trim(reaction)//trim(finalnuclide)//' '//trim(quantity)
          open (1, file = resfile(1:12), status = 'replace')
          call write_header(indent,topline,source,user,date,oformat)
          call write_target(indent)
          call write_reaction(indent,reaction,0.D0,0.D0,6,5)
          call write_residual(id2,iz,ia,finalnuclide)
          call write_quantity(id2,quantity)
          call write_datablock(id2,Ncol,Ninc,col,un)
          do k = 1, Ninc
            write(1, '(2es15.6)') en(k), xsnat(k)
          enddo
          close (unit = 1)
        endif
!
! Per ground state and isomer
!
        resline = ''
        Nk = 0
        do nex = 0, numlev
          do k = 1, Ninc
            xsnat(k) = 0.
            natbranch(k) = 0.
          enddo
          resexist = .false.
          absum = 0.
          do i = 1, isonum
            resfile = 'rp000000.L00'//natstring(i)
            write(resfile(3:8), '(2i3.3)') iz, ia
            write(resfile(11:12), '(i2.2)') nex
            inquire (file = resfile, exist = lexist)
            if (lexist) then
              resexist = .true.
              open (2, file = resfile, status = 'old', iostat = istat)
              if (istat /= 0) call read_error(resfile, istat)
              do 
                read(2,'(a)',iostat = istat) line
                if (istat == -1) exit
                key='residual:'
                keyix=index(line,trim(key))
                if (keyix > 0) then
                  k = 1
                  resline (k) = line
                  do 
                    read(2,'(a)',iostat = istat) line
                    if (istat == -1) exit
                    if (line(1:4) == '#   ') then
                      k = k + 1
                      resline (k) = line
                    else
                      Nk = k
                      exit
                    endif
                  enddo
                endif
                key='entries'
                keyix=index(line,trim(key))
                if (keyix > 0) then
                  read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Ninc
                  if (istat /= 0) call read_error(resfile, istat, eor = 'continue', eof = 'continue')
                  read(2,'(/)')
                  do k = 1, Ninc
                    read(2, '(3es15.6)', iostat = istat) en(k), xs, br
                    if (istat /= 0) call read_error(resfile, istat, eor = 'continue', eof = 'continue')
                    xsnat(k) = xsnat(k) + abun(i) * xs
                    natbranch(k) = natbranch(k) + abun(i) * br
                  enddo
                  exit
                endif
              enddo
              close (unit = 2)
              absum = absum + abun(i)
            endif
          enddo
          if (absum.gt.0.) then
            do k=1,Ninc
              natbranch(k)=natbranch(k)/absum
            enddo
          endif
          if (resexist) then
            Ncol=3
            topline=trim(targetnuclide)//trim(reaction)//trim(finalnuclide)//' '//trim(quantity)
            open (1, file = resfile(1:12), status = 'replace')
            call write_header(indent,topline,source,user,date,oformat)
            call write_target(indent)
            call write_reaction(indent,reaction,0.D0,0.D0,6,5)
            do k = 1, Nk
              write(1,'(a)') trim(resline(k))
            enddo
            call write_quantity(id2,quantity)
            call write_datablock(id2,Ncol,Ninc,col,un)
            do k = 1, Ninc
              write(1, '(3es15.6)') en(k), xsnat(k), natbranch(k)
            enddo
            close (unit = 1)
          endif
        enddo
    enddo
  enddo
  endif
!
! 6. Exclusive channel cross sections
!
  quantity='cross section'
  col(1)='E'
  col(2)='xs'
  Ncol=2
  MF = 3
  if (filechannels) then
    do npart = 0, maxchannel
      do ia = 0, numia
        do ih = 0, numih
          do it = 0, numit
            do id = 0, numid
              do ip = 0, numip
                do in = 0, numin
                  if (in + ip + id + it + ih + ia /= npart) cycle
                  xsnat = 0.
                  xsexist = .false.
                  do i = 1, isonum
                    xsfile = 'xs000000.tot'//natstring(i)
                    write(xsfile(3:8), '(6i1)') in, ip, id, it, ih, ia
                    inquire (file = xsfile, exist = lexist)
                    if (lexist) then
                      xsexist = .true.
                      open (2, status = 'old', file = xsfile, iostat = istat)
                      if (istat /= 0) call read_error(xsfile, istat)
                      do 
                        read(2,'(a)',iostat = istat) line
                        if (istat == -1) exit
                        key='type'
                        keyix=index(line,trim(key))
                        if (keyix > 0) read(line(keyix+len_trim(key)+2:80),'(a)', iostat = istat) reaction
                        do k= 1, 18
                          if (reaction(k:k) == '"') reaction(k:k) = ' '
                        enddo
                        reaction = trim(adjustl(reaction))
                        key='entries'
                        keyix=index(line,trim(key))
                        if (keyix > 0) then
                          read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Ninc
                          if (istat /= 0) call read_error(xsfile, istat, eor = 'continue', eof = 'continue')
                          read(2,'(/)')
                          do k = 1, Ninc
                            read(2, '(2es15.6)', iostat = istat) en(k), xs
                            if (istat /= 0) call read_error(xsfile, istat, eor = 'continue', eof = 'continue')
                            xsnat(k) = xsnat(k) + abun(i) * xs
                          enddo
                          exit
                        endif
                      enddo
                      close (unit = 2)
                    endif
                  enddo
                  if (xsexist) then
                    topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
                    MF = 0
                    ident = 100000 * in + 10000 * ip + 1000 * id + 100 * it + 10 * ih + ia
                    do i = 1,nummt
                      if (MTchannel(i) == ident) then
                        MT = i
                        exit
                      endif
                    enddo
                    open (1, status = 'replace', file = xsfile(1:12))
                    call write_header(indent,topline,source,user,date,oformat)
                    call write_target(indent)
                    call write_reaction(indent,reaction,0.D0,0.D0,MF,MT)
                    call write_quantity(id2,quantity)
                    call write_datablock(id2,Ncol,Ninc,col,un)
                    do k = 1, Ninc
                      write(1, '(2es15.6)') en(k), xsnat(k)
                    enddo
                    close (unit = 1)
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  endif
!
! 7. Fission cross sections
!
  if (filefission) then
    xsnat = 0.
    fissionexist = .false.
    do i = 1, isonum
      fisfile = 'fission.tot'//natstring(i)
      inquire (file = fisfile, exist = lexist)
      if (lexist) then
        fissionexist = .true.
        open (2, status = 'old', file = fisfile, iostat = istat)
        if (istat /= 0) call read_error(fisfile, istat)
        do 
          read(2,'(a)',iostat = istat) line
          if (istat == -1) exit
          key='entries'
          keyix=index(line,trim(key))
          if (keyix > 0) then
            read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Ninc
            if (istat /= 0) call read_error(fisfile, istat, eor = 'continue', eof = 'continue')
            read(2,'(/)')
            do k = 1, Ninc
              read(2, '(2es15.6)', iostat = istat) en(k), xs
              if (istat /= 0) call read_error(fisfile, istat, eor = 'continue', eof = 'continue')
              xsnat(k) = xsnat(k) + abun(i) * xs
            enddo
            exit
          endif
        enddo
        close (unit = 2)
      endif
    enddo
    if (fissionexist) then
      reaction='('//parsym(k0)//',f)'
      topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
      open (1, status = 'unknown', file = 'fission.tot')
      call write_header(indent,topline,source,user,date,oformat)
      call write_target(indent)
      call write_reaction(indent,reaction,0.D0,0.D0,3,18)
      call write_quantity(id2,quantity)
      call write_datablock(id2,Ncol,Ninc,col,un)
      do k = 1, Ninc
        write(1, '(2es15.6)') en(k), xsnat(k)
      enddo
      close (unit = 1)
    endif
  endif
!
! 8. Fission yields
!
! Excitation function per fission product
!
  if (flagmassdis) then
    do iz = 1, Ztarget
      do ia = 1, Atarget
        xs1nat = 0.
        xs2nat = 0.
        resexist = .false.
        do i = 1, isonum
          resfile = 'fp000000.tot'//natstring(i)
          write(resfile(3:8), '(2i3.3)') iz, ia
          inquire (file = resfile, exist = lexist)
          if (lexist) then
            resexist = .true.
            open (2, status = 'old', file = resfile, iostat = istat)
            if (istat /= 0) call read_error(resfile, istat)
            do 
              read(2,'(a)',iostat = istat) line
              if (istat == -1) exit
              key='entries'
              keyix=index(line,trim(key))
              if (keyix > 0) then
                read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Ninc
                if (istat /= 0) call read_error(resfile, istat, eor = 'continue', eof = 'continue')
                read(2,'(/)')
                do k = 1, Ninc
                  read(2, '(3es15.6)', iostat = istat) en(k), xs1, xs2
                  if (istat /= 0) call read_error(resfile, istat, eor = 'continue', eof = 'continue')
                  xs1nat(k) = xs1nat(k) + abun(i) * xs1
                  xs2nat(k) = xs2nat(k) + abun(i) * xs2
                enddo
                exit
              endif
            enddo
            close (unit = 2)
          endif
        enddo
        if (resexist) then
          open (1, status = 'unknown', file = resfile(1:12))
          write(3, '("# ", a1, " + nat-", a2, ": ff yield of ", i3, a2)') parsym(k0), Starget, ia, nuc(iz)
          write(3, '("# ")')
          write(3, '("# ")')
          write(3, '("# # energies =", i6)') Ninc
          write(3, '("#    E         xs")')
          if (flagffevap) then
            write(1, '("# E-incident   FF Yield   FP yield")')
            do nen = 1, Ninc
              write(1, '(3es15.6)') eninc(nen), xs1nat(nen), xs2nat(nen)
            enddo
          else
            write(1, '("# E-incident   FF Yield")')
            do nen = 1, Ninc
              write(1, '(2es15.6)') eninc(nen), xs1nat(nen)
            enddo
          endif
          close (unit = 1)
        endif
    enddo
  enddo
!
! Mass distribution per incident energy
!
    do k = 1, Ninc
      xs1nat = 0.
      xs2nat = 0.
      fyexist = .false.
      do i = 1, isonum
        fyfile = 'yield0000.000.fis'//natstring(i)
        write(fyfile(6:13), '(f8.3)') eninc(k)
        write(fyfile(6:9), '(i4.4)') int(eninc(k))
        inquire (file = fyfile, exist = lexist)
        if (lexist) then
          fyexist = .true.
          open (2, status = 'old', file = fyfile, iostat = istat)
          if (istat /= 0) call read_error(fyfile, istat)
          read(2, '(////)', iostat = istat)
          if (istat /= 0) call read_error(fyfile, istat, eor = 'continue', eof = 'continue')
          do ia = 1, isotope(i)
            read(2, '(3x, 2e15.4)', iostat = istat) xs1, xs2
            if (istat /= 0) call read_error(fyfile, istat, eor = 'continue', eof = 'continue')
            xs1nat(ia) = xs1nat(ia) + abun(i) * xs1
            xs2nat(ia) = xs2nat(ia) + abun(i) * xs2
          enddo
          close (unit = 2)
        endif
      enddo
      if (fyexist) then
        open (1, status = 'unknown', file = fyfile(1:16))
        write(3, '("# ", a1, " +  nat-", a2, ": mass yields")') parsym(k0), Starget
        write(3, '("# E-incident = ", f8.3)') eninc(k)
        write(3, '("# ")')
        write(3, '("# ")')
        write(3, '("# Mass    Yield   Corrected yield")')
        do ia = 1, Atarget
          write(1, '(i3, 3x, es12.4, 3x, es12.4)') ia, xs1nat(ia), xs2nat(ia)
        enddo
        close (unit = 1)
      endif
    enddo
  endif
!
! 9. Recoil spectra
!
  if (flagrecoil) then
    col(1)='E-out'
    col(2)='xs'
    un(2)='mb/MeV'
    Ncol=2
    quantity='recoil spectrum'
    reaction='('//parsym(k0)//',x)'
    zbeg = max(Ztarget - numZ - 2, 1)
    zend = Ztarget + 2
    abeg = max(isotope(1) - numN - 2, 1)
    aend = isotope(isonum) + 4
    do iz = zbeg, zend
      do ia = abeg, aend
        do k = 1, Ninc
          Estr=''
          write(Estr,'(es12.6)') eninc(k)
          specexist = .false.
          do k2 = 1, 7*numen2
            enspecnat(k2) = 0.
            xsspecnat(k2, 1) = 0.
          enddo
          do i = 1, isonum
            do k2 = 1, numen2
              enspec(i, k2) = 0.
              xsspec(i, k2, 1) = 0.
            enddo
!
! In general, the recoil spectra for the various isotopes are all different.
! Therefore, we first read the recoil energy grids and cross sections into memory.
!
            neniso(i) = 0
            recfile = 'rec000000spec0000.000.tot'//natstring(i)
            write(recfile(4:9), '(2i3.3)') iz, ia
            write(recfile(14:21), '(f8.3)') eninc(k)
            write(recfile(14:17), '(i4.4)') int(eninc(k))
            inquire (file = recfile, exist = lexist)
            if (lexist) then
              specexist = .true.
              open (2, status = 'old', file = recfile)
              if (istat /= 0) call read_error(recfile, istat, eor = 'continue', eof = 'continue')
              do 
                read(2,'(a)',iostat = istat) line
                if (istat == -1) exit
                key='entries'
                keyix=index(line,trim(key))
                if (keyix > 0) then
                  read(2,'(/)')
                  do k2 = 1, numen2
                    read(2, '(f8.3, es12.5)', iostat = istat) enspec(i, k2), xsspec(i, k2, 1)
                    if (istat /= 0) call read_error(recfile, istat, eor = 'continue', eof = 'continue')
                    neniso(i) = neniso(i) + 1
                  enddo
                endif
              enddo
              close (unit = 2)
            endif
          enddo
!
! Make one unifying energy grid by removing double energies and sorting the remaining energies.
!
          if (specexist) then
            nenen = 0
            do i = 1, isonum
Loop2:        do k2 = 1, neniso(i)
                E = enspec(i, k2)
                do nen = 1, nenen
                  if (E == enspecnat(nen)) cycle Loop2
                enddo
                nenen = nenen + 1
                enspecnat(nenen) = E
              enddo Loop2
            enddo
            do n1 = 1, nenen
              do n2 = n1, nenen
                if (enspecnat(n1) <= enspecnat(n2)) cycle
                entmp = enspecnat(n1)
                enspecnat(n1) = enspecnat(n2)
                enspecnat(n2) = entmp
              enddo
            enddo
!
! Interpolation and construction of natural spectra.
!
            do nen = 1, nenen
Loop3:        do i = 1, isonum
                do k2 = 1, neniso(i) - 1
                  Ea = enspec(i, k2)
                  Eb = enspec(i, k2 + 1)
                  if (enspecnat(nen) >= Ea .and. enspecnat(nen) < Eb) then
                    Efac = (enspecnat(nen) - Ea) / (Eb - Ea)
                    xs = xsspec(i, k2, 1) + Efac * (xsspec(i, k2 + 1, 1) - xsspec(i, k2, 1))
                    xsspecnat(nen, 1) = xsspecnat(nen, 1) + abun(i) * xs
                    cycle Loop3
                  endif
                enddo
              enddo Loop3
            enddo
            open (1, status = 'unknown', file = recfile(1:24))
            topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
            call write_header(indent,topline,source,user,date,oformat)
            call write_target(indent)
            call write_reaction(indent,reaction,0.D0,0.D0,6,5)
            call write_real(id2,'E-incident [MeV]',eninc(k))
            call write_quantity(id2,quantity)
            call write_datablock(id2,Ncol,nenen,col,un)
            do nen = 1, nenen
              write(1, '(2es15.6)') enspecnat(nen), xsspecnat(nen, 1)
            enddo
            close (unit = 1)
          endif
        enddo
      enddo
    enddo
  endif
!
! 10. Isotope production
!
! Merge output files in results for natural elements
!
  if (flagprod) then
    zbeg = max(Ztarget - numZ - 2, 1)
    zend = Ztarget + 2
    abeg = max(isotope(1) - numN - 2, 1)
    aend = isotope(isonum) + 4
    do iz = zbeg, zend
      do ia = abeg, aend
        do nex = - 1, min(numlev, 99)
          do k = 0, numtime
            actnat(k) = 0.
            Nisnat(k) = 0.
            Ynat(k) = 0.
            fracnat(k) = 0.
          enddo
          resexist = .false.
          do i = 1, isonum
            Yfile = 'Y000000.tot'//natstring(i)
            write(Yfile(2:7), '(2i3.3)') iz, ia
            if (nex >= 0) write(Yfile(9:11), '("L", i2.2)') nex
            inquire (file = Yfile, exist = lexist)
            if (lexist) then
              resexist = .true.
              open (2, status = 'old', file = Yfile, iostat = istat)
              if (istat /= 0) call read_error(Yfile, istat, eor = 'continue', eof = 'continue')
              do k = 1, 9
                read(2, '(a)',iostat = istat)
                if (istat /= 0) call read_error(Yfile, istat, eor = 'continue', eof = 'continue')
              enddo
              do k = 1, numtime
                read(2, '(f8.1, 3e15.5, f15.5)', iostat = istat) xs, act, Nis, Y, frac
                if (istat /= 0) call read_error(Yfile, istat, eor = 'continue', eof = 'continue')
                actnat(k) = actnat(k) + abun(i) * act
                Nisnat(k) = Nisnat(k) + abun(i) * Nis
                Ynat(k) = Ynat(k) + abun(i) * Y
                fracnat(k) = fracnat(k) + abun(i) * frac
              enddo
              close (unit = 2)
            endif
          enddo
          if (resexist) then
            open (3, status = 'unknown', file = Yfile(1:11))
            write(3, '("# Reaction: ", a1, " + nat", a2, a)') parsym(k0), nuc(Ztarget), trim(str(1)(22:132))
            do k = 2, 9
              write(3, '(a)') trim(str(k))
            enddo
            do k = 1, numtime
              write(3, '(f8.1, 3es15.6, f15.5)') Tgrid(k), actnat(k), Nisnat(k), Ynat(k), fracnat(k)
            enddo
          endif
        enddo
      enddo
    enddo
  endif
  return
end subroutine natural
! Copyright A.J. Koning 2021
