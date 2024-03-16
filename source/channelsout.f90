subroutine channelsout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of exclusive reaction channels
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
!   sgl               ! single precision kind
! All global variables
!   numia             ! maximum number of alphas in channel description
!   numid             ! maximum number of deuterons in channel description
!   numih             ! maximum number of helions in channel description
!   numin             ! maximum number of neutrons in channel description
!   numip             ! maximum number of protons in channel description
!   numit             ! maximum number of tritons in channel description
!   numlev            ! maximum number of discrete levels
! Variables for numerics
!   maxchannel        ! maximal number of outgoing particles in individual channel description
!   maxenrec          ! number of recoil energies
!   xseps             ! limit for cross sections
! Variables for output
!   filechannels      ! flag for exclusive channel cross sections on separate file
!   filerecoil        ! flag for recoil spectra on separate file
! Variables for basic reaction
!   flagendf          ! flag for information for ENDF - 6 file
!   flagrecoil        ! flag for calculation of recoils
! Variables for input energies
!   eninc             ! incident energy in MeV
!   flaginitpop       ! flag for initial population distribution
!   nin               ! counter for incident energy
!   Ninc              ! number of incident energies
! Variables for main input
!   Atarget           ! mass number of target nucleus
!   k0                ! index of incident particle
! Variables for output
!   flagblock         ! flag to block spectra, angle and gamma files
!   flagcheck         ! flag for output of numerical checks
!   flagcompo         ! flag for output of cross section components
!   filefission       ! flag for fission cross sections on separate file
!   flagspec          ! flag for output of spectra
! Variables for fission
!   flagfission       ! flag for fission
! Variables for incident channel
!   channelsum        ! sum over exclusive channel cross sections
!   xsabs             ! absorption cross section
!   xsgr              ! total smoothed giant resonance cross section
!   xsngnsum          ! sum over total (projectile, gamma - ejectile) cross section
!   xsparticle        ! total particle production cross section
!   xspreeq           ! preeq. cross section per particle type and outgoing energy
! Variables for energy grid
!   ebegin            ! first energy point of energy grid
!   egrid             ! outgoing energy grid
!   Einc              ! incident energy in MeV
! Variables for energies
!   eend              ! last energy point of energy grid
!   eendhigh          ! last energy point for energy grid for any particle
!   eninccm           ! center - of - mass incident energy in MeV
!   Ethrexcl          ! threshold incident energy for exclusive channel
!   idchannel         ! identifier for exclusive channel
!   Qexcl             ! Q - value for exclusive channel
!   reacstring        ! string for exclusive reaction channel
!   Ninclow         ! number of incident energies below Elow
! Variables for existence libraries
!   chanexist         ! flag for existence of exclusive cross section
!   chanfisexist      ! flag for existence of exclusive fission cros
!   chanisoexist      ! flag for existence of exclusive iso
!   chanopen          ! flag to open channel with first non - zero cross s
!   gamchanexist      ! flag for existence of exclusive discrete gam
!   recchanexist      ! flag for existence of exclusive recoil spectra
!   spchanexist       ! flag for existence of exclusive spectra
!   spfischanexist    ! flag for existence of exclusive fission spectra
! Variables for exclusive channels
!   Eavchannel        ! channel average energy
!   Especsum          ! total emission energy
!   exclbranch        ! exclusive channel yield per isomer
!   fisstring         ! string for exclusive fission reaction channel
!   gmult             ! continuum gamma multiplicity
!   idnum             ! counter for exclusive channel
!   xschancheck       ! integrated channel spectra
!   xschaniso         ! channel cross section per isomer
!   xschannel         ! channel cross section
!   xschannelsp       ! channel cross section spectra
!   xsfischancheck    ! integrated fission channel spectra
!   xsfischannel      ! fission channel cross section
!   xsfischannelsp    ! fission channel cross section spectra
!   xsgamchannel      ! gamma channel cross section
!   xsgamdischan      ! discrete gamma channel cross section
!   xsparcheck        ! total particle production cross section
!   xsratio           ! ratio of exclusive cross section over residual p
!   xsspeccheck       ! total particle production spectra
!   yieldchannel      ! relative yield
! Variables for multiple emission
!   Fcomp             ! compound population fraction per nucleus
!   Fdir              ! direct population fraction per nucleus
!   Fpreeq            ! preequilibrium population fraction per nucleus
!   xsinitpop         ! initial population cross section
!   xsmpreeq          ! multiple pre - equilibrium emission spectrum
!   xsngnspec         ! total (projectile, gamma - ejectile) spectrum
! Variables for binary reactions
!   xsdisctot         ! total cross section summed over discrete states
!   xsnonel           ! non - elastic cross
! Variables for binary emission spectra
!   xscomp            ! compound elastic cross section
! Variables for nuclides
!   parskip           ! logical to skip outgoing particle
! Constants
!   isochar           ! symbol of isomer
!   iso               ! counter for isotope
!   natstring         ! string extension for file names
!   parname           ! name of particle
!   parsym            ! symbol of particle
! Variables for levels
!   edis              ! energy of level
!   Liso              ! isomeric number of target
!   tau               ! lifetime of state in seconds
! Variables for level density
!   Nlast             ! last discrete level
! Variables for thermal cross sections
!   fexclbranch       ! exclusive channel yield per isomer
!   fxschaniso        ! channel cross section per isomer
!   fxschannel        ! channel cross section
!   fxsgamchannel     ! gamma channel cross section
!   fxsgamdischan     ! discrete gamma channel cross section
!   fxsratio          ! ratio of exclusive cross section over residual pr
! Variables for recoil
!   Erec              ! recoil energy
!   specrecoil        ! recoil spectrum
!
! *** Declaration of local data
!
  implicit none
  character(len=3)  :: isostring ! string to designate target isomer
  character(len=3)  :: massstring ! 
  character(len=6)  :: finalnuclide !
  character(len=8)  :: spstring  ! string
  character(len=13)  :: Estr
  character(len=12) :: isofile   ! file with isomeric cross section
  character(len=12) :: gamfile   ! giant resonance parameter file
  character(len=16) :: xsfile    ! file with channel cross sections
  character(len=21) :: spfile    ! file with spectra
  character(len=18) :: reaction   ! reaction
  character(len=15) :: col(8)    ! header
  character(len=15) :: un(8)    ! units
  character(len=80) :: quantity   ! quantity
  character(len=80) :: gamquant   ! quantity
  character(len=132) :: topline    ! topline
  integer           :: MF
  integer           :: MT
  integer           :: i
  integer           :: Z         ! 
  integer           :: A         ! 
  integer           :: kiso      ! 
  integer           :: i1        ! value
  integer           :: i2        ! value
  integer           :: ia        ! mass number from abundance table
  integer           :: id        ! counter for deuterons
  integer           :: idc       ! help variable
  integer           :: ident     ! exclusive channel identifier
  integer           :: ih        ! hole number
  integer           :: in        ! counter for neutrons
  integer           :: ip        ! particle number
  integer           :: it        ! counter for tritons
  integer           :: Ncomp     ! neutron number index for compound nucleus
  integer           :: nen       ! energy counter
  integer           :: nex0      ! excitation energy bin of compound nucleus
  integer           :: nex       ! excitation energy bin of compound nucleus
  integer           :: Ngam      ! counter for discrete gamma rays
  integer           :: NL        ! last discrete level
  integer           :: npart     ! number of particles in outgoing channel
  integer           :: type      ! particle type
  integer           :: Zcomp     ! proton number index for compound nucleus
  integer           :: Ncol      ! counter
  real(sgl)         :: emissum   ! integrated binary emission spectrum
  real(sgl)         :: xs        ! help variable
  real(sgl)         :: Egam      ! help variable
!
! ****************** Output of channel cross sections ******************
!
  Estr=''
  write(Estr,'(es13.6)') Einc
  un = 'mb'
  col(1)='E'
  un(1)='MeV'
  col(2)='xs'
  col(3)='gamma_xs'
  col(4)='xs/res.prod.xs'
  un(4) = ''
  col(5)='Direct'
  col(6)='Preequilibrium'
  col(7)='Compound'
  quantity='cross section'
  gamquant='gamma cross section and multiplicity'
  write(*, '(/" 6. Exclusive cross sections"/)')
  write(*, '(" 6a. Total exclusive cross sections "/)')
  write(*, '("     Emitted particles     cross section reaction", "         level    isomeric    isomeric    lifetime", &
 &  " relative yield")')
  write(*, '("    n   p   d   t   h   a", 40x, "cross section", "   ratio")')
  do npart = 0, maxchannel
    do ia = 0, numia
      do ih = 0, numih
        do it = 0, numit
          do id = 0, numid
            do ip = 0, numip
Loop1:        do in = 0, numin
                if (in + ip + id + it + ih + ia /= npart) cycle
                if ( .not. chanopen(in, ip, id, it, ih, ia)) cycle
                ident = 100000 * in + 10000 * ip + 1000 * id + 100 * it + 10 * ih + ia
                do idc = 0, idnum
                  if (idchannel(idc) == ident) then
                    if (xschannel(idc) < xseps) cycle Loop1
                    write(*, '(1x, 6i4, 3x, es12.5, 2x, a17, 40x, es12.5)') in, ip, id, it, ih, ia, xschannel(idc), &
 &                     reacstring(idc), yieldchannel(idc)
                    Zcomp = ip + id + it + 2 * ih + 2 * ia
                    Ncomp = in + id + 2 * it + ih + 2 * ia
                    NL = Nlast(Zcomp, Ncomp, 0)
                    do nex0 = 1, NL
                      if (tau(Zcomp, Ncomp, nex0) /= 0.) then
                        write(*, '(61x, "0    ", es12.5, f9.5)') xschaniso(idc, 0), exclbranch(idc, 0)
                        do nex = 1, NL
                          if (tau(Zcomp, Ncomp, nex) /= 0.) then
                            write(*, '(59x, i3, 4x, es12.5, f9.5, 2x, es12.5, " sec. ")')  levnum(Zcomp, Ncomp, nex), & 
 &                          xschaniso(idc, nex), exclbranch(idc, nex), tau(Zcomp, Ncomp, nex)
                          endif
                        enddo
                        cycle Loop1
                      endif
                    enddo
                  endif
                enddo
              enddo Loop1
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
!
! Write results on separate files
!
  isostring = '   '
  if (filechannels) then
    if (Liso >= 1) isostring = '('//isochar(Liso)//')'
    do npart = 0, maxchannel
    do ia = 0, numia
    do ih = 0, numih
    do it = 0, numit
    do id = 0, numid
    do ip = 0, numip
Loop2:    do in = 0, numin
      if (in + ip + id + it + ih + ia /= npart) cycle
      if ( .not. chanopen(in, ip, id, it, ih, ia)) cycle
      ident = 100000 * in + 10000 * ip + 1000 * id + 100 * it + 10 * ih + ia
      do idc = 0, idnum
        if (idchannel(idc) == ident) then
          if (xschannel(idc) < xseps .and. npart /= 0 .and. .not. chanexist(in, ip, id, it, ih, ia)) cycle Loop2
          Zcomp = ip + id + it + 2 * ih + 2 * ia
          Ncomp = in + id + 2 * it + ih + 2 * ia
          Z = ZZ(Zcomp, Ncomp, 0)
          A = AA(Zcomp, Ncomp, 0)
          massstring='   '
          write(massstring,'(i3)') A
          finalnuclide=trim(nuc(Z))//adjustl(massstring)
          reaction=reacstring(idc)
          un = 'mb'
          col(1)='E'
          un(1)='MeV'
          col(2)='xs'
          col(3)='gamma_xs'
          col(4)='xs/res.prod.xs'
          un(4) = ''
          MF = 3
          MT = 0
          do i = 1,nummt
            if (MTchannel(i) == ident) then
              MT = i
              exit
            endif
          enddo
!
! A. Total
!
          xsfile = 'xs000000.tot'//natstring(iso)
          write(xsfile(3:8), '(6i1)') in, ip, id, it, ih, ia
          if ( .not. chanexist(in, ip, id, it, ih, ia)) then
            chanexist(in, ip, id, it, ih, ia) = .true.
            open (unit = 1, file = xsfile, status = 'unknown')
            topline=trim(targetnuclide)//trim(reaction)//trim(finalnuclide)//' '//trim(quantity)
            if (flagcompo) then
              Ncol=7
            else
              Ncol=4
            endif
            call write_header(topline,source,user,date,oformat)
            call write_target
            call write_reaction(reaction,Qexcl(idc, 0),Ethrexcl(idc, 0),MF,MT)
            call write_residual(Z,A,finalnuclide)
            call write_datablock(quantity,Ncol,Ninc,col,un)
            if (flagcompo) then
              do nen = 1, Ninclow
                write(1, '(7es15.6)') eninc(nen), fxschannel(nen, idc), fxsgamchannel(nen, idc), &
 &                fxsratio(nen, idc), Fdir(Zcomp, Ncomp) * fxschannel(nen, idc), &
 &                Fpreeq(Zcomp, Ncomp) * fxschannel(nen, idc), Fcomp(Zcomp, Ncomp) * fxschannel(nen, idc)
                  fxschannel(nen, idc) = 0.
                  fxsgamchannel(nen, idc) = 0.
                  fxsratio(nen, idc) = 0.
              enddo
              do nen = Ninclow + 1, nin - 1
                write(1, '(7es15.6)') eninc(nen), 0., 0., 0., 0., 0., 0.
              enddo
            else
              do nen = 1, Ninclow
                write(1, '(4es15.6)') eninc(nen), fxschannel(nen, idc), fxsgamchannel(nen, idc), fxsratio(nen, idc)
                  fxschannel(nen, idc) = 0.
                  fxsgamchannel(nen, idc) = 0.
                  fxsratio(nen, idc) = 0.
              enddo
              do nen = Ninclow + 1, nin - 1
                write(1, '(4es15.6)') eninc(nen), 0., 0., 0.
              enddo
            endif
          else
            open (unit = 1, file = xsfile, status = 'old', position = 'append')
          endif
          if (flagcompo) then
            write(1, '(7es15.6)') Einc, xschannel(idc), xsgamchannel(idc), xsratio(idc), &
 &            Fdir(Zcomp, Ncomp) * xschannel(idc), Fpreeq(Zcomp, Ncomp) * xschannel(idc), Fcomp(Zcomp, Ncomp) * xschannel(idc)
          else
            write(1, '(4es15.6)') Einc, xschannel(idc), xsgamchannel(idc), xsratio(idc)
          endif
          close (unit = 1)
!
! B. Ground state and isomers
!
          NL = Nlast(Zcomp, Ncomp, 0)
          do nex0 = 1, NL
            if (tau(Zcomp, Ncomp, nex0) /= 0.) then
              do nex = 0, NL
                if (nex == 0 .or. tau(Zcomp, Ncomp, nex) /= 0.) then
                  isofile = 'xs000000.L00'
                  write(isofile(3:8), '(6i1)') in, ip, id, it, ih, ia
                  write(isofile(11:12), '(i2.2)') levnum(Zcomp, Ncomp, nex)
                  if ( .not. chanisoexist(in, ip, id, it, ih, ia, nex)) then
                    chanisoexist(in, ip, id, it, ih, ia, nex) = .true.
                    open (unit = 1, file = isofile, status = 'unknown')
                    if (nex == 0) then
                      kiso = 0
                    else
                      kiso = kiso + 1
                    endif
                    finalnuclide=trim(nuc(Z))//trim(adjustl(massstring))//isochar(kiso)
                    topline=trim(targetnuclide)//trim(reaction)//trim(finalnuclide)//' '//trim(quantity)
                    col(3)='Isomeric_ratio'
                    un(3)=''
                    Ncol=3
                    MF = 10
                    call write_header(topline,source,user,date,oformat)
                    call write_target
                    call write_reaction(reaction,Qexcl(idc, nex),Ethrexcl(idc, nex),MF,MT)
                    call write_residual(Z,A,finalnuclide)
                    call write_level(2,kiso,levnum(Zcomp, Ncomp, nex),edis(Zcomp, Ncomp, nex), &
 &                    jdis(Zcomp, Ncomp, nex),parlev(Zcomp, Ncomp, nex),tau(Zcomp, Ncomp, nex))
                    call write_datablock(quantity,Ncol,Ninc,col,un)
                    do nen = 1, Ninclow
                      write(1, '(3es15.6)') eninc(nen), fxschaniso(nen, idc, nex), fexclbranch(nen, idc, nex)
                        fxschaniso(nen, idc, nex) = 0.
                        fexclbranch(nen, idc, nex) = 0.
                    enddo
                    do nen = Ninclow + 1, nin - 1
                      write(1, '(3es15.6)') eninc(nen), 0., 0.
                    enddo
                  else
                    open (unit = 1, file = isofile, status = 'old', position = 'append')
                  endif
                  write(1, '(3es15.6)') Einc, xschaniso(idc, nex), exclbranch(idc, nex)
                  close (unit = 1)
                endif
              enddo
              exit
            endif
          enddo
!
! C. Discrete gamma-rays
!
          if (flagendf) then
            gamfile = 'xs000000.gam'
            write(gamfile(3:8), '(6i1)') in, ip, id, it, ih, ia
            topline=trim(targetnuclide)//trim(reaction)//trim(finalnuclide)//' '//trim(gamquant)
            un = ''
            col(1)='Parent_level'
            col(2)='Daughter_level'
            col(3)='xs'
            un(3)='mb'
            col(4)='Parent_energy'
            un(4)='MeV'
            col(5)='Daughter_energy'
            un(5)='MeV'
            col(6)='Gamma_energy'
            un(6)='MeV'
            Ncol=6
            MF = 13
            if ( .not. gamchanexist(in, ip, id, it, ih, ia)) then
              gamchanexist(in, ip, id, it, ih, ia) = .true.
              open (unit = 1, file = gamfile, status = 'unknown')
              do nen = 1, Ninclow
                Ngam = 0
                do i1 = 1, numlev
                  do i2 = 0, i1
                    if (fxsgamdischan(nen, idc, i1, i2) > 0.) Ngam = Ngam + 1
                  enddo
                enddo
                call write_header(topline,source,user,date,oformat)
                call write_target
                call write_reaction(reaction,Qexcl(idc, 0),Ethrexcl(idc, 0),MF,MT)
                call write_real(2,'E-incident [MeV]',eninc(nen))
                call write_residual(Z,A,finalnuclide)
                call write_datablock(gamquant,Ncol,Ngam,col,un)
                do i1 = 1, numlev
                  do i2 = 0, i1
                    if (fxsgamdischan(nen, idc, i1, i2) > 0.) then
                      Egam = edis(Zcomp, Ncomp, i1) - edis(Zcomp, Ncomp, i2)
                      write(1, '(2(i6, 9x), 4es15.6)') i1, i2, fxsgamdischan(nen, idc, i1, i2), &
 &                      edis(Zcomp, Ncomp, i1), edis(Zcomp, Ncomp, i2), Egam
                    endif
                  enddo
                enddo
              enddo
              do nen = Ninclow + 1, nin - 1
                Ngam = 0
                call write_header(topline,source,user,date,oformat)
                call write_target
                call write_reaction(reaction,Qexcl(idc, 0),Ethrexcl(idc, 0),MF,MT)
                call write_real(2,'E-incident [MeV]',eninc(nen))
                call write_residual(Z,A,finalnuclide)
                call write_datablock(gamquant,Ncol,Ngam,col,un)
              enddo
            else
              open (unit = 1, file = gamfile, status = 'old', position = 'append')
            endif
            Ngam = 0
            do i1 = 1, numlev
              do i2 = 0, i1
                if (xsgamdischan(idc, i1, i2) > 0.) Ngam = Ngam + 1
              enddo
            enddo
            call write_header(topline,source,user,date,oformat)
            call write_target
            call write_reaction(reaction,Qexcl(idc, 0),Ethrexcl(idc, 0),MF,MT)
            call write_real(2,'E-incident [MeV]',Einc)
            call write_residual(Z,A,finalnuclide)
            call write_datablock(gamquant,Ncol,Ngam,col,un)
            do i1 = 1, numlev
              do i2 = 0, i1
                if (xsgamdischan(idc, i1, i2) > 0.) then
                  Egam = edis(Zcomp, Ncomp, i1) - edis(Zcomp, Ncomp, i2)
                  write(1, '(2(i6, 9x), 4es15.6)') i1, i2, xsgamdischan(idc, i1, i2), edis(Zcomp, Ncomp, i1), &
 &                  edis(Zcomp, Ncomp, i2), Egam
                endif
              enddo
            enddo
            close (unit = 1)
          endif
        endif
      enddo
    enddo Loop2
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
  endif
!
! *************** Output of fission channel cross sections *************
!
  if (flagfission) then
    write(*, '(/" 6a2. Exclusive fission cross sections "/)')
    write(*, '("     Emitted particles     cross section reaction")')
    write(*, '("    n   p   d   t   h   a")')
    do npart = 0, maxchannel
      do ia = 0, numia
        do ih = 0, numih
          do it = 0, numit
            do id = 0, numid
              do ip = 0, numip
                do in = 0, numin
                  if (in + ip + id + it + ih + ia /= npart) cycle
                  if ( .not. chanopen(in, ip, id, it, ih, ia)) cycle
                  if (nin == Ninclow + 1) chanfisexist(in, ip, id, it, ih, ia) = .false.
                  ident = 100000 * in + 10000 * ip + 1000 * id + 100 * it + 10 * ih + ia
                  do idc = 0, idnum
                    if (idchannel(idc) == ident) goto 430
                  enddo
                  cycle
  430             if (xsfischannel(idc) < xseps) cycle
                  write(*, '(1x, 6i4, 3x, es12.5, 2x, a17)') in, ip, id, it, ih, ia, xsfischannel(idc), fisstring(idc)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  endif
  write(*, '(/" Absorption cross section                 :", f12.5)') xsabs
  write(*, '(/" Sum over exclusive channel cross sections:", f12.5)') channelsum
  write(*, '(" (n,gn) + (n,gp) +...(n,ga) cross sections:", f12.5)') xsngnsum
  write(*, '(" Total                                    :", f12.5)') channelsum + xsngnsum
  if (flaginitpop) then
    write(*, '(" Initial population cross section         :", f12.5)') xsinitpop
  else
    write(*, '(" Non-elastic cross section                :", f12.5)') xsnonel
  endif
!
! Write results on separate files
!
  MF = 3
  if (filefission) then
    do npart = 0, maxchannel
      do ia = 0, numia
        do ih = 0, numih
          do it = 0, numit
            do id = 0, numid
              do ip = 0, numip
Loop3:          do in = 0, numin
                  if (in + ip + id + it + ih + ia /= npart) cycle
                  if ( .not. chanopen(in, ip, id, it, ih, ia)) cycle
                  ident = 100000 * in + 10000 * ip + 1000 * id + 100 * it + 10 * ih + ia
                  do idc = 0, idnum
                    if (idchannel(idc) == ident) then
                      if (xsfischannel(idc) < xseps .and. npart /= 0 .and. .not. chanexist(in, ip, id, it, ih, ia)) cycle Loop3
                      Zcomp = ip + id + it + 2 * ih + 2 * ia
                      Ncomp = in + id + 2 * it + ih + 2 * ia
                      xsfile = 'xs000000.fis'
                      write(xsfile(3:8), '(6i1)') in, ip, id, it, ih, ia
                      reaction=fisstring(idc)
                      if ( .not. chanfisexist(in, ip, id, it, ih, ia)) then
                        chanfisexist(in, ip, id, it, ih, ia) = .true.
                        open (unit = 1, file = xsfile, status = 'unknown')
                        topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
                        Ncol=2
                        MT = 0
                        if (ident == 100000) MT = 19
                        if (ident == 200000) MT = 20
                        if (ident == 300000) MT = 21
                        if (ident == 400000) MT = 38
                        call write_header(topline,source,user,date,oformat)
                        call write_target
                        call write_reaction(reaction,Qexcl(idc, 0),Ethrexcl(idc, 0),MF,MT)
                        call write_datablock(quantity,Ncol,Ninc,col,un)
                        do nen = 1, nin - 1
                          write(1, '(2es15.6)') eninc(nen), 0.
                        enddo
                      else
                        open (unit = 1, file = xsfile, status = 'old', position = 'append')
                      endif
                      write(1, '(2es15.6)') Einc, xsfischannel(idc)
                      close (unit = 1)
                    endif
                  enddo
                enddo Loop3
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  endif
!
! ********** Check of total particle production cross section **********
!
  if (flagcheck) then
    write(*, '(/" Check of particle production cross sections")')
    write(*, '(" (Only applies if non-elastic cross section is", " exhausted by exclusive cross sections)"/)')
    do type = 1, 6
      if (parskip(type)) cycle
      write(*, '(1x, a8, "=", es12.5, "    Summed exclusive cross sections=", es12.5)') parname(type), xsparticle(type), &
 &      xsparcheck(type)
    enddo
  endif
!
! *************** Output of channel cross section spectra **************
!
  MF = 6
  if (flagspec) then
    col(1) = 'E-out'
    un(1) = 'MeV'
    do type = 0, 6
      col(type+2) = parname(type)
      un(type+2) = 'mb/MeV'
    enddo
    Ncol=8
    write(*, '(/" 6b. Exclusive spectra ")')
    do npart = 0, maxchannel
      do ia = 0, numia
        do ih = 0, numih
          do it = 0, numit
            do id = 0, numid
              do ip = 0, numip
                do in = 0, numin
                  if (in + ip + id + it + ih + ia /= npart) cycle
                  if ( .not. chanopen(in, ip, id, it, ih, ia)) cycle
                  ident = 100000 * in + 10000 * ip + 1000 * id + 100 * it + 10 * ih + ia
                  do idc = 0, idnum
                    if (idchannel(idc) == ident) goto 630
                  enddo
                  cycle
  630             if (xschannel(idc) < xseps) cycle
                  write(*, '(/"      Emitted particles     ", "cross section reaction      gamma cross section")')
                  write(*, '("    n   p   d   t   h   a")')
                  write(*, '(1x, 6i4, 3x, es12.5, 2x, a17, es12.5)') in, ip, id, it, ih, ia, xschannel(idc), reacstring(idc), &
 &                  xsgamchannel(idc)
                  write(*, '(/"  Outgoing spectra"/)')
                  write(*, '("   Energy  ", 7(a8, 4x)/)') (parname(type), type = 0, 6)
                  do nen = ebegin(0), eendhigh
                    write(*, '(1x, f8.3, 7es12.5)') egrid(nen), (xschannelsp(idc, type, nen), type = 0, 6)
                  enddo
                  if (filechannels) then
                    spstring = 'sp000000'
                    write(spstring(3:8), '(6i1)') in, ip, id, it, ih, ia
                    if (flagblock) then
                      spfile = spstring//'.tot'
                      if (.not. spchanexist(in,ip,id,it,ih,ia)) then
                        spchanexist(in,ip,id,it,ih,ia) = .true.
                        open (unit=1, file=spfile, status='unknown')
                      else
                        open (unit=1, file=spfile, status='unknown', position='append')
                      endif
                    else
                      spfile = spstring//'E0000.000.tot'
                      write(spfile(10:17), '(f8.3)') Einc
                      write(spfile(10:13), '(i4.4)') int(Einc)
                      open (unit=1, file=spfile, status='unknown')
                    endif
                    reaction=reacstring(idc)
                    quantity='emission spectrum'
                    topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
                    do i = 1,nummt
                      if (MTchannel(i) == ident) then
                        MT = i
                        exit
                      endif
                    enddo
                    call write_header(topline,source,user,date,oformat)
                    call write_target
                    call write_reaction(reaction,Qexcl(idc, 0),Ethrexcl(idc, 0),MF,MT)
                    call write_real(2,'E-incident [MeV]',Einc)
                    call write_datablock(quantity,Ncol,eendhigh-ebegin(0)+1,col,un)
                    if (npart == 0) then
                      do nen = ebegin(0), eendhigh
                        write(1, '(8es15.6)') egrid(nen), xschannelsp(idc, 0, nen), (xsngnspec(type, nen), type = 1, 6)
                      enddo
                    else
                      do nen = ebegin(0), eendhigh
                        write(1, '(8es15.6)') egrid(nen), (xschannelsp(idc, type, nen), type = 0, 6)
                      enddo
                    endif
                    close (unit = 1)
                  endif
                  if (flagcheck) then
                    emissum = xschancheck(idc)
                    xs = 0.
                    if (npart == 1 .and. in == 1) xs = xsdisctot(1)
                    if (npart == 1 .and. ip == 1) xs = xsdisctot(2)
                    if (npart == 1 .and. id == 1) xs = xsdisctot(3)
                    if (npart == 1 .and. it == 1) xs = xsdisctot(4)
                    if (npart == 1 .and. ih == 1) xs = xsdisctot(5)
                    if (npart == 1 .and. ia == 1) xs = xsdisctot(6)
                    emissum = emissum + xs
                    write(*, '(/"  E-av    ", 7(f7.3, 5x))') (Eavchannel(idc, type), type = 0, 6)
                    write(*, '("  multi   ", f7.3, 6(11x, i1))') gmult(idc), in, ip, id, it, ih, ia
                    write(*, '("  Total   ", 7(f7.3, 5x))') gmult(idc) * Eavchannel(idc, 0), &
 &                    in * Eavchannel(idc, 1), ip * Eavchannel(idc, 2), id * Eavchannel(idc, 3), it * Eavchannel(idc, 4), &
 &                    ih * Eavchannel(idc, 5), ia * Eavchannel(idc, 6)
                    write(*, '(/" Available energy:", f10.5)') Qexcl(idc, 0) + eninccm
                    write(*, '(" Emission energy :", f10.5)') Especsum(idc)
                    write(*, '(/" Check of integrated emission spectra:")')
                    if (npart == 0) then
                      write(*, '(" Cross section (x multiplicity)    =", es12.5)') gmult(idc)*xschannel(idc)
                    else
                      write(*, '(" Cross section (x multiplicity)    =", es12.5)') npart*xschannel(idc)
                    endif
                    if (npart == 1) then
                      write(*, '(" Integrated spectra + discrete c.s.=", es12.5)') emissum
                      write(*, '(" Integrated spectra                  =", es12.5)') xschancheck(idc)
                      write(*, '(" Discrete cross sections             =", es12.5)') xs
                    else
                      write(*, '(" Integrated spectra                =", es12.5)') emissum
                    endif
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
!
! *********** Output of fission channel cross section spectra **********
!
    if (flagfission) then
      write(*, '(/" 6b2. Exclusive fission spectra ")')
      do npart = 0, maxchannel
        do ia = 0, numia
          do ih = 0, numih
            do it = 0, numit
              do id = 0, numid
                do ip = 0, numip
                  do in = 0, numin
                    if (in + ip + id + it + ih + ia /= npart) cycle
                    if ( .not. chanopen(in, ip, id, it, ih, ia)) cycle
                    ident = 100000 * in + 10000 * ip + 1000 * id + 100 * it + 10 * ih + ia
                    do idc = 0, idnum
                      if (idchannel(idc) == ident) goto 730
                    enddo
                    cycle
  730               if (xsfischannel(idc) <= xseps) cycle
                    write(*, '(/"      Emitted particles     ", "cross section reaction")')
                    write(*, '("    n   p   d   t   h   a")')
                    write(*, '(1x, 6i4, 3x, es12.5, 2x, a17)') in, ip, id, it, ih, ia, xsfischannel(idc), fisstring(idc)
                    write(*, '(/"  Outgoing spectra"/)')
                    write(*, '("   Energy  ", 7(a8, 4x)/)') (parname(type), type = 0, 6)
                    do nen = ebegin(0), eend(0)
                      write(*, '(1x, f8.3, 7es12.5)') egrid(nen), (xsfischannelsp(idc, type, nen), type = 0, 6)
                    enddo
                    if (filechannels) then
                      spstring = 'sp000000'
                      write(spstring(3:8), '(6i1)') in, ip, id, it, ih, ia
                      if (flagblock) then
                        spfile = spstring//'.fis'
                        if (.not. spfischanexist(in,ip,id,it,ih,ia)) then
                          spfischanexist(in,ip,id,it,ih,ia) = .true.
                          open (unit=1, file=spfile, status='unknown')
                        else
                          open (unit=1, file=spfile, status='unknown', position='append')
                        endif
                      else
                        spfile = spstring//'E0000.000.fis'
                        write(spfile(10:17), '(f8.3)') Einc
                        write(spfile(10:13), '(i4.4)') int(Einc)
                        open (unit=1, file=spfile, status='unknown')
                      endif
                      reaction=fisstring(idc)
                      quantity='fission emission spectrum'
                      topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
                      call write_header(topline,source,user,date,oformat)
                      call write_target
                      call write_reaction(reaction,Qexcl(idc, 0),Ethrexcl(idc, 0),0,0)
                      call write_real(2,'E-incident [MeV]',Einc)
                      call write_datablock(quantity,Ncol,eendhigh-ebegin(0)+1,col,un)
                      do nen = ebegin(0), eendhigh
                        write(1, '(8es15.6)') egrid(nen), (xsfischannelsp(idc, type, nen), type = 0, 6)
                      enddo
                      close (unit = 1)
                    endif
                    if (flagcheck) then
                      write(*, '(/" Check of integrated emission spectra:")')
                      write(*, '(" Cross section (x multiplicity)=", es12.5)') npart * xsfischannel(idc)
                      write(*, '(" Integrated spectra            =", es12.5)') xsfischancheck(idc)
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
! ************* Check of total particle production spectra *************
!
    if (flagcheck) then
      write(*, '(/" Check of total particle production spectra"/)')
      write(*, '(" (Only applies if non-elastic cross section", " is exhausted by exclusive cross sections)")')
      do type = 0, 6
        if (parskip(type)) cycle
        write(*, '(/" Summed exclusive ", a8, " spectra"/)') parname(type)
        write(*, '("  Energy   Summed     Composite   Difference"/)')
        do nen = ebegin(type), eend(type)
          xs = xspreeq(type, nen) + xsmpreeq(type, nen) + xsgr(type, nen) + xscomp(type, nen)
          write(*, '(1x, f8.3, 2es12.5, 2x, es12.5)') egrid(nen), xsspeccheck(type, nen), xs, xsspeccheck(type, nen) - xs
        enddo
      enddo
    endif
  endif
!
! ***************** Output of channel recoil spectra *******************
!
  MF = 6
  if (flagrecoil) then
    write(*, '(/" 6c. Exclusive recoil spectra ")')
    do npart = 0, maxchannel
      do ia = 0, numia
        do ih = 0, numih
          do it = 0, numit
            do id = 0, numid
              do ip = 0, numip
                do in = 0, numin
                  if (in + ip + id + it + ih + ia /= npart) cycle
                  if ( .not. chanopen(in, ip, id, it, ih, ia)) cycle
                  ident = 100000 * in + 10000 * ip + 1000 * id + 100 * it + 10 * ih + ia
                  do idc = 0, idnum
                    if (idchannel(idc) == ident) goto 930
                  enddo
                  cycle
  930             if (xschannel(idc) < xseps) cycle
                  Zcomp = ip + id + it + 2 * ih + 2 * ia
                  Ncomp = in + id + 2 * it + ih + 2 * ia
                  Z = ZZ(Zcomp, Ncomp, 0)
                  A = AA(Zcomp, Ncomp, 0)
                  massstring='   '
                  write(massstring,'(i3)') A
                  finalnuclide=trim(nuc(Z))//adjustl(massstring)
                  write(*, '(/"      Emitted particles     ", "cross section reaction")')
                  write(*, '("    n   p   d   t   h   a")')
                  write(*, '(1x, 6i4, 3x, es12.5, 2x, a17)') in, ip, id, it, ih, ia, xschannel(idc), reacstring(idc)
                  write(*, '(/" Recoil spectrum"/)')
                  write(*, '("   Energy  Cross section"/)')
                  do nen = 0, maxenrec
                    write(*, '(1x, f8.3, es12.5)') Erec(Zcomp, Ncomp, nen), specrecoil(Zcomp, Ncomp, nen) * xsratio(idc)
                  enddo
                  if (filechannels .and. filerecoil) then
                    spstring = 'sp000000'
                    write(spstring(3:8), '(6i1)') in, ip, id, it, ih, ia
                    if (flagblock) then
                      spfile = spstring//'.rec'
                      if (.not. recchanexist(in,ip,id,it,ih,ia)) then
                        recchanexist(in,ip,id,it,ih,ia) = .true.
                        open (unit=1, file=spfile, status='unknown')
                      else
                        open (unit=1, file=spfile, status='unknown', position='append')
                      endif
                    else
                      spfile = spstring//'E0000.000.rec'
                      write(spfile(10:17), '(f8.3)') Einc
                      write(spfile(10:13), '(i4.4)') int(Einc)
                      open (unit=1, file=spfile, status='unknown')
                    endif
                    reaction=reacstring(idc)
                    quantity='recoil spectrum - exclusive'
                    topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
                    col(1) = 'E-out'
                    un(1) = 'MeV'
                    col(2) = 'xs'
                    un(2) = 'mb/MeV'
                    Ncol=2
                    do i = 1,nummt
                      if (MTchannel(i) == ident) then
                        MT = i
                        exit
                      endif
                    enddo
                    call write_header(topline,source,user,date,oformat)
                    call write_target
                    call write_reaction(reaction,Qexcl(idc, 0),Ethrexcl(idc, 0),MF,MT)
                    call write_real(2,'E-incident [MeV]',Einc)
                    call write_residual(Z,A,finalnuclide)
                    call write_datablock(quantity,Ncol,maxenrec+1,col,un)
                    do nen = 0, maxenrec
                      write(1, '(2es15.6)') Erec(Zcomp, Ncomp, nen), specrecoil(Zcomp, Ncomp, nen) * xsratio(idc)
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
  return
end subroutine channelsout
! Copyright A.J. Koning 2021
