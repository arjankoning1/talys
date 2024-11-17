subroutine astroout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of astrophysical reaction rates
!
! Author    : Stephane Hilaire and Stephane Goriely
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
!   numN              ! maximum number of neutrons from initial compound nucleus
!   numT              ! number of temperatures
!   numZ              ! maximum number of protons from initial compound nucleus
! Variables for astrophysics
!   astroE            ! energy, in MeV, for Maxwellian average
!   flagastroex       ! flag for calculation of astrophysics reaction rate to f
!   flagastrogs       ! flag for calculation of astrophysics reaction rate with
!   nonthermlev       ! non - thermalized level in the calculation of astrophysic
!   nTmax             ! effective number of temperatures for Maxwellian
! Variables for fission
!   flagfission       ! flag for fission
! Variables for main input
!   Atarget           ! mass number of target nucleus
!   k0                ! index of incident particle
!   Ltarget           ! excited level of target
!   Starget           ! symbol of target nucleus
!   Ztarget           ! charge number of target nucleus
! Variables for gamma rays
!   flagracap         ! flag for radiative capture model
! Variables for nuclides
!   AA                ! mass number of residual nucleus
!   T9                ! Temperature grid in 10 **9 K
!   targetP           ! parity of target
!   targetspin        ! spin of target
!   ZZ                ! charge number of residual nucleus
! Constants
!   nuc               ! symbol of nucleus
!   parsym            ! symbol of particle
! Variables for files
!   path              ! directory containing files to be read
! Variables for levels
!   edis              ! energy of level
!   tau               ! lifetime of state in seconds
! Variables for level density
!   Nlast             ! last discrete level
! Variables for astro
!   macsastro         ! Maxwellian - averaged thermonuclear reaction
!   macsastroex       ! thermonuclear reaction cross section to a
!   macsastrofis      ! thermonuclear reaction cross section for f
!   macsastroracap    ! Maxwellian - averaged thermonuc. reac. c.s.
!   maxNastro         ! maximal number of neutrons away from initi
!   maxZastro         ! maximal number of protons away from initia
!   partf             ! integrated partition function
!   rateastro         ! thermonuclear reaction rate
!   rateastroex       ! thermonuclear reaction rate to a given exc
!   rateastrofis      ! thermonuclear reaction rate factor for fis
!   rateastroracap    ! thermonuclear reaction rate for direct cap
!
! *** Declaration of local data
!
  implicit none
! integer, parameter :: numZ1=numZ+1                     ! numZ + 1
! integer, parameter :: numN1=numN+1                     ! numN + 1
  logical            :: lexist                           ! logical to determine existence
  character(len=1)   :: sym                              ! symbol
  character(len=3)   :: massstring
  character(len=6)   :: ZAstring
  character(len=6)   :: finalnuclide !
  character(len=7)   :: machar                           ! part of filename for MACS
  character(len=132) :: mafile                           ! file with MACS
  character(len=132) :: astrofile                        ! file with astro results
  character(len=1)   :: yesno                            ! y or n function
  character(len=18) :: reaction   ! reaction
  character(len=15) :: col(200)    ! header
  character(len=15) :: un(200)    ! units
  character(len=80) :: quantity   ! quantity
  character(len=132) :: topline    ! topline
  integer            :: A                                ! mass number of target nucleus
  integer            :: Acomp                            ! mass number index for compound nucleus
  integer            :: arespro((numZ+1)*(numN+1))       ! A of residual product
  integer            :: i                                ! counter
  integer            :: Ncol
  integer            :: type
  integer            :: j                                ! counter
  integer            :: ia                               ! mass number from abundance table
  integer            :: ires                             ! counter
  integer            :: iresprod                         ! counter
  integer            :: istat                            ! logical for file access
  integer            :: iwriterp                         ! counter
  integer            :: maxAastro                        ! maximal number of nucleons away from initial compound  nucleus for astrop
  integer            :: Ncomp                            ! neutron number index for compound nucleus
  integer            :: nex                              ! excitation energy bin of compound nucleus
  integer            :: Z                                ! charge number of target nucleus
  integer            :: Zix
  integer            :: Nix
  integer            :: Zcomp                            ! proton number index for compound nucleus
  integer            :: zrespro((numZ+1)*(numN+1))       ! Z of residual product
  real(sgl)          :: branch                           ! branching ratio to a given excited state
  real(sgl)          :: dmacs                            ! uncertainty of MACS
  real(sgl)          :: dxs                              ! uncertainty of experimental cross section
  real(sgl)          :: macs                             ! Maxwellian cross section
  real(sgl)          :: rateastroeps                     ! cutoff value
  real(sgl)          :: rateastrorp(numT,(numZ+1)*(numN+1))   ! rate for astro residual production
  real(sgl)          :: ratio                            ! if 0<ratio<1 then x is between xl1 and xl2
  real(sgl)          :: xs                               ! help variable
  real(dbl)          :: Qth
  real(dbl)          :: Eth
!
! ************************ Output of reaction rates ********************
!
! partf      : integrated partition function
!
  col=''
  col(1)='T'
  un(1)='10^9_K'
  col(2)='reaction_rate'
  un(2)='cm3/mol/s'
  col(3)='MACS'
  un(3)='mb'
  col(4)='G(T)'
  un(4)=''
  Ncol=4
  reaction='('//parsym(k0)//',x)'
  quantity='reaction rate'
  write(*, '(/" 8. Thermonuclear reaction rates",/)')
  maxAastro = maxZastro + maxNastro
  do Acomp = 0, maxAastro
    do Zcomp = 0, maxZastro
      Ncomp = Acomp - Zcomp
      if (Ncomp < 0 .or. Ncomp > maxNastro) cycle
      do j = 1, nTmax
        if (rateastro(Zcomp, Ncomp, j) > 0.) then
          Z = ZZ(Zcomp, Ncomp, 0)
          A = AA(Zcomp, Ncomp, 0)
          massstring = '   '
          write(massstring(1:3), '(i3)') A
          ZAstring = '000000'
          write(ZAstring(1:3), '(i3.3)') Z
          write(ZAstring(4:6), '(i3.3)') A
          finalnuclide=trim(nuc(Z))//adjustl(massstring)
          topline=trim(targetnuclide)//trim(reaction)//trim(finalnuclide)//' '//trim(quantity)
          astrofile = 'astrorate'//ZAstring//'.tot'
          open (unit = 1, file = astrofile, status = 'replace')
          call write_header(topline,source,user,date,oformat)
          call write_target
          call write_reaction(reaction,0.D0,0.D0,0,0)
          if (nTmax == 1) then
            call write_real(2,'E-average [MeV]',astroE)
            call write_char(2,'Excites states contribution',yesno(flagastrogs))
          endif
          call write_residual(Z,A,finalnuclide)
          call write_datablock(quantity,Ncol,nTmax,col,un)
!         if (nTmax == 1) then
!           write(*, '(/" Reaction rate for Z=", i3, " A=", i3, " (", i3, a2, ") at <E>=", f8.5, &
!&            " MeV (Excited States Contribution : ", a1, ")"/)') Z, A, A, nuc(Z), astroE, yesno(flagastrogs)
!         else
!           write(*, '(/" Reaction rate for Z=", i3, " A=", i3, " (", i3, a2, ")"/)') Z, A, A, nuc(Z)
!         endif
          do i = 1, nTmax
            write(1, '(4es15.6)') T9(i), rateastro(Zcomp, Ncomp, i), macsastro(Zcomp, Ncomp, i), partf(i)
          enddo
!         write(*, '("    T        G(T)        Rate       MACS "/)')
!         do i = 1, nTmax
!           write(*, '(1x, f8.4, 3es12.5)') T9(i), partf(i), rateastro(Zcomp, Ncomp, i), macsastro(Zcomp, Ncomp, i)
!         enddo
          close (unit = 1)
          call write_outfile(astrofile,flagoutall)
          if (flagastroex) then
            do nex = 0, Nlast(Zcomp, Ncomp, 0)
              if (nex > 0 .and. tau(Zcomp, Ncomp, nex) == 0.) cycle
              write(astrofile(17:19),'("L",i2.2)') nex
              open (unit = 1, file = astrofile, status = 'replace')
              call write_header(topline,source,user,date,oformat)
              call write_target
              call write_reaction(reaction,0.D0,0.D0,0,0)
              call write_residual(Z,A,finalnuclide)
              call write_level(2,-1,nex,edis(Zcomp, Ncomp, nex),jdis(Zcomp, Ncomp, nex),parlev(Zcomp, Ncomp, nex), &
 &              tau(Zcomp, Ncomp, nex))
              call write_datablock(quantity,Ncol,nTmax,col,un)
!             write(*, '(/" Reaction rate for Z=", i3, " A=", i3, " (", i3, a2, ") to the excited states L", i2.2, &
!&              " at E=", f12.5, " MeV", /)') Z, A, A, nuc(Z), nex, edis(Zcomp, Ncomp, nex)
!             write(*, '("    T       Rate         MACS      ", "Branching"/)')
              do i = 1, nTmax
                branch = 0.
                if (rateastro(Zcomp, Ncomp, i) > 0.) branch = rateastroex(Zcomp, Ncomp, i, nex) / rateastro(Zcomp, Ncomp, i)
!               write(*, '(1x, f8.4, 3es12.5)') T9(i), rateastroex(Zcomp, Ncomp, i, nex), macsastroex(Zcomp, Ncomp, i, nex), branch
                write(1, '(4es15.6)') T9(i), rateastroex(Zcomp, Ncomp, i, nex), macsastroex(Zcomp, Ncomp, i, nex), branch
              enddo
              close (unit = 1)
              call write_outfile(astrofile,flagoutall)
            enddo
          endif
          if (flagracap .and. Zcomp == 0 .and. Acomp == 0) then
            write(*, '(/"    T      Rate(Eq)    Rate(DC)    MACS(Eq)    MACS(DC)  "/)')
            do i = 1, nTmax
              write(*, '(1x, f8.4, 4es12.5)') T9(i), rateastro(Zcomp, Ncomp, i) - rateastroracap(i), &
 &              rateastroracap(i), macsastro(Zcomp, Ncomp, i) - macsastroracap(i), macsastroracap(i)
            enddo
          endif
          exit
        endif
      enddo
    enddo
  enddo
!
! Comparison with experimental MACS at 30 keV
!
  if (astroE >= 0.029 .and. astroE <= 0.031) then
    Z = Ztarget
    A = Atarget
    ratio = 0.
    macs = 0.
    dmacs = 0.
    machar = trim(nuc(Z)) // '.macs'
    mafile = trim(path)//'gamma/macs/'//machar
    inquire (file = mafile, exist = lexist)
    if (lexist) then
      open (unit = 2, file = mafile, status = 'old')
      do
        read(2, '(4x, i4, 2(es12.4))', iostat = istat) ia, xs, dxs
        if (istat == -1) exit
        if (A == ia) then
          macs = xs
          dmacs = dxs
          exit
        else
          cycle
        endif
      enddo
      close (unit = 2)
    endif
    if (macs > 0.) ratio = macsastro(0, 0, 1) / macs
    astrofile = 'macs.g'
    open (unit = 1, file = astrofile, status = 'replace')
    write(1, '("# Z   A     MACS(mb)    Exp(mb)     dEXP", "(mb)   MACS/Exp")')
    write(1, '(2i4, 4es12.5)') Z, A, macsastro(0, 0, 1), macs, dmacs, ratio
    close (unit = 1)
  endif
!
! Write results to separate files
!
  do type=1,4
    if (type == 1) then
      sym='g'
      Zix=0
      Nix=0
      Qth=Qres(0, 0, 0)
      Eth=Ethresh(0, 0, 0)
    endif
    if (type == 2) then
      sym='p'
      Zix=1
      Nix=0
      Qth=Qres(0, 0, 2)
      Eth=Ethresh(0, 0, 2)
    endif
    if (type == 3) then
      sym='a'
      Zix=2
      Nix=2
      Qth=Qres(0, 0, 6)
      Eth=Ethresh(0, 0, 6)
    endif
    if (type == 4) then
      if (.not.flagfission) exit
      sym='f'
      Qth=0.
      Eth=0.
    endif
    astrofile = 'astrorate.'//sym
    col(1)='T'
    un(1)='10^9_K'
    col(2)='reaction_rate'
    un(2)='cm3/mol/s'
    col(3)='MACS'
    un(3)='mb'
    col(4)='G(T)'
    un(4)=''
    Ncol=4
    reaction='('//parsym(k0)//','//sym//')'
    quantity='reaction rate'
    topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
    open (unit = 1, file = astrofile, status = 'replace')
    call write_header(topline,source,user,date,oformat)
    call write_target
    call write_reaction(reaction,Qth,Eth,0,0)
    if (nTmax == 1) then
      call write_real(2,'E-average [MeV]',astroE)
      call write_char(2,'Excites states contribution',yesno(flagastrogs))
    endif
    call write_datablock(quantity,Ncol,nTmax,col,un)
    if (type <= 3) then
      do i = 1, nTmax
        write(1, '(4es15.6)') T9(i), rateastro(Zix, Nix, i), macsastro(Zix, Nix, i), partf(i)
      enddo
    else
      do i = 1, nTmax
        write(1, '(4es15.6)') T9(i), rateastrofis(i), macsastrofis(i), partf(i)
      enddo
    endif
    close (unit = 1)
    call write_outfile(astrofile,flagoutall)
  enddo
!
! output partial rates(n,g) to given excited states in a specific file
!
  if (flagastroex) then
    col(4)='Branching_ratio'
    un(4)=''
    col(5)='G(T)'
    un(5)=''
    Ncol=5
    reaction='('//parsym(k0)//','//sym//')'
    quantity='reaction rate'
    do type=1,3
      if (type == 1) then
        sym='g'
        Zix=0
        Nix=0
      endif
      if (type == 2) then
        sym='p'
        Zix=1
        Nix=0
      endif
      if (type == 3) then
        sym='a'
        Zix=2
        Nix=2
      endif
      do nex = 1, Nlast(Zix, Nix, 0)
        if (tau(Zix, Nix, nex) /= 0.) goto 10
      enddo
      cycle 
   10 do nex = 0, Nlast(Zix, Nix, 0)
        if (nex == 0 .or. tau(Zix, Nix, nex) /= 0.) then
          astrofile = 'astrorate.'//sym//'.L00'
          write(astrofile(14:15), '(i2.2)') nex
          topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' - Level '//astrofile(14:15)
          open (unit = 1, file = astrofile, status = 'replace')
          call write_header(topline,source,user,date,oformat)
          call write_target
          call write_reaction(reaction,Qth,Eth,0,0)
          call write_level(2,-1,nex,edis(Zix, Nix, nex),jdis(Zix, Nix, nex),parlev(Zix, Nix, nex),tau(Zix, Nix, nex))
          if (nTmax == 1) then
            call write_real(2,'E-average [MeV]',astroE)
            call write_char(2,'Excites states contribution',yesno(flagastrogs))
          endif
          call write_datablock(quantity,Ncol,nTmax,col,un)
          do i = 1, nTmax
            branch = 0.
            if (rateastro(Zix, Nix, i) > 0.) branch = rateastroex(Zix, Nix, i, nex) / rateastro(Zix, Nix, i)
            write(1, '(5es15.6)') T9(i), rateastroex(Zix, Nix, i, nex), macsastroex(Zix, Nix, i, nex), branch, partf(i)
          enddo
          close (unit = 1)
          call write_outfile(astrofile,flagoutall)
        endif
      enddo
    enddo
  endif
  rateastroeps = 1.e-10
  iresprod = 0
  do Acomp = 0, maxAastro
    do Zcomp = 0, maxZastro
      Ncomp = Acomp - Zcomp
      if (Ncomp < 0 .or. Ncomp > maxNastro) cycle
      if (Zcomp == 0 .and. Ncomp == 0) cycle
      if (Zcomp == 0 .and. Ncomp == 1) cycle
      if (Zcomp == 1 .and. Ncomp == 0) cycle
      if (Zcomp == 2 .and. Ncomp == 2) cycle
      iwriterp = 0
      do i = 1, nTmax
        if (rateastro(Zcomp, Ncomp, i) > rateastroeps) iwriterp = 1
      enddo
      if (iwriterp == 1) then
        iresprod = iresprod + 1
        do i = 1, nTmax
          rateastrorp(i, iresprod) = rateastro(Zcomp, Ncomp, i)
        enddo
        zrespro(iresprod) = ZZ(Zcomp, Ncomp, 0)
        arespro(iresprod) = AA(Zcomp, Ncomp, 0)
      endif
    enddo
  enddo
  astrofile = 'astrorate.tot'
  Zix = parZ(k0)
  Nix = parN(k0)
  reaction='('//parsym(k0)//',x)'
  topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
  open (unit = 1, file = astrofile, status = 'replace')
  call write_header(topline,source,user,date,oformat)
  call write_target
  call write_level(2,-1,0,edis(Zix, Nix, 0),jdis(Zix, Nix, 0),parlev(Zix, Nix, 0),0.)
  call write_reaction(reaction,0.D0,0.D0,0,0)
  call write_integer(2,'reactions',iresprod+5)
  write(1,' ("#   astrogs: ",a1)') yesno(flagastrogs)
  if (nTmax == 1) then
    call write_real(2,'E-average [MeV]',astroE)
  else
    if ( .not. flagastrogs) then
      if (nonthermlev ==  - 1) then
        call write_char(2,'comment','fully thermalized target')
      endif
      if (nonthermlev > 0) then
        call write_integer(2,'comment: thermalized target on the GS excluding level ',nonthermlev)
      endif
      if (nonthermlev == 0) then
        call write_integer(2,'comment: thermalized target in level ',nonthermlev)
      endif
    else
      if (Ltarget == 0) then
        call write_char(2,'comment','non-thermalized target in its ground state')
      endif
      if (Ltarget > 0) then
        call write_integer(2,'comment: non-thermalized target in its excited state ',Ltarget)
      endif
    endif
  endif
  un='cm3/mol/s'
  col(1)='T'
  col(1)='10^9_K'
  col(2)='G(T)'
  un(2)=''
  massstring='   '
  write(massstring,'(i3)') AA(0, 0, 0)
  col(3)='('//parsym(k0)//',g)'//trim(nuc(ZZ(0,0,0)))//adjustl(trim(massstring))
  massstring='   '
  write(massstring,'(i3)') AA(0, 1, 0)
  col(4)='('//parsym(k0)//',n)'//trim(nuc(ZZ(0,1,0)))//adjustl(trim(massstring))
  massstring='   '
  write(massstring,'(i3)') AA(1, 0, 0)
  col(5)='('//parsym(k0)//',p)'//trim(nuc(ZZ(1,0,0)))//adjustl(trim(massstring))
  massstring='   '
  write(massstring,'(i3)') AA(2, 2, 0)
  col(6)='('//parsym(k0)//',a)'//trim(nuc(ZZ(2,2,0)))//adjustl(trim(massstring))
  col(7)='fission'
  do ires= 1, iresprod
    massstring='   '
    write(massstring,'(i3)') arespro(ires)
    col(7 + ires)=trim(nuc(zrespro(ires)))//trim(adjustl(massstring))
  enddo
  Ncol=iresprod + 7
  call write_datablock(quantity,Ncol,nTmax,col,un)
  do i = 1, nTmax
    write(1, '(200es15.6)') T9(i), partf(i), rateastro(0, 0, i), rateastro(0, 1, i), rateastro(1, 0, i), rateastro(2, 2, i), &
 &    rateastrofis(i), (rateastrorp(i, ires), ires = 1, iresprod)
  enddo
  close (unit = 1)
  call write_outfile(astrofile,flagoutall)
  return
end subroutine astroout
! Copyright A.J. Koning 2021
