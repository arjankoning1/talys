subroutine binary
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Binary reaction results
!
! Author    : Arjan Koning and Stephane Hilaire
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl              ! single precision kind
! All global variables
!   numJ             ! maximum J - value
! Variables for numerics
!   popeps           ! limit for population cross sections
!   xseps            ! limit for cross sections
! Variables for output
!   flagcheck        ! flag for output of numerical checks
!   flagpop          ! flag for output of population
!   flagspec         ! flag for output of spectra
! Variables for basic reaction
!   flagchannels     ! flag for exclusive channels calculation
!   flagrecoil       ! flag for calculation of recoils
! Variables for compound reactions
!   flagcomp         ! flag for compound angular distribution calculation
! Variables for input energies
!   flaginitpop      ! flag for initial population distribution
!   nin              ! counter for incident energy
!   nin0             ! counter for incident energy
! Variables for main input
!   k0               ! index of incident particle
!   Ltarget          ! excited level of target
! Variables for preequilibrium
!   pespinmodel      ! model for pre - equilibrium or compound spin distribution
! Variables for fission
!   flagfission      ! flag for fission
! Variables for gamma rays
!   flagracap        ! flag for radiative capture model
! Variables for OMP
!   flagomponly      ! flag to execute ONLY an optical model calculation
! Variables for energy grid
!   deltaE           ! energy bin around outgoing energies
!   ebegin           ! first energy point of energy grid
!   egrid            ! outgoing energy grid
!   Einc             ! incident energy in MeV
!   Etop             ! top of outgoing energy bin
! Variables for energies
!   eendhigh         ! last energy point for energy grid for any particle
!   Einc0            ! incident energy in MeV
!   eoutdis          ! outgoing energy of discrete state reaction
!   flagpreeq        ! flag for pre - equilibrium calculation
!   nendisc          ! last discrete bin
! Variables for excitation energy grid
!   deltaEx          ! excitation energy bin for population arrays
!   Ex               ! excitation energy
!   Exmax            ! maximum excitation energy for residual nucleus
!   maxex            ! maximum excitation energy bin for residual nucleus
! Variables for binary reactions
!   binemissum       ! integrated binary emission spectrum
!   Eaveragebin      ! average outgoing energy
!   feedbinary       ! feeding from first compound nucleus
!   sfactor          ! spin factor
!   xscompdisc       ! compound cross section for discrete state
!   xscompdisctot    ! compound cross section summed over discrete states
!   xscompel         ! compound elastic cross section
!   xscompnonel      ! total compound non - elastic cross section
!   xscompound       ! total compound cross section
!   xsconttot        ! total cross section for continuum
!   xsdircont        ! direct cross section for continuum
!   xsdirect         ! total direct cross section
!   xsdisc           ! total cross section for discrete state
!   xsdisctot        ! total cross section summed over discrete states
!   xselastot        ! total elastic cross section (shape + compound)
!   xsnonel          ! non - elastic cross
!   xspopdir         ! direct population cross section per nucleus
!   xspopex0         ! binary population cross section
! Variables for endf data
!   xscompel6        ! compound elastic cross section
!   xsnonel6         ! non - elastic cross section
! Variables for binary emission spectra
!   binnorm          ! normalization factor for binary spectra
!   xsbinemis        ! cross section for emission from first compound nucleus
! Variables for incident channel
!   partdecay        ! total decay per particle
!   popdecay         ! decay from population
!   preeqpop         ! pre - equilibrium population
!   preeqpopex       ! pre - equilibrium population
!   xscompcont       ! compound cross section for continuum
!   xsbinary         ! cross section from initial compound to residual nucleus
!   xsdirdisc        ! direct cross section for discrete state direct cross section
!   xsdirdiscsum     ! total direct cross section
!   xsdirdisctot     ! direct cross section summed over discrete states
!   xselasinc        ! total elastic cross section (neutrons only) for inc. channel
!   xsgrsum          ! sum over giant resonance cross sections
!   xsgrtot          ! total smoothed giant resonance cross section
!   xspop            ! population cross section
!   xspopex          ! population cross section summed over spin and parity
!   xspopexP         ! population cross section per parity
!   xspopnuc         ! population cross section per nucleus
!   xspopnucP        ! population cross section per nucleus per parity
!   xspreeqsum       ! total preequilibrium cross section summed over particles
!   xspreeqtot       ! preequilibrium cross section per particle type
!   xsreacinc        ! reaction cross section for incident channel
! Variables for nuclides
!   AA               ! mass number of residual nucleus
!   Nindex           ! neutron number index for residual nucleus
!   NN               ! neutron number of residual nucleus
!   parskip          ! logical to skip outgoing particle
!   targetP          ! parity of target
!   targetspin       ! spin of target
!   Zindex           ! charge number index for residual nucleus
!   ZZ               ! charge number of residual nucleus
! Constants
!   nuc              ! symbol of nucleus
!   pardis           ! parity distribution
!   parN             ! neutron number of particle
!   parname          ! name of particle
!   parZ             ! charge number of particle
! Variables for levels
!   jdis             ! spin of level
!   parlev           ! parity of level
! Variables for level density
!   Nlast            ! last discrete level
! Variables for direct capture
!   xsracape         ! direct radiative capture cross section
!   xsracapecont     ! direct radiative capture continuum cross section
!   xsracapedisc     ! direct radiative capture discrete cross section
!   xsracappop       ! population cross section for radiative capture
!   xsracappopex     ! population cross section for radiative capture
! Variables for preequilibrium initialization
!   maxJph           ! maximal spin for particle - hole states
!
!
! *** Declaration of local data
!
  implicit none
  character(len=18) :: reaction   ! reaction
  character(len=132) :: topline    ! topline
  character(len=21) :: binfile
  character(len=13) :: Estr
  character(len=15) :: col(200) 
  character(len=15) :: un(200)   
  character(len=80) :: quantity   ! quantity
  integer   :: A           ! mass number of target nucleus
  integer   :: J           ! spin of level
  integer   :: N           ! neutron number of residual nucleus
  integer   :: nen         ! energy counter
  integer   :: nex         ! excitation energy bin of compound nucleus
  integer   :: Nix         ! neutron number index for residual nucleus
  integer   :: Ncol
  integer   :: Nk
  integer   :: NL          ! last discrete level
  integer   :: odd         ! odd (1) or even (0) nucleus
  integer   :: parity      ! parity
  integer   :: type        ! particle type
  integer   :: Z           ! charge number of target nucleus
  integer   :: Zix         ! charge number index for residual nucleus
  real(sgl) :: ald         ! level density parameter
  real(sgl) :: Eaveragesum ! help variable
  real(sgl) :: Eex         ! excitation energy
  real(sgl) :: factor      ! multiplication factor
  real(sgl) :: frac        ! help variable
  real(sgl) :: ignatyuk    ! function for energy dependent level density parameter a
  real(sgl) :: popepsA     ! limit for population cross sections per energy
  real(sgl) :: spindis     ! Wigner spin distribution
  real(sgl) :: spincut     ! spin cutoff factor
  real(sgl) :: sc          ! spin cutoff factor
  real(sgl) :: term        ! help variable
  real(sgl) :: xscompall   ! total compound cross section summed over particles
  real(sgl) :: xsb         ! help variable
!
! * Add direct and pre-equilibrium cross sections to population arrays *
!
! No binary reaction for initial excitation energy population.
!
  xselastot = xselasinc
  if (flaginitpop) return
  if (flagomponly .and. .not. flagcomp) return
!
! Depending on whether the nucleus is even or odd, the quantum number J appearing in the array xspop represents either J or J+0.5.
!
  do type = 0, 6
    if (parskip(type)) cycle
    Zix = Zindex(0, 0, type)
    Nix = Nindex(0, 0, type)
    NL = Nlast(Zix, Nix, 0)
    do nex = 0, NL
      if (xsdirdisc(type, nex) /= 0.) then
        J = int(jdis(Zix, Nix, nex))
        parity = parlev(Zix, Nix, nex)
        term = xsdirdisc(type,nex)
        xspop(Zix, Nix, nex, J, parity) = xspop(Zix, Nix, nex, J, parity) + term
        xspopex(Zix, Nix, nex) = xspopex(Zix, Nix, nex) + term
        if (flagpop) then
          xspopnucP(Zix, Nix, parity) = xspopnucP(Zix, Nix, parity) + term
          xspopexP(Zix, Nix, nex, parity) = xspopexP(Zix, Nix, nex, parity) + term
          popdecay(type, nex, J, parity) = popdecay(type, nex, J, parity) + term
          partdecay(type, parity) = partdecay(type, parity) + term
        endif
      endif
      xspopex0(type, nex) = xspopex(Zix, Nix, nex)
    enddo
    xspopdir(Zix, Nix) = xsdirdisctot(type)
    xsbinary(type) = xsbinary(type) + xsdirdisctot(type)
    xspopnuc(Zix, Nix) = xspopnuc(Zix, Nix) + xsdirdisctot(type)
!
! ******* Assign spin distribution to pre-equilibrium population *******
!
! ignatyuk   : function for energy dependent level density parameter a
!
    if (flagpreeq) then
      if (pespinmodel <= 2) then
        popepsA = popeps / max(5 * maxex(Zix, Nix), 1)
        do nex = NL + 1, maxex(Zix, Nix)
          Eex = Ex(Zix, Nix, nex)
          ald = ignatyuk(Zix, Nix, Eex, 0)
          sc = spincut(Zix, Nix, ald, Eex, 0, 0)
          do parity = - 1, 1, 2
            do J = 0, maxJph
              if (xspopex(Zix, Nix, nex) > popepsA) sfactor(Zix, Nix, nex, J, parity) = &
                  xspop(Zix, Nix, nex, J, parity) / xspopex(Zix, Nix, nex)
              if (pespinmodel == 1 .and. sfactor(Zix, Nix, nex, J, parity) > 0.) then
                preeqpop(Zix, Nix, nex, J, parity) = sfactor(Zix, Nix, nex, J, parity) * preeqpopex(Zix, Nix, nex)
              else
                factor = spindis(sc,real(J)) * pardis
                preeqpop(Zix, Nix, nex, J, parity) = factor * preeqpopex(Zix, Nix, nex)
              endif
            enddo
          enddo
        enddo
      endif
      do nex = NL + 1, maxex(Zix, Nix)
        xspopex(Zix, Nix, nex) = xspopex(Zix, Nix, nex) + preeqpopex(Zix, Nix, nex)
        do parity = - 1, 1, 2
          do J = 0, maxJph
            term = preeqpop(Zix, Nix, nex, J, parity)
            xspop(Zix, Nix, nex, J, parity) = xspop(Zix, Nix, nex, J, parity) + term
            if (flagpop) then
              xspopnucP(Zix, Nix, parity) = xspopnucP(Zix, Nix, parity) + term
              xspopexP(Zix, Nix, nex, parity) = xspopexP(Zix, Nix, nex, parity) + term
              popdecay(type, nex, J, parity) = popdecay(type, nex, J, parity) + term
              partdecay(type, parity) = partdecay(type, parity) + term
            endif
          enddo
        enddo
      enddo
      xspopnuc(Zix, Nix) = xspopnuc(Zix, Nix) + xspreeqtot(type) + xsgrtot(type)
      xsbinary(type) = xsbinary(type) + xspreeqtot(type) + xsgrtot(type)
    endif
!
! ************* Other total binary cross sections **********************
!
    xscompdisctot(type) = 0.
    do nex = 0, NL
      if (type == k0 .and. nex == Ltarget) cycle
      xsdisc(type, nex) = xspopex0(type, nex)
      xscompdisc(type, nex) = xsdisc(type, nex) - xsdirdisc(type, nex)
      xscompdisctot(type) = xscompdisctot(type) + xscompdisc(type, nex)
    enddo
    xsdisctot(type) = xsdirdisctot(type) + xscompdisctot(type)
    xsdircont(type) = xspreeqtot(type) + xsgrtot(type)
    if (type == 0 .and. flagracap) then
      xsdisctot(type) = xsdisctot(type) + xsracapedisc
      xsdircont(type) = xsdircont(type) + xsracapecont
      xspopnuc(Zix, Nix) = xspopnuc(Zix, Nix) + xsracape
      xsbinary(type) = xsbinary(type) + xsracape
      do nex = 0, NL
        if (xsracappopex(nex) /= 0.) then
          J = int(jdis(Zix, Nix, nex))
          parity = parlev(Zix, Nix, nex)
          term=xsracappopex(nex)
          xspop(Zix, Nix, nex, J, parity) = xspop(Zix, Nix, nex, J, parity) + term
          xspopex(Zix, Nix, nex) = xspopex(Zix, Nix, nex) + term
          xspopex0(type, nex) = xspopex0(type, nex) + term
          if (flagpop) then
            xspopnucP(Zix, Nix, parity) = xspopnucP(Zix, Nix, parity) + term
            xspopexP(Zix, Nix, nex, parity) = xspopexP(Zix, Nix, nex, parity) + term
            popdecay(type, nex, J, parity) = popdecay(type, nex, J, parity) + term
            partdecay(type, parity) = partdecay(type, parity) + term
          endif
        endif
      enddo
      do nex = NL + 1, maxex(Zix, Nix)
        xspopex(Zix, Nix, nex) = xspopex(Zix, Nix, nex) + xsracappopex(nex)
        do parity = - 1, 1, 2
          do J = 0, numJ
            term = xsracappop(nex, J, parity)
            xspop(Zix, Nix, nex, J, parity) = xspop(Zix, Nix, nex, J, parity) + term
            if (flagpop) then
              xspopnucP(Zix, Nix, parity) = xspopnucP(Zix, Nix, parity) + term
              xspopexP(Zix, Nix, nex, parity) = xspopexP(Zix, Nix, nex, parity) + term
              popdecay(type, nex, J, parity) = popdecay(type, nex, J, parity) + term
              partdecay(type, parity) = partdecay(type, parity) + term
            endif
          enddo
        enddo
      enddo
    endif
    xsdirect(type) = xsdirdisctot(type) + xsdircont(type)
    if (xscompcont(type) < xseps) xscompcont(type) = 0.
    xsconttot(type) = xscompcont(type) + xsdircont(type)
    xscompound(type) = xscompdisctot(type) + xscompcont(type)
  enddo
  xscompel = xspopex0(k0, Ltarget)
  xscompel6(nin) = xscompel
  xselastot = xselasinc + xscompel
  xsnonel = max(xsreacinc - xscompel, 0.)
  xsnonel6(nin) = xsnonel
  xscompall = max(xsreacinc - xsdirdiscsum - xspreeqsum - xsgrsum - xsracape, 0.)
  xscompnonel = xscompall - xscompel
  xscompnonel = max(xscompnonel, 0.)
  if (k0 == 0) xstotinc = xselastot + xsnonel
!
! ***************** Create binary feeding channels *********************
!
! This is necessary for exclusive cross sections
!
  if (flagchannels) then
    do type = 0, 6
      if (parskip(type)) cycle
      Zix = Zindex(0, 0, type)
      Nix = Nindex(0, 0, type)
      do nex = 0, maxex(Zix, Nix)
        feedbinary(type, nex) = xspopex(Zix, Nix, nex)
      enddo
    enddo
    feedbinary(k0, Ltarget) = 0.
  endif
!
! *************** Interpolate decay on emission spectrum ***************
!
! binaryspectra: subroutine for creation of binary spectra
!
  if (flagspec .or. flagrecoil) call binaryspectra
!
! ********************** Average emission energy ***********************
!
  if (flagrecoil .or. flagspec) then
    do type = 0, 6
      if (parskip(type)) cycle
      binemissum(type) = 0.
      Eaveragesum = 0.
      do nen = ebegin(type), nendisc(type)
        binemissum(type) = binemissum(type) + xsbinemis(type, nen) * deltaE(nen)
        Eaveragesum = Eaveragesum + egrid(nen) * xsbinemis(type, nen) * deltaE(nen)
      enddo
      Zix = Zindex(0, 0, type)
      Nix = Nindex(0, 0, type)
      NL = Nlast(Zix, Nix, 0)
      if (eoutdis(type, NL) > 0.) then
        frac = Etop(nendisc(type)) - eoutdis(type, NL)
        binemissum(type) = binemissum(type) - xsbinemis(type, nendisc(type)) * frac
        Eaveragesum = Eaveragesum - egrid(nendisc(type)) * xsbinemis(type, nendisc(type)) * frac
      endif
      xsb = binemissum(type)
      do nex = 0, NL
        if (type == k0 .and. nex == Ltarget) cycle
        xsb = xsb + xsdisc(type, nex)
        Eaveragesum = Eaveragesum + xsdisc(type, nex) * eoutdis(type, nex)
      enddo
      if (xsb > 0.) then
        Eaveragebin(type) = Eaveragesum / xsb
      else
        Eaveragebin(type) = 0.
      endif
    enddo
  endif
!
! ************ Output of population after binary emission **************
!
  if (flagpop) then
    quantity = 'binary cross sections'
    reaction='('//parsym(k0)//',bin)'
    Estr=''
    write(Estr,'(es13.6)') Einc
    binfile = 'binE0000.000.out'
    write(binfile(5:12), '(f8.3)') Einc
    write(binfile(5:8), '(i4.4)') int(Einc)
    topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
    open (unit = 1, file = binfile, status = 'replace')
    call write_header(topline,source,user,date,oformat)
    call write_target
    call write_reaction(reaction,0.D0,0.D0,0,0)
    call write_real(2,'E-incident [MeV]',Einc)
    write(*, '(/" ########## BINARY CHANNELS ###########")')
    write(*, '(/" ++++++++++ BINARY CROSS SECTIONS ++++++++++"/)')
    if (flagfission) write(*, '(" fission  channel", 23x, ":", es12.5)') xsbinary( - 1)
    do type = 0, 6
      if (parskip(type)) cycle
      Z = ZZ(0, 0, type)
      N = NN(0, 0, type)
      A = AA(0, 0, type)
      write(*, '(1x, a8, " channel to Z=", i3, " N=", i3, " (", i3, a2, "):", es12.5)') parname(type), Z, N, A, nuc(Z), &
 &      xsbinary(type)
    enddo
    if (flagspec) then
      quantity = 'binary emission spectra'
      write(*, '(/" Binary emission spectra"/)')
      un = 'mb'
      col(1) = 'E-out'
      un(1) = 'MeV'
      do type = 0, 6
        col(2+type) = parname(type)
      enddo
      Ncol = 8
      Nk = eendhigh - ebegin(0) + 1
      call write_datablock(quantity,Ncol,Nk,col,un)
      do nen = ebegin(0), eendhigh
        write(1, '(8es15.6)') egrid(nen), (xsbinemis(type, nen), type = 0, 6)
      enddo
    endif
    if (flagspec .and. flagcheck) then
      write(*, '(/" ++++++++++ CHECK OF INTEGRATED ", "BINARY EMISSION SPECTRA ++++++++++"/)')
      write(*, '(13x, "Continuum cross section  Integrated", " spectrum  Compound normalization Average emission energy"/)')
      un = 'mb'
      col(1) = 'particle'
      un(1) = ''
      col(1) = 'particle'
      col(2) = 'Continuum'
      col(3) = 'Int. spectrum'
      col(4) = 'Comp. norm.'
      col(5) = 'Av. emission E.'
      un(5) = 'MeV'
      Ncol = 5
      Nk = 7
      quantity = 'check of binary emission spectra'
      call write_datablock(quantity,Ncol,Nk,col,un)
      do type = 0, 6
        if (parskip(type)) cycle
        write(1, '(3x, a8, 4x, 4es15.6)') parname(type), &
 &        xscompcont(type) + xspreeqtot(type) + xsgrtot(type), binemissum(type), binnorm(type), Eaveragebin(type)
      enddo
    endif
    write(*, '(/" ++++++++++ POPULATION AFTER BINARY EMISSION ++++++++++",/)')
    reaction='('//parsym(k0)//',x)' 
    quantity='post-binary population'
    do type = 0, 6
      if (parskip(type)) cycle
      Zix = Zindex(0, 0, type)
      Nix = Nindex(0, 0, type)
      NL = Nlast(Zix, Nix, 0)
      Z = ZZ(0, 0, type)
      N = NN(0, 0, type)
      A = AA(0, 0, type)
      if (xspopnuc(Zix, Nix) == 0.) cycle
      odd = mod(A, 2)
      call write_char(2,'ejectile',parname(type))
      call write_double(2,'post-binary population [mb]',xspopnuc(Zix, Nix))
      call write_real(2,'maximum excitation energy [MeV]',Exmax(Zix, Nix))
      if (maxex(Zix, Nix) > NL) then
        call write_integer(2,'number of discrete levels',NL)
        call write_integer(2,'number of continuum bins',maxex(Zix, Nix) - NL)
        call write_real(2,'continuum bin size [MeV]',deltaEx(Zix, Nix, maxex(Zix, Nix)))
      endif
      un = 'mb'
      col(1)='bin'
      col(2)='Ex'
      un(2) = ''
      col(3)='population'
      do J = 0, numJ
        col(3+2*J+1)='JP=    -'
        write(col(3+2*J+1)(4:7),'(f4.1)') J+0.5*odd
        col(3+2*J+2)='JP=    +'
        write(col(3+2*J+2)(4:7),'(f4.1)') J+0.5*odd
      enddo
      Ncol = 2*(numJ + 1) + 3
      Nk = maxex(Zix, Nix) + 1
      call write_datablock(quantity,Ncol,Nk,col,un)
      do nex = 0, maxex(Zix, Nix)
        write(1, '(3x, i6, 6x, 200es15.6)') nex, Ex(Zix, Nix, nex), &
 &        xspopex(Zix, Nix, nex), ((xspop(Zix, Nix, nex, J, parity), parity = - 1, 1, 2), J = 0, numJ)
      enddo
    enddo
  endif
  close (unit = 1)
  call write_outfile(binfile,flagoutall)
!
! Remove compound elastic scattering from population of target state.
!
  xspopex(parZ(k0), parN(k0), Ltarget) = 0.
  xspop(parZ(k0), parN(k0), Ltarget, int(targetspin), targetP) = 0.
  preeqpopex(parZ(k0), parN(k0), Ltarget) = 0.
  Einc0 = Einc
  nin0 = nin
!
! **************************** Recoils *********************************
!
! binaryrecoil: subroutine for recoil for binary reaction
!
  if (flagrecoil) call binaryrecoil
  return
end subroutine binary
! Copyright A.J. Koning 2021
