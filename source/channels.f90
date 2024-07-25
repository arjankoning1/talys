subroutine channels
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Exclusive reaction channels
!
! Author    : Marieke Duijvestijn and Arjan Koning
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
!   numchantot        ! maximum number of exclusive channels
!   numen             ! maximum number of outgoing energies
!   numex             ! maximum number of excitation energies
!   numia             ! maximum number of alphas in channel description
!   numid             ! maximum number of deuterons in channel description
!   numih             ! maximum number of helions in channel description
!   numin             ! maximum number of neutrons in channel description
!   numip             ! maximum number of protons in channel description
!   numit             ! maximum number of tritons in channel description
!   numlev            ! maximum number of discrete levels
!   numNchan          ! maximum number of outgoing neutron in individual channel description
!   numpar            ! number of particles
!   numZchan          ! maximum number of outgoing protons in individual channel description
! Variables for numerics
!   maxchannel        ! maximal number of outgoing particles in individual channel description
!   maxN              ! maximal number of neutrons away from initial compound nucleus
!   maxZ              ! maximal number of protons away from initial compound nucleus
!   xseps             ! limit for cross sections
! Variables for output
!   flagcheck         ! flag for output of numerical checks
!   flagspec          ! flag for output of spectra
! Variables for fission
!   flagfission       ! flag for fission
! Variables for input energies
!   flaginitpop       ! flag for initial population distribution
! Variables for main input
!   k0                ! index of incident particle
!   Ninit             ! neutron number of initial compound nucleus
!   Zinit             ! charge number of initial compound nucleus
! Variables for energy grid
!   deltaE            ! energy bin around outgoing energies
!   ebegin            ! first energy point of energy grid
!   egrid             ! outgoing energy grid
!   Etop              ! top of outgoing energy bin
! Variables for energies
!   eend              ! last energy point of energy grid
!   eoutdis           ! outgoing energy of discrete state reaction
!   Ethrexcl          ! threshold incident energy for exclusive channel
!   idchannel         ! identifier for exclusive channel
!   nendisc           ! last discrete bin
!   Qexcl             ! Q - value for exclusive channel
!   reacstring        ! string for exclusive reaction channel
! Variables for excitation energy grid
!   maxex             ! maximum excitation energy bin for residual nucleus
! Variables for existence libraries
!   chanopen          ! flag to open channel with first non - zero cross s
!   idnumfull         ! flag to designate maximum number of exclusive ch.
! Variables for exclusive channels
!   Eavchannel        ! channel average energy
!   Especsum          ! total emission energy
!   exclbranch        ! exclusive channel yield per isomer
!   fisstring         ! string for exclusive fission reaction channel
!   gamexcl           ! exclusive gamma cross section per excitation ener
!   gmult             ! continuum gamma multiplicity
!   idnum             ! counter for exclusive channel
!   specemis          ! exclusive emission contribution
!   xschancheck       ! integrated channel spectra
!   xschaniso         ! channel cross section per isomer
!   xschannel         ! channel cross section
!   xschannelsp       ! channel cross section spectra
!   xsexcl            ! exclusive cross section per excitation energy
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
!   feedexcl          ! feeding terms from compound excitation ene
!   fisfeedex         ! fission contribution from excitation energy bin
!   popexcl           ! population cross section of bin just before decay
!   xsinitpop         ! initial population cross section
! Variables for binary emission spectra
!   binemis           ! emission spectra from initial compound nucleus
! Variables for incident channel
!   channelsum        ! sum over exclusive channel cross sections
!   xsabs             ! absorption cross section
!   xspopnuc          ! population cross section per nucleus
!   xsreacinc         ! reaction cross section for incident channel
! Variables for nuclides
!   parinclude        ! logical to include outgoing particle
!   parskip           ! logical to skip outgoing particle
!   targetE           ! excitation energy of target
! Constants
!   parN              ! neutron number of particle
!   parsym            ! symbol of particle
!   parZ              ! charge number of particle
! Variables for levels
!   edis              ! energy of level
!   tau               ! lifetime of state in seconds
! Variables for level density
!   Nlast             ! last discrete level
! Variables for masses
!   S                 ! separation energy
!   specmass          ! specific mass for residual nucleus
!
! *** Declaration of local data
!
  implicit none
  character(len=1) :: apos                                                 ! apostrophe
  integer          :: Aix                                                  ! mass number index for residual nucleus
  integer          :: i                                                    ! counter
  integer          :: i2                                                   ! value
  integer          :: ia                                                   ! mass number from abundance table
  integer          :: iaend                                                ! end of particle summation
  integer          :: id                                                   ! counter for deuterons
  integer          :: idc                                                  ! help variable
  integer          :: idd                                                  ! identifier for channel
  integer          :: idend                                                ! end of particle summation
  integer          :: ident                                                ! exclusive channel identifier
  integer          :: identorg(0:numpar)                                   ! identifier for previous channel
  integer          :: idorg                                                ! identifier for previous channel
  integer          :: ih                                                   ! hole number
  integer          :: ihend                                                ! end of particle summation
  integer          :: in                                                   ! counter for neutrons
  integer          :: inend                                                ! help variable
  integer          :: ip                                                   ! particle number
  integer          :: ipend                                                ! end of particle summation
  integer          :: it                                                   ! counter for tritons
  integer          :: itend                                                ! end of particle summation
  integer          :: Ncomp                                                ! neutron number index for compound nucleus
  integer          :: nen                                                  ! energy counter
  integer          :: Nend                                                 ! maximal neutron number
  integer          :: nex                                                  ! excitation energy bin of compound nucleus
  integer          :: nexout                                               ! energy index for outgoing energy
  integer          :: Nix                                                  ! neutron number index for residual nucleus
  integer          :: NL                                                   ! last discrete level
  integer          :: npart                                                ! number of particles in outgoing channel
  integer          :: Ntot                                                 ! number of nucleon units in exit channel
  integer          :: type                                                 ! particle type
  integer          :: type2                                                ! particle type
  integer          :: Zcomp                                                ! proton number index for compound nucleus
  integer          :: Zend                                                 ! maximal charge number
  integer          :: Zix                                                  ! charge number index for residual nucleus
  integer          :: Ztot                                                 ! number of nucleon units in exit channel
  real(sgl)        :: Eaveragesum                                          ! help variable
  real(sgl)        :: emissum                                              ! integrated binary emission spectrum
  real(sgl)        :: fissum                                               ! help variable
  real(sgl)        :: frac                                                 ! help variable
  real(sgl)        :: specexcl(0:numchantot, 0:numpar, 0:numex+1, 0:numen) ! exclusive spectra per excitation energy
  real(sgl)        :: term                                                 ! help variable
  real(sgl)        :: term1                                                ! help variable
  real(sgl)        :: term2                                                ! help variable
  real(sgl)        :: term3                                                ! help variable
!
! ************************ Initialization ******************************
!
  channelsum = 0.
  xsabs = 0.
!
! Initially, all the flux is in the initial compound state.
!
  if (flaginitpop) then
    xsexcl(0, maxex(0, 0) + 1) = xsinitpop
  else
    xsexcl(0, maxex(0, 0) + 1) = xsreacinc
  endif
  gamexcl(0, maxex(0, 0) + 1) = 0.
!
! ********** Construction of exclusive channel cross sections **********
!
! 1. Loop over all residual nuclei, starting with the first residual nucleus, and then according to decreasing Z and N.
!
! Each idnum represents a different exclusive channel.
!
  idnum = -1
  Zend = min(numZchan, maxZ)
  Zend = min(Zend, Zinit)
  Nend = min(numNchan, maxN)
  Nend = min(Nend, Ninit)
  do Zix = 0, Zend
    do Nix = 0, Nend
!
! 2. To minimize the loops, the maximal possible number of each particle type is determined, given Zix and Nix.
!    E.g. for Zix=1, Nix=2 there can be at most one proton, 2 neutrons, one deuteron or 1 triton in the exit channel.
!    These maximal numbers are determined by simple formulae.
!
      Aix = Zix + Nix
      inend = 0
      ipend = 0
      idend = 0
      itend = 0
      ihend = 0
      iaend = 0
      if (Aix /= 0) then
        if (parinclude(1)) inend = min(Nix, maxchannel)
        if (parinclude(2)) ipend = min(Zix, maxchannel)
        if (parinclude(3)) idend = min(2 * Zix * Nix / Aix, maxchannel)
        if (parinclude(4)) itend = min(3 * Zix * Nix / (2 * Aix), maxchannel)
        if (parinclude(5)) ihend = min(3 * Zix * Nix / (2 * Aix), maxchannel)
        if (parinclude(6)) iaend = min(Zix * Nix / Aix, maxchannel)
      endif
      inend = min(inend, numin)
      ipend = min(ipend, numip)
      idend = min(idend, numid)
      itend = min(itend, numit)
      ihend = min(ihend, numih)
      iaend = min(iaend, numia)
!
! 3. Determine whether residual nucleus under consideration can be reached by particle combination.
!    Increase running index idnum and identifier ident for particle combination.
!
      do ih = 0, ihend
      do it = 0, itend
      do id = 0, idend
      do ia = 0, iaend
      do ip = 0, ipend
      do in = 0, inend
        if ( .not. chanopen(in, ip, id, it, ih, ia) .and. idnumfull) cycle
        if (idnum == numchantot) cycle
        Ztot = ip + id + it + 2 * ih + 2 * ia
        Ntot = in + id + 2 * it + ih + 2 * ia
        npart = in + ip + id + it + ih + ia
        if (npart > maxchannel) cycle
        if (Ztot /= Zix .or. Ntot /= Nix) cycle
        ident = 100000 * in + 10000 * ip + 1000 * id + 100 * it + 10 * ih + ia
        idnum = idnum + 1
        idchannel(idnum) = ident
!
! Initialization of arrays. Since the idnum counter may be reset at the end of the loop, in the case that the exclusive
! cross section is below a threshold, the various exclusive channel arrays need to be initialized to zero here.
!
        xschannel(idnum) = 0.
        yieldchannel(idnum) = 0.
        xsgamchannel(idnum) = 0.
        xsfischannel(idnum) = 0.
        xschancheck(idnum) = 0.
        xsfischancheck(idnum) = 0.
        xsratio(idnum) = 0.
        do i = 0, numlev
          Qexcl(idnum, i) = 0.
          Ethrexcl(idnum, i) = 0.
          xschaniso(idnum, i) = 0.
          exclbranch(idnum, i) = 0.
          do i2 = 0, numlev
            xsgamdischan(idnum, i, i2) = 0.
          enddo
        enddo
        Qexcl(0, 0) = S(0, 0, k0) + targetE
        do nexout = 0, numex + 1
          xsexcl(idnum, nexout) = 0.
          gamexcl(idnum, nexout) = 0.
        enddo
        if (flagspec) then
          do type = 0, 6
            Eavchannel(idnum, type) = 0.
            do nen = 0, numen
              xschannelsp(idnum, type, nen) = 0.
              xsfischannelsp(idnum, type, nen) = 0.
              do nexout = 0, numex + 1
                specexcl(idnum, type, nexout, nen) = 0.
              enddo
            enddo
          enddo
        endif
!
! 4. Determine source paths for exclusive channel.
!
! The exclusive channel under consideration may have been reached through different paths.
! Here the previous path is determined by looking at the type of the last emitted particle.
!
        do type = 0, 6
          if (parskip(type)) cycle
          identorg(type) = - 1
          if (type == 0) idd = ident
          if (type == 1) idd = ident - 100000
          if (type == 2) idd = ident - 10000
          if (type == 3) idd = ident - 1000
          if (type == 4) idd = ident - 100
          if (type == 5) idd = ident - 10
          if (type == 6) idd = ident - 1
          if (idd < 0) cycle
          do idorg = 0, idnum
            if (idchannel(idorg) == idd) then
              identorg(type) = idorg
!
! The Q-value for exclusive channels is determined, both for the ground state and isomers.
!
              Zcomp = Zix - parZ(type)
              Ncomp = Nix - parN(type)
              if (Qexcl(idnum, 0) == 0.) Qexcl(idnum, 0) = Qexcl(idorg, 0) - S(Zcomp, Ncomp, type)
              do nex = 0, Nlast(Zix, Nix, 0)
                Qexcl(idnum, nex) = Qexcl(idnum, 0) - edis(Zix, Nix, nex)
                if (Qexcl(idnum, nex) >= 0.) then
                  Ethrexcl(idnum, nex) = 0.d0
                else
                  Ethrexcl(idnum, nex) = - (Qexcl(idnum, nex) / specmass(parZ(k0), parN(k0), k0))
                endif
              enddo
              cycle
            endif
          enddo
        enddo
!
! Initialize fission channel.
!
        if (flagfission .and. Zix == 0 .and. Nix == 0) then
          nex = maxex(0, 0) + 1
          xsfischannel(idnum) = xsfischannel(idnum) + fisfeedex(0, 0, nex)
        endif
        do nexout = maxex(Zix, Nix), 0, - 1
          do type = 0, 6
            if (parskip(type)) cycle
            if (identorg(type) ==  -1) cycle
!
! Determine the compound nucleus belonging to the emitted particle type.
!
            idorg = identorg(type)
            Zcomp = Zix - parZ(type)
            Ncomp = Nix - parN(type)
!
! 5. Exclusive cross sections per excitation energy
!
! The inclusive cross section per excitation energy S is tracked in xsexcl.
! The inclusive spectrum per excitation energy dS/dEk' is tracked in specexcl.
!
! Special case for binary channel. The binary emission spectrum has already been determined in subroutine binemission.
!
            if (Zcomp == 0 .and. Ncomp == 0) then
              nex = maxex(Zcomp, Ncomp) + 1
              xsexcl(idnum, nexout) = xsexcl(idnum, nexout) + feedexcl(Zcomp, Ncomp, type, nex, nexout)
              if (type == 0) gamexcl(idnum, nexout) = gamexcl(idnum, nexout) + &
                feedexcl(Zcomp, Ncomp, 0, nex, nexout)
              if (flagspec) then
                do nen = ebegin(type), eend(type)
                  specexcl(idnum, type, nexout, nen) = specexcl(idnum, type, nexout, nen) + &
                    binemis(type, nexout, nen)
                enddo
              endif
            endif
!
! Calculation of inclusive cross section per excitation energy S
!
            do nex = maxex(Zcomp, Ncomp), 1, - 1
              if (feedexcl(Zcomp, Ncomp, type, nex, nexout) == 0.) cycle
              term1 = feedexcl(Zcomp, Ncomp, type, nex, nexout) / popexcl(Zcomp, Ncomp, nex)
              term2 = term1 * xsexcl(idorg, nex)
              xsexcl(idnum, nexout) = xsexcl(idnum, nexout) + term2
              if (type == 0) then
                gamexcl(idnum, nexout) = gamexcl(idnum, nexout) + term2
                if (nex <= Nlast(Zcomp, Ncomp, 0)) xsgamdischan(idnum, nex, nexout) = term2
              endif
                gamexcl(idnum, nexout) = gamexcl(idnum, nexout) + term1 * gamexcl(idorg, nex)
!   endif
!
! 6. Exclusive spectra per excitation energy
!
! Calculation of inclusive spectrum per excitation energy dS/dEk'
!
              if (flagspec) then
                do type2 = 0, 6
                  if (parskip(type2)) cycle
!
! Calculation of first term of exclusive spectrum.
!
                  do nen = ebegin(type2), eend(type2)
                    term3 = term1 * specexcl(idorg, type2, nex, nen)
                    specexcl(idnum, type2, nexout, nen) = specexcl(idnum, type2, nexout, nen) + term3
                  enddo
!
! Calculation of second term of exclusive spectrum.
! The spectrum of the last emitted particle is obtained by interpolating the feeding terms in subroutine specemission.
! The discrete gamma-ray production is excluded from the exclusive gamma spectra.
!
! specemission: subroutine for exclusive emission spectra
!
                  if (type2 == 0 .and. nex <= Nlast(Zcomp, Ncomp, 0)) cycle
                  if (type2 == type) then
                    call specemission(Zcomp, Ncomp, nex, idorg, type, nexout)
                    do nen = ebegin(type), eend(type)
                      specexcl(idnum, type, nexout, nen) = specexcl(idnum, type, nexout, nen) + &
                        specemis(nen)
                    enddo
                  endif
                enddo
              endif
            enddo
          enddo
!
! 7. Exclusive cross sections: total and per isomer
!
! At the end of the decay, the nucleus can end up at an isomer or the ground state.
! The exclusive cross sections for the particular channel can now be established.
! The same is true for the spectra.
! The spectra per isomer/ground state are always added.
!
          if (nexout > Nlast(Zix, Nix, 0)) cycle
          if (tau(Zix, Nix, nexout) /= 0..or.nexout == 0) then
            xschaniso(idnum, nexout) = xsexcl(idnum, nexout)
            xschannel(idnum) = xschannel(idnum) + xsexcl(idnum, nexout)
            xsgamchannel(idnum) = xsgamchannel(idnum) + gamexcl(idnum, nexout)
            if (flagspec) then
              do type = 0, 6
                if (parskip(type)) cycle
                do nen = ebegin(type), eend(type)
                  xschannelsp(idnum, type, nen) = xschannelsp(idnum, type, nen) + &
                    specexcl(idnum, type, nexout, nen)
                enddo
              enddo
            endif
          endif
        enddo
        if (xspopnuc(Zix, Nix) > 0.) xsratio(idnum) = xschannel(idnum) / xspopnuc(Zix, Nix)
        channelsum = channelsum + xschannel(idnum)
        if (in == 0) xsabs = xsabs + xschannel(idnum)
!
! For non-threshold reactions (positive Q-value) we always assign a minimum value to the exclusive cross section.
! (The transmission coefficients for these reactions might have been zero (from ECIS), but non-threshold reactions
! theoretically have a non-zero cross section.)
!
        if (Qexcl(idnum, 0) > 0..and.xschannel(idnum) <= xseps) then
          xschannel(idnum) = xseps
          xsgamchannel(idnum) = xseps
          xschaniso(idnum, 0) = xseps
          do i = 1, Nlast(Zix, Nix, 0)
            if (tau(Zix, Nix, i) /= 0.) then
              xschaniso(idnum, 0) = 0.5 * xseps
              xschaniso(idnum, i) = 0.5 * xseps
            endif
          enddo
        endif
!
! For each exclusive channel, the multiplicity or yield per isomer/ground state is determined.
!
        if (xschannel(idnum) /= 0.) then
          do i = 0, Nlast(Zix, Nix, 0)
            exclbranch(idnum, i) = xschaniso(idnum, i) / xschannel(idnum)
          enddo
        endif
!
! 8. Exclusive fission cross sections
!
! All exclusive fission cross sections and spectra are determined.
!
        if (flagfission) then
          do nex = maxex(Zix, Nix), 1, - 1
            if (popexcl(Zix, Nix, nex) /= 0.) then
              term = fisfeedex(Zix, Nix, nex) / popexcl(Zix, Nix, nex)
              xsfischannel(idnum) = xsfischannel(idnum) + term * xsexcl(idnum, nex)
              if (flagspec) then
                do type = 0, 6
                  if (parskip(type)) cycle
                  do nen = ebegin(type), eend(type)
                    xsfischannelsp(idnum, type, nen) = xsfischannelsp(idnum, type, nen) + &
                      term * specexcl(idnum, type, nex, nen)
                  enddo
                enddo
              endif
            endif
          enddo
          channelsum = channelsum + xsfischannel(idnum)
          xsabs = xsabs + xsfischannel(idnum)
        endif
!
! 9. Check.
!
! The summed exclusive cross sections should equal the particle production cross sections and analogously for the spectra.
! The result will appear in the output, if requested.
! Also the exclusive spectra are integrated so that they can be compared with the exclusive cross sections.
!
        if (flagcheck) then
          xsparcheck(1) = xsparcheck(1) + in * xschannel(idnum)
          xsparcheck(2) = xsparcheck(2) + ip * xschannel(idnum)
          xsparcheck(3) = xsparcheck(3) + id * xschannel(idnum)
          xsparcheck(4) = xsparcheck(4) + it * xschannel(idnum)
          xsparcheck(5) = xsparcheck(5) + ih * xschannel(idnum)
          xsparcheck(6) = xsparcheck(6) + ia * xschannel(idnum)
          if (flagspec) then
            fissum = 0.
            do type = 0, 6
              if (parskip(type)) cycle
              emissum = 0.
              Eaveragesum = 0.
              do nen = ebegin(type), eend(type)
                xsspeccheck(type, nen) = xsspeccheck(type, nen) + xschannelsp(idnum, type, nen)
                emissum = emissum + xschannelsp(idnum, type, nen) * deltaE(nen)
                Eaveragesum = Eaveragesum + egrid(nen) * xschannelsp(idnum, type, nen) * deltaE(nen)
                if (flagfission) fissum = fissum + xsfischannelsp(idnum, type, nen) * deltaE(nen)
              enddo
              if (type == 0) then
                do i = 1, numlev
                  do i2 = 0, i
                    if (xsgamdischan(idnum, i, i2) > 0.) then
                      emissum = emissum + xsgamdischan(idnum, i, i2)
                      Eaveragesum = Eaveragesum + xsgamdischan(idnum, i, i2) * &
                        (edis(Zcomp, Ncomp, i) - edis(Zcomp, Ncomp, i2))
                    endif
                  enddo
                enddo
              endif
              if (npart == 1) then
                NL = Nlast(Zix, Nix, 0)
                if (eoutdis(type, NL) > 0.) then
                  frac = Etop(nendisc(type)) - eoutdis(type, NL)
                  emissum = emissum - xschannelsp(idnum, type, nendisc(type)) * frac
                  Eaveragesum = Eaveragesum - egrid(nendisc(type)) * &
                    xschannelsp(idnum, type, nendisc(type)) * frac
                  if (flagfission) fissum = fissum - xsfischannelsp(idnum, type, nendisc(type)) * frac
                endif
              endif
              if (type > 0 .or. npart == 0) xschancheck(idnum) = xschancheck(idnum) + emissum
              if (emissum > 0.) Eavchannel(idnum, type) = Eaveragesum / emissum
            enddo
            if (xschannel(idnum) > 0.) then
              gmult(idnum) = xsgamchannel(idnum) / xschannel(idnum)
            else
              gmult(idnum) = 0.
            endif
            Especsum(idnum) = gmult(idnum) * Eavchannel(idnum, 0) + &
              in * Eavchannel(idnum, 1) + ip * Eavchannel(idnum, 2) + &
              id * Eavchannel(idnum, 3) + it * Eavchannel(idnum, 4) + &
              ih * Eavchannel(idnum, 5) + ia * Eavchannel(idnum, 6)
            if (flagfission) xsfischancheck(idnum) = fissum
          endif
        endif
!
! 10. Create reaction string for output
!
        apos = "'"
        reacstring(idnum) = '( ,            '
        reacstring(idnum)(2:2) = parsym(k0)
        i = 4
        if (npart == 0.) then
          reacstring(idnum)(4:4) = 'g'
          i = i + 1
        endif
        if (in == 1) then
          write(reacstring(idnum)(i:i), '(a1)') parsym(1)
          i = i + 1
        endif
        if (in > 1) then
          write(reacstring(idnum)(i:i+1), '(i1, a1)') in, parsym(1)
          i = i + 2
        endif
        if (ip == 1) then
          write(reacstring(idnum)(i:i), '(a1)') parsym(2)
          i = i + 1
        endif
        if (ip > 1) then
          write(reacstring(idnum)(i:i+1), '(i1, a1)') ip, parsym(2)
          i = i + 2
        endif
        if (id == 1) then
          write(reacstring(idnum)(i:i), '(a1)') parsym(3)
          i = i + 1
        endif
        if (id > 1) then
          write(reacstring(idnum)(i:i+1), '(i1, a1)') id, parsym(3)
          i = i + 2
        endif
        if (it == 1) then
          write(reacstring(idnum)(i:i), '(a1)') parsym(4)
          i = i + 1
        endif
        if (it > 1) then
          write(reacstring(idnum)(i:i+1), '(i1, a1)') it, parsym(4)
          i = i + 2
        endif
        if (ih == 1) then
          write(reacstring(idnum)(i:i), '(a1)') parsym(5)
          i = i + 1
        endif
        if (ih > 1) then
          write(reacstring(idnum)(i:i+1), '(i1, a1)') ih, parsym(5)
          i = i + 2
        endif
        if (ia == 1) then
          write(reacstring(idnum)(i:i), '(a1)') parsym(6)
          i = i + 1
        endif
        if (ia > 1) then
          write(reacstring(idnum)(i:i+1), '(i1, a1)') ia, parsym(6)
          i = i + 2
        endif
        if (npart == 1) then
          if (k0 == 1 .and. in == 1) then
            reacstring(idnum)(i:i) = apos
            i = i + 1
          endif
          if (k0 == 2 .and. ip == 1) then
            reacstring(idnum)(i:i) = apos
            i = i + 1
          endif
          if (k0 == 3 .and. id == 1) then
            reacstring(idnum)(i:i) = apos
            i = i + 1
          endif
          if (k0 == 4 .and. it == 1) then
            reacstring(idnum)(i:i) = apos
            i = i + 1
          endif
          if (k0 == 5 .and. ih == 1) then
            reacstring(idnum)(i:i) = apos
            i = i + 1
          endif
          if (k0 == 6 .and. ia == 1) then
            reacstring(idnum)(i:i) = apos
            i = i + 1
          endif
        endif
        reacstring(idnum)(i:i) = ')'
        if (flagfission) fisstring(idnum) = reacstring(idnum)(1:i-1)//'f)'
!
! Reset idnum counter in the case that the cross section is too small.
!
        if ((xschannel(idnum) >= xseps .and. .not. idnumfull) .or. npart == 0) then
          if ( .not. chanopen(in, ip, id, it, ih, ia)) opennum = opennum + 1
          chanopen(in, ip, id, it, ih, ia) = .true.
        endif
        if (xschannel(idnum) < xseps .and. npart > 1 .and. .not. chanopen(in, ip, id, it, ih, ia)) idnum = idnum - 1
        if (opennum == numchantot - 10) idnumfull = .true.
        if (idnum < 0) cycle
        if (xschannel(idnum) < 0.) xschannel(idnum) = xseps
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
    enddo
  enddo
!
! Calculate relative yield per channel
!
  if (channelsum > 0.) then
    do idc = 0, idnum
      yieldchannel(idc) = xschannel(idc) / channelsum
    enddo
  endif
!
! Set threshold energy for inelastic scattering to that of first excited state.
!
  do idc = 0, idnum
    if (k0 == 1 .and. idchannel(idc) == 100000) then
      Ethrexcl(idc, 0) = Ethrexcl(idc, 1)
      exit
    endif
  enddo
  return
end subroutine channels
! Copyright A.J. Koning 2021
