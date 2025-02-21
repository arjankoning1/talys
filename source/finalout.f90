    subroutine finalout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of final results
!
! Author    : Arjan Koning
!
! 2024-06-30: Original code
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
! Variables for basic reaction
!   flagchannels      ! flag for exclusive channels calculation
!   flagpartable      ! flag for output of model parameters on separate file
! Variables for numerics
!   maxchannel        ! maximal number of outgoing particles in individual channel description
!   maxN              ! maximal number of neutrons away from initial compound nucleus
!   maxZ              ! maximal number of protons away from initial compound nucleus
! Variables for output
!   flagexc           ! flag for output of excitation functions
!   flaggamdis        ! flag for output of discrete gamma - ray intensities
! Variables for compound reactions
!   flagcomp          ! flag for compound angular distribution calculation
!   xscaptherm        ! thermal capture cross section
! Variables for direct reactions
!   flagdisc          ! flag for output of discrete state cross sections
! Variables for input energies
!   Ninc            ! number of incident energies
! Variables for main input
!   Atarget           ! mass number of target nucleus
!   Ltarget           ! excited level of target
!   flagnatural       ! flag for calculation of natural element
!   k0                ! index of incident particle
!   Ztarget           ! charge number of target nucleus
! Variables for preequilibrium
!   Cbreak            ! adjustable parameter for break - up reactions
!   Cknock            ! adjustable parameter for knockout reactions
!   Cstrip            ! adjustable parameter for stripping / pick - up reactions
!   Kph               ! constant for single - particle level density par. (g = A / Kph)
!   M2constant        ! constant for matrix element in exciton model
!   M2limit           ! constant for asymptotic value for matrix element
!   M2shift           ! constant for energy shift for matrix element
!   phmodel           ! particle - hole state density model
!   Rgamma            ! adjustable parameter for pre - equilibrium gamma decay
!   Rnunu             ! ratio for two - component matrix element
!   Rnupi             ! ratio for two - component matrix element
!   Rpinu             ! ratio for two - component matrix element
!   Rpipi             ! ratio for two - component matrix element
!   Rspincutpreeq     ! adjustable constant (global) for preequilibrium spin cutoff factor
! Variables for level density
!   alev              ! level density parameter
!   alphald           ! alpha - constant for asymptotic level density parameter
!   betald            ! beta - constant for asymptotic level density parameter
!   cfermi            ! width of Fermi distribution for damping
!   cglobal           ! global constant to adjust tabulated level densities
!   gammashell1       ! gamma - constant for asymptotic level density parameter
!   gammashell2       ! gamma - constant for asymptotic level density parameter
!   ldmodel           ! level density model
!   pair              ! pairing energy
!   pairconstant      ! constant for pairing energy systematics
!   pglobal           ! global constant to adjust tabulated level densities
!   Pshiftconstant    ! global constant for pairing shift
!   Rspincut          ! adjustable constant (global) for spin cutoff factor
!   Ufermi            ! energy of Fermi distribution for damping
! Variables for OMP
!   adepthcor         ! adjustable parameter for depth of DF alpha potential
!   alphaomp          ! alpha optical model
!   aradialcor        ! adjustable parameter for shape of DF alpha potential
!   avadjust          ! adjustable factor for OMP (default 1.)
!   avdadjust         ! adjustable factor for OMP (default 1.)
!   avsoadjust        ! adjustable factor for OMP (default 1.)
!   awadjust          ! adjustable factor for OMP (default 1.)
!   awdadjust         ! adjustable factor for OMP (default 1.)
!   d1adjust          ! adjustable factor for OMP (default 1.)
!   d2adjust          ! adjustable factor for OMP (default 1.)
!   d3adjust          ! adjustable factor for OMP (default 1.)
!   flagjlm           ! flag for using semi - microscopic JLM OMP
!   flagomponly       ! flag to execute ONLY an optical model calculation
!   lv1adjust         ! adjustable parameter for JLM OMP
!   lvadjust          ! adjustable parameter for JLM OMP
!   lvsoadjust        ! adjustable parameter for JLM OMP
!   lw1adjust         ! adjustable parameter for JLM OMP
!   lwadjust          ! adjustable parameter for JLM OMP
!   lwsoadjust        ! adjustable parameter for JLM OMP
!   rcadjust          ! adjustable factor for OMP (default 1.)
!   rvadjust          ! adjustable factor for OMP (default 1.)
!   rvdadjust         ! adjustable factor for OMP (default 1.)
!   rvsoadjust        ! adjustable factor for OMP (default 1.)
!   rwadjust          ! adjustable factor for OMP (default 1.)
!   rwdadjust         ! adjustable factor for OMP (default 1.)
!   v1adjust          ! adjustable factor for OMP (default 1.)
!   v2adjust          ! adjustable factor for OMP (default 1.)
!   v3adjust          ! adjustable factor for OMP (default 1.)
!   v4adjust          ! adjustable factor for OMP (default 1.)
!   vso1adjust        ! adjustable factor for OMP (default 1.)
!   vso2adjust        ! adjustable factor for OMP (default 1.)
!   w1adjust          ! adjustable factor for OMP (default 1.)
!   w2adjust          ! adjustable factor for OMP (default 1.)
!   w3adjust          ! adjustable factor for OMP (default 1.)
!   w4adjust          ! adjustable factor for OMP (default 1.)
!   wso1adjust        ! adjustable factor for OMP (default 1.)
!   wso2adjust        ! adjustable factor for OMP (default 1.)
! Variables for fission
!   Cnubar1           ! adjustable parameter for nubar constant value
!   Cnubar2           ! adjustable parameter for nubar energy slope
!   Tmadjust          ! adjustable parameter for PFNS temperature
!   Fsadjust          ! adjustable parameter for PFNS scission fraction
!   flagfission       ! flag for fission
!   Rspincutff        ! adjustable parameter (global) for FF spin cutoff factor
! Variables for existence libraries
!   chanopen          ! flag to open channel with first non - zero cross section
!   fisexist          ! flag for existence of fission cross section
!   gamexist          ! flag for existence of gamma production cross section
!   rpexist           ! flag for existence of residual production cross section
!   rpisoexist        ! flag for existence of isomeric residual production cross section
! Variables for OMP
!   Rprime            ! potential scattering radius
!   Sstrength         ! s, p, d, etc - wave strength function
! Variables for incident channel
!   maxA              ! maximal number of nucleons away from initial compound nucleus
! Variables for energies
!   idchannel         ! identifier for exclusive channel
!   Ethresh           ! threshold incident energy for residual nucleus
!   Ethrexcl          ! threshold incident energy for exclusive channel
!   Qexcl             ! Q - value for exclusive channel
!   Qres              ! Q - value for residual nucleus
!   reacstring        ! string for exclusive reaction channel
! Variables for exclusive channels
!   idnum             ! counter for exclusive channel
! Variables for nuclides
!   AA                ! mass number of residual nucleus
!   Nindex            ! neutron number index for residual nucleus
!   parinclude        ! logical to include outgoing particle
!   parskip           ! logical to skip outgoing particle
!   Zindex            ! charge number index for residual nucleus
!   ZZ                ! charge number of residual nucleus
! Constants
!   cparity           ! parity (character)
!   nuc               ! symbol of nucleus
!   parname           ! name of particle
!   parsym            ! symbol of particle
! Variables for resonance parameters
!   D0theo            ! mean s - wave resonance spacing
! Variables for levels
!   edis              ! energy of level
!   jdis              ! spin of level
!   parlev            ! parity of level
!   tau               ! lifetime of state in seconds
! Variables for level density
!   Nlast             ! last discrete level
! Variables for masses
!   S                 ! separation energy
! Variables for preequilibrium
!   Esurf          ! well depth for surface interaction
!
! *** Declaration of local data
!
  implicit none
  character(len=1)   :: yesno                 ! y or n function
  character(len=1)   :: ejec                  ! symbol of ejectile 
  character(len=6)   :: discfile              ! file with elastic scattering angular distribution
  character(len=6)   :: contfile              ! file for continuum
  character(len=9)   :: totfile               ! file with total cross sections
  character(len=9)   :: reactstring(0:6)      ! reaction string
  character(len=10)  :: binfile               ! file for binary output
  character(len=11)  :: fisfile               ! fission file
  character(len=12)  :: rpfile                ! file with residual production cross sections
  character(len=12)  :: isofile               ! file with isomeric cross section
  character(len=12)  :: xsfile                ! file with channel cross sections
  character(len=19)  :: gamfile               ! giant resonance parameter file
  integer            :: A                     ! mass number of target nucleus
  integer            :: Acomp                 ! mass number index for compound nucleus
  integer            :: i                     ! counter
  integer            :: i1                    ! value
  integer            :: i2                    ! value
  integer            :: ia                    ! mass number from abundance table
  integer            :: id                    ! counter for deuterons
  integer            :: idc                   ! help variable
  integer            :: ident                 ! exclusive channel identifier
  integer            :: ih                    ! hole number
  integer            :: in                    ! counter for neutrons
  integer            :: ip                    ! particle number
  integer            :: it                    ! counter for tritons
  integer            :: Ncomp                 ! neutron number index for compound nucleus
  integer            :: nex                   ! excitation energy bin of compound nucleus
  integer            :: Nix                   ! neutron number index for residual nucleus
  integer            :: npart                 ! number of particles in outgoing channel
  integer            :: type                  ! particle type
  integer            :: Z                     ! charge number of target nucleus
  integer            :: Zcomp                 ! proton number index for compound nucleus
  integer            :: Zix                   ! charge number index for residual nucleus
!
! Write model parameters to separate file
!
  if (flagpartable) then
    write(51, '("## General parameters:")')
    write(51, '("##   level density")')
    write(51, '("ldmodel        ", i1)') ldmodel(0, 0)
    write(51, '("colenhance     ", a1)') yesno(flagcolall)
    if (ldmodel(0, 0) <= 3) then
      write(51, '("alphald        ", f10.5)') alphald(0, 0)
      write(51, '("betald         ", f10.5)') betald(0, 0)
      write(51, '("gammashell1    ", f10.5)') gammashell1(0, 0)
      write(51, '("gammashell2    ", f10.5)') gammashell2
      write(51, '("pairconstant   ", f10.5)') pairconstant
      write(51, '("pshiftconstant ", f10.5)') Pshiftconstant(0, 0)
      write(51, '("Rspincut       ", f10.5)') Rspincut
    endif
    write(51, '("Rspincutff     ", f10.5)') Rspincutff
    write(51, '("cglobal        ", es12.5)') cglobal
    write(51, '("pglobal        ", es12.5)') pglobal
    if (phmodel == 1) write(51, '("Kph            ", f10.5)') Kph
    write(51, '("## General parameters:")')
    write(51, '("##   photon strength function")')
    write(51, '("strength       ", i2)') strength
    write(51, '("strengthM1     ", i2)') strengthM1
    write(51, '("xscaptherm     ", es12.5)') xscaptherm(-1)
    write(51, '("## General parameters:")')
    write(51, '("##   pre-equilibrium")')
    write(51, '("M2constant     ", f10.5)') M2constant
    write(51, '("M2limit        ", f10.5)') M2limit
    write(51, '("M2shift        ", f10.5)') M2shift
    write(51, '("Rpipi          ", f10.5)') Rpipi
    write(51, '("Rnunu          ", f10.5)') Rnunu
    write(51, '("Rpinu          ", f10.5)') Rpinu
    write(51, '("Rnupi          ", f10.5)') Rnupi
    write(51, '("Rgamma         ", f10.5)') Rgamma
    write(51, '("Rspincutpreeq  ", f10.5)') Rspincutpreeq
    write(51, '("Esurf          ", f10.5)') Esurf
    do type = 1, 6
      write(51, '("Cstrip        ", a1, f10.5)') parsym(type), Cstrip(type)
      write(51, '("Cknock        ", a1, f10.5)') parsym(type), Cknock(type)
      write(51, '("Cbreak        ", a1, f10.5)') parsym(type), Cbreak(type)
    enddo
    write(51, '("## General parameters:")')
    write(51, '("##   fission")')
    write(51, '("fismodel       ", i1)') fismodel
    write(51, '("Cnubar1        ", f10.5)') Cnubar1
    write(51, '("Cnubar2        ", f10.5)') Cnubar2
    write(51, '("Tmadjust       ", f10.5)') Tmadjust
    write(51, '("Fsadjust       ", f10.5)') Fsadjust
    write(51, '("Cbarrier       ", f10.5)') Cbarrier
    write(51, '("## General parameters:")')
    write(51, '("##   optical model")')
    if (flagjlm) then
      write(51, '("lvadjust       ", f10.5)') lvadjust
      write(51, '("lwadjust       ", f10.5)') lwadjust
      write(51, '("lv1adjust      ", f10.5)') lv1adjust
      write(51, '("lw1adjust      ", f10.5)') lw1adjust
      write(51, '("lvsoadjust     ", f10.5)') lvsoadjust
      write(51, '("lwsoadjust     ", f10.5)') lwsoadjust
    endif
    write(51, '("deuteronomp    ", i1)') deuteronomp
    write(51, '("alphaomp       ", i1)') alphaomp
    if (alphaomp >= 3 .and. alphaomp <= 5) then
      write(51, '("aradialcor     ", f10.5)') aradialcor
      write(51, '("adepthcor      ", f10.5)') adepthcor
    endif
    do i = 1, 7
      if (i == 7) then
        type = 0
        ejec = 'e'
      else
        type = i
        ejec = parsym(type)
      endif
      if ( .not. flagjlm .or. type > 2) then
        write(51, '("v1adjust      ", a1, f10.5)') ejec, v1adjust(type)
        write(51, '("v2adjust      ", a1, f10.5)') ejec, v2adjust(type)
        write(51, '("v3adjust      ", a1, f10.5)') ejec, v3adjust(type)
        write(51, '("v4adjust      ", a1, f10.5)') ejec, v4adjust(type)
        write(51, '("rvadjust      ", a1, f10.5)') ejec, rvadjust(type)
        write(51, '("avadjust      ", a1, f10.5)') ejec, avadjust(type)
        write(51, '("w1adjust      ", a1, f10.5)') ejec, w1adjust(type)
        write(51, '("w2adjust      ", a1, f10.5)') ejec, w2adjust(type)
        write(51, '("w3adjust      ", a1, f10.5)') ejec, w3adjust(type)
        write(51, '("w4adjust      ", a1, f10.5)') ejec, w4adjust(type)
        write(51, '("rwadjust      ", a1, f10.5)') ejec, rwadjust(type)
        write(51, '("awadjust      ", a1, f10.5)') ejec, awadjust(type)
        write(51, '("rvdadjust     ", a1, f10.5)') ejec, rvdadjust(type)
        write(51, '("avdadjust     ", a1, f10.5)') ejec, avdadjust(type)
        write(51, '("d1adjust      ", a1, f10.5)') ejec, d1adjust(type)
        write(51, '("d2adjust      ", a1, f10.5)') ejec, d2adjust(type)
        write(51, '("d3adjust      ", a1, f10.5)') ejec, d3adjust(type)
        write(51, '("rwdadjust     ", a1, f10.5)') ejec, rwdadjust(type)
        write(51, '("awdadjust     ", a1, f10.5)') ejec, awdadjust(type)
        write(51, '("rvsoadjust    ", a1, f10.5)') ejec, rvsoadjust(type)
        write(51, '("avsoadjust    ", a1, f10.5)') ejec, avsoadjust(type)
        write(51, '("vso1adjust    ", a1, f10.5)') ejec, vso1adjust(type)
        write(51, '("vso2adjust    ", a1, f10.5)') ejec, vso2adjust(type)
        write(51, '("rwsoadjust    ", a1, f10.5)') ejec, rwsoadjust(type)
        write(51, '("awsoadjust    ", a1, f10.5)') ejec, awsoadjust(type)
        write(51, '("wso1adjust    ", a1, f10.5)') ejec, wso1adjust(type)
        write(51, '("wso2adjust    ", a1, f10.5)') ejec, wso2adjust(type)
        write(51, '("rcadjust      ", a1, f10.5)') ejec, rcadjust(type)
      endif
    enddo
    if (k0 == 1 .and. (parinclude(0) .or. flagcomp) .and. Rprime /= 0.) then
      write(51, '("## General parameters:")')
      write(51, '("##   resonance parameters")')
      write(51, '("## Z   A     S0        R      xs(therm)    D0", "         a         P        Sn")')
      write(51, '("##", 2i4, 7es10.3)') Ztarget, Atarget, Sstrength(0) * 1.e4, Rprime, xscaptherm(-1), D0theo(0, 0), &
 &      alev(0, 0), pair(0, 0), S(0, 0, 1)
    endif
  endif
  close (unit = 51)
!
! ****************** Integrated binary cross sections ******************
!
! flagexc    : flag for output of excitation functions
!
  if (flagomponly) return
  if (flagnatural .or. .not. flagexc) return
  write(*, '(/" ########## EXCITATION FUNCTIONS ###########"/)')
  write(*, '(" 1. Total (binary) cross sections"/)')
  totfile = 'all.tot'
  call write_outfile(totfile,flagoutall)
  if (.not. flagoutall) then
    write(*,'("file: total.tot")')
    write(*,'("file: elastic.tot")')
    write(*,'("file: nonelastic.tot")')
    write(*,'("file: reaction.tot")')
  endif
  write(*, '(/" 2. Binary non-elastic cross sections (non-exclusive)"/)')
  binfile = 'binary.tot'
  call write_outfile(binfile,flagoutall)
!
! ************** Total particle production cross sections **************
!
  write(*, '(/" 3. Total particle production cross sections"/)')
  do type = 0, 6
    if (parskip(type)) cycle
    totfile = ' prod.tot'
    write(totfile(1:1), '(a1)') parsym(type)
    call write_outfile(totfile,flagoutall)
  enddo
  if (flagfission) then
    write(*, '(/" 3b. Total fission cross sections "/)')
    fisfile = 'fission.tot'
    call write_outfile(fisfile,flagoutall)
  endif
!
! ******************** Residual production cross sections **************
!
  write(*, '(/" 4. Residual production cross sections"/)')
  do Acomp = 0, maxA
    do Zcomp = 0, maxZ
      Ncomp = Acomp - Zcomp
      if (Ncomp < 0 .or. Ncomp > maxN) cycle
      if ( .not. rpexist(Zcomp, Ncomp)) cycle
      Z = ZZ(Zcomp, Ncomp, 0)
      A = AA(Zcomp, Ncomp, 0)
      rpfile = 'rp000000.tot'
      write(rpfile(3:8), '(2i3.3)') Z, A
      call write_outfile(rpfile,flagoutall)
      do nex = 0, Nlast(Zcomp, Ncomp, 0)
        if ( .not. rpisoexist(Zcomp, Ncomp, nex)) cycle
        isofile = 'rp000000.L00'
        write(isofile(3:8), '(2i3.3)') Z, A
        write(isofile(11:12), '(i2.2)') levnum(Zcomp, Ncomp, nex)
        call write_outfile(isofile,flagoutall)
      enddo
    enddo
  enddo
!
! ********************** Fission cross sections ************************
!
  if (flagfission) then
    write(*, '(/" 4b. Fission cross sections per fissioning nuclide"/)')
    do Acomp = 0, maxA
      do Zcomp = 0, maxZ
        Ncomp = Acomp - Zcomp
        if (Ncomp < 0 .or. Ncomp > maxN) cycle
        if ( .not. fisexist(Zcomp, Ncomp)) cycle
        Z = ZZ(Zcomp, Ncomp, 0)
        A = AA(Zcomp, Ncomp, 0)
        rpfile = 'rp000000.fis'
        write(rpfile(3:8), '(2i3.3)') Z, A
        call write_outfile(rpfile,flagoutall)
      enddo
    enddo
  endif
!
! ******************** Reactions to discrete states ********************
!
  if (flagdisc) then
    do type = 0, 6
      if (type == k0) then
        reactstring(type) = 'Inelastic'
      else
        reactstring(type) = '  ( , )  '
        write(reactstring(type)(4:4), '(a1)') parsym(k0)
        write(reactstring(type)(6:6), '(a1)') parsym(type)
      endif
    enddo
    write(*, '(/" 5. Binary reactions to discrete levels and continuum"/)')
    Zcomp = 0
    Ncomp = 0
    do type = 0, 6
      if (parskip(type)) cycle
      Zix = Zindex(Zcomp, Ncomp, type)
      Nix = Nindex(Zcomp, Ncomp, type)
      Z = ZZ(Zcomp, Ncomp, type)
      A = AA(Zcomp, Ncomp, type)
      if (flagchannels) then
        contfile = '  .tot'
        write(contfile(1:2), '(2a1)') parsym(k0), parsym(type)
        call write_outfile(contfile,flagoutall)
      endif
      do nex = 0, Nlast(Zix, Nix, 0)
        if (type == k0 .and. nex == Ltarget) cycle
        discfile = '  .L00'
        write(discfile(1:2), '(2a1)') parsym(k0), parsym(type)
        write(discfile(5:6), '(i2.2)') nex
        call write_outfile(discfile,flagoutall)
      enddo
      if (flagchannels) then
        contfile = '  .con'
        write(contfile(1:2), '(2a1)') parsym(k0), parsym(type)
        call write_outfile(contfile,flagoutall)
      endif
    enddo
  endif
!
! ******************** Exclusive channels cross sections ***************
!
! 1. Exclusive cross sections
!
  if (flagchannels) then
    write(*, '(/" 6. Exclusive cross sections"/)')
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
                    if (idchannel(idc) == ident) then
                      xsfile = 'xs000000.tot'
                      write(xsfile(3:8), '(6i1)') in, ip, id, it, ih, ia
                      call write_outfile(xsfile,flagoutall)
                      Zcomp = ip + id + it + 2 * ih + 2 * ia
                      Ncomp = in + id + 2 * it + ih + 2 * ia
                      isofile = 'xs000000.tot'
                      write(isofile(3:8), '(6i1)') in, ip, id, it, ih, ia
                      do nex = 0, Nlast(Zcomp, Ncomp, 0)
                        write(isofile(11:12), '(i2.2)') levnum(Zcomp, Ncomp, nex)
                        call write_outfile(isofile,flagoutall)
                      enddo
                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  endif
!
! 2. Exclusive fission cross sections
!
    if (flagfission) then
      write(*, '(/" 6b. Exclusive fission cross sections"/)')
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
                      if (idchannel(idc) == ident) then
                        xsfile = 'xs000000.fis'
                        write(xsfile(3:8), '(6i1)') in, ip, id, it, ih, ia
                        call write_outfile(xsfile,flagoutall)
                      endif
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    endif
!
! ************************* Gamma-ray intensities **********************
!
  if (flaggamdis) then
    write(*, '(/" 7. Gamma-ray intensities"/)')
    do Zcomp = 0, maxZ
      do Ncomp = 0, maxN
        Z = ZZ(Zcomp, Ncomp, 0)
        A = AA(Zcomp, Ncomp, 0)
        do i1 = 0, numlev
          do i2 = 0, i1
            if ( .not. gamexist(Zcomp, Ncomp, i1, i2)) cycle
            gamfile = 'gam000000L00L00.tot'
            write(gamfile(4:9), '(2i3.3)') Z, A
            write(gamfile(11:12), '(i2.2)') i1
            write(gamfile(14:15), '(i2.2)') i2
            call write_outfile(gamfile,flagoutall)
          enddo
        enddo
      enddo
    enddo
  endif
  return
end subroutine finalout
! Copyright A.J. Koning 2021
