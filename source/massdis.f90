subroutine massdis
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Fission fragment yields
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
!   sgl             ! single precision kind
! All global variables
!   numelem         ! number of elements
!   numenin         ! number of incident energies
!   numJ            ! maximum J - value
!   numin           ! maximum number of neutrons in channel description
!   numN            ! maximum number of neutrons from initial compound nucleus
!   numneu          ! number of neutrons
!   numnu           ! number of neutrons from fission
!   numpop          ! number of population bins
!   numZ            ! maximum number of protons from initial compound nucleus
! Variables for fission
!   flagffspin      ! flag to use spin distribution in initial FF population
!   flagoutfy       ! flag for output detailed fission yield calculation
!   fymodel         ! fission yield model, 1: Brosa 2: GEF
!   ffmodel         ! fission fragment model, 1: GEF 2: HF3D (Okumura) 3: SPY 0: user
!   gefran          ! number of random events for GEF calculation
!   Rfiseps         ! ratio for limit for fission cross section per nucleus
!   yieldfile       ! fission yield file
! Variables for numerics
!   maxN            ! maximal number of neutrons away from initial compound nucleus
!   maxZ            ! maximal number of protons away from initial compound nucleus
! Variables for main input
!   Ainit           ! mass number of initial compound nucleus
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
!   Ztarget         ! charge number of target nucleus
! Constants
!   isochar         ! symbol of isomer
!   nuc             ! symbol of nucleus
!   parA            ! mass number of particle
!   parZ            ! charge number of particle
! Variables for total cross sections
!   xsfistot        ! total fission cross section
! Variables for energies
!   Etotal          ! total energy of compound system (target + projectile)
! Variables for excitation energy grid
!   Ex              ! excitation energy
!   fisfeedJP       ! fission contribution from excitation energy bin per J, P
!   maxex           ! maximum excitation energy bin for residual nucleus
!   maxJ            ! maximal J - value
! Variables for multiple emission
!   fisfeedex       ! fission contribution from excitation energy bin
!   xsfeed          ! cross section from compound to residual nucleus
! Variables for incident channel
!   xsbinary        ! cross section from initial compound to residual nucleus
! Variables for nuclides
!   AA              ! mass number of residual nucleus
!   Nindex          ! neutron number index for residual nucleus
!   parskip         ! logical to skip outgoing particle
!   Q               ! Q - value
!   Zindex          ! charge number index for residual nucleus
!   ZZ              ! charge number of residual nucleus
! Variables for files
!   path            ! directory containing files to be read
! Variables for mass distribution
!   disa            ! normalised fission fragment mass yield per excitation energy bin
!   disacor         ! normalised fission product mass yield per excitation energy bin
!   disaz           ! normalised fission fragment isotope yield per excitation energy
!   disazcor        ! normalised fission product isotope yield per excitation energy b
!   excfis          ! excitation energy at fission
!   EaverageA       ! average emission energy per A
!   EaverageZA      ! average emission energy per (Z,A)
!   Excff           ! excitation energy of fission fragment
!   dExcff          ! width of excitation energy of fission fragment
!   nubar           ! average nu
!   fpeps           ! ratio for limit for fission product cross section
!   fpratio         ! fission product isomeric ratio
!   nuA             ! nu per A
!   nuZA            ! nu per Z,A
!   Pdisnu          ! prompt fission neutrons distribution
!   TKE             ! total kinetic energy
!   xsApost         ! post - neutron emission corrected cross section
!   xsApre          ! pre - neutron emission cross section
!   xsfpex          ! excitation energy spectrum per fission fragment
!   xstotpost       ! post - neutron emission fission product cross section
!   xstotpre        ! pre - neutron emission fission product cross section
!   xsZApost        ! post - neutron emission corrected isotopic cross section
!   xsZApre         ! pre - neutron emission isotopic cross section
!   yieldApost      ! post - neutron emission corrected fission yield
!   yieldApre       ! pre - neutron emission fission yield
!   yieldfpex       ! fission yield per isomer
!   yieldtotpost    ! post - neutron emission fission product yield
!   yieldtotpre     ! pre - neutron emission fission product yield
!   yieldZApost     ! post - neutron emission corrected isotopic yield
!   yieldZApre      ! pre - neutron emission isotopic yield
!
! *** Declaration of local data
!
  implicit none
  include "gef.cmb"
  integer, parameter :: numZff=80                                   ! number of Z of fission fragments
  integer, parameter :: numNff=150                                  ! number of N of fission fragments
  integer, parameter :: numpair=2000                                !
  character(len=13)  :: fffile                                      ! fission fragment file
  character(len=132) :: gefpath                                     ! path for GEF files
  character(len=132) :: string                                      ! string
  character(len=2)   :: Sfile                                       !
  character(len=3)   :: massstring                                  ! string for mass number
  character(len=6)   :: nucstring                                   !
  character(len=8)   :: Estring                                     !
  character(len=132) :: Efile                                       !
  character(len=132) :: ffname                                      !
  character(len=132) :: ffpath                                      !
  character(len=132) :: Yfile(2)                                    ! file with production yields
  logical            :: lexist                                      ! logical to determine existence
  real(sgl)          :: beldm1(136, 203)                            ! binding energy from liquid drop model
  real(sgl)          :: Ebin(0:numpop)                              ! energy of bin
  real(sgl)          :: Etabtot(numZff, numNff, 1000)               ! tabulated energy
  real(sgl)          :: Exfis(1000)                                 ! excitation energy for fission
  real(sgl)          :: dEH                                         !
  real(sgl)          :: dEL                                         !
  real(sgl)          :: Efac                                        ! help variable
  real(sgl)          :: Eff(0:numenin)                              !
  real(sgl)          :: Efftab(2)                                   !
  real(sgl)          :: Extab(2)                                    !
  real(sgl)          :: Effrel                                      !
  real(sgl)          :: EH                                          !
  real(sgl)          :: ELL
  real(sgl)          :: Eheavy(2, numpair)                          !
  real(sgl)          :: Elight(2, numpair)                          !
  real(sgl)          :: fisepsA                                     ! fission tolerance
  real(sgl)          :: fisepsB                                     ! fission tolerance
  real(sgl)          :: Fmulti                                      ! factor for multi-chance fission
  real(sgl)          :: Jfis                                        ! spin of fissioning system
  real(sgl)          :: Jtabtot(numZff, numNff, 100)                ! total spin from GEF
  real(sgl)          :: partfisJ(0:numJ)                            ! partial fission spin distribution
  real(sgl)          :: partfisxs                                   ! partial fission cross section
  real(sgl)          :: popffEx(numZff, numNff, 0:numpop)           ! energy population of FF
  real(sgl)          :: popffJ(numZff, numNff, 0:numJ)              ! spin population of FF
  real(sgl)          :: sum                                         ! help variable
  real(sgl)          :: sumE                                        ! summed sensitivity x cross section
  real(sgl)          :: sumJ                                        ! sum over spin distribution
  real(sgl)          :: sumpost                                     ! sum over post-neutron FP's
  real(sgl)          :: sumpre                                      ! sum over pre-neutron FP's
  real(sgl)          :: sumxs                                       ! sum over emission channels
  real(sgl)          :: term                                        ! help variable
  real(sgl)          :: TKE0(2, numpair)                            !
  real(sgl)          :: TXE0(2, numpair)                            !
  real(sgl)          :: TK                                          !
  real(sgl)          :: TX                                          !
  real(sgl)          :: Wheavy(2, numpair)                          !
  real(sgl)          :: Wlight(2, numpair)                          !
  real(sgl)          :: Y(2, numpair)                               ! product yield (in ENDF-6 format)
  real(sgl)          :: Y0                                          !
  real(sgl)          :: ushell1(136, 203)                           ! shell correction
  real(sgl)          :: xsfis(1000)                                 ! fission cross section
  real(sgl)          :: xsfisFF                                     ! fission cross section per FF
  real(sgl)          :: xstabcomp(0:numZ, 0:numN, numZff, numNff)   ! Z, N cross section from GEF
  real(sgl)          :: xstabtot(numZff, numNff)                    ! total cross section from GEF
  real(sgl)          :: Ytabtot(numZff, numNff)                     ! yield from GEF
  integer            :: A                                           ! mass number of target nucleus
  integer            :: Afile                                       !
  integer            :: Aheavy(2, numpair)                          ! heaviest isotope per element
  integer            :: Alight(2, numpair)                          ! lightest isotope per element
  integer            :: gefwrite                                    ! integer for output detailed fission yield calculation
  integer            :: i                                           ! counter
  integer            :: istat                                       ! logical for file access
  integer            :: ia                                          ! mass number from abundance table
  integer            :: in                                          ! counter for neutrons
  integer            :: iskip                                       ! help variable
  integer            :: istep                                       ! help variable
  integer            :: iz                                          ! charge number of residual nucleus
  integer            :: inH                                         !
  integer            :: inL                                         !
  integer            :: iaH                                         !
  integer            :: iaL                                         !
  integer            :: iza                                         ! counter for Z,A combinations
  integer            :: izH                                         !
  integer            :: izL                                         !
  integer            :: j                                           ! counter
  integer            :: Jgef                                        ! counter
  integer            :: k                                           ! counter
  integer            :: nb                                          ! help variable
  integer            :: Ncomp                                       ! neutron number index for compound nucleus
  integer            :: nen                                         ! energy counter
  integer            :: nex                                         ! excitation energy bin of compound nucleus
  integer            :: nexend                                      ! last energy index
  integer            :: nexgef                                      ! counter
  integer            :: Nff                                         !
  integer            :: Nfftab                                      !
  integer            :: Ntotal(2)                                   !
  integer            :: Nix                                         ! neutron number index for residual nucleus
  integer            :: odd                                         ! odd (1) or even (0) nucleus
  integer            :: parity                                      ! parity
  integer            :: type                                        ! particle type
  integer            :: Z                                           ! charge number of target nucleus
  integer            :: Zcomp                                       ! proton number index for compound nucleus
  integer            :: Zheavy(2, numpair)                          !
  integer            :: Zlight(2, numpair)                          !
  integer            :: Zix                                         ! charge number index for residual nucleus
!
! ************************** Mass yields *******************************
!
! Initialization
!
  xsApre = 0.
  xsApost = 0.
  yieldApre = 0.
  yieldApost = 0.
  nuA = 0.
  nuZA = 0.
  EaverageZA = 0.
  EaverageA = 0.
  xsZApre = 0.
  xsZApost = 0.
  yieldZApre = 0.
  yieldZApost = 0.
  xsfpex = 0.
  yieldfpex = 0.
  fpratio = 0.
  Excff = 0.
  dExcff = 0.
  TKE = 0.
  Ntotal = 0
  Zlight = 0
  Alight = 0
  Zheavy = 0
  Aheavy = 0
  Elight = 0.
  Wlight = 0.
  Eheavy = 0.
  Wheavy = 0.
  Y = 0.
  TKE0 = 0.
  TXE0 = 0.
  yieldtotpre = 0.
  xstotpre = 0.
  yieldtotpost = 0.
  xstotpost = 0.
  Pdisnu = 0.
  nubar = 0.
  Eff = 0.
  xstabtot = 0.
  Etabtot = 0.
  Jtabtot = 0.
  popffEx = 0.
  popffJ = 0.
  xstabcomp = 0.
  fpeps = Rfiseps * xsfistot
  if (fpeps == 0.) return
  if (fymodel == 2 .or. fymodel ==3) then
    gefpath = trim(path)//'fission/gef/'
    if (flagoutfy) then
      gefwrite = 1
    else
      gefwrite = 0
    endif
!
! Read nuclear structure information for GEF
!
    open (unit = 4, file = trim(gefpath)//'beldm.dat', status = 'old')
    read(4, * ) beldm1
    close(4)
    do i = 1, 203
      do j = 1, 136
       beldm(i, j) = beldm1(j, i)
      end do
    end do
    open (unit = 4, file = trim(gefpath)//'ushell.dat', status = 'old')
    read(4, * ) ushell1
    close(4)
    do i = 1, 203
      do j = 1, 136
       ushel(i, j) = ushell1(j, i)
      end do
    end do
    open (unit = 4, file = trim(gefpath)//'nucprop.dat', status = 'old')
    do i = 1, 3885
      read(4, * ) (RNucTab(i, j), j = 1, 8)
    end do
    close(4)
  endif
!
! Loop over nuclides
!
! brosafy    : subroutine for fission fragment yields based on Brosa
!
  do Zcomp = 0, maxZ
    do Ncomp = 0, maxN
      Z = ZZ(Zcomp, Ncomp, 0)
      A = AA(Zcomp, Ncomp, 0)
      odd = mod(A, 2)
      Zix = Zindex(Zcomp, Ncomp, 0)
      Nix = Nindex(Zcomp, Ncomp, 0)
      if (xsfeed(Zcomp, Ncomp, - 1) <= fpeps) cycle
      if (Zcomp == 0 .and. Ncomp == 0) then
        nexend = maxex(Zcomp, Ncomp) + 1
      else
        nexend = maxex(Zcomp, Ncomp)
      endif
      fisepsA = fpeps / max(3 * maxex(Zcomp, Ncomp), 1)
      iskip = 0
      istep = 4
      if (fymodel == 2) then
        Exfis = 0.
        xsfis = 0.
        nen = 0
      endif
      do nex = nexend, 0, - 1
        if (Zcomp == 0 .and. Ncomp == 0 .and. nex == maxex(Zcomp, Ncomp) + 1) then
          excfis = Etotal
          partfisxs = xsbinary( - 1)
          do J = 0, numJ
            partfisJ(J) = 0.
            do parity = - 1, 1, 2
              partfisJ(J) = partfisJ(J) + fisfeedJP(Zcomp, Ncomp, nex, J, parity)
            enddo
          enddo
        else
          if (mod(iskip, istep) /= 0) then
            iskip = iskip + 1
            cycle
          endif
          if (nex - istep + 1 < 0) cycle
          if (Ex(Zcomp, Ncomp, nex - istep + 1) >= 30.) then
            partfisxs = 0.
            do J = 0, numJ
              partfisJ(J) = 0.
            enddo
            do i = 0, istep - 1
              partfisxs = partfisxs + fisfeedex(Zcomp, Ncomp, nex - i)
              do J = 0, numJ
                do parity = - 1, 1, 2
                  partfisJ(J) = partfisJ(J) + fisfeedJP(Zcomp, Ncomp, nex - i, J, parity)
                enddo
              enddo
            enddo
            if (partfisxs /= 0) then
              excfis = 0.
              do i = 0, istep - 1
                excfis = excfis + fisfeedex(Zcomp, Ncomp, nex - i) * Ex(Zcomp, Ncomp, nex - i)
              enddo
              excfis = excfis / partfisxs
            endif
            iskip = 1
          else
            excfis = Ex(Zcomp, Ncomp, nex)
            partfisxs = fisfeedex(Zcomp, Ncomp, nex)
            do J = 0, numJ
              partfisJ(J) = 0.
              do parity = - 1, 1, 2
                partfisJ(J) = partfisJ(J) + fisfeedJP(Zcomp, Ncomp, nex, J, parity)
              enddo
            enddo
          endif
        endif
        if (partfisxs > fisepsA) then
!
! Brosa
!
          if (fymodel == 1) then
            call brosafy(Zix, Nix)
            term = 0.5 * partfisxs
            do ia = 1, A
              xsApre(ia) = xsApre(ia) + term * disa(ia)
              xsApost(ia) = xsApost(ia) + term * disacor(ia)
              do iz = 1, Z
                in = ia - iz
                if (in < 1 .or. in > numneu) cycle
                xsZApre(iz, in) = xsZApre(iz, in) + term * disaz(ia, iz)
                xsZApost(iz, in) = xsZApost(iz, in) + term * disazcor(ia, iz)
              enddo
            enddo
          endif
!
! GEF
!
          if (fymodel == 2) then
            nen = nen + 1
            Exfis(nen) = excfis
            xsfis(nen) = partfisxs
          endif
!
! GEF + TALYS evaporation
!
          if (fymodel == 3 .and. A <= 350) then
            fisepsB = fisepsA / (5 * maxJ(Zcomp, Ncomp, nex)) * 0.5
            do J = 0, maxJ(Zcomp, Ncomp, nex)
              if (partfisJ(J) < fisepsB) cycle
              Jfis = real(J) + 0.5 * odd
              call gefsub(Z, A, excfis, Jfis)
              write( * , * ) " FF excitation for Z=", Z, " A=", A, " Ex= ", excfis, " J=", Jfis, " xs=", &
 &              partfisJ(J), " N_cases:", N_cases
              do iza = 1, N_cases
                iz = NZMkey(iza, 3)
                in = NZMkey(iza, 2)
                if (iz > numZff .or. in > numNff) cycle
                term = Ytab(iza) * partfisJ(J)
                xstabtot(iz, in) = xstabtot(iz, in) + term
                xstabcomp(Zcomp, Ncomp, iz, in) = xstabcomp(Zcomp, Ncomp, iz, in) + term
                do nexgef = 1, 1000
                  Etabtot(iz, in, nexgef) = Etabtot(iz, in, nexgef) + term * Etab(iza, nexgef)
                enddo
                do Jgef = 1, 100
                  Jtabtot(iz, in, Jgef) = Jtabtot(iz, in, Jgef) + term * Jtab(iza, Jgef)
                enddo
              enddo
            enddo
          endif
        endif
      enddo
!
! GEF
!
      if (fymodel == 2 .and. A <= 350) then
        call geftalys(real(Z), real(A), nen, Exfis, xsfis, gefwrite, gefran)
        do ia = 1, A
          xsApre(ia) = xsApre(ia) + 0.5 * ysum(ia)
          xsApost(ia) = xsApost(ia) + 0.5 * ysump(ia)
          if (ia <= 200) then
            do iz = 1, Z
              in = ia - iz
              if (in >= 1 .and. in <= numneu) then
                xsZApre(iz, in) = xsZApre(iz, in) + 0.5 * yAZ(ia, iz)
                yieldZApre(iz, in) = yAZ(ia, iz)
                xsZApost(iz, in) = xsZApost(iz, in) + 0.5 * yAZp(ia, iz)
              endif
            enddo
          endif
        enddo
        if (xsfistot > 0.) then
          Fmulti = Ncomp * xsfeed(Zcomp, Ncomp, - 1)
          do i = 1, numnu
            if (ann_sum(i) > 0.) Pdisnu(1, i) = Pdisnu(1, i) + (Fmulti + ann_sum(i)) / xsfistot
          enddo
          do ia = 1, A
            if (anpre_sum(ia) > 0.) nuA(1, ia) = nuA(1, ia) + (Fmulti + anpre_sum(ia)) / xsfistot
          enddo
          nubar(1) = nubar(1) + (Fmulti + anMean_sum) / xsfistot
        endif
      endif
    enddo
  enddo
  do type = 0, 6
    if (parskip(type)) cycle
    sum = 0.
    do i = 1, numin
      sum = sum + Pdisnu(type, i)
    enddo
    if (sum > 0.) then
      do i = 1, numin
        Pdisnu(type, i) = Pdisnu(type, i)/ sum
      enddo
    endif
  enddo
!
! GEF + TALYS evaporation (fymodel 3)
!  or
! yields + TALYS evaporation (fymodel 4)
!
  if (fymodel >= 3) then
    Ebin(0) = 0.
    do i = 1, numpop
      Ebin(i) = 0.1 * i
    enddo
    if (fymodel == 3) then
      sumxs = 0.
      do iz = 1, numZff
        do in = 1, numNff
          if (iz > numelem .or. in > numneu) cycle
          if (xstabtot(iz, in) == 0.) cycle
          sumxs = sumxs + xstabtot(iz, in)
          sumE = 0.
          do nexgef = 1, 1000
            sumE = sumE + Etabtot(iz, in, nexgef)
          enddo
          if (sumE > 0.) then
            do nexgef = 1, 1000
              Etabtot(iz, in, nexgef) = Etabtot(iz, in, nexgef) / sumE
            enddo
          endif
          sumJ = 0.
          do Jgef = 1, 100
            sumJ = sumJ + Jtabtot(iz, in, Jgef)
          enddo
          if (sumJ > 0.) then
            do Jgef = 1, 100
              Jtabtot(iz, in, Jgef) = Jtabtot(iz, in, Jgef) / sumJ
            enddo
          endif
        enddo
      enddo
    endif
    if (sumxs > 0.) then
      do iz = 1, numZff
        do in = 1, numNff
          if (iz > numelem .or. in > numneu) cycle
          Ytabtot(iz, in) = xstabtot(iz, in) / sumxs
        enddo
      enddo
    endif
    do iz = 1, numZff
      do in = 1, numNff
        if (iz > numelem .or. in > numneu) cycle
        ia = iz + in
        xsfisFF = xsfistot * Ytabtot(iz, in)
        xsApre(ia) = xsApre(ia) + xsfisFF
        xsZApre(iz, in) = xsZApre(iz, in) + xsfisFF
        yieldZApre(iz, in) = Ytabtot(iz, in)
        do nexgef = 1, 1000
          popffEx(iz, in, nexgef) = popffEx(iz, in, nexgef) + xsfisFF * Etabtot(iz, in, nexgef)
        enddo
        do J = 1, 30
          popffJ(iz, in, J) = Jtabtot(iz, in, J)
        enddo
      enddo
    enddo
    do iz = 1, Z
      do ia = iz + 1, A
        in = ia - iz
        if (in > numneu) cycle
        if (xsZApre(iz, in) == 0.) cycle
        sumJ = 0.
        do J = 1, 30
          sumJ = sumJ + popffJ(iz, in, J)
        enddo
        if (sumJ > 0.) then
          do J = 1, 30
            popffJ(iz, in, J) = popffJ(iz, in, J) / sumJ
          enddo
        endif
        nb = numpop
        do nex = 100, numpop
          if (popffEx(iz, in, nex) < 1.e-9) then
            nb = nex
            exit
          endif
        enddo
        nb = max(nb, 1)
        fffile = 'ff000000.ex'
        write(fffile(3:5), '(i3.3)') iz
        write(fffile(6:8), '(i3.3)') ia
        open (unit = 1, file = fffile, status = 'replace')
        if (flagffspin) then
          write(1, * ) nb, 30, 1, " Yff = ", yieldZApre(iz, in), " xs = ", xsZApre(iz, in)
          do nex = 0, nb
            write(1, '(f10.5, 30es12.5)') Ebin(nex), (0.5 * popffEx(iz, in, nex) * popffJ(iz, in, J), J = 1, 30)
          enddo
        else
          write(1, * ) nb, 0, 1, " Yff = ", yieldZApre(iz, in), " xs = ", xsZApre(iz, in)
          do nex = 0, nb
            write(1, '(f10.5, es12.5)') Ebin(nex), popffEx(iz, in, nex)
          enddo
        endif
        close(1)
      enddo
    enddo
!
! fymodel 4: Okumura model - read in yields and excitation energies
!
!            1: GEF (database by Ali Al-Adili and Fredrik Nordstroem)
!            2: HF3D (Okumura) (database by Toshihiko Kawano)
!            3: SPY (Okumura) (Jean-Francois Lemaitre)
!            4: Langevin-4D (Titech)
!
    if (fymodel == 4) then
      if (yieldfile(1:1) == ' ') then
        Effrel = Etotal
        ffpath = trim(path)//'fission/ff/'
        if (ffmodel == 0) ffname = 'user'
        if (ffmodel == 1) ffname = 'gef'
        if (ffmodel == 2) ffname = 'hf3d'
        if (ffmodel == 3) ffname = 'spy'
        if (ffmodel == 4) ffname = 'langevin4d'
        ffpath = trim(ffpath)//trim(ffname)//'/'
        massstring = '   '
        Afile = Ainit
        Sfile = nuc(Zinit)
        write(massstring, '(i3.3)') Afile
        nucstring = trim(Sfile) //massstring 
        ffpath = trim(ffpath)//trim(nucstring)//'/'
        Efile = trim(ffpath)//trim(nucstring)//'_'//trim(ffname)//'.E'
        inquire (file=Efile,exist=lexist)
        if (.not.lexist) then
          write(*,'(" TALYS-error: Non-existent FF file ",a)') trim(Efile)
          stop
        endif
        open (unit = 1, file = Efile, status = 'unknown')
        i = 1
        do
          read(1, * , iostat = istat) Eff(i)
          if (istat == -1) exit
          i = i + 1
        enddo
        close (1)
        Nff = i - 1
        if (Effrel <= Eff(1)) then
          Efftab(1) = Eff(1)
          Nfftab = 1
        else
         Nfftab = 2
         if (Effrel >= Eff(Nff)) then
            Efftab(1) = Eff(Nff)
            Efftab(2) = Emaxtalys
            Extab(2) = Efftab(1)
          else
            call locate(Eff, 0, Nff, Effrel, nen)
            Efftab(1) = Eff(nen)
            Efftab(2) = Eff(nen + 1)
            Extab(2) = Eff(nen + 1)
          endif
          Efac = (Effrel - Efftab(1)) / (Efftab(2) - Efftab(1))
        endif
        Extab(1)=Efftab(1)
        do k = 1, Nfftab
          Estring = '        '
          write(Estring, '(es8.2)') Extab(k)
          Estring(5:5) = 'e'
          Yfile(k) = trim(ffpath)//trim(nucstring)//'_'//Estring// 'MeV_'//trim(ffname)//'.ff'
        enddo
      else
        Nfftab = 1
        Yfile(1) = trim(yieldfile)
      endif
      do k = 1, Nfftab
        write(*, '(/, " Fission fragment yields read from ", a)') trim(Yfile(k))
        inquire (file=Yfile(k),exist=lexist)
        if (.not.lexist) then
          write(*,'(" TALYS-error: Non-existent FF file ",a)') trim(Yfile(k))
          stop
        endif
        open (unit = 1, file = Yfile(k), status = 'unknown')
        read(1, '(///13x, i6)') Ntotal(k)
        read(1, '(a)') string
        do i = 1, Ntotal(k)
          read(1, '(a)', iostat = istat) string
          if (istat == -1) exit
          read(string, * ) Zlight(k, i), Alight(k, i), Zheavy(k, i), Aheavy(k, i), Y(k, i), &
 &            TKE0(k, i), TXE0(k, i), Elight(k, i), Wlight(k, i), Eheavy(k, i), Wheavy(k, i)
        enddo
        close (1)
      enddo
      write(*, '(/, "  Fission fragment pairs")')
      write(*, '("  Zl  Al  Zh  Ah   Yield        TKE         TXE        ELight     dElight     Eheavy      dEheavy")')
      do i = 1, Ntotal(Nfftab)
        izL = Zlight(Nfftab, i)
        iaL = Alight(Nfftab, i)
        izH = Zheavy(Nfftab, i)
        iaH = Aheavy(Nfftab, i)
        inL = iaL - izL
        inH = iaH - izH
        ia = iaL + iaH
        if (inL < 1 .or. inL > numneu) cycle
        if (inH < 1 .or. inH > numneu) cycle
        if (ia /= Ainit) cycle
        ELL = Elight(Nfftab, i)
        dEL = Wlight(Nfftab, i)
        EH = Eheavy(Nfftab, i)
        dEH = Wheavy(Nfftab, i)
        Y0 = Y(Nfftab, i)
        TK = TKE0(Nfftab, i)
        TX = TXE0(Nfftab, i)
        if (Nfftab == 2) then
          do j = 1, Ntotal(1)
            if (Zlight(1, j) == izL .and. Alight(1, j) == iaL) then
              ELL = Elight(1, j) + Efac * (Elight(2, i) - Elight(1, j))
              dEL = Wlight(1, j) + Efac * (Wlight(2, i) - Wlight(1, j))
              Y0 = Y(1, j) + Efac * (Y(2, i) - Y(1, j))
              TK = TKE0(1, j) + Efac * (TKE0(2, i) - TKE0(1, j))
              TX = TXE0(1, j) + Efac * (TXE0(2, i) - TXE0(1, j))
            endif
            if (Zheavy(1, j) == izH .and. Aheavy(1, j) == iaH) then
              EH = Eheavy(1, j) + Efac * (Eheavy(2, i) - Eheavy(1, j))
              dEH = Wheavy(1, j) + Efac * (Wheavy(2, i) - Wheavy(1, j))
              Y0 = Y(1, j) + Efac * (Y(2, i) - Y(1, j))
              TK = TKE0(1, j) + Efac * (TKE0(2, i) - TKE0(1, j))
              TX = TXE0(1, j) + Efac * (TXE0(2, i) - TXE0(1, j))
            endif
          enddo
        endif
        xsfisFF = xsfistot * Y0
        xsApre(iaL) = xsApre(iaL) + xsfisFF
        xsApre(iaH) = xsApre(iaH) + xsfisFF
        xsZApre(izL, inL) = xsZApre(izL, inL) + xsfisFF
        xsZApre(izH, inH) = xsZApre(izH, inH) + xsfisFF
        Excff(izL, inL) = ELL
        Excff(izH, inH) = EH
        dExcff(izL, inL) = dEL
        dExcff(izH, inH) = dEH
        TKE(izL, inL) = TK
        TKE(izH, inH) = TK
        write(*, '(4i4, 7es12.5)') izL, iaL, izH, iaH, Y0, TK, TX, ELL, dEL, EH, dEH
      enddo
    endif
  endif
!
! fymodel 5: General model - read in full population per fission fragment
!
  if (fymodel == 5) then
    fffile = 'ff000000.ex'
    do iz = 1, Z
      do ia = iz + 1, A
        write(fffile(3:5), '(i3.3)') iz
        write(fffile(6:8), '(i3.3)') ia
        inquire (file = fffile , exist = lexist)
        if (lexist) then
          open (unit = 1, file = fffile, status = 'old')
          in = ia - iz
          if (in <= numneu) then
            read(1, * ) nb
            do nex = 0, nb
              read(1, '(f10.5, es12.5)') Ebin(nex), popffEx(iz, in, nex)
            enddo
          endif
          close (1)
        endif
      enddo
    enddo
  endif
!
! Normalization to fission yields (sum = 2)
!
  sumpre = 0.
  sumpost = 0.
  do ia = 1, Atarget
    sumpre = sumpre + xsApre(ia)
    sumpost = sumpost + xsApost(ia)
  enddo
  do iz = 1, Ztarget
    do ia = iz + 1, Atarget
      in = ia - iz
      if (iz > numelem .or. in > numneu) cycle
      if (xsZApre(iz, in) == 0.) cycle
      if (sumpre > 0.) yieldZApre(iz, in) = 2. * xsZApre(iz, in) / sumpre
      yieldApre(ia) = yieldApre(ia) + yieldZApre(iz, in)
      yieldtotpre = yieldtotpre + yieldZApre(iz, in)
      xstotpre = xstotpre + xsZApre(iz, in)
      if (fymodel <= 2) then
        if (sumpost > 0.) yieldZApost(iz, in) = 2. * xsZApost(iz, in) / sumpost
        yieldApost(ia) = yieldApost(ia) + yieldZApost(iz, in)
        yieldtotpost = yieldtotpost + yieldZApost(iz, in)
        xstotpost = xstotpost + xsZApost(iz, in)
      endif
    enddo
  enddo
  return
end subroutine massdis
! Copyright A.J. Koning 2021
