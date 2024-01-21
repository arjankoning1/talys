subroutine ffevap
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Evaporation of fission fragments
!
! Author    : Arjan Koning and Jean-Francois Lemaitre
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
!   numen2          ! maximum number of outgoing energies
!   numia           ! maximum number of alphas in channel description
!   numid           ! maximum number of deuterons in channel description
!   numih           ! maximum number of helions in channel description
!   numin           ! maximum number of neutrons in channel description
!   numip           ! maximum number of protons in channel description
!   numit           ! maximum number of tritons in channel description
!   numneu          ! number of neutrons
!   numnu           ! number of neutrons from fission
!   numpfns         ! number of energies for PFNS grid
!   numpop          ! number of population bins
!   numpar          ! number of particles
! Variables for output
!   flagmain        ! flag for main output
!   flagspec        ! flag for output of spectra
! Variables for compound reactions
!   flagcomp        ! flag for compound angular distribution calculation
! Variables for preequilibrium
!   epreeq          ! on - set incident energy for preequilibrium calculation
! Variables for numerics
!   maxN            ! maximal number of neutrons away from initial compound nucleus
!   maxZ            ! maximal number of protons away from initial compound nucleus
! Variables for fission
!   fymodel         ! fission yield model, 1: Brosa 2: GEF
!   Rfiseps         ! ratio for limit for fission cross section per nucleus
! Variables for gamma rays
!   flagracap       ! flag for radiative capture model
! Variables for input energies
!   enincmax        ! maximum incident energy
! Variables for main input
!   Atarget0        ! mass number of target nucleus
!   k0              ! index of incident particle
!   Ztarget0        ! charge number of target nucleus
! Variables for basic reaction
!   flagffruns      ! flag to designate subsequent evaporation of fission products
! Variables for OMP
!   flagompall      ! flag for new optical model calculation for all residual
! Variables for multiple emission
!   xsinitpop       ! initial population cross section
! Variables for total cross sections
!   xsfistot        ! total fission cross section
! Variables for incident channel
!   multiplicity    ! particle multiplicity
!   xspopex         ! population cross section summed over spin and parity
!   xspopnuc        ! population cross section per nucleus
! Variables for energies
!   idchannel       ! identifier for exclusive channel
! Variables for exclusive channels
!   idnum           ! counter for exclusive channel
!   xschannel       ! channel cross section
! Variables for spectra
!   Eaverage        ! average outgoing energy
!   espec           ! outgoing energy grid
!   xssumout        ! cross section summed over mechanisms
! Variables for nuclides
!   parinclude      ! logical to include outgoing particle
!   parskip         ! logical to skip outgoing particle
! Constants
!   onethird        ! 1 / 3
!   parA            ! mass number of particle
!   parZ            ! charge number of particle
!   sqrttwopi       ! sqrt(2. * pi)
! Variables for levels
!   tau             ! lifetime of state in seconds
! Variables for level density
!   Nlast           ! last discrete level
! Variables for existence libraries
!   fpexist         ! flag for existence of fission product
! Variables for mass distribution
!   Aff             ! mass number of fission fragment
!   dExcff          ! width of excitation energy of fission fragment
!   Excff           ! excitation energy of fission fragment
!   Eavpfns         ! average energy of prompt fission neutrons spectrum
!   EaverageA       ! average emission energy per A
!   EaverageZA      ! average emission energy per (Z,A)
!   Epfns           ! energy of PFNS
!   Epfnsaverage    ! average energy of PFNS
!   fpratio         ! fission product isomeric ratio
!   maxpfns         ! maximum energy of prompt fission neutrons spectrum
!   NEpfns          ! number of energies for PFNS
!   nuA             ! nu per A
!   nuZA            ! nu per Z,A
!   nubar           ! average nu
!   Pdisnu          ! prompt fission neutrons distribution
!   Pdisnuav        ! average prompt fission neutrons distribution
!   pfns            ! prompt fission neutron spectrum
!   pfnscm          ! prompt fission neutron spectrum in CM
!   TKE             ! total kinetic energy
!   xsApost         ! post - neutron emission corrected cross section
!   xsfpex          ! excitation energy spectrum per fission fragment
!   xstotpost       ! post - neutron emission fission product cross section
!   xsZApost        ! post - neutron emission corrected isotopic cross section
!   xsZApre         ! pre - neutron emission isotopic cross section
!   yieldApost      ! post - neutron emission corrected fission yield
!   yieldApre       ! pre - neutron emission fission yield
!   yieldfpex       ! fission yield per isomer
!   yieldtotpost    ! post - neutron emission fission product yield
!   yieldZApost     ! post - neutron emission corrected isotopic yield
!   yieldZApre      ! pre - neutron emission isotopic yield
!   Zff             ! charge number of fission fragment
!
! *** Declaration of local data
!
  implicit none
  character (len=8)  :: ffstring                     ! FF string
  character (len=14) :: fffile                       ! FF file
  integer   :: ACN                                   ! compound nucleus mass
  integer   :: i                                     ! counter
  integer   :: ia                                    ! mass number from abundance table
  integer   :: iaa                                   ! counter
  integer   :: iap                                   ! mass number
  integer   :: id                                    ! counter for deuterons
  integer   :: idc                                   ! help variable
  integer   :: ident                                 ! exclusive channel identifier
  integer   :: ih                                    ! hole number
  integer   :: in                                    ! counter for neutrons
  integer   :: inH                                   ! neutron number
  integer   :: izH                                   ! charge number
  integer   :: iaH                                   ! mass number
  integer   :: inn                                   ! neutron number
  integer   :: inp                                   ! unit for input
  integer   :: ip                                    ! particle number
  integer   :: it                                    ! counter for tritons
  integer   :: iz                                    ! charge number of residual nucleus
  integer   :: izp                                   ! chanrge number
  integer   :: nen                                   ! energy counter
  integer   :: nen2                                  ! energy counter
  integer   :: nex                                   ! excitation energy bin of compound nucleus
  integer   :: nex0                                  ! base number for discrete level
  integer   :: Nix                                   ! neutron number index for residual nucleus
  integer   :: npar                                  ! counter
  integer   :: type                                  ! particle type
  integer   :: ZCN                                   ! compound nucleus charge
  integer   :: Zix                                   ! charge number index for residual nucleus
  real(sgl) :: dE                                    ! help variable
  real(sgl) :: dEex                                  ! help variable
  real(sgl) :: E                                     ! incident energy
  real(sgl) :: Eav                                   ! average energy
  real(sgl) :: Eb                                    ! begin energy
  real(sgl) :: Ecm1                                  ! C.M. energy
  real(sgl) :: Ecm2                                  ! C.M. energy
  real(sgl) :: Ee                                    ! end energy
  real(sgl) :: Eex                                   ! excitation energy
  real(sgl) :: Ekinff                                ! kinetic energy of F.F.
  real(sgl) :: Ekintot                               ! total kinetic energy
  real(sgl) :: Emaxff                                ! maximum energy for population
  real(sgl) :: Esumpfns                              ! integrated PFNS
  real(sgl) :: fac                                   ! factor
  real(sgl) :: fac1                                  ! factor
  real(sgl) :: fac2                                  ! factor
  real(sgl) :: fiseps                                ! limit for fission cross section per excitation energy bin
  real(sgl) :: gau(0:numpop)                         !
  real(sgl) :: gauss                                 ! Gaussian contribution
  real(sgl) :: maxwell                               ! Maxwell distribution
  real(sgl) :: Pmultiff(numelem, numneu, 0:numpar, 0:numnu) !
  real(sgl) :: spec                                  ! spectrum
  real(sgl) :: speccm                                ! spectrum in CM
  real(sgl) :: sqrtE                                 ! square root of energy
  real(sgl) :: sqrtEkinff                            ! square root of kinetic energy of FF's
  real(sgl) :: sum                                   ! help variable
  real(sgl) :: summax                                ! integral over Maxwellian
  real(sgl) :: sumpfns                               ! integrated PFNS
  real(sgl) :: sumpfnscm                             ! integrated PFNS in CM
  real(sgl) :: sumpost                               ! sum over post-neutron FP's
  real(sgl) :: xs1                                   ! help variable
  real(sgl) :: xs2                                   ! help variable
  real(sgl) :: xsb                                   ! begin cross section
  real(sgl) :: xse                                   ! end cross section
  real(sgl) :: xsc                                   ! interpolated cross section
  real(sgl) :: xsexcpart(0:numpar, 0:numin)          ! partial cross section
  real(sgl) :: term                                  ! integral over Maxwellian
  real(sgl) :: yA                                    ! pre-neutron emission mass yield
  real(sgl) :: yZA                                   ! pre-neutron emission isotopic yield
!
! ********************** Loop over fission fragments *******************
!
! Do a full TALYS calculation for each fission fragment and incident energy
!
  flagffruns = .true.
  Epfnsaverage = 0.
  pfns = 0.
  pfnscm = 0.
  maxpfns = 0.
  fiseps = Rfiseps * xsfistot
  Pmultiff = 0.
  write(*, '(/" ########## Start of loop over fission fragments"/)')
!
! Use Viola systematics for first order guess of kinetic energy of FF
!
  ZCN = Ztarget0 + parZ(k0)
  ACN = Atarget0 + parA(k0)
  Ekintot = 0.1189 * Ztarget0 * Ztarget0 / ((Atarget0 + 1) **onethird) + 7.3
  do ia = 1, ACN
    do iz = 1, ZCN
      in = ia - iz
      if (in < 1 .or. in > numneu) cycle
      if (xsZApre(iz, in) <= fiseps) cycle
      Ekinff = real(ACN - ia) / real(ia) * Ekintot / ACN
!
! Okumura model
!
      ffstring = 'ff000000'
      write(ffstring(3:5), '(i3.3)') iz
      write(ffstring(6:8), '(i3.3)') ia
      if (fymodel == 4) then
        fffile=ffstring//'.ex'
        open (unit = 1, file = fffile, status = 'replace')
        write(1, '(i6,2i4,f12.5,2es12.5,f12.5,a)') numpop, 0, 1, Einc+S(0,0,1)+Q(1), yieldZApre(iz,in), xsZApre(iz,in), &
 &        TKE(iz,in)," (numpop numJ numparity Ex Yff xsff Ekin)"
        Eex = Excff(iz, in)
        dEex = dExcff(iz, in)
        Emaxff = min(Eex + 4. * dEex, 80.)
        fac1 = 1. / (dEex * sqrttwopi)
        fac2 = 1. / (2. * dEex **2)
        dE = Emaxff / numpop
        sum = 0.
        do nex = 0, numpop
          E = dE * nex
          gau(nex) = fac1 * exp( - ((E - Eex) **2) * fac2)
          sum = sum + dE * gau(nex)
        enddo
        do nex = 0, numpop
          E = dE * nex
          gauss = gau(nex) / sum
          write(1, '(f10.5, es12.5)') E, dE*gauss*xsZApre(iz, in)
        enddo
        close(1)
        Ekinff = real(ACN - ia) / real(ia) * TKE(iz, in) / ACN
      endif
      sqrtEkinff = sqrt(Ekinff)
      Aff = ia
      Zff = iz
      call evaptalys
!
! Add fission product cross sections
!
      do Zix = 0, maxZ
        do Nix = 0, maxN
          izp = iz - Zix
          if (izp < 1) cycle
          inp = in - Nix
          if (inp < 1) cycle
          iap = izp + inp
          xsZApost(izp, inp) = xsZApost(izp, inp) + xspopnuc(Zix, Nix)
          xsApost(iap) = xsApost(iap) + xspopnuc(Zix, Nix)
          do nex = 0, Nlast(Zix, Nix, 0)
            if (nex == 0 .or. tau(Zix, Nix, nex) /= 0.) then
              nex0 = min(nex, 1)
              xsfpex(izp, inp, nex0) = xsfpex(izp, inp, nex0) + xspopex(Zix, Nix, nex)
            endif
          enddo
        enddo
      enddo
!
! Add prompt fission particle and gamma production and spectra
!
      if (xsinitpop > 0.) then
        do npar = 0, numin
          do type = 0, 6
            xsexcpart(type, npar) = 0.
          enddo
          do iaa = 0, numia
            do ih = 0, numih
              do it = 0, numit
                do id = 0, numid
                  do ip = 0, numip
                    do inn = 0, numin
                      if (inn + ip + id + it + ih + iaa /= npar) cycle
                      ident = 100000 * inn + 10000 * ip + 1000 * id + 100 * it + 10 * ih + iaa
                      do idc = 0, idnum
                        if (idchannel(idc) == ident) then
                          xsc = xschannel(idc)
                          xsexcpart(1, inn) = xsexcpart(1, inn) + xsc
                          xsexcpart(2, ip) = xsexcpart(2, ip) + xsc
                          xsexcpart(3, id) = xsexcpart(3, id) + xsc
                          xsexcpart(4, it) = xsexcpart(4, it) + xsc
                          xsexcpart(5, ih) = xsexcpart(5, ih) + xsc
                          xsexcpart(6, iaa) = xsexcpart(6, iaa) + xsc
                        endif
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
        yA = yieldApre(ia)
        yZA = yieldZApre(iz, in)
        do type = 0, 6
          if (parskip(type)) cycle
          nubar(type) = nubar(type) + yZA * multiplicity(type)
          if (yA > 0.) then
            nuA(type, ia) = nuA(type, ia) + yZA / yA * multiplicity(type)
            nuZA(type, iz, in) = nuZA(type, iz, in) + yZA * multiplicity(type)
            do npar = 0, numin
              Pmultiff(iz, in, type, npar) = xsexcpart(type, npar) / xsinitpop
            enddo
            term = yZa * ACN / (ACN - ia)
            Epfnsaverage(type) = Epfnsaverage(type) + 0.5 * term * Eaverage(type)
            EaverageZA(type, iz, in) = Eaverage(type)
            EaverageA(type, ia) = EaverageA(type, ia) + yZA / yA * Eaverage(type)
            fffile = ffstring // '.' // parsym(type) // 'spec'
            open (unit=1, file=fffile, status='unknown')
            write(1, '("# Prompt fission ",a8," spectrum of ", i3, a2)') parname(type), ia, nuc(iz)
            write(1, '("# Number of energies:", i6)') NEpfns+1
            write(1, '("#    E-out spectrum_CM spectrum_lab")')
            do nen = 1, NEpfns
              E = Epfns(nen)
!
! C.M.
!
              speccm = 0.
              do nen2 = 0, eend(type) - 1
                Eb = espec(type, nen2)
                Ee = espec(type, nen2 + 1)
                xsb = xssumout(type, nen2)
                xse = xssumout(type, nen2 + 1)
                if (Eb == 0. .and. Ee == 0.) cycle
                if (E >= Eb .and. E <= Ee) then
                  call pol1(Eb, Ee, xsb, xse, E, speccm)
                  exit
                endif
              enddo
!
! Lab
!
              sum = 0.
              sqrtE = sqrt(E)
              Eb = (sqrtE - sqrtEkinff) **2
              Ee = (sqrtE + sqrtEkinff) **2
              fac = 1. / (4. * sqrtEkinff)
              if (abs(Ee - Eb) >= 1.e-5) then
                do nen2 = 1, eend(type) - 1
                  Ecm1 = espec(type, nen2)
                  Ecm2 = espec(type, nen2+1)
                  if (Ecm1 <= Eb .and. Ecm2 <= Eb) cycle
                  if (Ecm1 >= Ee .and. Ecm2 >= Ee) cycle
                  xs1 = 0.
                  xs2 = 0.
                  if (Ecm1 > 0.) xs1 = xssumout(type,nen2) / sqrt(Ecm1)
                  if (Ecm2 > 0.) xs2 = xssumout(type,nen2+1) / sqrt(Ecm2)
                  term = 0.5 * (xs1 + xs2)
                  if (Ecm1 < Eb.and.Ecm2 >= Eb) dE = Ecm2 - Eb
                  if (Ecm1 > Eb.and.Ecm2 < Ee) dE = Ecm2 - Ecm1
                  if (Ecm1 <= Ee.and.Ecm2 > Ee) dE = Ee - Ecm1
                  sum = sum + term * dE
                enddo
              endif
              spec = sum * fac
              if (type == 0) spec = speccm
              pfns(type, nen) = pfns(type, nen) + spec
              pfnscm(type, nen) = pfnscm(type, nen) + speccm
              write(1, '(f10.5, 2es12.5)') E, speccm, spec
            enddo
            close(1)
          endif
        enddo
      endif
    enddo
  enddo
  write(*, '(/" ########## End of loop over fission fragments"/)')
!
! Calculate P(nu)
!
  do ia = 1, ACN / 2
    do iz = 1, ZCN - 1
      in = ia - iz
      if (in < 1 .or. in > numneu) cycle
      if (xsZApre(iz, in) < fiseps .and. .not. fpexist(1, iz, in, -1)) cycle
      do type = 0, 6
        if (parskip(type)) cycle
        sum = 0.
        do npar = 0, numin
          sum = sum + Pmultiff(iz, in, type, npar)
        enddo
        Pmultiff(iz, in, type, npar) = Pmultiff(iz, in, type, npar) / sum
        izH = ZCN - iz
        iaH = ACN - ia
        inH = iaH - izH
        sum = 0.
        do npar = 0, numin
          sum = sum + Pmultiff(izH, inH, type, npar)
        enddo
        Pmultiff(izH, inH, type, npar) = Pmultiff(izH, inH, type, npar) / sum
        do npar = 0, numin
          do i = 0, npar
            Pdisnu(type, npar) = Pdisnu(type, npar) + yieldZApre(iz, in) * &
              Pmultiff(iz, in, type, i) * Pmultiff(izH, inH, type, npar - i)
          enddo
        enddo
      enddo
    enddo
  enddo
  do type = 0, 6
    if (parskip(type)) cycle
    sum = 0.
    do npar = 0, numin
      sum = sum + Pdisnu(type, npar)
    enddo
    Pdisnuav(type) = 0.
    if (sum > 0.) then
      do npar = 0, numin
        Pdisnu(type, npar) = Pdisnu(type, npar) / sum
        Pdisnuav(type) = Pdisnuav(type) + npar * Pdisnu(type, npar)
      enddo
    endif
  enddo
!
! Average energy, normalization and relation to Maxwellian
!
  do type = 0, 6
    if (parskip(type)) cycle
    sumpfns = 0.
    sumpfnscm = 0.
    Esumpfns = 0.
    Eavpfns(type) = 0.
    do nen = 1, NEpfns
      sumpfns = sumpfns + pfns(type, nen) * dEpfns(nen)
      sumpfnscm = sumpfnscm + pfnscm(type, nen) * dEpfns(nen)
    enddo
    if (sumpfns > 0.) then
      do nen = 1, NEpfns
        pfns(type, nen) = pfns(type, nen) / sumpfns
        if (sumpfnscm > 0.) pfnscm(type, nen) = pfnscm(type, nen) / sumpfnscm
        Esumpfns = Esumpfns + Epfns(nen) * pfns(type, nen) * dEpfns(nen)
        if (type == 0) then
          pfns(type, nen) = pfns(type, nen) * nubar(type)
          maxpfns(type, nen) = 1.
        else
          maxpfns(type, nen) = 0.
        endif
      enddo
      Eavpfns(type) = Esumpfns
      Eav = Eavpfns(type)
      if (type > 0) then
        if (Eav > 0.) then
          summax = 0.
          do nen = 1, NEpfns
            E = Epfns(nen)
            maxwell = sqrt(E) * exp( - E / (twothird*Eav))
            summax = summax + maxwell * dEpfns(nen)
          enddo
          do nen = 1, NEpfns
            E = Epfns(nen)
            maxwell = sqrt(E) * exp( - E / (twothird*Eav))
            if (maxwell > 0.) maxpfns(type, nen) = pfns(type, nen) / maxwell * summax
          enddo
        endif
      endif
    endif
  enddo
  sumpost = 0.
  do ia = 1, ACN
    sumpost = sumpost + xsApost(ia)
  enddo
  sumpost = 0.5 * sumpost
  if (sumpost > 0.) then
    xstotpost = 0.
    yieldtotpost = 0.
    do iz = 1, ZCN
      do ia = iz + 1, ACN
        in = ia - iz
        if (in > numneu) cycle
        if (xsZApost(iz, in) == 0.) cycle
        yieldZApost(iz, in) = xsZApost(iz, in) / sumpost
        yieldApost(ia) = yieldApost(ia) + yieldZApost(iz, in)
        xstotpost = xstotpost + xsZApost(iz, in)
        yieldtotpost = yieldtotpost + yieldZApost(iz, in)
        if (xsfpex(iz, in, 1) > 0.) then
          do nex = 0, 1
            yieldfpex(iz, in, nex) = xsfpex(iz, in, nex) / sumpost
            if (yieldZApost(iz, in) > 0.) fpratio(iz, in, nex) = yieldfpex(iz, in, nex) / yieldZApost(iz, in)
          enddo
        endif
      enddo
    enddo
  endif
!
! Reset variables to those of original target.
!
  flagffruns = .false.
  call talysinput
  flagmain = .false.
  call talysinitial
  flagmain = .true.
  if ( .not. flagompall) call basicxs(0, 0)
  if (parinclude(0)) call gamma(0, 0)
  if (enincmax >= epreeq .or. flagracap) then
    call preeqinit
    call excitoninit
  endif
  if (flagracap) call racapinit
  if (flagcomp) call compoundinit
!
! Output
!
  call massdisout
  call nubarout
  call nudisout
  if (flagspec) call pfnsout
  return
end subroutine ffevap
! Copyright A.J. Koning 2021
