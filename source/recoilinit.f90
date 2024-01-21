subroutine recoilinit
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialization of basic recoil information
!
! Author    : Stephane Hilaire
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl            ! single precision kind
! All global variables
!   numendisc      ! number of discrete outgoing energies
!   numenrec       ! maximum number of recoil energies
!   numN           ! maximum number of neutrons from initial compound nucleus
!   numpar         ! number of particles
!   numZ           ! maximum number of protons from initial compound nucleus
! Variables for numerics
!   nangle         ! number of angles
!   nanglecont     ! number of angles for continuum
!   nanglerec      ! number of recoil angles
!   maxenrec       ! number of recoil energies
!   maxN           ! maximal number of neutrons away from initial compound nucleus
!   maxZ           ! maximal number of protons away from initial compound nucleus
! Variables for energy grid
!   angle           ! angle in degrees
!   anglecont       ! angle in degrees for continuum
! Constants
!   deg2rad        ! conversion factor for degrees to radians
! Variables for recoil
!   areaejlab        ! Total surface of LAB ddx bins
!   areareclab       ! Total surface of LAB ddx bins
!   dEejlab          ! width of ejectile lab bin
!   ddxejlab         ! array containing the ddx spectrum of light part
!   ddxrec           ! array containing the lab double differential xs
!   ddxrectot        ! array containing the total recoil flux in a giv
!   Eejlab           ! center of ejectile lab bin
!   Eejlabmax        ! maximum energy of ejectile lab bin
!   Eejlabmin        ! minimum energy of ejectile lab bin
!   Erec             ! recoil energy
!   Erecmax          ! minimal energy limit of recoil bin
!   Erecmin          ! minimal energy limit of recoil bin
!   iejlab           ! number of ejectile lab bins
!   recoilint        ! total recoil integrated over spectrum
!   specrecoil       ! recoil spectrum
!   xsejlab          ! LAB ejectile spectrum
!   xsejlabint       ! LAB energy-integrated spectrum
! Variables for recoil initialization
!   angcontmax       ! maximum of angular bin
!   angcontmin       ! minimum of angular bin
!   angrecmax        ! array for upper bin angle values
!   angrecmin        ! array for lower bin angle values
!   cosangcontmax    ! cosine of maximum of angular bin
!   cosangcontmin    ! cosine of minimum of angular bin
!   cosangmax        ! cosine of maximum of angular bin
!   cosangmin        ! cosine of minimum of angular bin
!   cosrecmax        ! array for cosine of upper bin angle values
!   cosrecmin        ! array for cosine of lower bin angle values
!   dcosang          ! width of cosine bin width of cosine bin
!   dcosangcont      ! width of cosine bin
!   dcosangrec       ! array for deltacos
!   sinangcontmax    ! sine of maximum of angular bin
!   sinangcontmin    ! sine of minimum of angular bin
!   sinangmax        ! sine of maximum of angular bin
!   sinangmin        ! sine of minimum of angular bin
! Variables for direct reactions
!   elwidth        ! width of elastic peak in MeV
! Variables for main input
!   k0             ! index of incident particle
!   Ninit          ! neutron number of initial compound nucleus
!   Zinit          ! charge number of initial compound nucleus
! Variables for energy grid
!   egrid          ! outgoing energy grid
!   Einc           ! incident energy in MeV
!   maxen          ! total number of energies
! Variables for energies
!   eoutdis        ! outgoing energy of discrete state reaction
!   Etotal         ! total energy of compound system (target + projectile)
! Variables for excitation energy grid
!   maxex          ! maximum excitation energy bin for residual nucleus
! Variables for incident channel
!   xstotinc       ! total cross section (neutrons only) for incident channel
! Variables for nuclides
!   parskip        ! logical to skip outgoing particle
!   targetE        ! excitation energy of target
! Constants
!   amu            ! atomic mass unit in MeV
!   parmass        ! mass of particle in a.m.u.
!   parN           ! neutron number of particle
!   parZ           ! charge number of particle
!   twopi          ! 2 * pi
! Variables for level density
!   Nlast          ! last discrete level
! Variables for masses
!   nucmass        ! mass of nucleus
!   S              ! separation energy
! Variables for recoil initialization
!   dcosangcont      ! width of cosine bin
!   dcosangrec       ! array for deltacos
! Variables for recoil
!   areaejlab        ! Total surface of LAB ddx bins
!   areareclab       ! Total surface of LAB ddx bins
!   ddxrec           ! array containing the lab double differential xs
!   ddxrectot        ! array containing the total recoil flux in a giv
!   dEejlab          ! width of ejectile lab bin
!   ejectmass        ! Ejectile mass
!   Erec             ! recoil energy
!   Erecinit         ! first compound nucleus recoil energy
!   Erecmax          ! minimal energy limit of recoil bin
!   Erecmin          ! minimal energy limit of recoil bin
!   Eejlab           ! center of ejectile lab bin
!   Eejlabmax        ! maximum energy of ejectile lab bin
!   Eejlabmin        ! minimum energy of ejectile lab bin
!   iejlab           ! number of ejectile lab bins
!   irecinit         ! counter
!   recoilmass       ! Recoil mass
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: nbbins=50                    ! number of bins to calculate maximum recoil energy
  integer, parameter :: numres=5000                  ! number of residual bins
  integer            :: iang                         ! running variable for angle
  integer            :: iej                          ! counter
  integer            :: ierec                        ! counter
  integer            :: iex                          ! Loop over excitation energy bins counter
  integer            :: in                           ! counter for neutrons
  integer            :: iz                           ! charge number of residual nucleus
  integer            :: nl                           ! last discrete level
  integer            :: type                         ! particle type
  real(sgl)          :: angmax                       ! maximum of angular bin
  real(sgl)          :: angmin                       ! minimum of angular bin
  real(sgl)          :: angval                       ! angle
  real(sgl)          :: dang                         ! delta angle
  real(sgl)          :: dangrec                      ! help variable
  real(sgl)          :: Emaxinc                      ! incident lab energy
  real(sgl)          :: Emaxlab                      ! help variable
  real(sgl)          :: Enrjlabmax(0:numpar)         ! maximum ejectile energy in lab (classical case)
  real(sgl)          :: Erecmaxmax(0:numZ, 0:numN)   ! maximum estimated recoil energy of the nucleus
  integer            :: maxentype                    ! maximum energy per particle type
  real(sgl)          :: Elast                        ! help variable
  integer            :: nend                         ! help variable
  real(sgl)          :: ehigh                        ! highest energy
  integer            :: nhigh                        ! help variable
  integer            :: iejmax                       ! maximum energy index
  real(sgl)          :: wcos                         ! width of cosine bin
  real(sgl)          :: projectmass                  ! projectile mass
  real(sgl)          :: ekinprojlab                  ! incident energy
  real(sgl)          :: PtotCM                       ! compound initial momentum (or cm momentum)
  real(sgl)          :: compmass                     ! Composite system mass
  integer            :: numZN(numres, numpar)        ! number of ZN combinations
  real(sgl)          :: EexCMmax(numres)             ! maximum energy in C.M. frame
  integer            :: ires                         ! counter
  integer            :: iloop                        ! loop counter
  integer            :: io1                          ! help variable
  integer            :: io2                          ! help variable
  integer            :: io3                          ! help variable
  integer            :: io4                          ! help variable
  integer            :: io5                          ! help variable
  integer            :: io6                          ! help variable
  integer            :: if1                          ! help variable
  integer            :: if2                          ! help variable
  integer            :: if3                          ! help variable
  integer            :: if4                          ! help variable
  integer            :: if5                          ! help variable
  integer            :: if6                          ! help variable
  integer            :: numZcomp                     ! number of protons
  integer            :: numNcomp                     ! number of neutrons
  integer            :: numZres                      ! number of protons
  integer            :: numNres                      ! number of neutrons
  integer            :: numZk                        ! number of protons
  integer            :: numNk                        ! number of neutrons
  real(sgl)          :: EexCMloc                     ! maximum energy in C.M. frame
  integer            :: k                            ! designator for particle
  integer            :: numrestot                    ! number of residual bins
  real(sgl)          :: dEexCM                       ! width of energy bin
  real(sgl)          :: EexCM(numres, 0:nbbins)      ! energy in C.M. frame
  real(sgl)          :: PCM(numres, 0:nbbins)        ! momentum in C.M. frame
  integer            :: iexc                         ! counter
  real(sgl)          :: ECMbin                       ! energy of C.M. bin
  real(sgl)          :: vCMbin                       ! velocity of C.M. bin
  integer            :: ifinal                       ! help variable
  real(sgl)          :: recoilGSmass                 ! mass
  integer            :: iexfinal                     ! counter
  real(sgl)          :: ratio                        ! if 0<ratio<1 then x is between xl1 and xl2
  real(sgl)          :: EejCM                        ! ejectile C.M. energy
  real(sgl)          :: PrecCM                       ! C.M. recoil momentum
  real(sgl)          :: vrecCM                       ! C.M. recoil velocity
  real(sgl)          :: preclab                      ! momentum of recoil in LAB frame
  real(sgl)          :: vreclab                      ! velocity of recoil in LAB frame
  real(sgl)          :: vejeccm                      ! C.M. ejectile velocity
  real(sgl)          :: vejeclab                     ! LAB ejectile velocity
  real(sgl)          :: eejeclab                     ! LAB ejectile energy
  real(sgl)          :: Erecmaxloc(numres)           ! maximum recoil energy
  real(sgl)          :: Erecbin                      ! recoil energy in bin
  real(sgl)          :: derec                        ! help variable
  integer            :: irecmaxmax(0:numZ, 0:numN)   ! counter
  real(sgl)          :: wnrj                         ! width of recoil energy bin
!
! **** Initialization of arrays containing the results in the lab ******
!
  areareclab = 0.
  areaejlab = 0.
  ddxejlab = 0.
  ddxrec = 0.
  ddxrectot = 0.
  dEejlab = 0.
  Eejlab = 0.
  Eejlabmax = 0.
  Eejlabmin = 0.
  Erec = 0.
  Erecmax = 0.
  Erecmin = 0.
  iejlab = 0
  recoilint = 0.
  specrecoil = 0.
  xsejlab = 0.
  xsejlabint = 0.
!
! ******** Creation of the angular lab grid for the recoil nuclei ******
!
! Angular information necessary for recoil calculation
!
  dang = 180. / nangle
  do iang = 0, nangle
    angval = angle(iang)
    angmin = max(0., angval - 0.5 * dang)
    angmax = min(180., angval + 0.5 * dang)
    cosangmin(iang) = cos(angmin * deg2rad)
    cosangmax(iang) = cos(angmax * deg2rad)
    sinangmin(iang) = sin(angmin * deg2rad)
    sinangmax(iang) = sin(angmax * deg2rad)
    dcosang(iang) = abs(cosangmin(iang) - cosangmax(iang))
  enddo
  sinangmin(0) = 1.e-6
  sinangmax(nangle) = 1.e-6
  dang = 180. / nanglecont
  do iang = 0, nanglecont
    angval = anglecont(iang)
    angcontmin(iang) = max(0., angval - 0.5 * dang)
    angcontmax(iang) = min(180., angval + 0.5 * dang)
    cosangcontmin(iang) = cos(angcontmin(iang) * deg2rad)
    cosangcontmax(iang) = cos(angcontmax(iang) * deg2rad)
    sinangcontmin(iang) = sin(angcontmin(iang) * deg2rad)
    sinangcontmax(iang) = sin(angcontmax(iang) * deg2rad)
    dcosangcont(iang) = abs(cosangcontmin(iang) - cosangcontmax(iang))
  enddo
  do iang = nanglecont + 1, 2 * nanglecont + 1
    cosangcontmax(iang) = - cosangcontmax(iang - nanglecont - 1)
    cosangcontmin(iang) = - cosangcontmin(iang - nanglecont - 1)
    sinangcontmin(iang) = - sinangcontmin(iang - nanglecont - 1)
    sinangcontmax(iang) = - sinangcontmax(iang - nanglecont - 1)
    dcosangcont(iang) = dcosangcont(2 * nanglecont + 1 - iang)
  enddo
  sinangcontmin(0) = 1.e-6
  sinangcontmax(nanglecont) = 1.e-6
  sinangcontmax(2 * nanglecont + 1) = - 1.e-6
  sinangcontmin(nanglecont + 1) = - 1.e-6
  dangrec = 180.0 / nanglerec
  do iang = 0, nanglerec
    angval = dangrec * real(iang)
    angrecmin(iang) = max(0., angval - 0.5 * dangrec)
    angrecmax(iang) = min(180., angval + 0.5 * dangrec)
    cosrecmin(iang) = cos(angrecmin(iang) * deg2rad)
    cosrecmax(iang) = cos(angrecmax(iang) * deg2rad)
    dcosangrec(iang) = abs(cosrecmin(iang) - cosrecmax(iang))
  enddo
  do iang = nanglerec + 1, 2 * nanglerec + 1
    cosrecmax(iang) = - cosrecmax(iang - nanglerec - 1)
    cosrecmin(iang) = - cosrecmin(iang - nanglerec - 1)
    dcosangrec(iang) = abs(cosrecmin(iang) - cosrecmax(iang))
  enddo
!
! ********** Creation of the lab energy grid for the ejectiles *********
!
  Emaxinc = Einc + S(0, 0, k0) + targetE
  do type = 0, 6
    if (parskip(type)) cycle
    Emaxlab = Emaxinc - S(0, 0, type)
    Enrjlabmax(type) = Emaxlab
    iejlab(type) = 0
    if (Emaxlab < 0.) cycle
    maxentype = 0
    do iej = 0, maxen
      if (egrid(iej) >= Emaxlab) then
        maxentype = iej
        exit
      endif
    enddo
    if (maxen > 0 .and. egrid(maxen) < Emaxlab) maxentype = maxen
    if (maxentype > 0) maxentype = max(maxentype, 3)
    nl = Nlast(parZ(type), parN(type), 0)
    Elast = eoutdis(type, nl) - elwidth
    Elast = max(Elast, 0.)
    if (Elast > egrid(maxentype)) then
      nend = maxentype
    else
      call locate(egrid, 0, maxentype, Elast, nend)
    endif
    do iej = 1, nend
      Eejlab(type, iej) = egrid(iej)
    enddo
    ehigh = egrid(maxentype) - egrid(nend)
    nhigh = nint(10. * ehigh)
    nhigh = min(nhigh, numendisc)
    iejlab(type) = nend + nhigh
    do iej = nend + 1, nend + nhigh
      Eejlab(type, iej) = egrid(nend) + 0.1 * (iej - nend)
      if (Eejlab(type, iej) > Emaxlab) then
        iejlab(type) = min(iejlab(type), iej)
      endif
    enddo
    iejmax = iejlab(type)
    do iej = 2, iejmax - 1
      dEejlab(type, iej) = 0.5 * (Eejlab(type, iej + 1) - Eejlab(type, iej - 1))
      Eejlabmax(type, iej) = 0.5 * (Eejlab(type, iej) + Eejlab(type, iej + 1))
      Eejlabmin(type, iej) = 0.5 * (Eejlab(type, iej) + Eejlab(type, iej - 1))
    enddo
    dEejlab(type, 1) = Eejlab(type, 1) + 0.5 * (Eejlab(type, 2) - Eejlab(type, 1))
    Eejlabmin(type, 1) = 0.
    Eejlabmax(type, 1) = dEejlab(type, 1)
    Eejlabmin(type, iejmax) = Eejlabmax(type, max(0, iejmax - 1))
    Eejlabmax(type, iejmax) = 2. * Eejlab(type, iejmax) - Eejlabmin(type, iejmax)
    Eejlabmax(type, iejmax) = Emaxlab
    dEejlab(type, iejmax) = Eejlabmax(type, iejmax) - Eejlabmin(type, iejmax)
    Eejlab(type, iejmax) = 0.5 * (Emaxlab + Eejlabmin(type, iejmax))
  enddo
!
! ************ Creation of the area array for the ejectiles ************
!
  do type = 0, 6
    if (parskip(type)) cycle
    do iang = 0, 2 * nanglecont + 1
      wcos = dcosangcont(iang)
      do iej = 1, iejlab(type)
        areaejlab(type, iej, iang) = wcos * dEejlab(type, iej)
      enddo
      areaejlab(type, 0, iang) = 0.
    enddo
  enddo
!
! *************** Calculate maximum excitation energies ****************
! ***************   for all possible residual nuclei    ****************
!
  projectmass = parmass(k0) * amu
  ekinprojlab = Einc
  PtotCM = sqrt(ekinprojlab * (ekinprojlab + 2. * projectmass))
  compmass = nucmass(0, 0) * amu
  Erecinit = sqrt(PtotCM **2 + compmass **2) - compmass
  numZN(1, 1) = 0
  numZN(1, 2) = 0
  numZN(1, 3) = 0
  numZN(1, 4) = 0
  numZN(1, 5) = 0
  numZN(1, 6) = 0
  EexCMmax(1) = Etotal
  numrestot = 1
!
! We loop over all possible residual to see if further emission of light particles is possible.
! Starting from a given nucleus defined by (numZcomp,numNcomp) we look if the nucleus reached after emission of an ejectile
! (numZres,numNres) has already been obtained from a previous emission.
! If not the number a new residual is reached and we have to loop again.
! If yes, we compare the maximum excitation energy obtained from the previous emission with that of the current emission
! and keep the highest.
!
  do
    do ires = 1, numrestot
      iloop = 1
      if (ires == numrestot) iloop = 0
      numZcomp = numZN(ires, 2) + numZN(ires, 3) + numZN(ires, 4) + 2 * (numZN(ires, 5) + numZN(ires, 6))
      numNcomp = numZN(ires, 1) + numZN(ires, 3) + numZN(ires, 5) + 2 * (numZN(ires, 4) + numZN(ires, 6))
Loop1: do type = 1, 6
        if1 = numZN(ires, 1)
        if2 = numZN(ires, 2)
        if3 = numZN(ires, 3)
        if4 = numZN(ires, 4)
        if5 = numZN(ires, 5)
        if6 = numZN(ires, 6)
        if (type == 1) if1 = if1 + 1
        if (type == 2) if2 = if2 + 1
        if (type == 3) if3 = if3 + 1
        if (type == 4) if4 = if4 + 1
        if (type == 5) if5 = if5 + 1
        if (type == 6) if6 = if6 + 1
        numZres = numZcomp + parZ(type)
        if (numZres > maxZ + 2) cycle
        if (numZres > Zinit) cycle
        numNres = numNcomp + parN(type)
        if (numNres > maxN + 2) cycle
        if (numNres > Ninit) cycle
        EexCMloc = EexCMmax(ires) - S(numZcomp, numNcomp, type)
        if (EexCMloc <= 0.) cycle
        do k = 1, numrestot
          io1 = numZN(k, 1)
          io2 = numZN(k, 2)
          io3 = numZN(k, 3)
          io4 = numZN(k, 4)
          io5 = numZN(k, 5)
          io6 = numZN(k, 6)
          numZk = numZN(k, 2) + numZN(k, 3) + numZN(k, 4) + 2 * (numZN(k, 5) + numZN(k, 6))
          numNk = numZN(k, 1) + numZN(k, 3) + numZN(k, 5) + 2 * (numZN(k, 4) + numZN(k, 6))
          if ((numZk == numZres) .and. (numNk == numNres)) then
            EexCMmax(k) = max(EexCMmax(k), EexCMloc)
            cycle Loop1
          endif
        enddo
!
! A new residual is reached
!
        numrestot = numrestot + 1
        iloop = 1
        do k = 1, 6
          numZN(numrestot, k) = numZN(ires, k)
        enddo
        numZN(numrestot, type) = numZN(ires, type) + 1
        EexCMmax(numrestot) = EexCMloc
      enddo Loop1
    enddo
    if (iloop /= 1) exit
  enddo
!
! ******** Define excitation energy bins for all residual nuclei *******
!
  do ires = 1, numrestot
    dEexCM = EexCMmax(ires) / nbbins
    do iex = 0, nbbins
      EexCM(ires, iex) = min(iex * dEexCM, EexCMmax(ires))
      EexCM(ires, iex) = max(iex * dEexCM, 0.)
      PCM(ires, iex) = 0.
    enddo
  enddo
!
! ********      Loop over all reactions paths to determine      ********
! ******** The maximum recoil energy possible for each residual ********
!
  PCM(1, nbbins) = PtotCM
  do ires = 1, numrestot
    numZcomp = numZN(ires, 2) + numZN(ires, 3) + numZN(ires, 4) + 2 * (numZN(ires, 5) + numZN(ires, 6))
    numNcomp = numZN(ires, 1) + numZN(ires, 3) + numZN(ires, 5) + 2 * (numZN(ires, 4) + numZN(ires, 6))
    compmass = nucmass(numZcomp, numNcomp) * amu
    do iex = nbbins, 0, - 1
      ECMbin = EexCM(ires, iex)
      vCMbin = PCM(ires, iex) / (compmass + ECMbin)
      do type = 0, 6
        ejectmass = parmass(type) * amu
        if1 = numZN(ires, 1)
        if2 = numZN(ires, 2)
        if3 = numZN(ires, 3)
        if4 = numZN(ires, 4)
        if5 = numZN(ires, 5)
        if6 = numZN(ires, 6)
        if (type == 1) if1 = if1 + 1
        if (type == 2) if2 = if2 + 1
        if (type == 3) if3 = if3 + 1
        if (type == 4) if4 = if4 + 1
        if (type == 5) if5 = if5 + 1
        if (type == 6) if6 = if6 + 1
        do k = 1, numrestot
          io1 = numZN(k, 1)
          io2 = numZN(k, 2)
          io3 = numZN(k, 3)
          io4 = numZN(k, 4)
          io5 = numZN(k, 5)
          io6 = numZN(k, 6)
          if ((io1 == if1) .and. (io2 == if2) .and. (io3 == if3) .and. (io4 == if4) .and. (io5 == if5) .and. (io6 == if6)) then
            ifinal = k
            numZres = numZN(k, 2) + numZN(k, 3) + numZN(k, 4) + 2 * (numZN(k, 5) + numZN(k, 6))
            numNres = numZN(k, 1) + numZN(k, 3) + numZN(k, 5) + 2 * (numZN(k, 4) + numZN(k, 6))
!
! If the residual nucleus is not among the considered ones it means that it cannot be obtained because all possible residuals
! have been defined. We thus consider another ejectile.
!
! Loop over the various excitation energies bins of the residual nucleus to deduce the lab recoil energy for these bins.
!
            recoilGSmass = nucmass(numZres, numNres) * amu
            do iexfinal = nbbins, 0, - 1
              recoilmass = recoilGSmass + EexCM(ifinal, iexfinal)
              if (type /= 0) ratio = recoilmass / (ejectmass * (ejectmass + recoilmass))
!
! Determine ejectile energy to check if emission can occur
!
              EejCM = ECMbin - EexCM(ifinal, iexfinal)
              EejCM = EejCM - S(numZcomp, numNcomp, type)
              if (abs(EejCM) <= 1.0e-10) EejCM = 0.
              if (EejCM <= 0.) cycle
              if (type /= 0) then
                  PrecCM = ejectmass * sqrt(2 * ratio * EejCM)
                else
                  PrecCM = EejCM
              endif
              vrecCM = PrecCM / recoilmass
              vreclab = vCMbin + vrecCM
              preclab = recoilmass * vreclab
              PCM(ifinal, iexfinal) = max(PCM(ifinal, iexfinal), preclab)
              if (type /= 0) then
                  vejeccm = PrecCM / ejectmass
                  vejeclab = vCMbin + vejeccm
                  eejeclab = 0.5 * ejectmass * vejeclab **2
                  Enrjlabmax(type) = max(Enrjlabmax(type), eejeclab)
                else
                  eejeclab = EejCM
                  Enrjlabmax(type) = max(Enrjlabmax(type), eejeclab)
              endif
            enddo
          endif
        enddo
      enddo
    enddo
  enddo
!
! Maximum recoil energy for each residual
!
  do ires = 1, numrestot
    Erecmaxloc(ires) = 0.
    numZres = numZN(ires, 2) + numZN(ires, 3) + numZN(ires, 4) + 2 * (numZN(ires, 5) + numZN(ires, 6))
    numNres = numZN(ires, 1) + numZN(ires, 3) + numZN(ires, 5) + 2 * (numZN(ires, 4) + numZN(ires, 6))
    recoilGSmass = nucmass(numZres, numNres) * amu
    do iexc = 0, nbbins
      recoilmass = recoilGSmass + EexCM(ires, iexc)
      PrecCM = PCM(ires, iexc)
      Erecbin = PrecCM **2 / (2. * recoilmass)
      Erecmaxloc(ires) = max(Erecmaxloc(ires), Erecbin)
    enddo
  enddo
!
! Finally the recoils energy grids are defined and maximum excitation energies are stored
!
  do ires = 1, numrestot
    numZres = numZN(ires, 2) + numZN(ires, 3) + numZN(ires, 4) + 2 * (numZN(ires, 5) + numZN(ires, 6))
    numNres = numZN(ires, 1) + numZN(ires, 3) + numZN(ires, 5) + 2 * (numZN(ires, 4) + numZN(ires, 6))
    Erecmaxmax(numZres, numNres) = Erecmaxloc(ires)
    irecmaxmax(numZres, numNres) = ires
  enddo
  do iz = 0, maxZ + 2
    do in = 0, maxN + 2
      if (Erecmaxmax(iz, in) <= 0.) cycle
      derec = Erecmaxmax(iz, in) / (maxenrec + 1)
      ires = irecmaxmax(iz, in)
      if (ires == 0) cycle
      do iexc = 0, numenrec
        Erecmin(iz, in, iexc) = iexc * derec
        Erecmax(iz, in, iexc) = (iexc + 1) * derec
        Erec(iz, in, iexc) = iexc * derec + 0.5 * derec
      enddo
    enddo
  enddo
!
! ************* Creation of the area array for the recoils *************
!
  do iz = 0, maxZ + 2
    do in = 0, maxN + 2
      if (Erecmaxmax(iz, in) == 0.) cycle
      do ierec = 0, maxenrec
        wnrj = Erecmax(iz, in, ierec) - Erecmin(iz, in, ierec)
        do iang = 0, 2 * nanglerec + 1
          wcos = dcosangrec(iang)
          areareclab(iz, in, ierec, iang) = wnrj * wcos
        enddo
      enddo
    enddo
  enddo
!
! ******* Feeding of the first compound nucleus recoil ddx array  ******
!
  derec = Erecmaxmax(0, 0) / (maxenrec + 1)
  irecinit = min(int(Erecinit / derec), numenrec)
  if (areareclab(0, 0, irecinit, 0) /= 0.) ddxrec(0, 0, maxex(0, 0), irecinit, 0) = xstotinc / twopi / areareclab(0, 0, irecinit, 0)
  ddxrectot(0, 0, maxex(0, 0)) = xstotinc / twopi
  return
end subroutine recoilinit
! Copyright A.J. Koning 2021
