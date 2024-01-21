subroutine angdisrecoil
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Recoil angular distributions for discrete states
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
!   numangcont     ! maximum number of angles for continuum
!   numangrec      ! maximum number of recoil angles
!   numen2         ! maximum number of outgoing energies
!   numenrec       ! maximum number of recoil energies
!   numpar         ! number of particles
! Variables for numerics
!   maxenrec       ! number of recoil energies
!   nangle         ! number of angles
!   nanglecont     ! number of angles for continuum
!   nanglerec      ! number of recoil angles
! Variables for basic reaction
!   flaglabddx     ! flag for calculation of DDX in LAB system
! Variables for main input
!   k0             ! index of incident particle
! Variables for energies
!   Etotal         ! total energy of compound system (target + projectile)
! Variables for excitation energy grid
!   nexmax         ! maximum excitation energy bin for residual nucleus
! Variables for excitation energy grid
!   Ex             ! excitation energy
!   Exmax          ! maximum excitation energy for residual nucleus
!   maxex          ! maximum excitation energy bin for residual nucleus
! Variables for binary reactions
!   xsdisc         ! total cross section for discrete state
!   xselastot      ! total elastic cross section (shape + compound)
! Variables for incident channel
!   xselasinc      ! total elastic cross section (neutrons only) for inc. channel
! Variables for angular distributions
!   discad         ! discrete state angular distribution
! Variables for nuclides
!   Nindex         ! neutron number index for residual nucleus
!   parskip        ! logical to skip outgoing particle
!   Zindex         ! charge number index for residual nucleus
! Constants
!   amu            ! atomic mass unit in MeV
!   parmass        ! mass of particle in a.m.u.
!   twopi          ! 2 * pi
! Variables for level density
!   Nlast          ! last discrete level
! Variables for recoil initialization
!   cosangmax      ! cosine of maximum of angular bin
!   cosangmin      ! cosine of minimum of angular bin
!   dcosang        ! width of cosine bin width of cosine bin
!   sinangmax      ! sine of maximum of angular bin
!   sinangmin      ! sine of minimum of angular bin
! Variables for recoil
!   angcm          ! CM angle with respect to the LAB
!   areaejlab      ! Total surface of LAB ddx bins
!   areareclab     ! Total surface of LAB ddx bins
!   cosejcm1       ! CM recoil energy cosine corresponding to (E1, ang1)
!   cosejcm2       ! CM recoil energy cosine corresponding to (E1, ang2)
!   cosejlab11     ! LAB ejectile angle cosine corresponding to (E1, ang1)
!   cosejlab12     ! LAB ejectile angle cosine corresponding to (E1, ang2)
!   cosejlab21     ! LAB ejectile angle cosine corresponding to (E2, ang1)
!   cosejlab22     ! LAB ejectile angle cosine corresponding to (E2, ang2)
!   cosreclab11    ! LAB recoil angle cosine corresponding to (E1, ang1)
!   cosreclab12    ! LAB recoil angle cosine corresponding to (E1, ang2)
!   cosreclab21    ! LAB recoil angle cosine corresponding to (E2, ang1)
!   cosreclab22    ! LAB recoil angle cosine corresponding to (E2, ang2)
!   ddxejlab       ! array containing the ddx spectrum of light part
!   ddxrec         ! array containing the lab double differential xs
!   ddxrectot      ! array containing the total recoil flux in a giv
!   dEejlab        ! width of ejectile lab bin
!   Eejcm1         ! Lower limit of CM ejectile energy bin
!   Eejcm2         ! Upper limit of CM ejectile energy bin
!   Eejlab11       ! LAB ejectile energy corresponding to (E1, ang1)
!   Eejlab12       ! LAB ejectile energy corresponding to (E1, ang2)
!   Eejlab21       ! LAB ejectile energy corresponding to (E2, ang1)
!   Eejlab22       ! LAB ejectile energy corresponding to (E2, ang2)
!   Eejlabmax      ! maximum energy of ejectile lab bin
!   Eejlabmin      ! minimum energy of ejectile lab bin
!   ejectmass      ! Ejectile mass
!   Erecinit       ! first compound nucleus recoil energy
!   Ereclab11      ! LAB recoil energy corresponding to (E1, ang1)
!   Ereclab12      ! LAB recoil energy corresponding to (E1, ang2)
!   Ereclab21      ! LAB recoil energy corresponding to (E2, ang1)
!   Ereclab22      ! LAB recoil energy corresponding to (E2, ang2)
!   Erecmax        ! minimal energy limit of recoil bin
!   Erecmin        ! minimal energy limit of recoil bin
!   iejlab         ! number of ejectile lab bins
!   recoilmass     ! Recoil mass
!   sinejcm1       ! CM recoil energy cosine corresponding to (E2, ang1)
!   sinejcm2       ! CM recoil energy cosine corresponding to (E2, ang2)
!   sinejlab11     ! LAB ejectile angle sine corresponding to (E1, ang1)
!   sinejlab12     ! LAB ejectile angle sine corresponding to (E1, ang2)
!   sinejlab21     ! LAB ejectile angle sine corresponding to (E2, ang1)
!   sinejlab22     ! LAB ejectile angle sine corresponding to (E2, ang2)
!   sinreclab11    ! LAB recoil angle sine corresponding to (E1, ang1)
!   sinreclab12    ! LAB recoil angle sine corresponding to (E1, ang2)
!   sinreclab21    ! LAB recoil angle sine corresponding to (E2, ang1)
!   sinreclab22    ! LAB recoil angle sine corresponding to (E2, ang2)
!   vcm            ! Compound nucleus velocity
!   vejcm1         ! velocity corresponding to Eejcm1
!   vejcm2         ! velocity corresponding to Eejcm2
!   vreccm1        ! Recoil velocity corresponding to Eejcm1
!   vreccm2        ! Recoil velocity corresponding to Eejcm2
! Variables for masses
!   nucmass        ! mass of nucleus
!   S              ! separation energy
!
! *** Declaration of local data
!
  implicit none
  integer   :: iang                                          ! running variable for angle
  integer   :: Ncomp                                         ! neutron number index for compound nucleus
  integer   :: nex                                           ! excitation energy bin of compound nucleus
  integer   :: Nix                                           ! neutron number index for residual nucleus
  integer   :: type                                          ! particle type
  integer   :: Zcomp                                         ! proton number index for compound nucleus
  integer   :: Zix                                           ! charge number index for residual nucleus
  real(sgl) :: angnorm                                       ! angular normalisation
  real(sgl) :: compmass                                      ! Composite system mass
  real(sgl) :: ddxCM                                         ! Double differential cross section for a given CM  energy-a
  real(sgl) :: Exrec1                                        ! Recoil Excitation energy corresponding to Eejcm1
  real(sgl) :: Exrec2                                        ! Recoil Excitation energy corresponding to Eejcm
  real(sgl) :: fluxCM                                        ! Total flux that must be distributed in the LAB
  real(sgl) :: inorm                                         ! counter
  real(sgl) :: pejcm1                                        ! Impulsion corresponding to Eejcm1
  real(sgl) :: pejcm2                                        ! Impulsion corresponding to Eejcm2
  real(sgl) :: ratio                                         ! if 0<ratio<1 then x is between xl1 and xl2
  real(sgl) :: SS                                            ! separation energy
  real(sgl) :: sumang                                        ! integral over angles
  real(sgl) :: xsCMlev                                       ! discrete cross section
  real(sgl) :: xsCMlev0                                      ! elastic cross section
  integer   :: irenorm                                       ! index indicating if lab array renormalisation is needed
  real(sgl) :: Emaxi                                         ! help variable
  real(sgl) :: dEori                                         ! help variable
  real(sgl) :: dEnew                                         ! help variable
  real(sgl) :: renorm                                        ! renormalisation factor
  real(sgl) :: labsurf                                       ! Total covered surface in the LAB
  real(sgl) :: labsurf1                                      ! total surface covered in the LAB by the ejectile image  in
  real(sgl) :: labsurf2                                      ! total surface covered in the LAB by the ejectile image  in
  real(sgl) :: scovej1(1:numen2, 0:2*numangcont+1)           ! surface covered for each LAB energy-angle bin
  real(sgl) :: scovej2(1:numen2, 0:2*numangcont+1)           ! surface covered for each LAB energy-angle bin
  real(sgl) :: scovrec1(0:numenrec, 0:2*numangrec+1)         ! surface covered for each LAB energy-angle bin
  real(sgl) :: scovrec2(0:numenrec, 0:2*numangrec+1)         ! surface covered for each LAB energy-angle bin
  integer   :: xlim1(2)                                      ! limits in the LAB ejectile energy grid
  integer   :: ylim1(2)                                      ! limits in the LAB ejectile angular grid
  integer   :: xlim2(2)                                      ! limits in the LAB ejectile energy grid
  integer   :: ylim2(2)                                      ! limits in the LAB ejectile angular grid
  real(sgl) :: ddxLAB                                        ! Double differential cross section for the LAB  energy-angu
  integer   :: iymax                                         ! maximum y-loop index
  integer   :: iymaxp1                                       ! maximum y-loop index
  integer   :: ix                                            ! help variable
  integer   :: iy                                            ! help variable
  integer   :: iymod                                         ! help variable
  integer   :: iys                                           ! help variable
  real(sgl) :: ddxejlabdis(0:numpar, 0:numen2, 0:numangcont) ! DDX in lab system for discrete states
  real(sgl) :: ddxrecadd                                     ! help variable
  real(sgl) :: sumddxrec                                     ! help variable
  real(sgl) :: fluxadd                                       ! flux added to recoil bin
  real(sgl) :: surfbin                                       ! LAB DDX area contribution
  integer   :: iarec2                                        ! counter
  integer   :: inex                                          ! counter
!
! Recoil local variables
!
! Initialisations (angcm=0 since first compound system decay)
!
  Zcomp = 0
  Ncomp = 0
  compmass = nucmass(Zcomp, Ncomp) * amu
  vcm = sqrt(2. * Erecinit / compmass)
  angcm = 0.
!
! Loop over ejectile type
!
  do type = 0, 6
    if (parskip(type)) cycle
    Zix = Zindex(Zcomp, Ncomp, type)
    Nix = Nindex(Zcomp, Ncomp, type)
    SS = S(Zcomp, Ncomp, type)
    ejectmass = parmass(type) * amu
    recoilmass = nucmass(Zix, Nix) * amu
    if (type /= 0) then
      ratio = recoilmass / (ejectmass * (ejectmass + recoilmass))
    endif
!
! Loop over residual nucleus discrete states
!
    if (type == k0) then
        xsCMlev0 = xselastot
    else
    xsCMlev0 = 0.
    endif
    do nex = 0, Nlast(Zix, Nix, 0)
      sumang = 0.
      if (nex /= 0) xsCMlev0 = 0.
      xsCMlev = xsCMlev0 + xsdisc(type, nex)
      if (nex > nexmax(type)) cycle
!
! Initialise ddxejlabdis
!
      if (flaglabddx) then
        do ix = 0, numen2
          do iy = 0, numangcont
            ddxejlabdis(type, ix, iy) = 0.
          enddo
        enddo
      endif
!
! Determine the residual nucleus bin associated with nex
!
      do
        if (nex == 0) then
          Exrec1 = Ex(Zix, Nix, 0)
          Exrec2 = 0.5 * (Ex(Zix, Nix, 1) + Ex(Zix, Nix, 0))
          exit
        endif
        if (nex == nexmax(type)) then
          Exrec2 = Exmax(Zix, Nix)
          Exrec1 = 0.5 * (Exmax(Zix, Nix) + Ex(Zix, Nix, nex - 1))
          exit
        endif
        Exrec1 = 0.5 * (Ex(Zix, Nix, nex) + Ex(Zix, Nix, nex - 1))
        Exrec2 = 0.5 * (Ex(Zix, Nix, nex) + Ex(Zix, Nix, nex + 1))
        exit
      enddo
!
! Determination of ejectile energies in the CM frame
!
      Eejcm1 = Etotal - SS - Exrec1
      Eejcm2 = Etotal - SS - Exrec2
      if (Eejcm1 < 0..or.Eejcm2 < 0.) cycle
!
! Determine recoil and ejectile momentum (or velocities)
!
      if (type /= 0) then
        pejcm1 = ejectmass * sqrt(2. * ratio * Eejcm1)
        pejcm2 = ejectmass * sqrt(2. * ratio * Eejcm2)
        vejcm1 = pejcm1 / ejectmass
        vejcm2 = pejcm2 / ejectmass
      else
        pejcm1 = Eejcm1
        pejcm2 = Eejcm2
      endif
      vreccm1 = pejcm1 / recoilmass
      vreccm2 = pejcm2 / recoilmass
!
! Loop over ejectile angles in the CM frame
!
      inorm = 0
      do
        do iang = 0, nangle
          cosejcm1 = cosangmin(iang)
          cosejcm2 = cosangmax(iang)
          sinejcm1 = sinangmin(iang)
          sinejcm2 = sinangmax(iang)
!
! Total flux that must be spread in the LAB frame
!
          ddxCM = discad(type, nex, iang)
          if (inorm == 0) sumang = sumang + ddxCM * abs(cosejcm1 - cosejcm2) * twopi
!
! Renormalisation of angular distribution if inorm=1
!
          angnorm = 1.
          if (inorm == 0) cycle
          if (iang == 0) then
            if (type /= k0) then
              if (sumang > 1.e-14) then
                angnorm = xsdisc(type, nex) / sumang
              endif
            else
              if (nex == 0 .and. sumang > 1.e-14) then
                angnorm = (xsdisc(type, nex) + xselasinc) / sumang
              else
                if (sumang > 1.e-14) then
                  angnorm = xsdisc(type, nex) / sumang
                endif
              endif
            endif
          endif
          fluxCM = ddxCM * dcosang(iang) * angnorm
          if (fluxCM == 0.) cycle
!
! EJECTILE and RECOIL TREATMENT
!
! Ejectile and recoil LAB angles and energies deduced from CM points
!
          call cm2lab(type)
!
! We check if ejectile energies are greater than the maximum value.
! If it is so we accordingly modify the maximum ejectile energy and renormalise the ejectile area array as well as
! the ejectile spectrum array in the lab
!
          if (flaglabddx) then
            irenorm = 0
            Emaxi = Eejlabmax(type, iejlab(type))
            dEori = Emaxi - Eejlabmin(type, iejlab(type))
            if (Eejlab11 > Emaxi) then
              Emaxi = Eejlab11
              irenorm = 1
            endif
            if (Eejlab12 > Emaxi) then
              Emaxi = Eejlab12
              irenorm = 1
            endif
            if (Eejlab21 > Emaxi) then
              Emaxi = Eejlab21
              irenorm = 1
            endif
            if (Eejlab22 > Emaxi) then
              Emaxi = Eejlab22
              irenorm = 1
            endif
            if (irenorm == 1) then
              Eejlabmax(type, iejlab(type)) = Emaxi
              dEnew = Emaxi - Eejlabmin(type, iejlab(type))
              dEejlab(type, iejlab(type)) = dEnew
              renorm = dEori / dEnew
              do iarec2 = 0, nanglecont
                areaejlab(type, iejlab(type), iarec2) = areaejlab(type, iejlab(type), iarec2) / renorm
                ddxejlab(type, iejlab(type), iarec2) = ddxejlab(type, iejlab(type), iarec2) * renorm
              enddo
            endif
!
! Calculate lab bins occupation for the ejectile in the LAB.
! We assume the image in the LAB of a CM triangle is a triangle.
! This is theoretically wrong but should be a good approximation.
!
            call labsurface(Zcomp, Ncomp, type, 1, Eejlab11, cosejlab11, sinejlab11, Eejlab12, cosejlab12, sinejlab12, Eejlab21, &
 &            cosejlab21, sinejlab21, labsurf1, xlim1, ylim1, scovej1, scovrec1)
            call labsurface(Zcomp, Ncomp, type, 1, Eejlab21, cosejlab21, sinejlab21, Eejlab22, cosejlab22, sinejlab22, Eejlab12, &
 &            cosejlab12, sinejlab12, labsurf2, xlim2, ylim2, scovej2, scovrec2)
            labsurf = labsurf1 + labsurf2
            if (labsurf == 0.) cycle
            ddxLAB = fluxCM / labsurf
!
! store ejectile flux in lab ddx array
!
            iymax = 2 * nanglecont + 1
            iymaxp1 = iymax + 1
            do ix = xlim1(1), xlim1(2)
              do iy = ylim1(1), ylim1(2)
                iymod = mod(iy, iymaxp1)
                if (iymod <= nanglecont) then
                  iys = iymod
                else
                  iys = iymax - iymod
                endif
                if (areaejlab(type, ix, iymod) > 0.) ddxejlabdis(type, ix, iys) = ddxejlabdis(type, ix, iys) + &
 &                ddxLAB / areaejlab(type, ix, iymod) * scovej1(ix, iymod)
              enddo
            enddo
            do  ix = xlim2(1), xlim2(2)
              do iy = ylim2(1), ylim2(2)
                iymod = mod(iy, iymaxp1)
                if (iymod <= nanglecont) then
                    iys = iymod
                  else
                    iys = iymax - iymod
                endif
                if (areaejlab(type, ix, iymod) > 0.) ddxejlabdis(type, ix, iys) = ddxejlabdis(type, ix, iys) + &
 &                ddxLAB / areaejlab(type, ix, iymod) * scovej2(ix, iymod)
              enddo
            enddo
          endif
!
! We check if recoil energies are greater than the maximum value calculated.
! If it is so we accordingly modify the maximum recoil energy, the recoil area array as well as the recoil spectrum array
! in the lab.
!
          irenorm = 0
          Emaxi = Erecmax(Zix, Nix, maxenrec)
          dEori = Emaxi - Erecmin(Zix, Nix, maxenrec)
          if (Ereclab11 > Emaxi) then
            Emaxi = Ereclab11
            irenorm = 1
          endif
          if (Ereclab12 > Emaxi) then
            Emaxi = Ereclab12
            irenorm = 1
          endif
          if (Ereclab21 > Emaxi) then
            Emaxi = Ereclab21
            irenorm = 1
          endif
          if (Ereclab22 > Emaxi) then
            Emaxi = Ereclab22
            irenorm = 1
          endif
          if (irenorm == 1) then
            Erecmax(Zix, Nix, maxenrec) = Emaxi
            dEnew = Emaxi - Erecmin(Zix, Nix, maxenrec)
            renorm = dEori / dEnew
            do iarec2 = 0, nanglerec
              areareclab(Zix, Nix, maxenrec, iarec2) = areareclab(Zix, Nix, maxenrec, iarec2) / renorm
              do inex = 0, maxex(Zix, Nix)
                ddxrec(Zix, Nix, inex, maxenrec, iarec2) = ddxrec(Zix, Nix, inex, maxenrec, iarec2) * renorm
              enddo
            enddo
          endif
!
! Calculate lab bins occupation for the recoil in the LAB.
! We assume the image in the LAB of a CM triangle is a triangle.
! This is theoretically wrong but should be a good approximation.
!
          call labsurface(Zcomp, Ncomp, type, 0, Ereclab11, cosreclab11, sinreclab11, Ereclab12, cosreclab12, sinreclab12, &
 &          Ereclab21, cosreclab21, sinreclab21, labsurf1, xlim1, ylim1, scovej1, scovrec1)
          call labsurface(Zcomp, Ncomp, type, 0, Ereclab21, cosreclab21, sinreclab21, Ereclab22, cosreclab22, sinreclab22, &
 &          Ereclab12, cosreclab12, sinreclab12, labsurf2, xlim2, ylim2, scovej2, scovrec2)
          labsurf = labsurf1 + labsurf2
          if (labsurf == 0.) cycle
          ddxLAB = fluxCM / labsurf
!
! Store recoil flux in lab ddx array
!
          sumddxrec = 0.
          iymax = 2 * nanglerec + 1
          iymaxp1 = iymax + 1
          do ix = xlim1(1), xlim1(2)
            do iy = ylim1(1), ylim1(2)
              iymod = mod(iy, iymaxp1)
              surfbin = areareclab(Zix, Nix, ix, iymod)
              if (surfbin /= 0.) then
                fluxadd = ddxLAB * scovrec1(ix, iymod)
                ddxrecadd = fluxadd / surfbin
                if (iymod <= nanglerec) then
                  iys = iymod
                else
                  iys = iymax - iymod
                endif
                ddxrec(Zix, Nix, nex, ix, iys) = ddxrec(Zix, Nix, nex, ix, iys) + ddxrecadd
                sumddxrec = sumddxrec + fluxadd
              endif
            enddo
          enddo
          do  ix = xlim2(1), xlim2(2)
            do  iy = ylim2(1), ylim2(2)
              iymod = mod(iy, iymaxp1)
              surfbin = areareclab(Zix, Nix, ix, iymod)
              if (surfbin /= 0.) then
                fluxadd = ddxLAB * scovrec2(ix, iymod)
                ddxrecadd = fluxadd / surfbin
                if (iymod <= nanglerec) then
                  iys = iymod
                else
                  iys = iymax - iymod
                endif
                ddxrec(Zix, Nix, nex, ix, iys) = ddxrec(Zix, Nix, nex, ix, iys) + ddxrecadd
                sumddxrec = sumddxrec + fluxadd
              endif
            enddo
          enddo
          ddxrectot(Zix, Nix, nex) = ddxrectot(Zix, Nix, nex) + sumddxrec
        enddo
!
! The renormalisation factor is calculated => we redo loop 160
!
        if (inorm == 0) then
          inorm = 1
          cycle
        else
          exit
        endif
      enddo
!
! We add the renormalized discrete component to the total ddxejlab array
!
      if (flaglabddx) then
        renorm = 1.0
        if (sumang /= xsCMlev .and. sumang >= 1.e-30) renorm = xsCMlev / sumang
        do ix = 1, iejlab(type)
          do iy = 0, nanglecont
            ddxejlab(type, ix, iy) = ddxejlab(type, ix, iy) + renorm * ddxejlabdis(type, ix, iy)
          enddo
        enddo
      endif
    enddo
  enddo
  return
end subroutine angdisrecoil
! Copyright A.J. Koning 2021
