subroutine comprecoil(Zcomp, Ncomp, nex, type, nexout, nenbeg, nenend)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Recoils from compound decay
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
!   sgl              ! single precision kind
! All global variables
!   numangcont       ! maximum number of angles for continuum
!   numangrec        ! maximum number of recoil angles
!   numen            ! maximum number of outgoing energies
!   numen2           ! maximum number of outgoing energies
!   numenrec         ! maximum number of recoil energies
! Variables for numerics
!   maxenrec         ! number of recoil energies
!   nanglecont       ! number of angles for continuum
!   nanglerec        ! number of recoil angles
! Variables for spectra
!   compspect        ! compound part of spectrum
!   preeqspect       ! multiple pre - equilibrium part of spectrum
! Variables for basic reaction
!   flaglabddx       ! flag for calculation of DDX in LAB system
!   flagrecoilav     ! flag for average velocity in recoil calculation
! Variables for energy grid
!   Ebottom          ! bottom of outgoing energy bin
!   egrid            ! outgoing energy grid
!   Einc             ! incident energy in MeV
!   Etop             ! top of outgoing energy bin
! Variables for nuclides
!   Nindex           ! neutron number index for residual nucleus
!   Zindex           ! charge number index for residual nucleus
! Constants
!   amu              ! atomic mass unit in MeV
!   deg2rad          ! conversion factor for degrees to radians
!   fourpi           ! 4. * pi
!   parmass          ! mass of particle in a.m.u.
! Variables for recoil initialization
!   angcontmax       ! maximum of angular bin
!   angcontmin       ! minimum of angular bin
!   cosangcontmax    ! cosine of maximum of angular bin
!   cosangcontmin    ! cosine of minimum of angular bin
!   dcosangcont      ! width of cosine bin
!   sinangcontmax    ! sine of maximum of angular bin
!   sinangcontmin    ! sine of minimum of angular bin
! Variables for recoil
!   angcm            ! CM angle with respect to the LAB
!   areaejlab        ! Total surface of LAB ddx bins
!   areareclab       ! Total surface of LAB ddx bins
!   cosejcm1         ! CM recoil energy cosine corresponding to (E1, ang1)
!   cosejcm2         ! CM recoil energy cosine corresponding to (E1, ang2)
!   cosejlab11       ! LAB ejectile angle cosine corresponding to (E1, ang1)
!   cosejlab12       ! LAB ejectile angle cosine corresponding to (E1, ang2)
!   cosejlab21       ! LAB ejectile angle cosine corresponding to (E2, ang1)
!   cosejlab22       ! LAB ejectile angle cosine corresponding to (E2, ang2)
!   cosreclab11      ! LAB recoil angle cosine corresponding to (E1, ang1)
!   cosreclab12      ! LAB recoil angle cosine corresponding to (E1, ang2)
!   cosreclab21      ! LAB recoil angle cosine corresponding to (E2, ang1)
!   cosreclab22      ! LAB recoil angle cosine corresponding to (E2, ang2)
!   ddxejlab         ! array containing the ddx spectrum of light part
!   ddxrec           ! array containing the lab double differential xs
!   ddxrectot        ! array containing the total recoil flux in a giv
!   Eejcm1           ! Lower limit of CM ejectile energy bin
!   Eejcm2           ! Upper limit of CM ejectile energy bin
!   Eejlab11         ! LAB ejectile energy corresponding to (E1, ang1)
!   Eejlab12         ! LAB ejectile energy corresponding to (E1, ang2)
!   Eejlab21         ! LAB ejectile energy corresponding to (E2, ang1)
!   Eejlab22         ! LAB ejectile energy corresponding to (E2, ang2)
!   ejectmass        ! Ejectile mass
!   Erec             ! recoil energy
!   Ereclab11        ! LAB recoil energy corresponding to (E1, ang1)
!   Ereclab12        ! LAB recoil energy corresponding to (E1, ang2)
!   Ereclab21        ! LAB recoil energy corresponding to (E2, ang1)
!   Ereclab22        ! LAB recoil energy corresponding to (E2, ang2)
!   recoilmass       ! Recoil mass
!   sinejcm1         ! CM recoil energy cosine corresponding to (E2, ang1)
!   sinejcm2         ! CM recoil energy cosine corresponding to (E2, ang2)
!   sinejlab11       ! LAB ejectile angle sine corresponding to (E1, ang1)
!   sinejlab12       ! LAB ejectile angle sine corresponding to (E1, ang2)
!   sinejlab21       ! LAB ejectile angle sine corresponding to (E2, ang1)
!   sinejlab22       ! LAB ejectile angle sine corresponding to (E2, ang2)
!   sinreclab11      ! LAB recoil angle sine corresponding to (E1, ang1)
!   sinreclab12      ! LAB recoil angle sine corresponding to (E1, ang2)
!   sinreclab21      ! LAB recoil angle sine corresponding to (E2, ang1)
!   sinreclab22      ! LAB recoil angle sine corresponding to (E2, ang2)
!   vcm              ! Compound nucleus velocity
!   vejcm1           ! velocity corresponding to Eejcm1
!   vejcm2           ! velocity corresponding to Eejcm2
!   vreccm1          ! Recoil velocity corresponding to Eejcm1
!   vreccm2          ! Recoil velocity corresponding to Eejcm2
! Variables for masses
!   nucmass          ! mass of nucleus
!
! *** Declaration of local data
!
  implicit none
  integer   :: Ncomp                                 ! neutron number index for compound nucleus
  integer   :: nen                                   ! energy counter
  integer   :: nenbeg                                ! help variable
  integer   :: nenend                                ! help variable
  integer   :: nex                                   ! excitation energy bin of compound nucleus
  integer   :: nexout                                ! energy index for outgoing energy
  integer   :: Nix                                   ! neutron number index for residual nucleus
  integer   :: type                                  ! particle type
  integer   :: Zcomp                                 ! proton number index for compound nucleus
  integer   :: Zix                                   ! charge number index for residual nucleus
  real(sgl) :: ang                                   ! angle
  real(sgl) :: Eout                                  ! outgoing energy
  real(sgl) :: xsgridad(numen, 0:numangcont)         ! angular distribution
  integer   :: iymax                                 ! maximum y-loop index
  integer   :: iymaxp1                               ! maximum y-loop index
  integer   :: ix                                    ! help variable
  integer   :: iy                                    ! help variable
  integer   :: ierec                                 ! counter
  integer   :: iarec                                 ! counter
  integer   :: ierecbeg                              ! begin of energy count
  integer   :: ierecend                              ! begin of energy count
  real(sgl) :: compmass                              ! Composite system mass
  real(sgl) :: ratio                                 ! if 0<ratio<1 then x is between xl1 and xl2
  real(sgl) :: fracCM                                ! fraction
  real(sgl) :: fracCMloc(0:numenrec)                 ! fraction
  real(sgl) :: vcmloc(0:numenrec)                    ! C.M. velocity
  real(sgl) :: dEejcm                                ! CM width of ejectile emission energy bin
  real(sgl) :: pejcm1                                ! Impulsion corresponding to Eejcm1
  real(sgl) :: pejcm2                                ! Impulsion corresponding to Eejcm2
  integer   :: iang                                  ! running variable for angle
  real(sgl) :: fluxCM                                ! Total flux that must be distributed in the LAB
  real(sgl) :: kalbach                               ! Kalbach function
  real(sgl) :: labs1                                 ! recoil energy
  real(sgl) :: labs2                                 ! recoil energy
  real(sgl) :: labsurf                               ! Total covered surface in the LAB
  real(sgl) :: scovej1(1:numen2, 0:2*numangcont+1)   ! surface covered for each LAB energy-angle bin
  real(sgl) :: scovej2(1:numen2, 0:2*numangcont+1)   ! surface covered for each LAB energy-angle bin
  real(sgl) :: scovrec1(0:numenrec, 0:2*numangrec+1) ! surface covered for each LAB energy-angle bin
  real(sgl) :: scovrec2(0:numenrec, 0:2*numangrec+1) ! surface covered for each LAB energy-angle bin
  integer   :: xlim1(2)                              ! limits in the LAB ejectile energy grid
  integer   :: ylim1(2)                              ! limits in the LAB ejectile angular grid
  integer   :: xlim2(2)                              ! limits in the LAB ejectile energy grid
  integer   :: ylim2(2)                              ! limits in the LAB ejectile angular grid
  real(sgl) :: ddxLAB                                ! Double differential cross section for the LAB  energy-angular bin
  integer   :: iymod                                 ! help variable
  integer   :: iys                                   ! help variable
  real(sgl) :: sumddxr                               ! help variable
  real(sgl) :: ddxrecadd                             ! help variable
  real(sgl) :: ddxejadd                              ! addition to double differential cross section for the LAB
  real(sgl) :: surfbin                               ! LAB DDX area contribution
  real(sgl) :: fluxadd                               ! flux added to recoil bin
  real(sgl) :: vcmav                                 ! mean center of mass velocity for a bin
!
! *** Calculate recoil and light particle spectra in the LAB frame *****
!
! We are decaying from a compound nucleus bin with index nex and excitation energies characteristics with a given ejectile
! (i.e. we are still within the type loop) in a residual nucleus bin nexout characterised by an excitation energy
! Exout=Ex(Zix,Nix,nexout) and the bin limits Ex1min (lower boundary of residual bin) and Ex1plus (upper boundary of residual bin).
! The corresponding ejectile emission energies have been deduced and are given by the egrid(ieject) values were ieject is between
! nenbeg and nenend emax (maximal emission energy) and emin (minimal emission energy) and Eout (emission energy)
!
! The bin that decays is also described by a distribution of kinetic energies and angles in the LAB frame given by
! the combination of the elements of the array ddxrec(Zcomp,Ncomp,nex,ierec,iarec) which gives the spectrum as function of
! the recoil energies ierec and lab angle iarec and the elements of the array ddxrectot(Zcomp,Ncomp,nex) which gives the
! integrated cross section distributed over the energy-angular bin.
! In other words, ddxrectot(Zcomp,Ncomp,nex) is the sum of the elements of ddxrectot(Zcomp,Ncomp,nex,ierec,iarec)
! times the attached lab surface element areareclab(Zcomp,Ncomp,ierec,iarec).
! Therefore, for each decaying bin, the fraction of the decay which corresponds to a nucleus moving with kinetic energy
! E=Erec(Zcomp,Ncomp,ierec) and the angular direction in the LAB is the ratio of
! ddxrec(Zcomp,Ncomp,nex,ierec,iarec)*areareclab(Zcomp,Ncomp,ierec,iarec) with ddxrectot(Zcomp,Ncomp,nex,ierec,iarec)
!
! Loop over the compound nucleus kinetic energies
!
  do nen = nenbeg, nenend
    Eout = egrid(nen)
    do iang = 0, nanglecont
      ang = 0.5 * (angcontmin(iang) + angcontmax(iang)) * deg2rad
      xsgridad(nen, iang) = compspect(nen) / fourpi + preeqspect(nen) * kalbach(type, Einc, Eout, ang)
    enddo
  enddo
  Zix = Zindex(Zcomp, Ncomp, type)
  Nix = Nindex(Zcomp, Ncomp, type)
  ejectmass = parmass(type) * amu
  recoilmass = nucmass(Zix, Nix) * amu
  if (type /= 0) then
    ratio = recoilmass / (ejectmass * (ejectmass + recoilmass))
  endif
  compmass = nucmass(Zcomp, Ncomp) * amu
  sumddxr = 0.
!
! Loop over recoil angles and velocities to calculate the fraction of the bin population with a given angle-velocity direction
!
  angcm = 0.
  vcmav = 0.
  do ierec = 0, maxenrec
    vcmloc(ierec) = sqrt(2. * Erec(Zcomp, Ncomp, ierec) / compmass)
    fracCMloc(ierec) = 0.
    do iarec = 0, nanglerec
!
! Calculate fraction of the decay having this CM velocity characteristics (i.e. fracCMloc)
!
      if (ddxrectot(Zcomp, Ncomp, nex) > 0.) fracCMloc(ierec) = fracCMloc(ierec) + &
        ddxrec(Zcomp, Ncomp, nex, ierec, iarec) * areareclab(Zcomp, Ncomp, ierec, iarec) / &
        ddxrectot(Zcomp, Ncomp, nex)
    enddo
    if (flagrecoilav) vcmav = vcmav + fracCMloc(ierec) * vcmloc(ierec)
  enddo
  if (flagrecoilav .and. type > 0) then
    ierecbeg = 1
    ierecend = 1
    vcmloc(1) = vcmav
    fracCMloc(1) = 1.
  else
    ierecbeg = 0
    ierecend = maxenrec
  endif
!
! Determine CM outgoing energies corresponding to the decay from nex to nexout
!
  do ierec = ierecbeg, ierecend
    vcm = vcmloc(ierec)
    fracCM = fracCMloc(ierec)
    do nen = nenbeg, nenend
      if (compspect(nen) + preeqspect(nen) == 0.) cycle
      Eejcm1 = Ebottom(nen)
      if (Eejcm1 == 0.) Eejcm1 = 0.5 * egrid(1)
      Eejcm2 = Etop(nen)
      dEejcm = Eejcm2 - Eejcm1
      if (type /= 0) then
        pejcm1 = ejectmass * sqrt(2 * ratio * Eejcm1)
        pejcm2 = ejectmass * sqrt(2 * ratio * Eejcm2)
        vejcm1 = pejcm1 / ejectmass
        vejcm2 = pejcm2 / ejectmass
      else
        pejcm1 = Eejcm1
        pejcm2 = Eejcm2
      endif
!
! Deduce recoil excitation energies and velocities corresponding to the ejectile energies (not perfect but should be good enough)
!
      vreccm1 = pejcm1 / recoilmass
      vreccm2 = pejcm2 / recoilmass
!
! loop over ejectile angles in the CM frame
!
      do iang = 0, nanglecont
        cosejcm1 = cosangcontmin(iang)
        cosejcm2 = cosangcontmax(iang)
        sinejcm1 = sinangcontmin(iang)
        sinejcm2 = sinangcontmax(iang)
!
! Total flux that must be spread in the LAB frame
!
        fluxCM = xsgridad(nen, iang) * dEejcm * dcosangcont(iang) * fracCM
        if (fluxCM == 0.) cycle
!
! EJECTILE and RECOIL TREATMENT
!
! The ejectile as well as recoil LAB angles and energies corresponding to the CM points are deduced
!
        call cm2lab(type)
!
! Calculate lab bins occupation for the ejectile in the LAB.
! We assume the image in the LAB of a CM triangle is a triangle.
! This is theoretically wrong but should be a good approximation.
!
        if (flaglabddx) then
         call labsurface(Zcomp, Ncomp, type, 1, Eejlab11, cosejlab11, sinejlab11, Eejlab12, cosejlab12, sinejlab12, Eejlab21, &
 &         cosejlab21, sinejlab21, labs1, xlim1, ylim1, scovej1, scovrec1)
         call labsurface(Zcomp, Ncomp, type, 1, Eejlab21, cosejlab21, sinejlab21, Eejlab22, cosejlab22, sinejlab22, Eejlab12, &
 &         cosejlab12, sinejlab12, labs2, xlim2, ylim2, scovej2, scovrec2)
         labsurf = labs1 + labs2
         if (labsurf == 0.) cycle
         ddxLAB = fluxCM / labsurf
!
! store ejectile flux in lab ddx array
!
! Deduce the real ejectile angles by coupling with compound nucleus moving directions in the LAB
!
          iymax = 2 * nanglecont + 1
          iymaxp1 = iymax + 1
          do ix = xlim1(1), xlim1(2)
            do iy = ylim1(1), ylim1(2)
              iymod = mod(iy, iymaxp1)
              surfbin = areaejlab(type, ix, iymod)
              if (surfbin /= 0.) then
                fluxadd = ddxLAB * scovej1(ix, iymod)
                ddxejadd = fluxadd / surfbin
                if (iymod <= nanglecont) then
                  iys = iymod
                else
                  iys = iymax - iymod
                endif
                ddxejlab(type, ix, iys) = ddxejlab(type, ix, iys) + ddxejadd
              endif
            enddo
          enddo
          do ix = xlim2(1), xlim2(2)
            do iy = ylim2(1), ylim2(2)
              iymod = mod(iy, iymaxp1)
              surfbin = areaejlab(type, ix, iymod)
              if (surfbin /= 0.) then
                fluxadd = ddxLAB * scovej2(ix, iymod)
                ddxejadd = fluxadd / surfbin
                if (iymod <= nanglecont) then
                  iys = iymod
                else
                  iys = iymax - iymod
                endif
                ddxejlab(type, ix, iys) = ddxejlab(type, ix, iys) + ddxejadd
              endif
            enddo
          enddo
        endif
!
! Calculate lab bins occupation for the recoil in the LAB
! We assume the image in the LAB of a CM triangle is a triangle
! This is theoretically wrong but should be a good approximation
!
        call labsurface(Zcomp, Ncomp, type, 0, Ereclab11, cosreclab11, sinreclab11, Ereclab12, cosreclab12, sinreclab12, &
 &        Ereclab21, cosreclab21, sinreclab21, labs1, xlim1, ylim1, scovej1, scovrec1)
        call labsurface(Zcomp, Ncomp, type, 0, Ereclab21, cosreclab21, sinreclab21, Ereclab22, cosreclab22, sinreclab22, &
 &        Ereclab12, cosreclab12, sinreclab12, labs2, xlim2, ylim2, scovej2, scovrec2)
        labsurf = labs1 + labs2
        if (labsurf == 0.) cycle
        ddxLAB = fluxCM / labsurf
!
! store recoil flux in lab ddx array
!
! Deduce the real ejectile angles by coupling with compound nucleus moving directions in the LAB
!
        iymax = 2 * nanglerec + 1
        iymaxp1 = iymax + 1
        do ix = xlim1(1), xlim1(2)
          do iy = ylim1(1), ylim1(2)
            iymod = mod(iy, iymaxp1)
            surfbin = areareclab(Zix, Nix, ix, iymod)
            if (surfbin /= 0.) then
              fluxadd = ddxLAB * scovrec1(ix, iymod)
              ddxrecadd = fluxadd / surfbin
              sumddxr = sumddxr + fluxadd
              if (iymod <= nanglerec) then
                iys = iymod
              else
                iys = iymax - iymod
              endif
              ddxrec(Zix, Nix, nexout, ix, iys) = ddxrec(Zix, Nix, nexout, ix, iys) + ddxrecadd
            endif
          enddo
        enddo
        do ix = xlim2(1), xlim2(2)
          do iy = ylim2(1), ylim2(2)
            iymod = mod(iy, iymaxp1)
            surfbin = areareclab(Zix, Nix, ix, iymod)
            if (surfbin /= 0.) then
              fluxadd = ddxLAB * scovrec2(ix, iymod)
              ddxrecadd = fluxadd / surfbin
              sumddxr = sumddxr + fluxadd
              if (iymod <= nanglerec) then
                iys = iymod
              else
                iys = iymax - iymod
              endif
              ddxrec(Zix, Nix, nexout, ix, iys) = ddxrec(Zix, Nix, nexout, ix, iys) + ddxrecadd
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
  ddxrectot(Zix, Nix, nexout) = ddxrectot(Zix, Nix, nexout) + sumddxr
  return
end subroutine comprecoil
! Copyright A.J. Koning 2021
