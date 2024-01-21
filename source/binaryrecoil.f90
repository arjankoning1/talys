subroutine binaryrecoil
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Recoil for binary reaction
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
!   numen2           ! maximum number of outgoing energies
!   numenrec         ! maximum number of recoil energies
! Variables for numerics
!   maxenrec         ! number of recoil energies
!   nanglecont       ! number of angles for continuum
!   nanglerec        ! number of recoil angles
! Variables for output
!   flagspec         ! flag for output of spectra
! Variables for basic reaction
!   flaglabddx       ! flag for calculation of DDX in LAB system
! Variables for energy grid
!   ebegin           ! first energy point of energy grid
!   Ebottom          ! bottom of outgoing energy bin
!   egrid            ! outgoing energy grid
!   Etop             ! top of outgoing energy bin
! Variables for energies
!   eend             ! last energy point of energy grid
!   Etotal           ! total energy of compound system (target + projectile)
! Variables for excitation energy grid
!   Ex               ! excitation energy
!   maxex            ! maximum excitation energy bin for residual nucleus
! Variables for binary reactions
!   binemissum       ! integrated binary emission spectrum
! Variables for binary emission spectra
!   xsbinemis        ! cross section for emission from first compound nucleus
!   xsbinemisad      ! angular distribution for emission from first compound nucleus
! Variables for compound nucleus from target
!   Exinc            ! excitation energy of entrance bin
! Variables for incident channel
!   xstotinc         ! total cross section (neutrons only) for incident channel
! Variables for nuclides
!   Nindex           ! neutron number index for residual nucleus
!   parskip          ! logical to skip outgoing particle
!   Zindex           ! charge number index for residual nucleus
! Constants
!   amu              ! atomic mass unit in MeV
!   parmass          ! mass of particle in a.m.u.
!   twopi            ! 2 * pi
! Variables for recoil initialization
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
!   dEejlab          ! width of ejectile lab bin
!   Eejcm1           ! Lower limit of CM ejectile energy bin
!   Eejcm2           ! Upper limit of CM ejectile energy bin
!   Eejlab11         ! LAB ejectile energy corresponding to (E1, ang1)
!   Eejlab12         ! LAB ejectile energy corresponding to (E1, ang2)
!   Eejlab21         ! LAB ejectile energy corresponding to (E2, ang1)
!   Eejlab22         ! LAB ejectile energy corresponding to (E2, ang2)
!   Eejlabmax        ! maximum energy of ejectile lab bin
!   Eejlabmin        ! minimum energy of ejectile lab bin
!   ejectmass        ! Ejectile mass
!   Erecinit         ! first compound nucleus recoil energy
!   Ereclab11        ! LAB recoil energy corresponding to (E1, ang1)
!   Ereclab12        ! LAB recoil energy corresponding to (E1, ang2)
!   Ereclab21        ! LAB recoil energy corresponding to (E2, ang1)
!   Ereclab22        ! LAB recoil energy corresponding to (E2, ang2)
!   Erecmax          ! minimal energy limit of recoil bin
!   Erecmin          ! minimal energy limit of recoil bin
!   iejlab           ! number of ejectile lab bins
!   irecinit         ! counter
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
!   S                ! separation energy
!
! *** Declaration of local data
!
  implicit none
  integer   :: nen                                   ! energy counter
  integer   :: Nix                                   ! neutron number index for residual nucleus
  integer   :: type                                  ! particle type
  integer   :: Zix                                   ! charge number index for residual nucleus
  real(sgl) :: compmass                              ! Composite system mass
  real(sgl) :: dEejcm                                ! CM width of ejectile emission energy bin
  real(sgl) :: Exrec1                                ! Recoil Excitation energy corresponding to Eejcm1
  real(sgl) :: Exrec2                                ! Recoil Excitation energy corresponding to Eejcm
  real(sgl) :: pejcm1                                ! Impulsion corresponding to Eejcm1
  real(sgl) :: pejcm2                                ! Impulsion corresponding to Eejcm2
  real(sgl) :: ratio                                 ! if 0<ratio<1 then x is between xl1 and xl2
  real(sgl) :: SS                                    ! separation energy
  real(sgl) :: sumejbinlab(0:6)                      ! Integrated ejectile binary spectrum in the LAB
  real(sgl) :: sumrecbinlab(0:6)                     ! Integrated recoil binary spectrum in the LAB
  integer   :: iex                                   ! Loop over excitation energy bins counter
  integer   :: iex1                                  ! Recoil excitation energy index corresponding to Eejcm1
  integer   :: iex2                                  ! Recoil excitation energy index corresponding to Eejcm1
  integer   :: iexmin                                ! minimum of iex1 and iex2
  integer   :: iexmax                                ! maximum of iex1 and iex2
  integer   :: numbinrec                             ! number of excitation energy bin covered
  integer   :: iang                                  ! running variable for angle
  real(sgl) :: ddxCM                                 ! Double differential cross section for a given CM  energy-angular b
  real(sgl) :: fluxCM                                ! Total flux that must be distributed in the LAB
  integer   :: irenorm                               ! index indicating if lab array renormalisation is needed
  real(sgl) :: Emaxi                                 ! help variable
  real(sgl) :: dEori                                 ! help variable
  real(sgl) :: dEnew                                 ! help variable
  real(sgl) :: renorm                                ! renormalisation factor
  real(sgl) :: labsurf1                              ! total surface covered in the LAB by the ejectile image  in the LAB
  real(sgl) :: labsurf2                              ! total surface covered in the LAB by the ejectile image  in the LAB
  real(sgl) :: labsurf                               ! Total covered surface in the LAB
  real(sgl) :: scovej1(numen2, 0:2*numangcont+1)     ! surface covered for each LAB energy-angle bin
  real(sgl) :: scovej2(numen2, 0:2*numangcont+1)     ! surface covered for each LAB energy-angle bin
  real(sgl) :: scovrec1(0:numenrec, 0:2*numangrec+1) ! surface covered for each LAB energy-angle bin
  real(sgl) :: scovrec2(0:numenrec, 0:2*numangrec+1) ! surface covered for each LAB energy-angle bin
  integer   :: xlim1(2)                              ! limits in the LAB ejectile energy grid
  integer   :: ylim1(2)                              ! limits in the LAB ejectile angular grid
  integer   :: xlim2(2)                              ! limits in the LAB ejectile energy grid
  integer   :: ylim2(2)                              ! limits in the LAB ejectile angular grid
  real(sgl) :: ddxLAB                                ! Double differential cross section for the LAB  energy-angular bin
  integer   :: iymax                                 ! maximum y-loop index
  integer   :: iymaxp1                               ! maximum y-loop index
  integer   :: ix                                    ! help variable
  integer   :: iy                                    ! help variable
  integer   :: iymod                                 ! help variable
  integer   :: iys                                   ! help variable
  real(sgl) :: ddxrecadd                             ! help variable
  real(sgl) :: surfbin                               ! LAB DDX area contribution
  real(sgl) :: fluxadd                               ! flux added to recoil bin
  real(sgl) :: sum                                   ! help variable
  real(sgl) :: sumenrec                              ! sum over energies
  real(sgl) :: sumenej                               ! sum over energies
  integer   :: iarec2                                ! counter
!
! ************************** Local variables ***************************
!
  if (flagspec) then
    do type = 0, 6
      if (parskip(type)) cycle
!
! *** Calculate recoil and light particle spectra in the LAB frame *****
!
! Initialisations (angcm=0 since first compound system decay)
!
      Zix = Zindex(0, 0, type)
      Nix = Nindex(0, 0, type)
      SS = S(0, 0, type)
      ejectmass = parmass(type) * amu
      recoilmass = nucmass(Zix, Nix) * amu
      compmass = nucmass(0, 0) * amu
      if (type /= 0) then
        ratio = recoilmass / (ejectmass * (ejectmass + recoilmass))
      endif
      vcm = sqrt(2. * Erecinit / compmass)
      angcm = 0.
      sumejbinlab(type) = 0.
      sumrecbinlab(type) = 0.
!
! Loop over ejectile CM outgoing energy bins (i.e. egrid(nen)) corresponding to the light particle spectrum component in the
! CM frame as calculated in subroutine binaryspectra.f which have been stored in xsbinemis(type,nen) and xsbinemisad(type,nen,iang)
!
      do nen = ebegin(type), eend(type)
        if (xsbinemis(type, nen) == 0.) cycle
        Eejcm1 = Ebottom(nen)
        if (Eejcm1 == 0.) Eejcm1 = 0.5 * egrid(1)
        if (Eejcm1 > Etotal - S(0, 0, type)) cycle
        Eejcm2 = Etop(nen)
        dEejcm = Eejcm2 - Eejcm1
        if (type /= 0) then
          pejcm1 = ejectmass * sqrt(2. * ratio * Eejcm1)
          pejcm2 = ejectmass * sqrt(2. * ratio * Eejcm2)
          vejcm1 = pejcm1 / ejectmass
          vejcm2 = pejcm2 / ejectmass
        else
          pejcm1 = Eejcm1
          pejcm2 = Eejcm2
        endif
!
! Deduce recoil excitation energies and velocities and find the excitation energy bin in which Exrec1 and Exrec2 are located
!
        Exrec1 = Exinc - SS - Eejcm1
        Exrec2 = Exinc - SS - Eejcm2
        vreccm1 = pejcm1 / recoilmass
        vreccm2 = pejcm2 / recoilmass
        iex1 = 0
        iex2 = 0
        do iex = 0, maxex(Zix, Nix)
          if (Exrec1 >= Ex(Zix, Nix, iex)) iex1 = iex
          if (Exrec2 >= Ex(Zix, Nix, iex)) iex2 = iex
        enddo
        iexmin = min(iex1, iex2)
        iexmax = max(iex1, iex2)
        numbinrec = iexmax - iexmin + 1
!
! A specific treatment should be done we cover more than one bin right now we do not go into the details
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
          ddxCM = xsbinemisad(type, nen, iang)
          fluxCM = ddxCM * dEejcm * dcosangcont(iang)
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
            call labsurface(0, 0, type, 1, Eejlab11, cosejlab11, sinejlab11, Eejlab12, cosejlab12, sinejlab12, Eejlab21, &
 &            cosejlab21, sinejlab21, labsurf1, xlim1, ylim1, scovej1, scovrec1)
            call labsurface(0, 0, type, 1, Eejlab21, cosejlab21, sinejlab21, Eejlab22, cosejlab22, sinejlab22, Eejlab12, &
 &            cosejlab12, sinejlab12, labsurf2, xlim2, ylim2, scovej2, scovrec2)
            labsurf = labsurf1 + labsurf2
            if (labsurf == 0.) cycle
            ddxLAB = fluxCM / labsurf
!
! Store ejectile flux in lab ddx array
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
                if (areaejlab(type, ix, iymod) > 0.) ddxejlab(type, ix, iys) = ddxejlab(type, ix, iys) + &
                  ddxLAB / areaejlab(type, ix, iymod) * scovej1(ix, iymod)
              enddo
            enddo
            do ix = xlim2(1), xlim2(2)
              do iy = ylim2(1), ylim2(2)
                iymod = mod(iy, iymaxp1)
                if (iymod <= nanglecont) then
                    iys = iymod
                  else
                    iys = iymax - iymod
                endif
                if (areaejlab(type, ix, iymod) > 0.) ddxejlab(type, ix, iys) = ddxejlab(type, ix, iys) + &
                  ddxLAB / areaejlab(type, ix, iymod) * scovej2(ix, iymod)
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
              do iex = 0, maxex(Zix, Nix)
                ddxrec(Zix, Nix, iex, maxenrec, iarec2) = ddxrec(Zix, Nix, iex, maxenrec, iarec2) * renorm
              enddo
            enddo
          endif
!
! Calculate lab bins occupation for the recoil in the LAB.
! We assume the image in the LAB of a CM triangle is a triangle.
! This is theoretically wrong but should be a good approximation.
!
          call labsurface(0, 0, type, 0, Ereclab11, cosreclab11, sinreclab11, Ereclab12, cosreclab12, sinreclab12, Ereclab21, &
 &          cosreclab21, sinreclab21, labsurf1, xlim1, ylim1, scovej1, scovrec1)
          call labsurface(0, 0, type, 0, Ereclab21, cosreclab21, sinreclab21, Ereclab22, cosreclab22, sinreclab22, Ereclab12, &
 &          cosreclab12, sinreclab12, labsurf2, xlim2, ylim2, scovej2, scovrec2)
          labsurf = labsurf1 + labsurf2
          if (labsurf == 0.) cycle
          ddxLAB = fluxCM / labsurf / numbinrec
!
! Store recoil flux in lab ddx array
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
                if (iymod <= nanglerec) then
                  iys = iymod
                else
                  iys = iymax - iymod
                endif
                do iex = iexmin, iexmax
                  ddxrec(Zix, Nix, iex, ix, iys) = ddxrec(Zix, Nix, iex, ix, iys) + ddxrecadd
                  ddxrectot(Zix, Nix, iex) = ddxrectot(Zix, Nix, iex) + fluxadd
                enddo
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
                if (iymod <= nanglerec) then
                  iys = iymod
                else
                  iys = iymax - iymod
                endif
                do iex = iexmin, iexmax
                  ddxrec(Zix, Nix, iex, ix, iys) = ddxrec(Zix, Nix, iex, ix, iys) + ddxrecadd
                  ddxrectot(Zix, Nix, iex) = ddxrectot(Zix, Nix, iex) + fluxadd
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
!
! Integrate the recoil and ejectile spectra in the LAB
!
     if (flaglabddx) then
       sumenej = 0.
       do nen = 1, iejlab(type)
         sum = 0.
         do iang = 0, nanglecont
           sum = sum + ddxejlab(type, nen, iang) * dcosangcont(iang)
         enddo
         if (nen <= iejlab(type)) sumenej = sumenej + sum * twopi * dEejlab(type, nen)
       enddo
       sumejbinlab(type) = sumenej
     endif
     sumenrec = 0.
     do iex = 0, maxex(Zix, Nix)
       sum = 0.
       do ix = 0, maxenrec
         do iy = 0, nanglerec
           sum = sum + ddxrec(Zix, Nix, iex, ix, iy) * areareclab(Zix, Nix, ix, iy)
         enddo
       enddo
       sumenrec = sumenrec + sum * twopi
     enddo
     if ((Zix == 0) .and. (Nix == 0)) then
       sumrecbinlab(type) = sumenrec - xstotinc
     else
       sumrecbinlab(type) = sumenrec
     endif
    enddo
  endif
!
! Renormalise the recoil and ejectile spectra
!
  do type = 0, 6
    if (parskip(type)) cycle
    Zix = Zindex(0, 0, type)
    Nix = Nindex(0, 0, type)
    if (flaglabddx) then
      if (sumejbinlab(type) <= 1.e-14) then
        renorm = 1.
      else
        renorm = binemissum(type) / sumejbinlab(type)
      endif
      do nen = 1, iejlab(type)
        do iang = 0, nanglecont
          ddxejlab(type, nen, iang) = ddxejlab(type, nen, iang) * renorm
        enddo
      enddo
    endif
    if (sumrecbinlab(type) <= 1.e-14) then
      renorm = 1.
    else
      renorm = binemissum(type) / sumrecbinlab(type)
    endif
    do iex = 0, maxex(Zix, Nix)
      if ((type == 0) .and. (iex == maxex(0, 0))) then
        ddxrectot(Zix, Nix, iex) = ddxrectot(Zix, Nix, iex)
      else
        ddxrectot(Zix, Nix, iex) = ddxrectot(Zix, Nix, iex) * renorm
      endif
      do ix = 0, maxenrec
        do iy = 0, nanglerec
          if ((type == 0) .and. (iex == maxex(0, 0)) .and. (ix == irecinit) .and. (iy == 0)) then
            ddxrec(Zix, Nix, iex, ix, iy) = ddxrec(Zix, Nix, iex, ix, iy)
          else
            ddxrec(Zix, Nix, iex, ix, iy) = ddxrec(Zix, Nix, iex, ix, iy) * renorm
          endif
        enddo
      enddo
    enddo
  enddo
  return
end subroutine binaryrecoil
! Copyright A.J. Koning 2021
