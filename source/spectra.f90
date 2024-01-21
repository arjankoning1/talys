subroutine spectra
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Creation of spectra
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
!   sgl              ! single precision kind
! All global variables
!   numendisc        ! number of discrete outgoing energies
! Variables for direct reactions
!   elwidth          ! width of elastic peak in MeV
! Variables for output
!   flagddx          ! flag for output of double - differential cross sections
! Variables for basic reaction
!   flagEchannel     ! flag for channel energy for emission spectrum
!   flagendf         ! flag for information for ENDF - 6 file
! Variables for numerics
!   maxN             ! maximal number of neutrons away from initial compound nucleus
!   maxZ             ! maximal number of protons away from initial compound nucleus
!   nangle           ! number of angles
!   nanglecont       ! number of angles for continuum
!   segment          ! help array for storing segment intersection points
! Variables for main input
!   Ainit            ! mass number of initial compound nucleus
!   flagnatural      ! flag for calculation of natural element
!   k0               ! index of incident particle
!   Ltarget          ! excited level of target
! Variables for incident channel
!   xsbinary         ! cross section from initial compound to residual nucleus
!   xsgr             ! total smoothed giant resonance cross section
!   xsparticle       ! total particle production cross section
!   xspreeq          ! preeq. cross section per particle typ and outgoing energye
! Variables for energy grid
!   ebegin           ! first energy point of energy grid
!   egrid            ! outgoing energy grid
!   Einc             ! incident energy in MeV
! Variables for energies
!   eend             ! last energy point of energy grid
!   eninccm          ! center - of - mass incident energy in MeV
!   eoutdis          ! outgoing energy of discrete state reaction
!   flagadd          ! flag for addition of discrete states to spectra flag
!   flagaddel        ! flag for addition of elastic peak to spectra
! Variables for multiple emission
!   Eaveragemul      ! average outgoing energy
!   xsfeed           ! cross section from compound to residual nucleus
!   xsmpreeq         ! multiple pre - equilibrium emission spectrum
!   xsmpreeqad       ! multiple preequilibrium angular distribution
! Variables for binary reactions
!   Eaveragebin      ! average outgoing energy
!   xsdisc           ! total cross section for discrete state
!   xselastot        ! total elastic cross section (shape + compound)
! Variables for binary emission spectra
!   xscomp           ! compound elastic cross section
!   xscompad         ! compound emission angular distribution
! Variables for spectra
!   buratio          ! break - up ratio
!   Eaverage         ! average outgoing energy
!   eendout          ! last energy point of energy grid
!   espec            ! outgoing energy grid
!   preeqratio       ! pre - equilibrium ratio
!   xsdiscout        ! smoothed angular distribution for discrete state
!   xsdiscoutad      ! smoothed angular distribution for discrete state
!   xscompout        ! compound emission angular distribution
!   xscompoutad      ! compound emission angular distribution
!   xsmpreeqout      ! multiple preequilibrium angular distribution
!   xsmpreeqoutad    ! multiple preequilibrium angular distribution
!   xspreeqbuout     ! preequilibrium cross section for breakup
!   xspreeqkiout     ! preequilibrium cross section for knockout and inelastic
!   xspreeqout       ! preequilibrium angular distribution per particle type an
!   xspreeqoutad     ! preequilibrium angular distribution per particle type
!   xspreeqpsout     ! preequilibrium cross section for pickup and stripping
!   xssumout         ! cross section summed over mechanisms
!   xssumoutad       ! angular distribution summed over mechanisms
! Variables for angular distributions
!   discad           ! discrete state angular distribution
! Variables for nuclides
!   parskip          ! logical to skip outgoing particle
! Variables for giant resonances
!   xsgrad           ! smoothed giant resonance angular distribution
! Constants
!   parA             ! mass number of particle
!   parN             ! neutron number of particle
!   parZ             ! charge number of particle
!   sqrttwopi        ! sqrt(2. * pi)
! Variables for level density
!   Nlast            ! last discrete level
! Variables for preequilibrium
!   xspreeqad        ! preequilibrium angular distribution per particle type and outg
!   xspreeqbu        ! preequilibrium cross section per particle type and outgoing energy for brea
!   xspreeqki        ! preequilibrium cross section per particle type and outgoing energy for knoc
!   xspreeqps        ! preequilibrium cross section per particle type and outgoing energy for pick
!
! *** Declaration of local data
!
  implicit none
  integer   :: fine        ! refinement factor for high-energy tail of spectrum
  integer   :: i           ! counter
  integer   :: iang        ! running variable for angle
  integer   :: iangdisc    ! counter for discrete angle
  integer   :: nen         ! energy counter
  integer   :: nend        ! help variable
  integer   :: nenout      ! counter for outgoing energy
  integer   :: nhigh       ! help variable
  integer   :: NL          ! last discrete level
  integer   :: type        ! particle type
  integer   :: Zcomp       ! proton number index for compound nucleus
  integer   :: Ncomp       ! neutron number index for compound nucleus
  real(sgl) :: Ares        ! mass number of residual nucleus
  real(sgl) :: convfac1    ! conversion factor for reference system
  real(sgl) :: convfac2    ! conversion factor for reference system
  real(sgl) :: convfac3    ! conversion factor for reference system
  real(sgl) :: diswidth    ! width of discrete level peak
  real(sgl) :: Ea          ! start energy of local adjustment
  real(sgl) :: Eb          ! end energy of local adjustment
  real(sgl) :: Ehigh       ! help variable
  real(sgl) :: Elast       ! help variable
  real(sgl) :: Eout        ! outgoing energy
  real(sgl) :: fac1        ! help variable
  real(sgl) :: fac2        ! help variable
  real(sgl) :: factor      ! multiplication factor
  real(sgl) :: gauss       ! Gaussian contribution
  real(sgl) :: xssum       ! help variable
  real(sgl) :: Eaveragesum ! help variable
!
! ********** Add smoothed discrete cross sections to spectra ***********
!
! locate            : subroutine to find value in ordered table
!
  Eaverage = 0.
  xsdisc(k0, Ltarget) = xselastot
  do type = 0, 6
    if (parskip(type)) cycle
    if (xsparticle(type) == 0.) cycle
    NL = Nlast(parZ(type), parN(type), 0)
    if (flagadd .or. type >= 2) then
      if (flagnatural) then
        Elast = Einc - 6. - elwidth
      else
        Elast = eoutdis(type, NL) - elwidth
      endif
      Elast = max(Elast, 0.)
      call locate(egrid, ebegin(type), eend(type), Elast, nend)
    else
      nend = eend(type)
      eendout(type) = eend(type)
    endif
    convfac1 = 1.
    convfac2 = 0.
    convfac3 = 0.
    if (flagEchannel) then
      Ares = real(Ainit - parA(type))
      convfac1 = Ares / Ainit
      convfac2 = real(parA(type)) / Ainit * real(parA(k0)) / Ainit * Einc
      if (convfac1 > 0..and.convfac2 > 0.) convfac3 = 2. * sqrt(convfac1 * convfac2)
    endif
    do nen = ebegin(type), nend
      espec(type, nen) = convfac1 * egrid(nen) + convfac2 + convfac3 * sqrt(egrid(nen))
      xsdiscout(type, nen) = xsgr(type, nen)
      xspreeqout(type, nen) = xspreeq(type, nen)
      xspreeqpsout(type, nen) = xspreeqps(type, nen)
      xspreeqkiout(type, nen) = xspreeqki(type, nen)
      xspreeqbuout(type, nen) = xspreeqbu(type, nen)
      xsmpreeqout(type, nen) = xsmpreeq(type, nen)
      xscompout(type, nen) = xscomp(type, nen)
      if (flagddx) then
        do iang = 0, nanglecont
          xsdiscoutad(type, nen, iang) = xsgrad(type, nen, iang)
          xspreeqoutad(type, nen, iang) = xspreeqad(type, nen, iang)
          xsmpreeqoutad(type, nen, iang) = xsmpreeqad(type, nen, iang)
          xscompoutad(type, nen, iang) = xscompad(type, nen, iang)
        enddo
      endif
    enddo
    if (flagadd .or. type >= 2) then
      if (flagendf) then
        fine = 2
      else
        fine = 10
      endif
      Ehigh = egrid(eend(type)) - egrid(nend)
      nhigh = int(fine * segment * Ehigh)
      nhigh = min(nhigh, numendisc - 8 * fine)
      eendout(type) = nend + nhigh + 8 * fine
      do nenout = nend + 1, eendout(type)
        espec(type, nenout) = egrid(nend) + (nenout - nend) / (segment * real(fine))
        Eout = espec(type, nenout)
        espec(type, nenout) = convfac1 * Eout + convfac2 + convfac3 * sqrt(Eout)
        call locate(egrid, nend, eend(type), Eout, nen)
        nen = min(nen, numen-1)
        Ea = Eout - egrid(nen)
        Eb = egrid(nen + 1) - egrid(nen)
        if (Ea < Eb) then
          factor = Ea / Eb
        else
          factor = 1.
        endif
        xsdiscout(type, nenout) = xsgr(type, nen) + factor * (xsgr(type, nen + 1) - xsgr(type, nen))
        xspreeqout(type, nenout) = xspreeq(type, nen) + factor * (xspreeq(type, nen + 1) - xspreeq(type, nen))
        xspreeqpsout(type, nenout) = xspreeqps(type, nen) + factor * (xspreeqps(type, nen + 1) - xspreeqps(type, nen))
        xspreeqkiout(type, nenout) = xspreeqki(type, nen) + factor * (xspreeqki(type, nen + 1) - xspreeqki(type, nen))
        xspreeqbuout(type, nenout) = xspreeqbu(type, nen) + factor * (xspreeqbu(type, nen + 1) - xspreeqbu(type, nen))
        xsmpreeqout(type, nenout) = xsmpreeq(type, nen) + factor * (xsmpreeq(type, nen + 1) - xsmpreeq(type, nen))
        xscompout(type, nenout) = xscomp(type, nen) + factor * (xscomp(type, nen + 1) - xscomp(type, nen))
        if (flagddx) then
          do iang = 0, nanglecont
            xsdiscoutad(type, nenout, iang) = xsgrad(type, nen, iang) + &
              factor * (xsgrad(type, nen + 1, iang) - xsgrad(type, nen, iang))
            xspreeqoutad(type, nenout, iang) = xspreeqad(type, nen, iang) + factor * (xspreeqad(type, nen + 1, iang) - &
 &            xspreeqad(type, nen, iang))
            xsmpreeqoutad(type, nenout, iang) = xsmpreeqad(type, nen, iang) + &
 &            factor * (xsmpreeqad(type, nen + 1, iang) - xsmpreeqad(type, nen, iang))
            xscompoutad(type, nenout, iang) = xscompad(type, nen, iang) + factor * (xscompad(type, nen + 1, iang) - &
 &            xscompad(type, nen, iang))
          enddo
        endif
        do i = 0, NL
          if (i == Ltarget .and. type == k0 .and. k0 > 1) cycle
          if (i == Ltarget .and. .not. flagaddel) cycle
          if (eoutdis(type, i) <= 0.) cycle
          diswidth = elwidth * (eoutdis(type, i) / eninccm) **1.5
          fac1 = 1. / (diswidth * sqrttwopi)
          fac2 = 1. / (2. * diswidth **2)
          gauss = fac1 * exp( - (Eout - eoutdis(type, i)) **2 * fac2)
          xsdiscout(type, nenout) = xsdiscout(type, nenout) + gauss * xsdisc(type, i)
          if (flagddx) then
            do iang = 0, nanglecont
              iangdisc = int( real(iang * nangle) / nanglecont)
              xsdiscoutad(type, nenout, iang) = xsdiscoutad(type, nenout, iang) + &
                gauss * discad(type, i, iangdisc)
            enddo
          endif
        enddo
      enddo
    endif
!
! ************************ Create total spectra ************************
!
    do nen = ebegin(type), eendout(type)
      xssumout(type, nen) = xspreeqout(type, nen) + xsmpreeqout(type, nen) + xscompout(type, nen)
      if (xssumout(type, nen) /= 0.) then
        preeqratio(type, nen) = (xspreeqout(type, nen) + xsmpreeqout(type, nen)) / xssumout(type, nen)
        if (k0 > 2 .or. type > 2) then
          buratio(type, nen) = xspreeqbuout(type, nen) / xssumout(type, nen)
        else
          buratio(type, nen) = 0.
        endif
      else
        preeqratio(type, nen) = 0.
        buratio(type, nen) = 0.
      endif
      xssumout(type, nen) = xssumout(type, nen) + xsdiscout(type, nen)
      if (flagddx) then
        do iang = 0, nanglecont
          xssumoutad(type, nen, iang) = xsdiscoutad(type, nen, iang) + xspreeqoutad(type, nen, iang) + &
            xsmpreeqoutad(type, nen, iang) + xscompoutad(type, nen, iang)
        enddo
      endif
    enddo
!
! ************************ Average emission energy *********************
!
    xssum = xsbinary(type)
    Eaveragesum = Eaveragebin(type) * xsbinary(type)
    do Zcomp = 0, maxZ
      do Ncomp = 0, maxN
        if (.not.flaginitpop .and. Zcomp == 0 .and. Ncomp == 0) cycle
        xssum = xssum + xsfeed(Zcomp, Ncomp, type)
        Eaveragesum = Eaveragesum + Eaveragemul(Zcomp, Ncomp, type) * xsfeed(Zcomp, Ncomp, type)
      enddo
    enddo
    if (xssum > 0.) Eaverage(type) = Eaveragesum / xssum
  enddo
  return
end subroutine spectra
! Copyright A.J. Koning 2021
