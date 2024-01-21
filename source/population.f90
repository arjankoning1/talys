subroutine population
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Processing of pre-equilibrium spectra into population bins
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
! 2022-11-14: Current version
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl            ! single precision kind
! All global variables
!   numpar         ! number of particles
! Variables for output
!   flagcheck      ! flag for output of numerical checks
! Variables for preequilibrium
!   flag2comp      ! flag for two - component pre - equilibrium model
!   pespinmodel    ! model for pre - equilibrium or compound spin distribution
! Variables for OMP
!   flagomponly    ! flag to execute ONLY an optical model calculation
! Variables for energy grid
!   ebegin         ! first energy point of energy grid
!   egrid          ! outgoing energy grid
! Variables for energies
!   eend           ! last energy point of energy grid
!   Etotal         ! total energy of compound system (target + projectile)
!   flaggiant      ! flag for collective contribution from giant resonances
!   flagmulpre     ! flag for multiple pre - equilibrium calculation
!   mulpreZN       ! logical for multiple pre - equilibrium per nucleus
! Variables for excitation energy grid
!   deltaEx        ! excitation energy bin for population arrays
!   Ex             ! excitation energy
!   maxex          ! maximum excitation energy bin for residual nucleus
! Variables for nuclides
!   Nindex         ! neutron number index for residual nucleus
!   parskip        ! logical to skip outgoing particle
!   Zindex         ! charge number index for residual nucleus
! Variables for giant resonances
!   xsgrstate      ! smoothed giant resonance cross section per state
! Constants
!   parA           ! mass number of particle
!   parN           ! neutron number of particle
!   parname        ! name of particle
!   parZ           ! charge number of particle
! Variables for level density
!   Nlast          ! last discrete level
! Variables for masses
!   S              ! separation energy
! Variables for incident channel
!   preeqpop       ! pre - equilibrium population
!   preeqpopex     ! pre - equilibrium population
!   xsgr           ! total smoothed giant resonance cross section
!   xsgrtot        ! total smoothed giant resonance cross section
!   xspreeq        ! preeq. cross section per particle type and outgoing energy
!   xspreeqtot     ! preequilibrium cross section per particle type
! Variables for preequilibrium initialization
!   maxJph         ! maximal spin for particle - hole states
!   maxpar         ! maximal particle number
! Variables for preequilibrium
!   p0             ! initial particle number
!   pnu0           ! initial neutron number
!   ppi0           ! initial proton number
!   xspopph        ! population cross section p
!   xspopph2       ! population cross section p
!   xspreeqJP      ! preeq. cross section per particle type, outgoing energy, J, P
!   xsstep         ! preeq. cross section per particle type, stage and outgoing E
!   xsstep2        ! two - component preequilibrium cross section
!
! *** Declaration of local data
!
  implicit none
  integer   :: h                 ! help variable
  integer   :: hnu               ! neutron hole number
  integer   :: hpi               ! proton hole number
  integer   :: J                 ! spin of level
  integer   :: na1               ! help variable
  integer   :: na2               ! help variable
  integer   :: nb1               ! help variable
  integer   :: nb2               ! help variable
  integer   :: nen1              ! energy counter
  integer   :: nen2              ! energy counter
  integer   :: nex               ! excitation energy bin of compound nucleus
  integer   :: nexout            ! energy index for outgoing energy
  integer   :: Nix               ! neutron number index for residual nucleus
  integer   :: NL                ! last discrete level
  integer   :: p                 ! particle number
  integer   :: parity            ! parity
  integer   :: pc                ! Brosa parameter
  integer   :: pcnu              ! composite neutron particle number
  integer   :: pcpi              ! composite proton particle number
  integer   :: pnu               ! neutron particle number
  integer   :: ppi               ! proton particle number
  integer   :: type              ! particle type
  integer   :: Zix               ! charge number index for residual nucleus
  real(sgl) :: Ea1               ! help variable
  real(sgl) :: Ea2               ! start energy of local adjustment
  real(sgl) :: Eb1               ! help variable
  real(sgl) :: Eb2               ! end energy of local adjustment
  real(sgl) :: Ehigh             ! help variable
  real(sgl) :: Elow              ! help variable
  real(sgl) :: Eout              ! outgoing energy
  real(sgl) :: norm              ! normalization factor
  real(sgl) :: SS                ! separation energy
  real(sgl) :: xs                ! help variable
  real(sgl) :: xsa               ! help variable
  real(sgl) :: xsb               ! help variable
  real(sgl) :: xscheck(0:numpar) ! help variable to check total population
  real(sgl) :: xshigh            ! high energy cross section
  real(sgl) :: xslow             ! help variable
!
! ****** Process pre-equilibrium spectra into population arrays ********
!
! locate     : subroutine to find value in ordered table outgoing energy
!
! The pre-equilibrium cross sections have been calculated on the emission energy grid.
! They are interpolated on the excitation energy grids of the level populations (both for the total and the
! spin/parity-dependent cases) to enable futher decay of the residual nuclides.
!
  if (flagomponly) return
  do type = 0, 6
    if (parskip(type)) cycle
    Zix = Zindex(0, 0, type)
    Nix = Nindex(0, 0, type)
    SS = S(0, 0, type)
    NL = Nlast(Zix, Nix, 0)
    if (ebegin(type) >= eend(type)) cycle
!
! We include multiple pre-equilibrium emission only for neutrons and protons.
!
    if (flagmulpre .and. (type == 1 .or. type == 2)) mulpreZN(Zix, Nix) = .true.
!
! Loop over excitation energies.
! Determine the emission energy that corresponds with the excitation energy.
!
    do nexout = NL + 1, maxex(Zix, Nix)
      Eout = Etotal - SS - Ex(Zix, Nix, nexout)
      if (Eout < egrid(ebegin(type))) cycle
      Elow = Eout - 0.5 * deltaEx(Zix, Nix, nexout)
      call locate(egrid, ebegin(type), eend(type), Elow, nen1)
      na1 = nen1
      nb1 = nen1 + 1
      Ea1 = egrid(na1)
      Eb1 = min(egrid(nb1), Etotal - SS)
      Ehigh = Eout + 0.5 * deltaEx(Zix, Nix, nexout)
      call locate(egrid, ebegin(type), eend(type), Ehigh, nen2)
      na2 = nen2
      nb2 = nen2 + 1
      Ea2 = egrid(na2)
      Eb2 = min(egrid(nb2), Etotal - SS)
      xsa = xspreeq(type, na1)
      xsb = xspreeq(type, nb1)
!
! Add contribution from giant resonances
!
! pol1     : subroutine for polynomial interpolation of first order
!
      if (flaggiant) then
        xsa = xsa + xsgr(type, na1)
        xsb = xsb + xsgr(type, nb1)
      endif
      call pol1(Ea1, Eb1, xsa, xsb, Elow, xslow)
      xsa = xspreeq(type, na2)
      xsb = xspreeq(type, nb2)
      if (flaggiant) then
        xsa = xsa + xsgr(type, na2)
        xsb = xsb + xsgr(type, nb2)
      endif
!
! Determine interpolated value.
!
      call pol1(Ea2, Eb2, xsa, xsb, Ehigh, xshigh)
      xs = 0.5 * (xslow + xshigh) * deltaEx(Zix, Nix, nexout)
      if (xs < 1.e-30) xs = 0.
      preeqpopex(Zix, Nix, nexout) = xs
!
! If the pre-equilibrium spin distribution is chosen, the spectrum is interpolated on the spin/parity dependent population.
!
      if (pespinmodel >= 3) then
        do parity = - 1, 1, 2
          do J = 0, maxJph
            xsa = xspreeqJP(type, na1, J, parity) + xscollcontJP(type, J, parity, na1)
            xsb = xspreeqJP(type, nb1, J, parity) + xscollcontJP(type, J, parity, nb1)
            if (flaggiant .and. J <= 3) then
              xsa = xsa + xsgrstate(type, J, 1, na1) + xsgrstate(type, J, 2, na1)
              xsb = xsb + xsgrstate(type, J, 1, nb1) + xsgrstate(type, J, 2, nb1)
            endif
            call pol1(Ea1, Eb1, xsa, xsb, Elow, xslow)
            xsa = xspreeqJP(type, na2, J, parity) + xscollcontJP(type, J, parity, na2)
            xsb = xspreeqJP(type, nb2, J, parity) + xscollcontJP(type, J, parity, nb2)
            if (flaggiant .and. J <= 3) then
              xsa = xsa + xsgrstate(type, J, 1, na2) + xsgrstate(type, J, 2, na2)
              xsb = xsb + xsgrstate(type, J, 1, nb2) + xsgrstate(type, J, 2, nb2)
            endif
            call pol1(Ea2, Eb2, xsa, xsb, Ehigh, xshigh)
            xs = 0.5 * (xslow + xshigh) * deltaEx(Zix, Nix, nexout)
            if (xs < 1.e-30) xs = 0.
            preeqpop(Zix, Nix, nexout, J, parity) = xs
          enddo
        enddo
      endif
!
! A similar interpolation is done for multiple pre-equilibrium emission.
!
      if (mulpreZN(Zix, Nix)) then
!
! 1. One-component model
!
        if ( .not. flag2comp) then
          do pc = p0, maxpar
            p = pc - parA(type)
            h = pc - p0
            if (p < 0 .or. h < 0) cycle
            xsa = xsstep(type, pc, na1)
            xsb = xsstep(type, pc, nb1)
            call pol1(Ea1, Eb1, xsa, xsb, Elow, xslow)
            xsa = xsstep(type, pc, na2)
            xsb = xsstep(type, pc, nb2)
            call pol1(Ea2, Eb2, xsa, xsb, Ehigh, xshigh)
            xs = 0.5 * (xslow + xshigh) * deltaEx(Zix, Nix, nexout)
            if (xs < 1.e-30) xs = 0.
            xspopph(Zix, Nix, nexout, p, h) = xs
          enddo
        else
!
! 2. Two-component model
!
          do pcpi = ppi0, maxpar
            ppi = pcpi - parZ(type)
            hpi = pcpi - ppi0
            if (ppi < 0 .or. hpi < 0) cycle
            do pcnu = pnu0, maxpar
              pnu = pcnu - parN(type)
              hnu = pcnu - pnu0
              if (pnu < 0 .or. hnu < 0) cycle
              xsa = xsstep2(type, pcpi, pcnu, na1)
              xsb = xsstep2(type, pcpi, pcnu, nb1)
              call pol1(Ea1, Eb1, xsa, xsb, Elow, xslow)
              xsa = xsstep2(type, pcpi, pcnu, na2)
              xsb = xsstep2(type, pcpi, pcnu, nb2)
              call pol1(Ea2, Eb2, xsa, xsb, Ehigh, xshigh)
              xs = 0.5 * (xslow + xshigh) * deltaEx(Zix, Nix, nexout)
              if (xs < 1.e-30) xs = 0.
              xspopph2(Zix, Nix, nexout, ppi, hpi, pnu, hnu) = xs
            enddo
          enddo
        endif
      endif
    enddo
  enddo
!
! ***************** Correct for interpolation errors *******************
!
! Normalization of the pre-equilibrium population cross section.
! Due to interpolation, the part of the continuum population that comes from pre-equilibrium is not exactly equal to
! the total pre-equilibrium cross section.
! The normalization is done over the whole excitation energy range.
!
  do type = 0, 6
    if (parskip(type)) cycle
    Zix = Zindex(0, 0, type)
    Nix = Nindex(0, 0, type)
    NL = Nlast(Zix, Nix, 0)
    xscheck(type) = 0.
    do nex = NL + 1, maxex(Zix, Nix)
      xscheck(type) = xscheck(type) + preeqpopex(Zix, Nix, nex)
    enddo
    norm = 1.
    if (xscheck(type) /= 0.) norm = (xspreeqtot(type) + xsgrtot(type)) / xscheck(type)
    do nex = NL + 1, maxex(Zix, Nix)
      preeqpopex(Zix, Nix, nex) = preeqpopex(Zix, Nix, nex) * norm
      if (pespinmodel >= 3) then
        do parity = - 1, 1, 2
          do J = 0, maxJph
            preeqpop(Zix, Nix, nex, J, parity) = preeqpop(Zix, Nix, nex, J, parity) * norm
          enddo
        enddo
      endif
    enddo
  enddo
!
! ************************* Check of population ************************
!
! The difference that was already corrected above is illustrated by the following table.
!
  if (flagcheck) then
    write(*, '(/" ########## POPULATION CHECK ##########"/)')
    write(*, '(" Particle Pre-equilibrium Population"/)')
    do type = 0, 6
      if (parskip(type)) cycle
      xs = xspreeqtot(type)
      if (flaggiant) xs = xs + xsgrtot(type)
      write(*, '(1x, a8, 2(f12.5))') parname(type), xs, xscheck(type)
    enddo
  endif
  return
end subroutine population
! Copyright A.J. Koning 2021
