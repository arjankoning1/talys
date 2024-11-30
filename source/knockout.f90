subroutine knockout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Contribution of knockout and complex inelastic reactions
!
! Author    : Arjan Koning and Vivian Dimitriou
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl          ! single precision kind
! Variables for preequilibrium
!   Cknock       ! adjustable parameter for knockout reactions
!   pairmodel    ! model for preequilibrium pairing energy
! Variables for main input
!   Atarget      ! mass number of target nucleus
!   k0           ! index of incident particle
!   Ntarget      ! neutron number of target nucleus
!   Ztarget      ! charge number of target nucleus
! Variables for energy grid
!   deltaE       ! energy bin around outgoing energies
! Variables for energy grid
!   ebegin       ! first energy point of energy grid
!   egrid        ! outgoing energy grid
!   Einc         ! incident energy in MeV
! Variables for energies
!   eend         ! last energy point of energy grid
!   eninccm      ! center - of - mass incident energy in MeV
!   Etotal       ! total energy of compound system (target + projectile)
! Variables for incident channel
!   xsreacinc    ! reaction cross section for incident channel
! Variables for inverse channel data
!   xsreac       ! reaction cross section
! Variables for nuclides
!   coulbar      ! Coulomb barrier
!   parskip      ! logical to skip outgoing particle
!   Q            ! Q - value
! Constants
!   parA         ! mass number of particle
!   parmass      ! mass of particle in a.m.u.
!   parN         ! neutron number of particle
!   parspin      ! spin of particle
!   parZ         ! charge number of particle
! Variables for ECIS
!   projmass     ! mass of projectile
! Variables for masses
!   S            ! separation energy
! Variables for preequilibrium
!   xspreeqki    ! preequilibrium cross section per particle type and outgoing energy for knock-out
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: key           ! keyword
  logical           :: flaginel      ! flag for inelastic scattering
  logical           :: flagknock     ! flag for knockout
  integer           :: ndelta        ! number of transferred particles
  integer           :: nen           ! energy counter
  integer           :: type          ! particle type
  integer           :: type2         ! particle type
  real(sgl)         :: AKO           ! Pauli correction factor
  real(sgl)         :: Ccl           ! adjustment factor
  real(sgl)         :: Ckn           ! adjustable parameter for knockout reactions
  real(sgl)         :: dE            ! help variable
  real(sgl)         :: denom         ! help variable
  real(sgl)         :: denomki(6)    ! denominator for knockout formula
  real(sgl)         :: ejec2sp1      ! 2*spin +1 of ejectile
  real(sgl)         :: ejecmass      ! mass of ejectile
  real(sgl)         :: emax          ! maximal emission energy within bin decay
  real(sgl)         :: Eout          ! outgoing energy
  real(sgl)         :: Eres          ! total energy of residual system
  real(sgl)         :: factor        ! multiplication factor
  real(sgl)         :: factor1       ! help variable
  real(sgl)         :: gscomp(6)     ! single-particle level density parameter
  real(sgl)         :: P             ! pairing energy
  real(sgl)         :: phi           ! help variable
  real(sgl)         :: Pn(6)         ! help variable
  real(sgl)         :: preeqpair     ! pre-equilibrium pairing energy
  real(sgl)         :: proj2sp1      ! 2*spin +1 of projectile
  real(sgl)         :: sigav         ! average cross section for emission channel
  real(sgl)         :: term1         ! help variable
  real(sgl)         :: termin0       ! help variable
  real(sgl)         :: termk0        ! help variable
  real(sgl)         :: termki        ! term for knockout and inelastic
  real(sgl)         :: total         ! help variable
  real(sgl)         :: U             ! excitation energy minus pairing energy
!
! ************************** Kalbach model *****************************
!
! The knockout model is described in C. Kalbach, "Preequilibrium reactions with complex channels", Phys. Rev. C71, 034606 (2005).
!
  projmass = parmass(k0)
  proj2sp1 = 2. * parspin(k0) + 1.
  do type = 1, 6
    if (parskip(type)) cycle
    if (k0 <= 2 .and. type <= 2) cycle
!
! Factors for knockout and inelastic processes.
!
! Calculation of terms independent of emission energy.
!
    ejecmass = parmass(type)
    ejec2sp1 = 2. * parspin(type) + 1.
    ndelta = abs(parA(k0) - parA(type))
    flagknock = ((k0 == 1 .or. k0 == 2) .and. type == 6)
    flaginel = (k0 == type)
    if ( .not. (flagknock .or. flaginel)) cycle
    termki = 0.
    AKO = 0.
    Ccl = 1. / 14.
    phi = 0.08
    if (Ntarget > 116 .and. Ntarget < 126) phi = 0.02 + 0.06 * (126 - Ntarget) / 10.
    if (Ntarget >= 126 .and. Ntarget < 129) phi = 0.02 + 0.06 * (Ntarget - 126) / 3.
    denom = Atarget - 2. * phi * Ztarget + 0.5 * phi * Ztarget
    Pn(1) = (Ntarget - phi * Ztarget) / denom
    Pn(2) = (Ztarget - phi * Ztarget) / denom
    Pn(6) = 0.5 * phi * Ztarget / denom
    gscomp(1) = Atarget / 13.
    gscomp(2) = Atarget / 13.
    gscomp(3) = Atarget / 52.
    gscomp(4) = Atarget / 156.
    gscomp(5) = Atarget / 156.
    gscomp(6) = Atarget / 208.
    if (flagknock) AKO = 1. / (2. * (gscomp(k0) **2)) + 1. / (2. * (gscomp(6) **2))
!
! Denominator of cluster emission formula
!
    do type2 = 1, 6
      denomki(type2) = 0.
      if (parskip(type2)) cycle
      emax = eninccm + Q(type2)
      dE = emax - coulbar(type2)
      if (dE > 2.) then
        total = 0.
        do nen = ebegin(type2), eend(type2)
          Eout = egrid(nen)
          if (Eout < coulbar(type2)) cycle
          total = total + xsreac(type2, nen) * deltaE(nen)
        enddo
        sigav = total / dE
      else
        sigav = xsreac(type2, eend(type2))
      endif
      denomki(type2) = (2. * parspin(type2) + 1.) * sigav * (emax + 2. * coulbar(type2)) * max(dE **2, 100.)
    enddo
!
! Knockout and inelastic terms
!
! adjust        : subroutine for energy-dependent parameter adjustment
!
    key = 'cknock'
    call adjust(Einc, key, 0, 0, type, 0, factor)
    Ckn = factor * Cknock(type)
    if (flagknock) then
      termk0 = projmass * denomki(k0) * gscomp(k0) * gscomp(6) **2 / (6. * gscomp(k0)) + &
        ejecmass * denomki(6) * gscomp(k0) * gscomp(6) **2 / (6. * gscomp(6))
      if (termk0 /= 0.) then
        term1 = Ccl * xsreacinc * ejec2sp1 * ejecmass
        termki = Ckn * term1 * Pn(6) * gscomp(k0) * gscomp(6) / termk0
      endif
    else
      do type2 = 1, 6
        if (type2 >= 3 .and. type2 <= 5) cycle
        termin0 = projmass * denomki(k0) * gscomp(k0) * gscomp(type2) **2 / (6. * gscomp(k0)) + &
          projmass * denomki(type2) * gscomp(k0) * gscomp(type2) **2 / (6. * gscomp(type2))
        if (termin0 /= 0.) then
          term1 = Ccl * xsreacinc * proj2sp1 * projmass
          termki = termki + Ckn * term1 * Pn(type2) * gscomp(type2) * gscomp(type2) / termin0
        endif
      enddo
    endif
!
! Calculation of knockout or inelastic spectra.
!
    do nen = ebegin(type), eend(type)
      Eout = egrid(nen)
      factor1 = xsreac(type, nen) * Eout
      P = preeqpair(parZ(type), parN(type), ndelta, Etotal, pairmodel)
      Eres = Etotal - S(0, 0, type) - Eout - P
!
! Check if outgoing energy exceeds maximal possible energy
!
      if (Eres < 0.) cycle
!
! Knockout term that depends on emission energy.
!
      if (flagknock .or. flaginel) then
        U = max(Eres - AKO, 0.)
        xspreeqki(type, nen) = termki * factor1 * U
      endif
    enddo
  enddo
  return
end subroutine knockout
! Copyright A.J. Koning 2021
