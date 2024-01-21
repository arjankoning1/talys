subroutine exciton2
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Two-component exciton model
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
!   sgl            ! single precision kind
! Variables for preequilibrium
!   flagpeout      ! flag for output of pre - equilibrium results
!   pespinmodel    ! model for pre - equilibrium or compound spin distribution
!   preeqmode      ! designator for pre - equilibrium model
! Variables for energy grid
!   ebegin         ! first energy point of energy grid
! Variables for energies
!   eend           ! last energy point of energy grid
!   Etotal         ! total energy of compound system (target + projectile)
! Variables for nuclides
!   parskip        ! logical to skip outgoing particle
! Variables for preequilibrium initialization
!   maxJph         ! maximal spin for particle - hole states
!   maxpar         ! maximal particle number
!   RnJ            ! spin distribution for particle - hole states
!   RnJsum         ! (2J + 1) * sum over spin distributions
! Variables for incident channel
!   xspreeq        ! preeq. cross section per particle type and outgoing energy
! Variables for preequilibrium
!   Ecomp          ! total energy of composite system
!   p0             ! initial particle number
!   pnu0           ! initial neutron number
!   ppi0           ! initial proton number
!   Spre           ! time - integrated strength of two - component exciton state
!   wemission2     ! two - component emission rate
!   xsflux         ! cross section flux
!   xspreeqJP      ! preeq. cross section per particle type, outgoing energy, J, P
!   xsstep         ! preeq. cross section per particle type, stage and outgoing E
!   xsstep2        ! two - component preequilibrium cross section
!
! *** Declaration of local data
!
  implicit none
  integer   :: h      ! help variable
  integer   :: hnu    ! neutron hole number
  integer   :: hpi    ! proton hole number
  integer   :: J      ! spin of level
  integer   :: n      ! exciton number
  integer   :: Ncomp  ! neutron number index for compound nucleus
  integer   :: nen    ! energy counter
  integer   :: p      ! particle number
  integer   :: parity ! parity
  integer   :: pnu    ! neutron particle number
  integer   :: ppi    ! proton particle number
  integer   :: type   ! particle type
  integer   :: Zcomp  ! proton number index for compound nucleus
  real(sgl) :: factor ! multiplication factor
  real(sgl) :: xs     ! help variable
!
! ******************** Never-come-back solution ************************
!
! exchange2: subroutine for calculation of two-component exchange terms
! lifetime2: subroutine for calculation of lifetime of two-component exciton state
!
! For each exciton number, we subsequently calculate the emission rates and exchange terms (both in subroutine exchange2)
! and the lifetime of the exciton state.
!
  Zcomp = 0
  Ncomp = 0
  Ecomp = Etotal
  call exchange2(Zcomp, Ncomp)
  do p = p0, maxpar
    do ppi = ppi0, maxpar
      hpi = ppi - ppi0
      do pnu = pnu0, maxpar
        hnu = pnu - pnu0
        if (ppi + pnu == p) call lifetime2(ppi, hpi, pnu, hnu)
      enddo
    enddo
  enddo
!
! ************* Calculation of pre-equilibrium cross section ***********
!
! The lifetimes and emission rates are processed into pre-equilibrium cross sections.
!
Loop1:  do ppi = ppi0, maxpar
    hpi = ppi - ppi0
    do pnu = pnu0, maxpar
      hnu = pnu - pnu0
      p = ppi + pnu
      if (p > maxpar) cycle Loop1
      h = hpi + hnu
      n = p + h
      factor = xsflux * Spre(ppi, hpi, pnu, hnu)
      do type = 0, 6
        if (parskip(type)) cycle
        if (preeqmode == 4 .and. (type == 1 .or. type == 2)) cycle
        do nen = ebegin(type), eend(type)
          xs = factor * wemission2(type, ppi, hpi, pnu, hnu, nen)
          xsstep(type, p, nen) = xsstep(type, p, nen) + xs
          xsstep2(type, ppi, pnu, nen) = xs
          xspreeq(type, nen) = xspreeq(type, nen) + xs
!
! Create J-dependent pre-equilibrium cross sections using the spin distribution.
! The result is normalized with the sum of RnJ over J.
! As an alternative, the Hauser-Feshbach spin distribution is adopted.
!
          if (pespinmodel >= 3) then
            do parity = - 1, 1, 2
              do J = 0, maxJph
                xspreeqJP(type, nen, J, parity) = xspreeqJP(type, nen, J, parity) + &
                  xs * 0.5 * (2. * J + 1.) * RnJ(n, J) / RnJsum(n)
              enddo
            enddo
          endif
        enddo
      enddo
    enddo
  enddo Loop1
!
! exciton2out    : subroutine for output of two-component exciton model parameters
!
  if (flagpeout) call exciton2out
  return
end subroutine exciton2
! Copyright A.J. Koning 2021
