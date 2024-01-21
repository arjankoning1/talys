subroutine exciton
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Exciton model
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
! Variables for exciton model
!   tauexc         ! lifetime of exciton state
! Variables for incident channel
!   xspreeq        ! preeq. cross section per particle type and outgoing energy
! Variables for preequilibrium
!   Ecomp          ! total energy of composite system
!   p0             ! initial particle number
!   wemission      ! emission rate per particle, exciton number and energy
!   xsflux         ! cross section flux
!   xspreeqJP      ! preeq. cross section per particle type, outgoing energy, J, P
!   xsstep         ! preeq. cross section per particle type, stage and outgoing E
!
! *** Declaration of local data
!
  implicit none
  integer   :: h      ! help variable
  integer   :: J      ! spin of level
  integer   :: n      ! exciton number
  integer   :: Ncomp  ! neutron number index for compound nucleus
  integer   :: nen    ! energy counter
  integer   :: p      ! particle number
  integer   :: parity ! parity
  integer   :: type   ! particle type
  integer   :: Zcomp  ! proton number index for compound nucleus
  real(sgl) :: factor ! multiplication factor
  real(sgl) :: xs     ! help variable
!
! ******************** Never-come-back solution ************************
!
! emissionrate: subroutine for emission rate
! lifetime    : subroutine for calculation of lifetime of exciton state
!
! For each exciton number, we subsequently calculate the emission rates and the lifetime of the exciton state according to the
! never-come-back approximation.
!
  Zcomp = 0
  Ncomp = 0
  Ecomp = Etotal
  do p = p0, maxpar
    h = p - p0
    call emissionrate(Zcomp, Ncomp, p, h)
    call lifetime(Zcomp, Ncomp, p, h)
  enddo
!
! ************* Calculation of pre-equilibrium cross section ***********
!
! The lifetimes and emission rates are processed into pre-equilibrium cross sections.
!
  do p = p0, maxpar
    h = p - p0
    n = p + h
    factor = xsflux * tauexc(p, h)
    do type = 0, 6
      if (parskip(type)) cycle
!
! Neutron and proton emission can also be calculated with MSD/MSC, if requested.
!
      if (preeqmode == 4 .and. (type == 1 .or. type == 2)) cycle
      do nen = ebegin(type), eend(type)
        xs = factor * wemission(type, p, h, nen)
        xsstep(type, p, nen) = xs
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
!
! excitonout     : subroutine for output of two-component exciton model parameters
!
  if (flagpeout) call excitonout
  return
end subroutine exciton
! Copyright A.J. Koning 2021
