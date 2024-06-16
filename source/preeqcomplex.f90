subroutine preeqcomplex
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Pre-equilibrium complex particle emission
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
!   sgl             ! single precision kind
! Variables for preequilibrium
!   breakupmodel    ! model for break - up reaction: 1. Kalbach 2. Avrigeanu
!   pespinmodel     ! model for pre - equilibrium or compound spin distribution
! Variables for main input
!   k0              ! index of incident particle
! Variables for energy grid
!   deltaE          ! energy bin around outgoing energies
!   ebegin          ! first energy point of energy grid
! Variables for energies
!   eend            ! last energy point of energy grid
! Variables for preequilibrium initialization
!   maxJph          ! maximal spin for particle - hole states
! Variables for incident channel
!   xspreeq         ! preeq. cross section per particle type and outgoing energy
! Variables for preequilibrium
!   xsflux          ! cross section flux
!   xspreeqbu       ! preequilibrium cross section per particle type and outgoing energy for break-up
!   xspreeqJP       ! preeq. cross section per particle type, outgoing energy, J, P
!   xspreeqki       ! preequilibrium cross section per particle type and outgoing energy for knockout
!   xspreeqps       ! preequilibrium cross section per particle type and outgoing energy for pickup and stripping
!
! *** Declaration of local data
!
  implicit none
  integer   :: J         ! spin of level
  integer   :: nen       ! energy counter
  integer   :: parity    ! parity
  integer   :: type      ! particle type
  real(sgl) :: factor    ! multiplication factor
  real(sgl) :: damper    ! fermi damping function
  real(sgl) :: expo      ! exponent
  real(sgl) :: c1        ! help variable
  real(sgl) :: c2        ! help variable
  real(sgl) :: pecompsum ! help variable
  real(sgl) :: xspecomp  ! pre-equilibrium complex particle cross section
!
! ********************** Various components ****************************
!
! stripping   : subroutine for contribution of stripping and pickup reactions
! knockout    : subroutine for contribution of knockout reactions
! breakup     : subroutine for contribution of breakup reactions
! breakupAVR  : subroutine for contribution of breakup reactions, Avrigeanu model
!
  call stripping
  call knockout
  if (k0 > 2) then
    if (breakupmodel == 1 .or. k0 /= 3) then
      call breakup
    else
      call breakupAVR
    endif
  endif
!
! *************************** Corrections ******************************
!
! Prevent complex particle pre-equilibrium to cause fluctuations at low energies 
! Prevent complex particle pre-equilibrium to exceed the reaction cross section.
!
  damper = 1. 
  c1 = 0.5 * coulbar(k0)
  c2 = 0.1 * coulbar(k0)
  expo = (Einc - c1) / c2
  if (expo <= 80.) damper = 1. - 1. / (1. + exp(expo))
  pecompsum = 0.
  do type = 1, 6
    do nen = ebegin(type), eend(type)
      xspreeqps(type, nen) = xspreeqps(type, nen) * damper 
      xspreeqki(type, nen) = xspreeqki(type, nen) * damper
      xspreeqbu(type, nen) = xspreeqbu(type, nen) * damper
      xspecomp = xspreeqps(type, nen) + xspreeqki(type, nen) + xspreeqbu(type, nen)
      pecompsum = pecompsum + xspecomp * deltaE(nen)
    enddo
  enddo
  if (pecompsum > xsflux) then
    factor = xsflux / pecompsum
    do type = 1, 6
      do nen = ebegin(type), eend(type)
        xspreeqps(type, nen) = xspreeqps(type, nen) * factor 
        xspreeqki(type, nen) = xspreeqki(type, nen) * factor
        xspreeqbu(type, nen) = xspreeqbu(type, nen) * factor
      enddo
    enddo
  endif
!
! If the pre-equilibrium spin distribution is chosen, we assume that the spin distribution for pickup,
! stripping and knockout is the same as in the exciton model.
!
  do type = 1, 6
    do nen = ebegin(type), eend(type)
      xspecomp = xspreeqps(type, nen) + xspreeqki(type, nen) + xspreeqbu(type, nen)
      if (pespinmodel >= 3 .and. xspreeq(type, nen) /= 0.) then
        do parity = - 1, 1, 2
          do J = 0, maxJph
            factor = xspreeqJP(type, nen, J, parity) / xspreeq(type, nen)
            xspreeqJP(type, nen, J, parity) = xspreeqJP(type, nen, J, parity) + factor * xspecomp
          enddo
        enddo
      endif
      xspreeq(type, nen) = xspreeq(type, nen) + xspecomp
    enddo
  enddo
  return
end subroutine preeqcomplex
! Copyright A.J. Koning 2021
