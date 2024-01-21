subroutine preeq
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Pre-equilibrium reactions
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
! Variables for basic reaction
!   flagang         ! flag for output of angular distributions
!   flagrecoil      ! flag for calculation of recoils
! Variables for output
!   flagddx         ! flag for output of double - differential cross sections
! Variables for preequilibrium
!   Esurf0          ! well depth for surface interaction
!   flag2comp       ! flag for two - component pre - equilibrium model
!   flagpecomp      ! flag for Kalbach complex particle emission model
!   flagpeout       ! flag for output of pre - equilibrium results
!   flagsurface     ! flag for surface effects in exciton model
!   preeqmode       ! designator for pre - equilibrium model
! Variables for main input
!   k0              ! index of incident particle
! Variables for OMP
!   flagomponly     ! flag to execute ONLY an optical model calculation
! Variables for energy grid
!   Einc            ! incident energy in MeV
! Variables for incident channel
!   xsdirdiscsum    ! total direct cross section
!   xsgrsum         ! sum over giant resonance cross sections
!   xsreacinc       ! reaction cross section for incident channel
! Constants
!   parA            ! mass number of particle
!   parN            ! neutron number of particle
!   parZ            ! charge number of particle
! Variables for preequilibrium initialization
!   Efermi          ! depth of Fermi well
! Variables for preequilibrium
!   Esurf           ! well depth for surface interaction
!   h0              ! initial hole number
!   hnu0            ! initial neutron hole number
!   hpi0            ! initial proton hole number
!   p0              ! initial particle number
!   pnu0            ! initial neutron number
!   ppi0            ! initial proton number
!   xsflux          ! cross section flux
!
! *** Declaration of local data
!
  implicit none
  real(sgl) :: surface            ! well depth for first hole
!
! ********************** General initializations ***********************
!
  if (flagomponly) return
  p0 = parA(k0)
  h0 = 0
  ppi0 = parZ(k0)
  hpi0 = 0
  pnu0 = parN(k0)
  hnu0 = 0
!
! Special case for photonuclear reactions
!
  if (k0 == 0) then
    p0 = 1
    h0 = 0
    pnu0 = 1
    hnu0 = 0
    ppi0 = 0
    hpi0 = 0
  endif
!
! Initialization for surface interaction
!
  Esurf = Efermi
  if (flagsurface) then
    if(Esurf0 ==  - 1.) then
      if (k0 == 1) Esurf = surface(1, Einc)
      if (k0 == 2) Esurf = surface(2, Einc)
    else
      Esurf = Esurf0
    endif
  endif
!
! *********************** Pre-equilibrium model ************************
!
! exciton     : subroutine for exciton model
! exciton2    : subroutine for two-component exciton model
!
! Pre-equilibrium models:
! preeqmode= 1: Exciton model: Analytical solution - matrix element
! preeqmode= 2: Exciton model: Numerical solution - matrix element
! preeqmode= 3: Exciton model: Numerical solution - optical model in transition rates
! preeqmode= 4: Multi-step direct/Multi-step compound
!
  if (flagpeout) write(*, '(/, " ########## PRE-EQUILIBRIUM ##########")')
!
! Correct reaction cross section for direct and giant resonance effects
!
  xsflux = xsreacinc - xsdirdiscsum - xsgrsum
  xsflux = max(xsflux, 0.)
!
! Choice between one-component and two-component exciton model.
!
  if ( .not. flag2comp) then
    call exciton
  else
    call exciton2
  endif
!
! Quantum-mechanical pre-equilibrium models.
!
! msd         : subroutine for multi-step direct model
! preeqcomplex: subroutine for pre-equilibrium complex particle emission
! preeqcorrect: subroutine to correct pre-equilibrium cross sections for direct effects
!
  if (preeqmode == 4) call msd
!
! For pre-equilibrium reactions involving complex particles, pickup, stripping and knockout contributions are added.
!
  if (flagpecomp) call preeqcomplex
!
! Direct contributions are subtracted from the pre-equilibrium cross sections and pre-equilibrium cross sections are collapsed onto
! discrete states.
!
  if (k0 > 0) call preeqcorrect
!
! The total pre-equilibrium cross sections are assembled.
!
! preeqtotal: subroutine for total pre-equilibrium cross sections
! preeqang  : subroutine for pre-equilibrium angular distributions
! preeqout  : subroutine for output of pre-equilibrium cross sections
!
  call preeqtotal
  if (flagang .or. flagddx .or. flagrecoil) call preeqang
  if (flagpeout) call preeqout
  return
end subroutine preeq
! Copyright A.J. Koning 2021
