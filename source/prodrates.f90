subroutine prodrates
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculate reaction rates
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
!   sgl          ! single precision kind
! All global variables
!   numenrp      ! number of incident energies for residual products
! Variables for medical isotope production
!   Area         ! target area in cm^2
!   Eback        ! lower end of energy range in MeV for isotope
!   Ebeam        ! incident energy in MeV for isotope production
!   Ibeam        ! beam current in mA for isotope production
!   rhotarget    ! target material density
! Variables for numerics
!   maxN         ! maximal number of neutrons away from initial compound nucleus
!   maxZ         ! maximal number of protons away from initial compound nucleus
! Variables for main input
!   k0           ! index of incident particle
! Constants
!   parZ         ! charge number of particle
!   qelem        ! elementary charge in C
! Variables for levels
!   Nisomer      ! number of isomers for this nuclide
! Variables for existence libraries
!   prodexist    ! logical to determine existence of residual production
! Variables for isotope production
!   Erp          ! incident energy
!   heat         ! produced heat
!   Mtar         ! active target mass
!   Nenrp        ! number of incident energies for residual production cross
!   prate        ! production rate per isotope
!   projnum      ! number of incident particles [s^ - 1]
!   targetdx     ! effective thickness of target
!   Vtar         ! active target volume
!   xsrp         ! residual production cross section in mb
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: numint=100           ! number of integration points
  integer            :: is                   ! isotope counter: -1=total, 0=ground state 1=isomer
  integer            :: N                    ! neutron number of residual nucleus
  integer            :: nE                   ! number of energies
  integer            :: nen                  ! energy counter
  integer            :: Ninte                ! number of integration points
  integer            :: Nix                  ! neutron number index for residual nucleus
  integer            :: Zix                  ! charge number index for residual nucleus
  real(sgl)          :: dE                   ! help variable
  real(sgl)          :: dEdx                 ! stopping power
  real(sgl)          :: dxdE(numint)         ! 1/stopping power
  real(sgl)          :: dxdEsum              ! integration sum
  real(sgl)          :: E                    ! incident energy
  real(sgl)          :: Ea                   ! start energy of local adjustment
  real(sgl)          :: Eb                   ! end energy of local adjustment
  real(sgl)          :: Eint(numint)         ! energy on integration grid
  real(sgl)          :: Erpgrid(0:numenrp)   ! incident energy for cross section in MeV
  real(sgl)          :: ratesum              ! integration sum
  real(sgl)          :: xs                   ! help variable
  real(sgl)          :: xsa                  ! help variable
  real(sgl)          :: xsb                  ! help variable
!
! *********** Determine integration grid and stopping power ************
!
! stoppingpower: subroutine to calculate stopping power
!
  Ninte = 100
  dE = (Ebeam - Eback) / Ninte
  dxdEsum = 0.
  do nE = 1, Ninte
    Eint(nE) = Eback + (nE - 0.5) * dE
    call stoppingpower(Eint(nE), dEdx)
    dxdE(nE) = 1. / dEdx
    dxdEsum = dxdEsum + dxdE(nE)
  enddo
  targetdx = dxdEsum * dE
  Vtar = Area * targetdx
  Mtar = rhotarget * Vtar
  heat = Ibeam * (Ebeam - Eback)
!
! ********************* Calculate reaction rates ***********************
!
! locate   : subroutine to find value in ordered table
! pol1     : subroutine for interpolation of first order
!
  projnum = Ibeam / (1000. * parZ(k0) * qelem)
  do Zix = - 1, maxZ
    do Nix = - 1, maxN
      do is = - 1, Nisomer(Zix, Nix)
        prate(Zix, Nix, is) = 0.
        ratesum = 0.
        if ( .not. prodexist(Zix, Nix, is)) cycle
        N = Nenrp(Zix, Nix, is)
        if (N == 0) cycle
        do nen = 1, N
          Erpgrid(nen) = Erp(Zix, Nix, is, nen)
        enddo
        do nE = 1, Ninte
          E = Eint(nE)
          if (E < Erpgrid(1)) cycle
          if (E > Erpgrid(N)) exit
          call locate(Erpgrid, 1, N, E, nen)
          if (nen == 0) cycle
          Ea = Erpgrid(nen)
          Eb = Erpgrid(nen + 1)
          xsa = xsrp(Zix, Nix, is, nen)
          xsb = xsrp(Zix, Nix, is, nen + 1)
          call pol1(Ea, Eb, xsa, xsb, E, xs)
          ratesum = ratesum + dxdE(nE) * xs
        enddo
        prate(Zix, Nix, is) = projnum / Vtar * ratesum * dE * 1.e-27
      enddo
    enddo
  enddo
  return
end subroutine prodrates
! Copyright A.J. Koning 2021
