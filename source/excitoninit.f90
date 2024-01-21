subroutine excitoninit
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialization of exciton model parameters
!
! Author    : Marieke Duijvestijn and Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl           ! single precision kind
! Variables for preequilibrium
!   flag2comp     ! flag for two - component pre - equilibrium model
!   preeqmode     ! designator for pre - equilibrium model
! Variables for main input
!   Ainit         ! mass number of initial compound nucleus
!   k0            ! index of incident particle
!   Ninit         ! neutron number of initial compound nucleus
!   Zinit         ! charge number of initial compound nucleus
! Variables for nuclides
!   Nindex        ! neutron number index for residual nucleus
!   parskip       ! logical to skip outgoing particle
!   Zindex        ! charge number index for residual nucleus
! Constants
!   amupi2h3c2    ! amu / (pi * pi * clight * clight * hbar **3) in mb ** - 1.MeV ** - 2.s ** - 1
!   parA          ! mass number of particle
!   parN          ! neutron number of particle
!   parspin       ! spin of particle
!   parZ          ! charge number of particle
!   pi2h3c2       ! 1 / (pi * pi * clight * clight * hbar **3) in mb ** - 1.MeV ** - 3.s ** - 1
! Variables for masses
!   redumass      ! reduced mass
! Variables for preequilibrium initialization
!   maxpar        ! maximal particle number
!   nfac          ! n!
!   numparx       ! maximum number of particles
! Variables for exciton model initialization
!   Qfactor       ! Q - factor for neutron / proton distinction
!   wfac          ! factor for emission rate
! Variables for exciton model
!   depletion     ! depletion factor at each stage
!   Gnupi         ! two - component branching ratio
!   Gnuplus       ! two - component branching ratio
!   Gpinu         ! two - component branching ratio
!   Gpiplus       ! two - component branching ratio
!   Lexc          ! exchange term
!   tauexc        ! lifetime of exciton state
!   tauexc2       ! lifetime of two - component exciton state
!   Wompfac       ! adjustable constant for OMP based transition rates

!
! *** Declaration of local data
!
  implicit none
  integer   :: aejec   ! mass number of leading particle
  integer   :: hnu     ! neutron hole number
  integer   :: hpi     ! proton hole number
  integer   :: Ncomp   ! neutron number index for compound nucleus
  integer   :: nejec   ! neutron number of leading particle
  integer   :: Nix     ! neutron number index for residual nucleus
  integer   :: nnu     ! neutron exciton number
  integer   :: npi     ! proton exciton number
  integer   :: nproj   ! neutron number of particle
  integer   :: p       ! particle number
  integer   :: pnu     ! neutron particle number
  integer   :: ppi     ! proton particle number
  integer   :: type    ! particle type
  integer   :: Zcomp   ! proton number index for compound nucleus
  integer   :: zejec   ! charge number of leading particle
  integer   :: Zix     ! charge number index for residual nucleus
  integer   :: zproj   ! charge number of particle
  real(sgl) :: factor1 ! help variable
  real(sgl) :: factor2 ! help variable
  real(sgl) :: factor3 ! help variable
  real(sgl) :: factor4 ! help variable
  real(sgl) :: factor5 ! help variable
  real(sgl) :: factor6 ! help variable
  real(sgl) :: factor7 ! help variable
  real(sgl) :: NoverA  ! help variable
  real(sgl) :: Qdenom  ! help variable
  real(sgl) :: Qenum   ! help variable
  real(sgl) :: ZoverA  ! help variable
!
! Factors for the emission rate.
!
  wfac(0) = pi2h3c2
  Zcomp = 0
  Ncomp = 0
  do type = 1, 6
    if (parskip(type)) cycle
    Zix = Zindex(Zcomp, Ncomp, type)
    Nix = Nindex(Zcomp, Ncomp, type)
    wfac(type) = amupi2h3c2 * (2. * parspin(type) + 1.) * redumass(Zix, Nix, type)
  enddo
!
! ******************* Q-factors for emission rate **********************
!
  if (.not. flag2comp) then
    ZoverA = real(Zinit) / real(Ainit)
    NoverA = real(Ninit) / real(Ainit)
    zproj = parZ(k0)
    nproj = parN(k0)
    do type = 1, 6
      do p = 1, numparx
        Qfactor(type, p) = 1.
      enddo
      if (parskip(type)) cycle
      zejec = parZ(type)
      nejec = parN(type)
      aejec = parA(type)
      do p = 1, maxpar
        if (p < aejec) cycle
        Qenum = 0.
        Qdenom = 0.
        do ppi = zproj, p - nproj
          hpi = ppi - zproj
          npi = ppi + hpi
          pnu = p - ppi
          hnu = p - nproj - ppi
          nnu = pnu + hnu
          factor1 = ZoverA **npi
          factor2 = NoverA **nnu
          factor3 = 1. / real(nfac(ppi) * nfac(pnu))
          factor4 = real(nfac(hpi) * nfac(hnu))
          Qdenom = Qdenom + factor1 * factor2 * factor3 / factor4
          if (ppi < zejec .or. pnu < nejec) cycle
          factor5 = ZoverA **(npi - zejec)
          factor6 = NoverA **(nnu - nejec)
          factor7 = 1. / real(nfac(ppi - zejec) * nfac(pnu - nejec))
          Qenum = Qenum + factor5 * factor6 * factor7 / factor4
        enddo
        if (Qdenom > 0.) Qfactor(type, p) = (Qenum / Qdenom) * real(nfac(p - aejec)) / real(nfac(p))
      enddo
    enddo
  endif
!
! ************* Potential for optical model exciton model **************
!
! bonetti: subroutine for determination of effective absorption optical potential
!
  if (preeqmode == 3) call bonetti
!
! Initialization
!
  depletion = 0.
  Gnupi = 0.
  Gnuplus = 0.
  Gpinu = 0.
  Gpiplus = 0.
  Lexc = 0.
  tauexc = 0.
  tauexc2 = 0.
  Wompfac = 0.
  return
end subroutine excitoninit
! Copyright A.J. Koning 2021
