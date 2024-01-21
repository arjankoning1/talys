subroutine preeqinit
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialization of general pre-equilibrium parameters
!
! Author    : Arjan Koning
!
! 2023-08-15: Original code
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
!   g             ! single - particle level density parameter
!   gn            ! single - particle neutron level density parameter
!   gp            ! single - particle proton level density parameter
! Variables for main input
!   Atarget       ! mass number of target nucleus
!   k0            ! index of incident particle
! Variables for energies
!   flagmulpre    ! flag for multiple pre - equilibrium calculation
! Variables for preequilibrium initialization
!   Apauli        ! two - component Pauli blocking correction
!   Apauli2       ! two - component Pauli blocking correction
!   Efermi        ! depth of Fermi well
!   maxexc        ! maximal exciton number
!   maxpar        ! maximal particle number
!   ncomb         ! n! / (n!(n - k)!)
!   nfac          ! n!
!   numexc        ! maximum exciton number
!   numparx       ! maximum number of particles
!   Rblann        ! Blann's factor
!
! *** Declaration of local data
!
  implicit none
  integer   :: h        ! hole number
  integer   :: hnu      ! neutron hole number
  integer   :: hpi      ! proton hole number
  integer   :: i        ! counter
  integer   :: J        ! spin
  integer   :: k        ! counter
  integer   :: n        ! exciton number
  integer   :: p        ! particle number
  integer   :: pnu      ! neutron particle number
  integer   :: ppi      ! proton particle number
  real(sgl) :: Epnu     ! Pauli energy
  real(sgl) :: Epp      ! Pauli energy
  real(sgl) :: Eppi     ! Pauli energy
  real(sgl) :: factor   ! help factor
  real(sgl) :: factorn  ! help variable
  real(sgl) :: factorp  ! help variable
  real(sgl) :: gs       ! single-particle level density parameter
  real(sgl) :: gsn      ! single-particle neutron level density parameter
  real(sgl) :: gsp      ! single-particle proton level density parameter
!
! ********************** General initializations ***********************
!
  maxexc = numexc
  maxpar = maxexc / 2
  Efermi = 38.
!
! ******************* Factorials and combinatorials ********************
!
  nfac(0) = 1.
  ncomb(0, 0) = 1.
  do n = 1, numexc
    nfac(n) = real(n) * nfac(n - 1)
    ncomb(n, 0) = 1.
    do k = 1, n
      ncomb(n, k) = nfac(n) / (nfac(k) * nfac(n - k))
    enddo
  enddo
!
! ******************** Pauli correction factors ************************
!
! 1. One component
!
  gs = g(0, 0)
  do p = - 1, numparx + 1
    do h = - 1, numparx + 1
      if (p ==  - 1 .or. h ==  - 1) then
        Apauli(p, h) = 0.
      else
        Epp = max(p, h) **2 / gs
        factor = (p * p + h * h + p + h) / (4. * gs)
        Apauli(p, h) = Epp - factor
      endif
    enddo
  enddo
!
! 2. Two components (Kalbach, Phys. Rev. C33, 818 (1986)).
!
  gsp = gp(0, 0)
  gsn = gn(0, 0)
  do ppi = - 1, numparx + 1
    do hpi = - 1, numparx + 1
      Eppi = max(ppi, hpi) **2 / gsp
      factorp = (ppi * ppi + hpi * hpi + ppi + hpi) / (4. * gsp)
      do pnu = - 1, numparx + 1
        do hnu = - 1, numparx + 1
          if (ppi ==  - 1 .or. hpi ==  - 1 .or. pnu ==  - 1 .or. hnu ==  - 1) then
            Apauli2(ppi, hpi, pnu, hnu) = 0.
          else
            Epnu = max(pnu, hnu) **2 / gsn
            factorn = (pnu * pnu + hnu * hnu + pnu + hnu) / (4. * gsn)
            Apauli2(ppi, hpi, pnu, hnu) = Eppi + Epnu - factorp - factorn
          endif
        enddo
      enddo
    enddo
  enddo
!
! ************* Blann's R-factor for multiple pre-equilibrium **********
!
! The R-factor takes into account the neutron-proton distinction in multiple preequilibrium.
! The indices for the R-factor are Rblann(itype,type,nstep) where type is the second emitted particle,
! itype is the first emitted particle and nstep is the stage of primary pre-equilibrium emission.
! For this factor it is assumed that N=Z and that an effective interaction ratio of  p-n  to p-p or n-n of a
! factor 3 applies, see Blann and Vonach, Phys. Rev. C28, 1475 (1983).
! Eqs. (9) and (10) and the paragraph below. The numbers are taken from GNASH.
!
  if (flagmulpre .and. .not. flag2comp) then
!
! 1. Incident neutron
!
    if (k0 == 1) then
      Rblann(2, 2, 1) = 0.
      Rblann(2, 1, 1) = 1.
      Rblann(2, 2, 2) = 0.25
      Rblann(2, 1, 2) = 0.75
      Rblann(2, 2, 3) = 0.375
      Rblann(2, 1, 3) = 0.625
      Rblann(2, 2, 4) = 0.438
      Rblann(2, 1, 4) = 0.562
      Rblann(2, 2, 5) = 0.469
      Rblann(2, 1, 5) = 0.531
      Rblann(1, 1, 1) = 0.25
      Rblann(1, 2, 1) = 0.75
      Rblann(1, 1, 2) = 0.375
      Rblann(1, 2, 2) = 0.625
      Rblann(1, 1, 3) = 0.438
      Rblann(1, 2, 3) = 0.562
      Rblann(1, 1, 4) = 0.469
      Rblann(1, 2, 4) = 0.531
      Rblann(1, 1, 5) = 0.484
      Rblann(1, 2, 5) = 0.516
    else
!
! 1. Incident proton
!
      Rblann(1, 1, 1) = 0.
      Rblann(1, 2, 1) = 1.
      Rblann(1, 1, 2) = 0.25
      Rblann(1, 2, 2) = 0.75
      Rblann(1, 1, 3) = 0.375
      Rblann(1, 2, 3) = 0.625
      Rblann(1, 1, 4) = 0.438
      Rblann(1, 2, 4) = 0.562
      Rblann(1, 1, 5) = 0.469
      Rblann(1, 2, 5) = 0.531
      Rblann(2, 1, 1) = 0.75
      Rblann(2, 2, 1) = 0.25
      Rblann(2, 1, 2) = 0.625
      Rblann(2, 2, 2) = 0.375
      Rblann(2, 1, 3) = 0.562
      Rblann(2, 2, 3) = 0.438
      Rblann(2, 1, 4) = 0.531
      Rblann(2, 2, 4) = 0.469
      Rblann(2, 1, 5) = 0.516
      Rblann(2, 2, 5) = 0.484
    endif
    do i = 1, 2
      do J = 1, 2
        do k = 6, numparx
          Rblann(i, J, k) = 0.5
        enddo
      enddo
    enddo
  endif
!
! excitoninit : subroutine for initialization of exciton model parameters
!
  call excitoninit
  return
end subroutine preeqinit
! Copyright A.J. Koning 2021
