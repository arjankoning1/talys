subroutine partfunc
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculate partition function
!
! Author    : Stephane Hilaire, Stephane Goriely, Arjan Koning
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
!   dbl            ! double precision kind
! All global variables
!   numex          ! maximum number of excitation energies
!   numJ           ! maximum J - value
! Variables for astrophysics
!   flagastrogs    ! flag for calculation of astrophysics react
!   nonthermlev    ! non - thermalized level in the calculation o
!   nTmax          ! effective number of temperatures for Maxwellian
! Variables for discrete levels
!   nlev           ! number of levels for nucleus
! Variables for main input
!   Atarget        ! mass number of target nucleus
!   k0             ! index of incident particle
!   Ltarget        ! excited level of target
! Variables for level density
!   ldmodel        ! level density model
! Variables for nuclides
!   T9             ! Temperature grid in 10 **9 K
! Variables for excitation energy grid
!   deltaEx        ! excitation energy bin for population arrays
!   Ex             ! excitation energy
!   maxex          ! maximum excitation energy bin for residual nucle
! Constants
!   parN           ! neutron number of particle
!   parZ           ! charge number of particle
! Variables for levels
!   edis           ! energy of level
!   jdis           ! spin of level
! Variables for level density
!   Nlast          ! last discrete level
! Variables for astro
!   partf          ! integrated partition function
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: numdiv=100              ! number of subdivisions considered for each continuum bin  in target nucleus
  integer            :: i                       ! level
  integer            :: idiv                    ! counter
  integer            :: Ir                      ! residual spin
  integer            :: ldmod                   ! level density model
  integer            :: nex                     ! excitation energy bin of compound nucleus
  integer            :: nexout                  ! energy index for outgoing energy
  integer            :: Nix                     ! neutron number index for residual nucleus
  integer            :: NL                      ! last discrete level
  integer            :: Pprime                  ! parity
  integer            :: Zix                     ! charge number index for residual nucleus
  real(sgl)          :: dex                     ! energy bin
  real(sgl)          :: Elev                    ! excitation energy in the target nucleus
  real(sgl)          :: MeVkT                   ! E/kT expressed in MeV/T9
  real(sgl)          :: Rspin                   ! residual spin
  real(sgl)          :: spindeg                 ! (2J+1) degeneracy of target level
  real(dbl)          :: density                 ! level density
  real(dbl)          :: fex                     ! help variable
  real(dbl)          :: rho                     ! integrated level density
  real(dbl)          :: sum                     ! help variable
  real(dbl)          :: sumjp(0:numex*numdiv)   ! help variable
!
! ************************ Partition Function **************************
!
  if (flagastrogs) then
    do i = 1, nTmax
      partf(i) = 1.
    enddo
    return
  endif
  call exgrid(0,0)
  MeVkT = 11.605
  Zix = parZ(k0)
  Nix = parN(k0)
  NL = Nlast(Zix, Nix, 0)
  ldmod = ldmodel(Zix, Nix)
  do i = 1, nTmax
    do nex = 0, numdiv * numex
      sumjp(nex) = 0.
    enddo
    sum = 0.
    nex = 0
    do nexout = 0, max(maxex(Zix, Nix) - 1, nlev(Zix, Nix))
!
! For discrete states, the begin and end points of the target
! spin/parity summation are both set equal to the target discrete
! level spin/parity.
!
      if (nexout == nonthermlev) cycle
      dex = deltaEx(Zix, Nix, nexout) / real(numdiv)
      Elev = edis(Zix, Nix, nexout) - edis(Zix, Nix, Ltarget)
      if (nexout <= NL) then
        spindeg = 2. * jdis(Zix, Nix, nexout) + 1
        fex = MeVkT * Elev / T9(i)
        if (fex > 80.) cycle
        sum = sum + spindeg * exp( - fex)
        sumjp(0) = spindeg * exp( - fex)
      else
!
! For decay to the continuum we use a spin and parity dependent level
! density.
! See subroutine exgrid for change in case of non-equiprobable parities.
!
! partf      : integrated partition function
!
        do idiv = 1, numdiv
          Elev = Ex(Zix, Nix, nexout) - Ex(Zix, Nix, Ltarget) + (real(idiv) - 0.5 * numdiv) * dex
          fex = MeVkT * Elev / T9(i)
          if (fex > 80.) cycle
          nex = nex + 1
          do Pprime = - 1, 1, 2
            do Ir = 0, numJ
              Rspin = Ir + mod(Atarget, 2) / 2.
              spindeg = 2. * (Ir + mod(Atarget, 2) / 2.) + 1.
              rho = density(Zix, Nix, Elev, Rspin, Pprime, 0, ldmod)
              sumjp(nex) = sumjp(nex) + spindeg * rho * exp( - fex) * dex
            enddo
          enddo
        enddo
      endif
    enddo
    do idiv = 1, nex
      sum = sum + (sumjp(idiv - 1) + sumjp(idiv)) / 2.
    enddo
    partf(i) = sum / (2. * jdis(Zix, Nix, Ltarget) + 1.)
  enddo
  return
end subroutine partfunc
! Copyright A.J. Koning 2021
