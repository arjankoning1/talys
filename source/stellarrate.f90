subroutine stellarrate
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculate reaction rate for a Maxwell-Boltzmann
!
! Author    : Stephane Goriely and Stephane Hilaire
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl               ! single precision kind
!   dbl               ! double precision kind
! Variables for fission
!   flagfission       ! flag for fission
! Variables for astrophysics
!   flagastrogs       ! flag for calculation of astrophysics reaction rate with
!   nTmax             ! effective number of temperatures for Maxwellian
! Variables for gamma rays
!   flagracap         ! flag for radiative capture model
! Variables for input energies
!   eninc             ! incident energy in MeV
!   Ninc            ! number of incident energies
! Variables for main input
!   k0                ! index of incident particle
! Variables for nuclides
!   T9                ! Temperature grid in 10 **9 K
! Constants
!   parN              ! neutron number of particle
!   parZ              ! charge number of particle
! Variables for levels
!   tau               ! lifetime of state in seconds
! Variables for level density
!   Nlast             ! last discrete level
! Variables for direct capture initialization
!   xsracap           ! direct radiative capture cross section
! Variables for masses
!   redumass          ! reduced mass
!   specmass          ! specific mass for residual nucleus
! Variables for astro
!   macsastro         ! Maxwellian - averaged thermonuclear reaction
!   macsastroex       ! thermonuclear reaction cross section to a
!   macsastrofis      ! thermonuclear reaction cross section for f
!   macsastroracap    ! Maxwellian - averaged thermonuc. reac. c.s.
!   maxNastro         ! maximal number of neutrons away from initi
!   maxZastro         ! maximal number of protons away from initia
!   partf             ! integrated partition function
!   rateastro         ! thermonuclear reaction rate
!   rateastroex       ! thermonuclear reaction rate to a given exc
!   rateastrofis      ! thermonuclear reaction rate factor for fis
!   rateastroracap    ! thermonuclear reaction rate for direct cap
!   xsastro           ! cross section for astrophysical calculatio
!   xsastroex         ! cross section for astrophysical calculati
!   xsastrofis        ! astrophysical fission cross section
!
! *** Declaration of local data
!
  implicit none
  integer   :: i       ! counter
  integer   :: iloop   ! loop counter
  integer   :: Ncomp   ! neutron number index for compound nucleus
  integer   :: nen     ! energy counter
  integer   :: nex     ! excitation energy bin of compound nucleus
  integer   :: nexloop ! counter
  integer   :: Nix     ! neutron number index for residual nucleus
  integer   :: nloop   ! help variable
  integer   :: Zcomp   ! proton number index for compound nucleus
  integer   :: Zix     ! charge number index for residual nucleus
  real(sgl) :: am      ! reduced mass
  real(sgl) :: convfac ! conversion factor for the average velocity at  temperature T9
  real(sgl) :: e       ! energy
  real(sgl) :: ee0     ! incident energy
  real(sgl) :: ee1     ! help variable
  real(sgl) :: fact0   ! reaction rate factor for photons
  real(sgl) :: fex     ! help variable
  real(sgl) :: smacs   ! total conversion factor for the average velocity at  temperature T9
  real(dbl) :: dE      ! help variable
  real(dbl) :: fact    ! reaction rate factor for particles
  real(dbl) :: rate    ! Maxwellian reaction rate
  real(dbl) :: sum     ! help variable
  real(dbl) :: term    ! help variable
  real(dbl) :: xs      ! help variable
!
! ************************ Partition Function **************************
!
! partfunc   : subroutine to calculate partition function
!
! Initialization
!
  fact0 = 3.951776e+17
  Zix = parZ(k0)
  Nix = parN(k0)
  am = redumass(Zix, Nix, k0)
!
! Calculation of the T-dependent partition function
!
  call partfunc
!
! Calculation of thermonuclear reaction rate factor
!
! partf       : integrated partition function
!
  convfac = 2.45484e+5 / sqrt(redumass(parZ(k0), parN(k0), k0))
  do Ncomp = 0, maxNastro
    do Zcomp = 0, maxZastro
      nloop = 1
      if (flagfission .and. Ncomp == 0 .and. Zcomp == 0) nloop = 2
      if (flagracap .and. Ncomp == 0 .and. Zcomp == 0) nloop = 3
      do iloop = 0, nloop
        nexloop = 0
        if (iloop == 0) nexloop = Nlast(Zcomp, Ncomp, 0)
        do nex = 0, nexloop
          if (nex > 0 .and. tau(Zcomp, Ncomp, nex) == 0.) cycle
          do i = 1, nTmax
            fact = 3.7335e+10 / (sqrt(am * (T9(i) **3)))
            sum = 0.
            do nen = 1, Ninc
              e = real(eninc(nen) * specmass(Zix, Nix, k0))
              if (nen < Ninc) then
                ee1 = real(eninc(nen + 1) * specmass(Zix, Nix, k0))
              else
                ee1 = e
              endif
              if (nen > 1) then
                ee0 = real(eninc(nen - 1) * specmass(Zix, Nix, k0))
              else
                ee0 = e
              endif
              if (iloop == 0) xs = xsastroex(Zcomp, Ncomp, nen, nex) / 1000.
              if (iloop == 1) xs = xsastro(Zcomp, Ncomp, nen) / 1000.
              if (iloop == 2) xs = xsastrofis(nen) / 1000.
              if (iloop == 3) xs = xsracap(nen) / 1000.
              fex = 11.605 * e / T9(i)
              if (fex > 80.) cycle
              if (k0 /= 0) then
                term = fact * e * xs * exp( - fex)
              else
                term = fact0 * e **2 * xs / (exp(fex) - 1.)
              endif
              dE = (ee1 - ee0) / 2.
              sum = sum + term * dE
            enddo
            if (flagastrogs .or. iloop == 3) then
              rate = sum
            else
              rate = sum / partf(i)
            endif
            smacs = convfac * sqrt(T9(i))
            if (iloop == 0) then
              rateastroex(Zcomp, Ncomp, i, nex) = rate
              macsastroex(Zcomp, Ncomp, i, nex) = rate / smacs
            endif
            if (iloop == 1) then
              rateastro(Zcomp, Ncomp, i) = rate
              macsastro(Zcomp, Ncomp, i) = rate / smacs
            endif
            if (iloop == 2) then
              rateastrofis(i) = rate
              macsastrofis(i) = rate / smacs
            endif
            if (iloop > 2) then
              rateastroracap(i) = rate
              macsastroracap(i) = rate / smacs
              rateastro(Zcomp, Ncomp, i) = rateastro(Zcomp, Ncomp, i) + rate
              macsastro(Zcomp, Ncomp, i) = macsastro(Zcomp, Ncomp, i) + rate / smacs
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
  return
end subroutine stellarrate
! Copyright A.J. Koning 2021
