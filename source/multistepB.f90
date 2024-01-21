subroutine multistepB
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Multi-step direct cross sections on outgoing energy grid
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
!   sgl           ! single precision kind
! Variables for output
!   flagddx       ! flag for output of double - differential cross sections
! Variables for numerics
!   nanglecont    ! number of angles for continuum
! Variables for energy grid
!   ebegin        ! first energy point of energy grid
!   egrid         ! outgoing energy grid
! Variables for energies
!   eend          ! last energy point of energy grid
!   eninccm       ! center - of - mass incident energy in MeV
! Variables for nuclides
!   parskip       ! logical to skip outgoing particle
! Variables for MSD
!   Emsd          ! minimal outgoing energy for MSD calculation
!   maxmsd        ! number of MSD steps
!   msdbins2      ! number of energy points for MSD calculation
!   msdstep       ! continuum n - step direct cross section
!   msdstep0      ! n - step cross section for MSD
!   msdstepad     ! continuum n - step direct angular distribution
!   msdstepad0    ! n - step angular distribution for MSD
!
! *** Declaration of local data
!
  implicit none
  integer   :: iang ! running variable for angle
  integer   :: na   ! help variable
  integer   :: nb   ! help variable
  integer   :: nc   ! counter
  integer   :: nen  ! energy counter
  integer   :: nen2 ! energy counter
  integer   :: ns   ! parameter for fission
  integer   :: type ! particle type
  real(sgl) :: Ea   ! start energy of local adjustment
  real(sgl) :: Eb   ! end energy of local adjustment
  real(sgl) :: Ec   ! help variable for energy
  real(sgl) :: Eout ! outgoing energy
  real(sgl) :: xs   ! help variable
  real(sgl) :: xsa  ! help variable
  real(sgl) :: xsb  ! help variable
  real(sgl) :: xsc  ! interpolated cross section
!
! ********** Interpolate continuum multi-step cross sections ***********
!
! locate     : subroutine to find value in ordered table
! pol2       : subroutine for polynomial interpolation of second order
!
  do type = 1, 2
    if (parskip(type)) cycle
    do nen = ebegin(type), eend(type)
      Eout = egrid(nen)
      if (Eout > eninccm) cycle
      call locate(Emsd, 0, msdbins2, Eout, nen2)
      if (nen2 > 1) then
        na = nen2 - 1
        nb = nen2
        nc = nen2 + 1
      else
        na = nen2
        nb = nen2 + 1
        nc = nen2 + 2
      endif
      Ea = Emsd(na)
      Eb = Emsd(nb)
      Ec = Emsd(nc)
      do ns = 2, maxmsd
        xsa = max(msdstep0(type, ns, na), 1.e-30)
        xsb = max(msdstep0(type, ns, nb), 1.e-30)
        xsc = max(msdstep0(type, ns, nc), 1.e-30)
        call pol2(Ea, Eb, Ec, xsa, xsb, xsc, Eout, xs)
        if (xs < 1.e-30) xs = 0.
        msdstep(type, ns, nen) = xs
        if (flagddx) then
          do iang = 0, nanglecont
            xsa = max(msdstepad0(type, ns, na, iang), 1.e-30)
            xsb = max(msdstepad0(type, ns, nb, iang), 1.e-30)
            xsc = max(msdstepad0(type, ns, nc, iang), 1.e-30)
            call pol2(Ea, Eb, Ec, xsa, xsb, xsc, Eout, xs)
            if (xs < 1.e-30) xs = 0.
            msdstepad(type, ns, nen, iang) = xs
          enddo
        endif
      enddo
    enddo
  enddo
  return
end subroutine multistepB
! Copyright A.J. Koning 2021
