subroutine multistepA
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Multi-step direct cross sections for MSD
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
! Variables for output
!   flagddx        ! flag for output of double - differential cross sections
! Variables for numerics
!   nanglecont     ! number of angles for continuum
! Variables for main input
!   k0             ! index of incident particle
! Variables for energy grid
!   anglecont      ! angle in degrees for continuum
! Variables for nuclides
!   parskip        ! logical to skip outgoing particle
! Constants
!   amu4pi2h2c2    ! amu / (4 * pi * pi * clight * clight * hbar **2) in mb ** - 1.MeV ** - 1
!   deg2rad        ! conversion factor for degrees to radians
!   parmass        ! mass of particle in a.m.u.
!   parN           ! neutron number of particle
!   parZ           ! charge number of particle
!   pi             ! pi
! Variables for masses
!   specmass       ! specific mass for residual nucleus
! Variables for MSD
!   dEmsd          ! energy bin for MSD
!   Emsd           ! minimal outgoing energy for MSD calculation
!   maxmsd         ! number of MSD steps
!   msdbins2       ! number of energy points for MSD calculation
!   msdstep0       ! n - step cross section for MSD
!   msdstepad0     ! n - step angular distribution for MSD
!   nangleint      ! number of possibilities to link intermediate angle to final angle
!   xscont         ! continuum one - step direct cross section
!   xscontad       ! continuum one - step direct angular distribution for MSD
!
! *** Declaration of local data
!
  implicit none
  integer   :: iang    ! running variable for angle
  integer   :: iangint ! intermediate angle index
  integer   :: iangout ! outgoing angle index
  integer   :: nenint  ! energy index
  integer   :: nenout  ! counter for outgoing energy
  integer   :: ns      ! counter
  integer   :: type    ! particle type
  real(sgl) :: angint  ! intermediate angle
  real(sgl) :: Cmulti  ! constant for second and higher steps
  real(sgl) :: dang    ! delta angle
  real(sgl) :: fac1    ! help variable
  real(sgl) :: fac2    ! help variable
  real(sgl) :: term1   ! help variable
  real(sgl) :: term2   ! help variable
  real(sgl) :: tot1    ! help variable
  real(sgl) :: total   ! help variable
!
! ************** Angle-integrated multi-step cross sections ************
!
  Cmulti = real(specmass(parZ(k0), parN(k0), k0)*amu4pi2h2c2)
  do type = 1, 2
    if (parskip(type)) cycle
    do ns = 2, maxmsd
      do nenout = 1, msdbins2
        msdstep0(type, ns, nenout) = 0.
        total = 0.
        do nenint = 1, nenout - 1
          term1 = msdstep0(type, ns - 1, nenout) * Emsd(nenint) * xscont(type, type, nenint, nenout)
          term2 = msdstep0(k0, ns - 1, nenout) * Emsd(nenint) * xscont(k0, type, nenint, nenout)
          total = total + 0.5 * (term1 + term2)
        enddo
        msdstep0(type, ns, nenout) = total * real(parmass(type)) * Cmulti * dEmsd
      enddo
    enddo
  enddo
!
! ******************* Multi-step angular distributions *****************
!
  if ( .not. flagddx) return
  dang = pi / nanglecont
  do type = 1, 2
    if (parskip(type)) cycle
    do ns = 2, maxmsd
      do nenout = 1, msdbins2
        do iang = 0, nanglecont
          msdstepad0(type, ns, nenout, iang) = 0.
        enddo
        do iangout = 0, nanglecont
          total = 0.
          do nenint = 1, nenout - 1
            do iangint = 0, nanglecont
              tot1 = 0.
              fac1 = msdstepad0(type, ns - 1, nenint, iangint) * Emsd(nenint)
              fac2 = msdstepad0(k0, ns - 1, nenint, iangint) * Emsd(nenint)
              do iang = 0, nanglecont
                term1 = fac1 * nangleint(iangout, iangint, iang) * xscontad(type, type, nenint, nenout, iang)
                term2 = fac2 * nangleint(iangout, iangint, iang) * xscontad(k0, type, nenint, nenout, iang)
                tot1 = tot1 + 0.5 * (term1 + term2)
              enddo
              angint = anglecont(iangint) * deg2rad
              total = total + tot1 * dang * dang * sin(angint)
            enddo
            msdstepad0(type, ns, nenout, iangout) = total * real(parmass(type)) * Cmulti * dEmsd
          enddo
        enddo
      enddo
    enddo
  enddo
  return
end subroutine multistepA
! Copyright A.J. Koning 2021
