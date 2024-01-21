subroutine onecontinuumA(itype, type)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Unnormalized one-step direct cross sections for MSD
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
! Variables for preequilibrium
!   flaggshell    ! flag for energy dependence of single particle level density parameter g
!   g             ! single - particle level density parameter
! Variables for output
!   flagddx       ! flag for output of double - differential cross sections
! Variables for numerics
!   nanglecont    ! number of angles for continuum
! Variables for level density
!   alev          ! level density parameter
! Variables for nuclides
!   Nindex        ! neutron number index for residual nucleus
!   Zindex        ! charge number index for residual nucleus
! Variables for masses
!   S             ! separation energy
!   specmass      ! specific mass for residual nucleus
! Variables for MSD
!   Emsd          ! minimal outgoing energy for MSD calculation
!   Emsdin        ! incident MSD energy
!   Emsdout       ! outgoing MSD energy
!   Exmsd         ! excitation energy for MSD energy grid
!   maxJmsd       ! maximal spin for MSD calculation
!   msdbins2      ! number of energy points for MSD calculation
!   numJmsd       ! maximum spin for MSD
!   xscont1       ! continuum one - step direct cross section
!   xscontad1     ! continuum one - step direct angular distribution for MSD (unnormalized)
!   xsdw          ! DWBA angular distribution as a function of incident energy, outgoing energy, ang. mom. and angle
!   xsdwin        ! DWBA cross section as a function of incident energy, outgoing energy and angular momentum
!
! *** Declaration of local data
!
  implicit none
  integer   :: h                 ! hole numberle
  integer   :: iang              ! running variable for angle
  integer   :: itype             ! help variable
  integer   :: J                 ! spin of level
  integer   :: nen1              ! energy counter
  integer   :: nen2              ! energy counter
  integer   :: Nix               ! neutron number index for residual nucleus
  integer   :: p                 ! particle number
  integer   :: type              ! particle type
  integer   :: Zix               ! charge number index for residual nucleus
  real(sgl) :: gs                ! single-particle level density parameter
  real(sgl) :: ignatyuk          ! function for energy dependent level density parameter a
  real(sgl) :: omega             ! particle-hole state density
  real(sgl) :: omegaJ(0:numJmsd) ! help variable
  real(sgl) :: QQ                ! Q-value
  real(sgl) :: rJ                ! help variable
  real(sgl) :: total             ! help variable
  real(sgl) :: xs                ! help variable
!
! ************* Calculate continuum one-step cross sections ************
!
! ignatyuk  : function for energy dependent level density parameter a
!
  h = 1
  p = 1
  Zix = Zindex(0, 0, type)
  Nix = Nindex(0, 0, type)
  gs = g(Zix, Nix)
  QQ = S(0, 0, itype) - S(0, 0, type)
  do nen1 = 0, msdbins2
    Emsdin = specmass(Zix, Nix, itype) * Emsd(nen1)
    do nen2 = nen1, msdbins2
      Emsdout = Emsd(nen2)
      Exmsd = Emsdin - Emsdout + QQ
      if (Exmsd < 0..or.(Exmsd + 0.1) >= Emsdin) cycle
      if (flaggshell) gs = g(Zix, Nix) * ignatyuk(Zix, Nix, Exmsd, 0) / alev(Zix, Nix)
      do J = 0, maxJmsd
        rJ = real(J)
        omegaJ(J) = omega(Zix, Nix, p, h, gs, Exmsd, rJ)
      enddo
      total = 0.
      do J = 0, maxJmsd
        xs = xsdwin(nen1, nen2, J, 0)
        total = total + omegaJ(J) * xs
      enddo
      xscont1(itype, type, nen1, nen2) = total
      if (flagddx) then
        do iang = 0, nanglecont
          total = 0.
          do J = 0, maxJmsd
            xs = xsdw(nen1, nen2, J, iang, 0)
            total = total + omegaJ(J) * xs
          enddo
          xscontad1(itype, type, nen1, nen2, iang) = total
        enddo
      endif
    enddo
  enddo
  return
end subroutine onecontinuumA
! Copyright A.J. Koning 2021
