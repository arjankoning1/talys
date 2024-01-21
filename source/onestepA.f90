subroutine onestepA(type)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Unnormalized one-step direct cross sections for outgoing energy grid
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
! Variables for energy grid
!   ebegin        ! first energy point of energy grid
!   egrid         ! outgoing energy grid
! Variables for energies
!   eend          ! last energy point of energy grid
!   eninccm       ! center - of - mass incident energy in MeV
! Variables for nuclides
!   Nindex        ! neutron number index for residual nucleus
!   Q             ! Q - value
!   Zindex        ! charge number index for residual nucleus
! Variables for MSD
!   Emsd          ! minimal outgoing energy for MSD calculation
!   Exmsd         ! excitation energy for MSD energy grid
!   maxJmsd       ! maximal spin for MSD calculation
!   msdbins2      ! number of energy points for MSD calculation
!   msdstep1      ! continuum one - step direct cross section (unnormalized)
!   msdstepad1    ! continuum one - step direct angular distribution (unnormalized)
!   numJmsd       ! maximum spin for MSD
!   xsdw          ! DWBA angular distribution as a function of incident energy, outgoing energy, ang. mom. and angle
!   xsdwin        ! DWBA cross section as a function of incident energy, outgoing energy and angular momentum
!
! *** Declaration of local data
!
  implicit none
  integer   :: h                 ! hole number
  integer   :: iang              ! running variable for angle
  integer   :: J                 ! spin of level
  integer   :: na                ! help variable
  integer   :: nb                ! help variable
  integer   :: nc                ! help variable
  integer   :: nen               ! energy counter
  integer   :: nen2              ! energy counter
  integer   :: Nix               ! neutron number index for residual nucleus
  integer   :: p                 ! particle number
  integer   :: type              ! particle type
  integer   :: Zix               ! charge number index for residual nucleus
  real(sgl) :: Ea                ! start energy
  real(sgl) :: Eb                ! end energy
  real(sgl) :: Ec                ! help variable for energy
  real(sgl) :: Eout              ! outgoing energy
  real(sgl) :: gs                ! single-particle level density parameter
  real(sgl) :: ignatyuk          ! function for energy dependent level density parameter a
  real(sgl) :: omega             ! particle-hole state density
  real(sgl) :: omegaJ(0:numJmsd) ! help variable
  real(sgl) :: rJ                ! help variable
  real(sgl) :: total             ! help variable
  real(sgl) :: xs                ! help variable
  real(sgl) :: xsa               ! help variable
  real(sgl) :: xsb               ! help variable
  real(sgl) :: xsc               ! interpolated cross section
  real(sgl) :: xsi               ! help variable
!
! ************* Calculate continuum one-step cross sections ************
!
! ignatyuk   : function for energy dependent level density parameter a
! locate     : subroutine to find value in ordered table
! pol2       : subroutine for polynomial interpolation of second order
!
  Zix = Zindex(0, 0, type)
  Nix = Nindex(0, 0, type)
  h = 1
  p = 1
  gs = g(Zix, Nix)
  do nen = ebegin(type), eend(type)
    Eout = egrid(nen)
    if (Eout > eninccm) cycle
    Exmsd = eninccm - Eout + Q(type)
    if (flaggshell) gs = g(Zix, Nix) * ignatyuk(Zix, Nix, Exmsd, 0) / alev(Zix, Nix)
    do J = 0, maxJmsd
      rJ = real(J)
      omegaJ(J) = omega(Zix, Nix, p, h, gs, Exmsd, rJ)
    enddo
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
    total = 0.
    do J = 0, maxJmsd
      xsa = log(max(xsdwin(0, na, J, 0), 1.e-30))
      xsb = log(max(xsdwin(0, nb, J, 0), 1.e-30))
      xsc = log(max(xsdwin(0, nc, J, 0), 1.e-30))
      call pol2(Ea, Eb, Ec, xsa, xsb, xsc, Eout, xsi)
      xs = exp(xsi)
      if (xs < 1.e-30) xs = 0.
      total = total + omegaJ(J) * xs
    enddo
    msdstep1(type, nen) = total
    if (flagddx) then
      do iang = 0, nanglecont
        total = 0.
        do J = 0, maxJmsd
          xsa = log(max(xsdw(0, na, J, iang, 0), 1.e-30))
          xsb = log(max(xsdw(0, nb, J, iang, 0), 1.e-30))
          xsc = log(max(xsdw(0, nc, J, iang, 0), 1.e-30))
          call pol2(Ea, Eb, Ec, xsa, xsb, xsc, Eout, xsi)
          xs = exp(xsi)
          if (xs < 1.e-30) xs = 0.
          total = total + omegaJ(J) * xs
        enddo
        msdstepad1(type, nen, iang) = total
      enddo
    endif
  enddo
  return
end subroutine onestepA
! Copyright A.J. Koning 2021
