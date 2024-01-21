subroutine eciscompound
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Create ECIS input file for compound cross section
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
! Variables for gamma rays
!   egr           ! energy of GR
!   ggr           ! width of GR
! Variables for main input
!   k0            ! index of incident particle
!   Zinit         ! charge number of initial compound nucleus
! Variables for energy grid
!   Einc          ! incident energy in MeV
! Variables for nuclides
!   Nindex        ! neutron number index for residual nucleus
!   parinclude    ! logical to include outgoing particle
!   parskip       ! logical to skip outgoing particle
!   targetspin    ! spin of target
!   Zindex        ! charge number index for residual nucleus
! Constants
!   parN          ! neutron number of particle
!   parZ          ! charge number of particle
! Variables for ECIS
!   angbeg        ! first angle
!   angend        ! last angle
!   anginc        ! angle increment
!   ecis1         ! 50 input flags ('T' or 'F') for ECIS
!   ecis2         ! 50 input flags ('T' or 'F') for ECIS
!   iterm         ! number of iterations
!   ncoll         ! number of nuclear states
!   njmax         ! maximal number of j - values in ECIS
!   npp           ! number of optical potentials
!   prodZ         ! product of charges of projectile and target nucleus
!   projmass      ! mass of projectile
!   resmass       ! mass of residual nucleus
!   rmatch        ! matching radius
!   spin          ! spin of incident particle
!   tarparity     ! parity of target nucleus
!   title         ! title of ECIS input file
! Variables for masses
!   S             ! separation energy
!   specmass      ! specific mass for residual nucleus
! Variables for optical model
!   av            ! real volume diffuseness
!   avd           ! real surface diffuseness
!   avso          ! real spin - orbit diffuseness
!   aw            ! imaginary volume diffuseness
!   awd           ! imaginary surface diffuseness
!   awso          ! imaginary spin - orbit diffuseness
!   d2            ! parameter for imaginary surface OMP
!   d3            ! parameter for imaginary surface OMP
!   disp          ! flag for dispersive optical model
!   ef            ! Fermi energy
!   rc            ! Coulomb radius
!   rv            ! real volume radius
!   rvd           ! real surface radius
!   rvso          ! real spin - orbit radius
!   rw            ! imaginary volume radius
!   rwd           ! imaginary surface radius
!   rwso          ! imaginary spin - orbit radius
!   v             ! real volume depth
!   vd            ! real surface depth
!   vso           ! real spin - orbit depth
!   w             ! imaginary volume depth
!   w2            ! parameter for imaginary volume OMP
!   wd            ! imaginary surface depth
!   wso           ! imaginary spin - orbit depth
! Variables for ECIS calculation of compound cross sections (reference only)
!   aldcomp       ! level density parameter with indices (Z, N)
!   bz1           ! elastic enhancement factor
!   E0comp        ! constant of temperature formula
!   ejeccomp      ! mass of projectile
!   elevelcomp    ! energy of level
!   Excomp        ! Matching Ex
!   jcomp         ! spin of level
!   masscomp      ! mass of nucleus with indices (Z, N)
!   ncont         ! number of continua
!   nsp1          ! number of uncoupled states and continua
!   nsp2          ! number of uncoupled states with angular distribution
!   pcomp         ! parity of level
!   prodZcomp     ! product of charges
!   spincomp      ! spin of incident particle
!   tempcomp      ! nuclear temperature
!   tgo           ! slow s - wave neutron gamma width / spacing
!   typecomp      ! particle type
!   Umcomp        ! matching point for U (excitation energy - pairing)
!
! *** Declaration of local data
!
  implicit none
  integer   :: i               ! counter
  integer   :: kopt            ! optical model index for fast particle
  integer   :: nex             ! excitation energy bin of compound nucleus
  integer   :: Nix             ! neutron number index for residual nucleus
  integer   :: Zix             ! charge number index for residual nucleus
  real(sgl) :: eopt            ! incident energy
!
! *********************** Write standard input *************************
!
! optical      : subroutine for determination of optical potential
!
  write(1, '(a72)') title
  write(1, '(a50)') ecis1
  write(1, '(a50)') ecis2
  write(1, '(4i5)') ncoll, njmax, iterm, npp
  write(1, '(10x, f10.5, 10x, 3("    1.e-10"))') rmatch
  write(1, '()')
  write(1, '(2i5, 10x, i5)') nsp1, nsp2, ncont
  if (Einc >= 0.01) then
    write(1, '(f5.2, 2i2, a1, 5f10.5)') targetspin, 0, 1, tarparity, Einc, spin, projmass, resmass, prodZ
  else
    write(1, '(f5.2, 2i2, a1, es10.3, 4f10.5)') targetspin, 0, 1, tarparity, Einc, spin, projmass, resmass, prodZ
  endif
  do nex = 1, nsp1
    write(1, '(f5.2, 2i2, a1, 5f10.5)') jcomp(nex), 0, nex+1, pcomp(nex), elevelcomp(nex), spincomp(nex), ejeccomp(nex), &
 &    masscomp(nex), prodZcomp(nex)
  enddo
  do nex = 0, nsp1
    Zix = parZ(typecomp(nex))
    Nix = parN(typecomp(nex))
    eopt = Einc - real(elevelcomp(nex) / specmass(Zix, Nix, typecomp(nex)))
    eopt = max(eopt, 0.001)
    kopt = typecomp(nex)
    call optical(Zix, Nix, kopt, eopt)
    write(1, '(3f10.5)') v, rv, av
    write(1, '(3f10.5)') w, rw, aw
    write(1, '(3f10.5)') vd, rvd, avd
    write(1, '(3f10.5)') wd, rwd, awd
    write(1, '(3f10.5)') vso, rvso, avso
    write(1, '(3f10.5)') wso, rwso, awso
    write(1, '(3f10.5)') rc, 0., 0.
    write(1, '(3f10.5)') 0., 0., 0.
  enddo
  write(1, '(3f10.5)') angbeg, anginc, angend
  write(1, '(f10.5)') bz1
  if (parinclude(0)) write(1, '(es10.3, 4f10.5)') tgo, S(0, 0, 1), 0., egr(0, 0, 1, 1, 1), ggr(0, 0, 1, 1, 1)
  do nex = 0, ncont
    if (parskip(0) .and. nex == 0) cycle
    write(1, '(7es10.3)') real(Zinit), aldcomp(nex), Umcomp(nex), tempcomp(nex), 0., E0comp(nex), Excomp(nex)
  enddo
  write(1, '(3i5)') 1, 1, 0
  Zix = Zindex(0, 0, k0)
  Nix = Nindex(0, 0, k0)
  if (disp(Zix, Nix, k0)) then
    do i = 1, npp
      write(1, '(10x, 2i5)') 2, 2
      write(1, '(10x, f10.5, 40x, f10.5)') ef(Zix, Nix, k0), w2(Zix, Nix, k0)
      write(1, '(20x, 2f10.5)') d3(Zix, Nix, k0), d2(Zix, Nix, k0)
    enddo
  endif
  return
end subroutine eciscompound
! Copyright A.J. Koning 2021
