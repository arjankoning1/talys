subroutine ecisdwbamac(itype, type)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Create ECIS input file for macroscopic DWBA calculation for
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
!   sgl        ! single precision kind
! Variables for main input
!   Zinit      ! charge number of initial compound nucleus
! Variables for nuclides
!   Nindex     ! neutron number index for residual nucleus
!   Zindex     ! charge number index for residual nucleus
!   ZZ         ! charge number of residual nucleus
! Constants
!   parmass    ! mass of particle in a.m.u.
!   parN       ! neutron number of particle
!   parspin    ! spin of particle
!   parZ       ! charge number of particle
! Variables for ECIS
!   angbeg     ! first angle
!   angend     ! last angle
!   anginc     ! angle increment
!   ecis1      ! 50 input flags ('T' or 'F') for ECIS
!   ecis2      ! 50 input flags ('T' or 'F') for ECIS
!   iterm      ! number of iterations
!   ncoll      ! number of nuclear states
!   njmax      ! maximal number of j - values in ECIS
!   npp        ! number of optical potentials
!   rmatch     ! matching radius
!   title      ! title of ECIS input file
! Variables for masses
!   nucmass    ! mass of nucleus
! Variables for optical model
!   av         ! real volume diffuseness
!   avd        ! real surface diffuseness
!   avso       ! real spin - orbit diffuseness
!   aw         ! imaginary volume diffuseness
!   awd        ! imaginary surface diffuseness
!   awso       ! imaginary spin - orbit diffuseness
!   d2         ! parameter for imaginary surface OMP
!   d3         ! parameter for imaginary surface OMP
!   disp       ! flag for dispersive optical model
!   ef         ! Fermi energy
!   rc         ! Coulomb radius
!   rv         ! real volume radius
!   rvd        ! real surface radius
!   rvso       ! real spin - orbit radius
!   rw         ! imaginary volume radius
!   rwd        ! imaginary surface radius
!   rwso       ! imaginary spin - orbit radius
!   v          ! real volume depth
!   vd         ! real surface depth
!   vso        ! real spin - orbit depth
!   w          ! imaginary volume depth
!   w2         ! parameter for imaginary volume OMP
!   wd         ! imaginary surface depth
!   wso        ! imaginary spin - orbit depth
! Variables for MSD
!   betamsd    ! deformation parameter
!   Emsdin     ! incident MSD energy
!   Emsdout    ! outgoing MSD energy
!   Exmsd      ! excitation energy for MSD energy grid
!   maxJmsd    ! maximal spin for MSD calculation
!
! *** Declaration of local data
!
  implicit none
  character(len=1) :: ch        ! character
  integer          :: i         ! counter
  integer          :: itype     ! help variable
  integer          :: J         ! spin of level
  integer          :: kopt      ! optical model index for fast particle
  integer          :: Nix       ! neutron number index for residual nucleus
  integer          :: Nixi      ! help variable
  integer          :: type      ! particle type
  integer          :: Z         ! charge number of target nucleus
  integer          :: Zix       ! charge number index for residual nucleus
  integer          :: Zixi      ! help variable
  real(sgl)        :: eopt      ! incident energy
!
! **************************** Write input *****************************
!
  Zix = Zindex(0, 0, type)
  Nix = Nindex(0, 0, type)
  Z = ZZ(0, 0, type)
  if (disp(Zix, Nix, itype)) then
    ecis1(10:10) = 'T'
  else
    ecis1(10:10) = 'F'
  endif
  write(9, '(a72)') title
  write(9, '(a50)') ecis1
  write(9, '(a50)') ecis2
  write(9, '(4i5)') ncoll, njmax, iterm, npp
  write(9, '(10x, f10.5, 10x, 3("    1.e-10"))') rmatch
  write(9, '(f5.2, 2i2, a1, 5f10.5)') 0., 1, 2, '+', Emsdin, &
 &  parspin(itype), parmass(itype), nucmass(parZ(itype), parN(itype)), real(Zinit - parZ(itype)) * parZ(itype)
  do J = 0, maxJmsd
    if (mod(J, 2) == 0) then
      ch = '+'
    else
      ch = '-'
    endif
    write(9, '(f5.2, 2i2, a1, 5f10.5)') real(J), 0, 3, ch, Exmsd, parspin(type), parmass(type), nucmass(Zix, Nix), &
 &    real(Z * parZ(type))
    write(9, '(2i5)') 1, J+1
  enddo
  do J = 0, maxJmsd
    write(9, '(i5, 5x, f10.5)') J, betamsd
  enddo
  do i = 1, 3
!
! Transition potential and potential for incident channel
!
    if (i <= 2) then
      Zixi = parZ(itype)
      Nixi = parN(itype)
      kopt = itype
      if (i == 1) then
        eopt = 0.5 * (Emsdin + Emsdout)
      else
        eopt = Emsdin
      endif
    endif
!
! Potential for outgoing channel
!
! optical : subroutine for determination of optical potential
!
    if (i == 3) then
      Zixi = parZ(type)
      Nixi = parN(type)
      kopt = type
      eopt = Emsdout
    endif
    call optical(Zixi, Nixi, kopt, eopt)
    write(9, '(3f10.5)') v, rv, av
    write(9, '(3f10.5)') w, rw, aw
    write(9, '(3f10.5)') vd, rvd, avd
    write(9, '(3f10.5)') wd, rwd, awd
    write(9, '(3f10.5)') vso, rvso, avso
    write(9, '(3f10.5)') wso, rwso, awso
    write(9, '(3f10.5)') rc, 0., 0.
    write(9, '(3f10.5)') 0., 0., 0.
  enddo
  write(9, '(3f10.5)') angbeg, anginc, angend
  if (disp(Zixi, Nixi, itype)) then
    do i = 1, 3
      write(9, '(10x, 2i5)') 2, 2
      write(9, '(10x, f10.5, 40x, f10.5)') ef(Zixi, Nixi, itype), w2(Zixi, Nixi, itype)
      write(9, '(20x, 2f10.5)') d3(Zixi, Nixi, itype), d2(Zixi, Nixi, itype)
    enddo
  endif
  return
end subroutine ecisdwbamac
! Copyright A.J. Koning 2021
