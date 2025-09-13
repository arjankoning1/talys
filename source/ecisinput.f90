subroutine ecisinput(Zix, Nix, kopt, e, rotational, vibrational, jlmloc)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Create ECIS input file
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
!   sgl          ! single precision kind
! Variables for ECIS
!   angbeg       ! first angle
!   angend       ! last angle
!   anginc       ! angle increment
!   d2disp       ! constant for imaginary potential
!   d3disp       ! constant for imaginary potential
!   ecis1        ! 50 input flags ('T' or 'F') for ECIS
!   ecis2        ! 50 input flags ('T' or 'F') for ECIS
!   efer         ! Fermi energy
!   Elevel       ! energy of level
!   ecisstep     ! integration step size for ECIS OMP calculation
!   iband        ! band number of level
!   idvib        ! identifier for existence of vibrational state inside rotational
!   iph          ! help variable
!   iqm          ! largest order of deformation
!   iqmax        ! maximum l - value of multipole expansion
!   iterm        ! number of iterations
!   Jband        ! angular momentum
!   Jlevel       ! spin of level
!   Kmag         ! magnetic quantum number
!   legendre     ! logical for output of Legendre coefficients
!   Nband        ! number of vibrational bands
!   ncoll        ! number of nuclear states
!   njmax        ! maximal number of j - values in ECIS
!   npp          ! number of optical potentials
!   nrad         ! number of radial points
!   Nrotbeta     ! number of deformation parameters for rotational nucleus
!   rotbeta      ! deformation parameters for rotational nucleus
!   Plevel       ! parity of level
!   prodZ        ! product of charges of projectile and target nucleus
!   projmass     ! mass of projectile
!   resmass      ! mass of residual nucleus
!   rmatch       ! matching radius
!   spin         ! spin of incident particle
!   tarparity    ! parity of target nucleus
!   tarspin      ! spin of target nucleus
!   title        ! title of ECIS input file
!   vibbeta      ! vibrational deformation parameter
!   w2disp       ! constant for imaginary potential
! Variables for JLM
!   normjlm      ! JLM potential normalization factors
!   potjlm       ! JLM potential depth values
!   radjlm       ! radial points for JLM potential
! Variables for optical model
!   av           ! real volume diffuseness
!   avd          ! real surface diffuseness
!   avso         ! real spin - orbit diffuseness
!   aw           ! imaginary volume diffuseness
!   awd          ! imaginary surface diffuseness
!   awso         ! imaginary spin - orbit diffuseness
!   disp         ! flag for dispersive optical model
!   rc           ! Coulomb radius
!   rv           ! real volume radius
!   rvd          ! real surface radius
!   rvso         ! real spin - orbit radius
!   rw           ! imaginary volume radius
!   rwd          ! imaginary surface radius
!   rwso         ! imaginary spin - orbit radius
!   v            ! real volume depth
!   vd           ! real surface depth
!   vso          ! real spin - orbit depth
!   w            ! imaginary volume depth
!   wd           ! imaginary surface depth
!   wso          ! imaginary spin - orbit depth
!
! *** Declaration of local data
!
  implicit none
  logical           :: jlmloc         ! flag for JLM OMP
  logical           :: rotational     ! flag for rotational input
  logical           :: vibrational    ! flag for vibrational input
  character(len=32) :: Eformat        ! format string
  integer           :: i              ! level
  integer           :: k              ! designator for particle
  integer           :: kopt           ! optical model index for fast particle
  integer           :: Nix            ! neutron number index for residual nucleus
  integer           :: Zix            ! charge number index for residual nucleus
  real(sgl)         :: e              ! energy
  real(sgl)         :: eopt           ! incident energy
!
! *********************** Write standard input *************************
!
!   rotational model
!
  write(9, '(a72)') title
  write(9, '(a50)') ecis1
  write(9, '(a50)') ecis2
  write(9, '(4i5)') ncoll, njmax, iterm, npp
  write(9, '(2f10.5, 10x, "    1.e-10    1.e-10    1.e-30")') ecisstep, rmatch
  if (legendre) write(9, '()')
  if (e >= 0.01) then
    Eformat='(f5.2,2i2,a1,5f10.5)'
  else
    Eformat='(f5.2,2i2,a1,es10.3,4f10.5)'
  endif
  write(9, fmt = Eformat) tarspin, idvib(1), 1, tarparity, e, spin, projmass, resmass, prodZ
  if (vibrational) write(9, '()')
  do i = 2, ncoll
    write(9, '(f5.2, 2i2, a1, 5f10.5)') Jlevel(i), idvib(i), min(i, npp), Plevel(i), Elevel(i)
    if (vibrational) then
      write(9, '(3i5)') iph(i), iband(i), 1
    else
      if (idvib(i) >= 1) write(9, '(2i5)') iph(i), iband(i)
      if (ecis1(3:3) == 'T') write(9, '()')
    endif
  enddo
  do i = 1, Nband
    write(9, '(2i5, f10.5)') Jband(i), Kmag(i), vibbeta(i)
  enddo
  if (rotational) then
    write(9, '(2i5, f10.1)') iqm, iqmax, tarspin
    write(9, '(7f10.5)') (rotbeta(i), i = 1, Nrotbeta)
  endif
!
! The optical model parameters can be calculated at the ground state
! and every excited state.
!
! 1. Phenomenological OMP
!
! optical      : subroutine for determination of optical potential
!
  if ( .not. jlmloc) then
    do i = 1, npp
      eopt = e - real(Elevel(i) * (resmass + projmass) / resmass)
      call optical(Zix, Nix, kopt, eopt)
      if (abs(v) >= 1000.) then
        write(9, '(es10.3, 2f10.5)') v, rv, av
      else
        write(9, '(3f10.5)') v, rv, av
      endif
      if (abs(w) >= 1000.) then
        write(9, '(es10.3, 2f10.5)') w, rw, aw
      else
        write(9, '(3f10.5)') w, rw, aw
      endif
      write(9, '(3f10.5)') vd, rvd, avd
      if (abs(wd) >= 1000.) then
        write(9, '(es10.3, 2f10.5)') wd, rwd, awd
      else
        write(9, '(3f10.5)') wd, rwd, awd
      endif
      write(9, '(3f10.5)') vso, rvso, avso
      write(9, '(3f10.5)') wso, rwso, awso
      write(9, '(3f10.5)') rc, 0., 0.
      write(9, '(3f10.5)') 0., 0., 0.
    enddo
    write(9, '(3f10.5)') angbeg, anginc, angend
    if (disp(Zix, Nix, kopt)) then
      do i = 1, npp
        write(9, '(10x, 2i5)') 2, 2
        write(9, '(10x, f10.5, 40x, f10.5)') efer, w2disp
        write(9, '(20x, 2f10.5)') d3disp, d2disp
      enddo
    endif
  endif
!
! 2. JLM OMP
!
! foldalpha: subroutine for double folding potential
! mom      : subroutine for microscopic optical model (Eric Bauge)
!
! Write the potentials in ECIS external input format
!
  if (jlmloc) then
    write(9, '(3f10.5)') angbeg, anginc, angend
    if (kopt == 6) then
      call foldalpha(Zix, Nix, e)
    else
      call mom(Zix, Nix, dble(prodZ), dble(e))
    endif
    write(9, '(2i5)') 1, 1
!
! 1.  real central
! 2.  Imaginary central
!
    do k = 1, 2
      write(9, '(9i5)') 1, 1, 0, k, 0, 0, 0, 1, -1
      write(9, '(f10.5)') normjlm(Zix, Nix, k)
      write(9, '(2(f10.5, es20.6))') (radjlm(Zix, Nix, i), - potjlm(Zix, Nix, i, k), i = 1, nrad - 2)
      write(9, '(2(f10.5, es20.6), a4)') (radjlm(Zix, Nix, i), - potjlm(Zix, Nix, i, k), i = nrad - 1, nrad), "last"
    enddo
!
! 3.  real surface (not used)
! 4.  Imaginary surface (not used)
!
    do k = 3, 4
      write(9, '(8i5)') 1, 1, 0, k, 0, 0, 0, -1
      write(9, '()')
    enddo
!
! 5.  real spin-orbit
! 6.  Imaginary spin-orbit
!
    if (kopt /= 6) then
      do k = 5, 6
        write(9, '(9i5)') 1, 1, 0, k, 0, 0, 0, 1, -1
        write(9, '(f10.5)') normjlm(Zix, Nix, k)
        write(9, '(2(f10.5, es20.6))') (radjlm(Zix, Nix, i), - 0.5 * potjlm(Zix, Nix, i, k), i = 1, nrad - 2)
        write(9, '(2(f10.5, es20.6), a4)') (radjlm(Zix, Nix, i), &
 &        - 0.5 * potjlm(Zix, Nix, i, k), i = nrad - 1, nrad), "last"
      enddo
    endif
!
! 7.  Coulomb
! 8.  Coulomb spin-orbit (not used)
!
    write(9, '(9i5)') 1, 1, 0, 7, 0, 0, 0, -1, -1
    write(9, '(2f10.5)') prodZ, rc
    if (kopt /= 6) then
      write(9, '(9i5)') 1, 1, 0, 8, 0, 0, 0, -1, -1
      write(9, '(3f10.5)') 0.0, 1.12, 0.55
    endif
  endif
  return
end subroutine ecisinput
! Copyright A.J. Koning 2021
