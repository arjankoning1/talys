subroutine incidentread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read ECIS results for incident energy
!
! Author    : Arjan Koning, Eric Bauge and Pascal Romain
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl              ! single precision kind
!   dbl              ! double precision kind
! All global variables
!   numl             ! number of l values
! Variables for numerics
!   nangle           ! number of angles
!   transeps         ! absolute limit for transmission coefficient
! Variables for direct reactions
!   flaginccalc      ! flag for new ECIS calculation for incident channel
!   flagspher        ! flag to force spherical optical model
! Variables for input energies
!   nin              ! counter for incident energy
!   Ninc           ! number of incident energies
! Variables for main input
!   k0               ! index of incident particle
!   Ltarget          ! excited level of target
! Variables for energy grid
!   ecisstatus       ! status of ECIS file
!   translimit       ! limit for transmission coefficient
! Variables for incident channel
!   directad         ! direct angular distribution
!   dleg             ! direct reaction Legendre coefficient
!   dorigin          ! origin of direct cross section (Direct or Preeq)
!   lmaxinc          ! maximal l - value for transm. coeff. for incident channel
!   ruth             ! elastic / Rutherford ratio
!   Tjlinc           ! transm. coeff. as a function of spin and l for inc. channel
!   Tlinc            ! transm. coeff. as a function of l for incident channel
!   xscollconttot    ! total collective cross section in the continuum
!   xscoupled        ! inelastic cross section from coupled channels
!   xsdirdisc        ! direct cross section for discrete state direct cross section
!   xsdirdiscsum     ! total direct cross section
!   xsdirdisctot     ! direct cross section summed over discrete states
!   xselasinc        ! total elastic cross section (neutrons only) for inc. channel
!   xsoptinc         ! optical model reaction cross section for incident channel
!   xsreacinc        ! reaction cross section for incident channel
!   xstotinc         ! total cross section (neutrons only) for incident channel
! Variables for nuclides
!   Nindex           ! neutron number index for residual nucleus
!   Zindex           ! charge number index for residual nucleus
! Constants
!   parspin          ! spin of particle
! Variables for deformation parameters
!   colltype         ! type of collectivity (D, V or R)
!   indexcc          ! level index for coupled channel
! Variables for levels
!   jdis             ! spin of level
! Variables for level density
!   Nlast            ! last discrete level
!
! *** Declaration of local data
!
  implicit none
  character(len=72) :: line           ! input line
  integer           :: i              ! counter
  integer           :: iang           ! running variable for angle
  integer           :: ii             ! counter
  integer           :: infileang      ! file with angular distributions
  integer           :: infilecs       ! file with total cross sections
  integer           :: infilein       ! file with inelastic cross sections
  integer           :: infileleg      ! file with Legendre coefficients
  integer           :: infiletr       ! file with transmission coefficients
  integer           :: ispin          ! spin index
  integer           :: ist            ! counter
  integer           :: istat          ! logical for file access
  integer           :: itype          ! help variable
  integer           :: k              ! designator for particle
  integer           :: l              ! multipolarity
  integer           :: lev            ! level number
  integer           :: Nix            ! neutron number index for residual nucleus
  integer           :: nJ             ! number of total J values for transmission coefficients
  integer           :: nL             ! number of Legendre coefficients
  integer           :: nS             ! number of states
  integer           :: nSt            ! number of states
  integer           :: Zix            ! charge number index for residual nucleus
  real(sgl)         :: eps            ! help variable
  real(sgl)         :: factor         ! multiplication factor
  real(sgl)         :: groundspin2    ! 2 * spin of ground state
  real(sgl)         :: jres           ! j-value
  real(sgl)         :: rj             ! help variable
  real(sgl)         :: xsr            ! help variable
  real(dbl)         :: ddl            ! direct reaction Legendre coefficient
  real(dbl)         :: Tcoef          ! transmission coefficients as a function of spin and  l-value
  real(dbl)         :: teps           ! help variable
  real(dbl)         :: xs             ! help variable
!
! ************ Read total, reaction and elastic cross section **********
!
! If the ECIS calculation has already been done in a previous run, we can read from existing files.
!
  if (flaginccalc) then
    open (unit = 3, status = 'unknown', file = 'ecis.inccs')
    open (unit = 7, status = 'unknown', file = 'ecis.inctr')
    open (unit = 8, status = 'unknown', file = 'ecis.incang')
    open (unit = 9, status = 'unknown', file = 'ecis.incleg')
    open (unit = 10, status = 'unknown', file = 'ecis.incin')
    infilecs = 3
    infiletr = 7
    infileang = 8
    infileleg = 9
    infilein = 10
    open (unit = 13, status = 'unknown', file = 'incident.cs')
    open (unit = 17, status = 'unknown', file = 'incident.tr')
    open (unit = 18, status = 'unknown', file = 'incident.ang')
    open (unit = 19, status = 'unknown', file = 'incident.leg')
    open (unit = 20, status = 'unknown', file = 'incident.in')
    do
      read(3, '(a72)', iostat = istat) line
      if(istat == -1) exit
      write(13, '(a72)') line
    enddo
    rewind 3
    do
      read(7, '(a72)', iostat = istat) line
      if(istat == -1) exit
      write(17, '(a72)') line
    enddo
    rewind 7
    do
      read(8, '(a72)', iostat = istat) line
      if(istat == -1) exit
      write(18, '(a72)') line
    enddo
    rewind 8
    do
      read(9, '(a72)', iostat = istat) line
      if(istat == -1) exit
      write(19, '(a72)') line
    enddo
    rewind 9
    do
      read(10, '(a72)', iostat = istat) line
      if(istat == -1) exit
      write(20, '(a72)') line
    enddo
    rewind 10
  else
    infilecs = 13
    infiletr = 17
    infileang = 18
    infileleg = 19
    infilein = 20
  endif
  read(infilecs, '()')
  if (k0 == 1) then
    read(infilecs, *, iostat=istat ) xs
    if (istat /= 0) then
      write(*, '(" TALYS-warning: Unreadable cross section in file ", i3, " set to zero")') infilecs
      xs = 0.
    endif
    xstotinc = real(xs)
  endif
  read(infilecs, *, iostat=istat ) xs
  if (istat /= 0) then
    write(*, '(" TALYS-warning: Unreadable cross section in file ", i3, " set to zero")') infilecs
    xs = 0.
  endif
  xsreacinc = real(xs)
  xsoptinc = real(xs)
  if (k0 == 1) then
    read(infilecs, *, iostat=istat ) xs
    if (istat /= 0) then
      write(*, '(" TALYS-warning: Unreadable cross section in file ", i3, " set to zero")') infilecs
      xs = 0.
    endif
    xselasinc = real(xs)
  endif
  eps = -1.e-3
  if (xselasinc < eps .or. xsreacinc < eps .or. xstotinc < eps) then
    write(*, '(" TALYS-warning: Negative OMP cross section")')
    write(*, '(" Elastic : ", es12.5)') xselasinc
    write(*, '(" Reaction: ", es12.5)') xsreacinc
    write(*, '(" Total   : ", es12.5)') xstotinc
    if (xselasinc < eps) then
      xselasinc = 0.
      write(*, '(" --> Elastic : ", es12.5)') xselasinc
    endif
    if (xsreacinc < eps) then
      xsreacinc = 0.
      write(*, '(" --> Reaction: ", es12.5)') xsreacinc
    endif
    if (xstotinc < eps) then
      xstotinc = xselasinc + xsreacinc
      write(*, '(" --> Total   : ", es12.5)') xstotinc
    endif
    xsoptinc = xsreacinc
  endif
  xselasinc = max(xselasinc, 0.)
  xsreacinc = max(xsreacinc, 0.)
  xstotinc = max(xstotinc, 0.)

!
! ******************* Read transmission coefficients *******************
!
! For rotational nuclei, the rotational transmission coefficients are transformed into into spherical equivalents.
!
  Zix = Zindex(0, 0, k0)
  Nix = Nindex(0, 0, k0)
  groundspin2 = int(2. * jdis(Zix, Nix, 0))
  read(infiletr, '(55x, i5)') nJ
  do i = 1, nJ
    read(infiletr, '(f10.1, 5x, i5)') rj, nS
    do k = 1, nS
      read(infiletr, '(i3, i6, f9.1, e20.8)') lev, l, jres, Tcoef
      if (l > numl) cycle
      if (lev == 1) then
        if (colltype(Zix, Nix) /= 'S' .and. .not.flagspher) then
          factor = (2. * rj + 1.) / (2. * jres + 1.) / (groundspin2 + 1.)
        else
          factor = 1.
        endif
        if (parspin(k0) == 0.5) then
          ispin = int(2. * (jres - real(l)))
        else
          ispin = int(jres - real(l))
        endif
        Tjlinc(ispin, l) = Tjlinc(ispin, l) + factor * max(real(Tcoef), 0.)
      endif
    enddo
  enddo
!
! *************** Processing of transmission coefficients **************
!
! Transmission coefficients averaged over spin and determination of maximal l-value.
! ECIS stops its output of transmission coefficients somewhat too early.
! For the highest l values the transmission coefficient for (l+spin) is not written in the output.
! Since these are small numbers we put them equal to the value for (l-spin).
!
! 1. Spin 1/2 particles: Neutrons, protons, tritons and Helium-3
!
  if (k0 /= 3 .and. k0 /= 6) then
    do l = 0, numl
      if (Tjlinc( - 1, l) /= 0 .and. Tjlinc(1, l) == 0) Tjlinc(1, l) = Tjlinc( - 1, l)
      if (Tjlinc( - 1, l) == 0 .and. Tjlinc(1, l) /= 0 .and. l > 0) Tjlinc( - 1, l) = Tjlinc(1, l)
      Tlinc(l) = ((l + 1) * Tjlinc(1, l) + l * Tjlinc( - 1, l)) / (2 * l + 1)
      teps = Tlinc(0) * translimit / (2 * l + 1)
      teps = max(teps, transeps)
      if (Tjlinc( - 1, l) < teps .and. Tjlinc(1, l) < teps) then
        lmaxinc = l - 1
        exit
      endif
      lmaxinc = l
    enddo
  endif
!
! 2. Spin 1 particles: Deuterons
!
  if (k0 == 3) then
    do l = 0, numl
      if (Tjlinc( - 1, l) /= 0 .and. Tjlinc(0, l) == 0) Tjlinc(0, l) = Tjlinc( - 1, l)
      if (Tjlinc( - 1, l) /= 0 .and. Tjlinc(1, l) == 0) Tjlinc(1, l) = Tjlinc( - 1, l)
      if (Tjlinc( - 1, l) == 0 .and. Tjlinc(1, l) /= 0 .and. l > 0) Tjlinc( - 1, l) = Tjlinc(1, l)
      Tlinc(l) = ((2 * l + 3) * Tjlinc(1, l) + (2 * l + 1) * Tjlinc(0, l) + &
        (2 * l - 1) * Tjlinc( - 1, l)) / (3 * (2 * l + 1))
      teps = Tlinc(0) * translimit / (2 * l + 1)
      teps = max(teps, transeps)
      if (Tjlinc( - 1, l) < teps .and. Tjlinc(0, l) < teps .and. Tjlinc(1, l) < teps) then
        lmaxinc = l - 1
        exit
      endif
      lmaxinc = l
    enddo
  endif
!
! 3. Spin 0 particles: Alpha-particles
!
  if (k0 == 6) then
    do l = 0, numl
      Tlinc(l) = Tjlinc(0, l)
      teps = Tlinc(0) * translimit / (2 * l + 1)
      teps = max(teps, transeps)
      if (Tlinc(l) < teps) then
        lmaxinc = l - 1
        exit
      endif
      lmaxinc = l
    enddo
  endif
!
! ******************* Direct reaction Legendre coefficients ************
!
! We read the Legendre coefficients for the direct component of the reaction only.
! The compound nucleus coefficients are calculated by TALYS later on.
! For coupled-channels reactions, the inelastic Legendre coefficients are also read.
!
  read(infileleg, '(55x, i5)') nSt
  do ist = 1, nSt
    read(infileleg, '(5x, i5)') nL
    do k = 1, nL
      read(infileleg, '(2i5, e20.10)') i, l, ddl
      if (l > 3 * numl) cycle
      ii = i - 1
      if (i /= 1) ii = indexcc(Zix, Nix, i)
      dleg(k0, ii, l) = real(ddl)
    enddo
  enddo
!
! ******************* Read elastic angular distributions ***************
!
! For charged particles, we also read the elastic/rutherford ratio.
!
  read(infileang, '(55x, i5)') nSt
  read(infileang, '(12x, i3)') nS
  do iang = 0, nangle
    do k = 1, nS
      read(infileang, '(i3, 12x, e12.5)', iostat = istat) itype, xs
      if (istat /= 0) cycle
      xsr = 1.e38
      if (xs <= 1.e38) xsr = real(xs)
      if (itype == 0) directad(k0, Ltarget, iang) = xsr
      if (itype == 1) ruth(iang) = xsr
    enddo
  enddo
!
! For coupled-channels calculations, we also read the discrete inelastic angular distributions and cross sections.
!
  if (colltype(Zix, Nix) /= 'S' .and. .not.flagspher) then
    xscoupled = 0.
    do ist = 2, nSt
      read(infileang, '(i5, 7x, i3)') i, nS
      do iang = 0, nangle
        do k = 1, nS
          read(infileang, '(i3, 12x, e12.5)', iostat = istat) itype, xs
          if (istat /= 0) cycle
          ii = indexcc(Zix, Nix, i)
          if (itype == 0) directad(k0, ii, iang) = real(xs)
        enddo
      enddo
    enddo
    read(infilein, '(55x, i5)') nSt
    do ist = 1, nSt
      read(infilein, '(e12.5, i3)') xs, i
      ii = indexcc(Zix, Nix, i + 1)
      xsdirdisc(k0, ii) = real(xs)
      if (ii <= Nlast(Zix, Nix, 0)) then
        xsdirdisctot(k0) = xsdirdisctot(k0) + xsdirdisc(k0, ii)
      else
        xscollconttot(k0) = xscollconttot(k0) + xsdirdisc(k0, ii)
      endif
      xscoupled = xscoupled + xsdirdisc(k0, ii)
      dorigin(k0, ii) = 'Direct'
    enddo
    xsdirdiscsum = xsdirdisctot(k0)
  endif
!
! Close files
!
  close (unit = 3, status = ecisstatus)
  close (unit = 7, status = ecisstatus)
  close (unit = 8, status = ecisstatus)
  close (unit = 9, status = ecisstatus)
  close (unit = 10, status = ecisstatus)
  if (nin == Ninc) then
    close (unit = 13, status = ecisstatus)
    close (unit = 17, status = ecisstatus)
    close (unit = 18, status = ecisstatus)
    close (unit = 19, status = ecisstatus)
    close (unit = 20, status = ecisstatus)
  endif
  return
end subroutine incidentread
! Copyright A.J. Koning 2021
