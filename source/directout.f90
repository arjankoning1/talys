subroutine directout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of direct reaction cross sections
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
!   sgl              ! single precision kind
! All global variables
!   numang           ! maximum number of angles
!   numlev2          ! maximum number of levels
! Variables for output
!   flagddx          ! flag for output of double - differential cross sections
!   flagspec         ! flag for output of spectra
! Variables for basic reaction
!   flagang          ! flag for output of angular distributions
! Variables for numerics
!   nangle           ! number of angles
!   nanglecont       ! number of angles for continuum
! Variables for main input
!   k0               ! index of incident particle
! Variables for discrete levels
!   nlev             ! number of levels for nucleus
! Variables for energy grid
!   angle            ! angle in degrees
!   anglecont        ! angle in degrees for continuum
!   ebegin           ! first energy point of energy grid
!   egrid            ! outgoing energy grid
! Variables for energies
!   eend             ! last energy point of energy grid
!   eoutdis          ! outgoing energy of discrete state reaction
!   flaggiant        ! flag for collective contribution from giant resonances
! Variables for incident channel
!   directad         ! direct angular distribution
!   xscollconttot    ! total collective cross section in the continuum
!   xsdirdisc        ! direct cross section for discrete state direct cross section
!   xsdirdisctot     ! direct cross section summed over discrete states
!   xsgr             ! total smoothed giant resonance cross section
!   xsgrtot          ! total smoothed giant resonance cross section
! Variables for nuclides
!   Nindex           ! neutron number index for residual nucleus
!   Zindex           ! charge number index for residual nucleus
! Variables for giant resonances
!   eoutgr           ! emission energy
!   grcollad         ! giant resonance angular distribution
!   xscollcont       ! total collective cross section in the continuum
!   xsgrcoll         ! giant resonance cross section
!   xsgrstate        ! smoothed giant resonance cross section per state
! Constants
!   cparity          ! parity (character)
! Variables for deformation parameters
!   betagr           ! deformation parameter for giant resonance
!   deform           ! deformation parameter
!   deftype          ! deformation length (D) or parameter (B)
!   Egrcoll          ! energy of giant resonance
!   Ggrcoll          ! width of giant resonance
! Variables for levels
!   edis             ! energy of level
!   jdis             ! spin of level
!   parlev           ! parity of level
!
! *** Declaration of local data
!
  implicit none
  integer   :: i                                 ! counter
  integer   :: iang                              ! running variable for angle
  integer   :: ilev                              ! counter for discrete levels
  integer   :: n1                                ! number of coordinate grid points
  integer   :: n2                                ! counter
  integer   :: nen                               ! energy counter
  integer   :: Nix                               ! neutron number index for residual nucleus
  integer   :: nrest                             ! help variable
  integer   :: nset                              ! help variable
  integer   :: plev(numlev2)                     ! parity of level
  integer   :: Zix                               ! charge number index for residual nucleus
  real(sgl) :: dad(numlev2, 0:numang)            ! direct angular distribution
  real(sgl) :: elev(numlev2)                     ! energy of level
  real(sgl) :: jlev(numlev2)                     ! spin of level
!
! *********************** Inelastic cross sections *********************
!
!
  write(*, '(/" ++++++++++ DIRECT CROSS SECTIONS ++++++++++"/)')
  write(*, '(" Direct inelastic cross sections"/)')
  write(*, '(" Level  Energy   E-out      J/P   Cross section", "  Def. par."/)')
  Zix = Zindex(0, 0, k0)
  Nix = Nindex(0, 0, k0)
  do i = 1, numlev2
    if (xsdirdisc(k0, i) /= 0.) write(*, '(1x, i3, 2f10.5, f7.1, a1, f12.5, 5x, a1, f9.5)') i, edis(Zix, Nix, i), eoutdis(k0, i), &
 &    jdis(Zix, Nix, i), cparity(parlev(Zix, Nix, i)), xsdirdisc(k0, i), deftype(Zix, Nix), deform(Zix, Nix, i)
  enddo
  write(*, '(/" Discrete direct inelastic cross section:", f12.5, "   Level 1-", i3)') xsdirdisctot(k0), nlev(Zix, Nix)
  write(*, '(" Collective cross section in continuum  :", f12.5)') xscollconttot(k0)
  if (flagang) then
    write(*, '(/" Direct inelastic angular distributions")')
    ilev = 0
    do i = 1, numlev2
      if (xsdirdisc(k0, i) /= 0.) then
        ilev = ilev + 1
        elev(ilev) = edis(Zix, Nix, i)
        jlev(ilev) = jdis(Zix, Nix, i)
        plev(ilev) = parlev(Zix, Nix, i)
        do iang = 0, nangle
          dad(ilev, iang) = directad(k0, i, iang)
        enddo
      endif
    enddo
    nset = ilev / 10
    nrest = mod(ilev, 10)
    do n1 = 1, nset
      n2 = 10 * (n1 - 1)
      write(*, '(/" Angle", 10(" Ex=", f6.3, "  "))') (elev(i), i = n2 + 1, n2 + 10)
      write(*, '(4x, 10("    JP=", f4.1, a1)/)') (jlev(i), cparity(plev(i)), i = n2 + 1, n2 + 10)
      do iang = 0, nangle
        write(*, '(1x, f5.1, 10es12.5)') angle(iang), (dad(i, iang), i = n2 + 1, n2 + 10)
      enddo
    enddo
    if (nrest > 0) then
      write(*, '(/" Angle  ", 10("Ex=", f6.3, "   "))') (elev(i), i = 10 * nset + 1, 10 * nset + nrest)
      write(*, '(4x, 10("    JP=", f4.1, a1)/)') (jlev(i), cparity(plev(i)), i = 10 * nset + 1, 10 * nset + nrest)
      do iang = 0, nangle
        write(*, '(1x, f5.1, 10es12.5)') angle(iang), (dad(i, iang), i = nset + 1, nset + nrest)
      enddo
    endif
  endif
!
! *********************** Giant resonances *****************************
!
  if ( .not. flaggiant) return
  write(*, '(/" ++++++++++ GIANT RESONANCES ++++++++++"/)')
  write(*, '("      Cross section   Exc. energy Emis. energy", "   Width    Deform. par."/)')
  write(*, '(" GMR  :", 5f12.5)') xsgrcoll(k0, 0, 1), Egrcoll(0, 1), eoutgr(k0, 0, 1), Ggrcoll(0, 1), betagr(0, 1)
  write(*, '(" GQR  :", 5f12.5)') xsgrcoll(k0, 2, 1), Egrcoll(2, 1), eoutgr(k0, 2, 1), Ggrcoll(2, 1), betagr(2, 1)
  write(*, '(" LEOR :", 5f12.5)') xsgrcoll(k0, 3, 1), Egrcoll(3, 1), eoutgr(k0, 3, 1), Ggrcoll(3, 1), betagr(3, 1)
  write(*, '(" HEOR :", 5f12.5)') xsgrcoll(k0, 3, 2), Egrcoll(3, 2), eoutgr(k0, 3, 2), Ggrcoll(3, 2), betagr(3, 2)
  write(*, '(/" Total:", f12.5/)') xsgrtot(k0)-xscollconttot(k0)
  if (flagddx) then
    write(*, '(" Average angular distributions", /)')
    write(*, '(" Angle    GMR         GQR         LEOR      HEOR"/)')
    do iang = 0, nanglecont
      write(*, '(1x, f5.1, 4es12.5)') anglecont(iang), grcollad(k0, 0, 1, iang), grcollad(k0, 2, 1, iang), &
 &      grcollad(k0, 3, 1, iang), grcollad(k0, 3, 2, iang)
    enddo
  endif
  if (flagspec) then
    write(*, '(/" Giant resonance spectra", /)')
    write(*, '("   Energy   Total       GMR        GQR       ", "LEOR       HEOR     Collective"/)')
    do nen = ebegin(k0), eend(k0)
      write(*, '(1x, f8.3, 6es11.4)') egrid(nen), xsgr(k0, nen), xsgrstate(k0, 0, 1, nen), xsgrstate(k0, 2, 1, nen), &
 &      xsgrstate(k0, 3, 1, nen), xsgrstate(k0, 3, 2, nen), xscollcont(k0, nen)
    enddo
  endif
  return
end subroutine directout
! Copyright A.J. Koning 2021
