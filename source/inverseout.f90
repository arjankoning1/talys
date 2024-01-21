subroutine inverseout(Zcomp, Ncomp)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Reaction output for outgoing channels
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
! Variables for direct reactions
!   flagtransen    ! flag for output of transmission coefficients per energy
! Variables for energy grid
!   ebegin         ! first energy point of energy grid
!   egrid          ! outgoing energy grid
! Variables for energies
!   eend           ! last energy point of energy grid
! Variables for inverse channel data
!   Tjl            ! transmission coefficient per particle, energy, spin and l - value
!   Tl             ! transmission coefficients per particle, energy and l - value
!   xselas         ! total elastic cross section (shape + compound)
!   xsopt          ! optical model reaction cross section
!   xsreac         ! reaction cross section
!   xstot          ! total cross section (neutrons only)
!  Variables for gamma-ray strength functions
!   lmax        ! maximal l - value for transmission coefficients
! Variables for nuclides
!   Nindex         ! neutron number index for residual nucleus
!   parskip        ! logical to skip outgoing particle
!   Zindex         ! charge number index for residual nucleus
! Constants
!   parname        ! name of particle
! Variables for masses
!   specmass       ! specific mass for residual nucleus
!
! *** Declaration of local data
!
  implicit none
  integer   :: l                ! multipolarity
  integer   :: Ncomp            ! neutron number index for compound nucleus
  integer   :: nen              ! energy counter
  integer   :: Nix              ! neutron number index for residual nucleus
  integer   :: type             ! particle type
  integer   :: Zcomp            ! proton number index for compound nucleus
  integer   :: Zix              ! charge number index for residual nucleus
  real(sgl) :: e                ! energy
!
! **************** Transmission coefficients per energy ****************
!
  write(*, '(/" ########## TRANSMISSION COEFFICIENTS AND", " INVERSE REACTION CROSS SECTIONS ##########")')
!
! For each energy, the whole set of transmission coefficients is given as a function of the l-value and spin value.
!
  if (flagtransen) then
    do type = 1, 6
      if (parskip(type)) cycle
      Zix = Zindex(Zcomp, Ncomp, type)
      Nix = Nindex(Zcomp, Ncomp, type)
      do nen = ebegin(type), eend(type)
        e = real(egrid(nen) / specmass(Zix, Nix, type))
        write(*, '(/" Transmission coefficients for incident ", a8, " at ", f10.5, " MeV"/)') parname(type), e
!
! 1. Spin 1/2 particles: Neutrons, protons, tritons and Helium-3
!
        if (type /= 3 .and. type /= 6 .and. lmax(type, nen) >= 0) then
          write(*, '("   L   T(L-1/2,L)   T(L+1/2,L)    Tav(L)"/)')
          do l = 0, lmax(type, nen)
            write(*, '(1x, i3, 3es13.5)') l, Tjl(type, nen, -1, l), Tjl(type, nen, 1, l), Tl(type, nen, l)
          enddo
        endif
!
! 2. Spin 1 particles: Deuterons
!
        if (type == 3 .and. lmax(type, nen) >= 0) then
          write(*, '("   L    T(L-1,L)     T(L,L)       ", "T(L+1,L)     Tav(L)"/)')
          do l = 0, lmax(type, nen)
            write(*, '(1x, i3, 4es13.5)') l, Tjl(type, nen, -1, l), Tjl(type, nen, 0, l), Tjl(type, nen, 1, l), Tl(type, nen, l)
          enddo
        endif
!
! 3. Spin 0 particles: Alpha-particles
!
        if (type == 6 .and. lmax(type, nen) >= 0) then
          write(*, '("   L     T(L)"/)')
          do l = 0, lmax(type, nen)
            write(*, '(1x, i3, es13.5)') l, Tjl(type, nen, 0, l)
          enddo
        endif
      enddo
    enddo
!
! ************ Transmission coefficients per angular momentum **********
!
  else
!
! For each l-value, the whole set of transmission coefficients is given as a function of the energy and spin value.
!
    do type = 1, 6
      if (parskip(type)) cycle
      Zix = Zindex(Zcomp, Ncomp, type)
      Nix = Nindex(Zcomp, Ncomp, type)
      if (ebegin(type) >= eend(type)) cycle
      do l = 0, lmax(type, eend(type))
        write(*, '(/" Transmission coefficients for incident ", a8, " and l= ", i2/)') parname(type), l
!
! 1. Spin 1/2 particles: Neutrons, protons, tritons and Helium-3
!
        if (type /= 3 .and. type /= 6) then
          write(*, '("    Energy    T(L-1/2,L)   T(L+1/2,L)     ", "Tav(L)"/)')
          do nen = ebegin(type), eend(type)
            e = real(egrid(nen) / specmass(Zix, Nix, type))
            write(*, '(1x, f10.5, 3es13.5)') e, Tjl(type, nen, - 1, l), Tjl(type, nen, 1, l), Tl(type, nen, l)
          enddo
        endif
!
! 2. Spin 1 particles: Deuterons
!
        if (type == 3) then
          write(*, '("    Energy     T(L-1,L)     T(L,L)       ", "T(L+1,L)     Tav(L)"/)')
          do nen = ebegin(type), eend(type)
            e = real(egrid(nen) / specmass(Zix, Nix, type))
            write(*, '(1x, f10.5, 4es13.5)') e, Tjl(type, nen, - 1, l), Tjl(type, nen, 0, l), Tjl(type, nen, 1, l), Tl(type, nen, l)
          enddo
        endif
!
! 3. Spin 0 particles: Alpha-particles
!
        if (type == 6) then
          write(*, '("    Energy      T(L)"/)')
          do nen = ebegin(type), eend(type)
            e = real(egrid(nen) / specmass(Zix, Nix, type))
            write(*, '(1x, f10.5, es13.5)') e, Tjl(type, nen, 0, l)
          enddo
        endif
      enddo
    enddo
  endif
!
! **************** Cross sections for inverse channels *****************
!
  do type = 1, 6
    if (parskip(type)) cycle
    if (ebegin(type) >= eend(type)) cycle
    Zix = Zindex(Zcomp, Ncomp, type)
    Nix = Nindex(Zcomp, Ncomp, type)
    if (type == 1) then
      write(*, '(/" Total cross sections for ", a8/)') parname(1)
      write(*, '("      E        total      reaction    elastic", "   OMP reaction"/)')
      do nen = ebegin(1), eend(type)
        e = real(egrid(nen) / specmass(Zix, Nix, type))
        write(*, '(1x, f10.5, 4es12.4)') e, xstot(1, nen), xsreac(1, nen), xselas(1, nen), xsopt(1, nen)
      enddo
    else
      write(*, '(/" Total cross sections for ", a8/)') parname(type)
      write(*, '("      E       reaction  OMP reaction"/)')
      do nen = ebegin(type), eend(type)
        e = real(egrid(nen) / specmass(Zix, Nix, type))
        write(*, '(1x, f10.5, 2es12.4)') e, xsreac(type, nen), xsopt(type, nen)
      enddo
    endif
  enddo
  return
end subroutine inverseout
! Copyright A.J. Koning 2021
