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
  character(len=3)   :: massstring !
  character(len=6)   :: finalnuclide !
  character(len=8)   :: tjlfile !
  character(len=7)   :: crossfile !
  character(len=15)  :: col(5)    ! header
  character(len=15)  :: un(5)    ! header
  character(len=80)  :: quantity   ! quantity
  character(len=132) :: topline   ! topline
  integer   :: l                ! multipolarity
  integer   :: Ncomp            ! neutron number index for compound nucleus
  integer   :: nen              ! energy counter
  integer   :: Nix              ! neutron number index for residual nucleus
  integer   :: Ncol
  integer   :: type             ! particle type
  integer   :: Zcomp            ! proton number index for compound nucleus
  integer   :: Zix              ! charge number index for residual nucleus
  integer   :: Z
  integer   :: A
  real(sgl) :: e                ! energy
!
! **************** Transmission coefficients per energy ****************
!
  write(*, '(/" ########## TRANSMISSION COEFFICIENTS AND INVERSE REACTION CROSS SECTIONS ##########")')
!
! For each energy, the whole set of transmission coefficients is given as a function of the l-value and spin value.
!
  do type = 1, 6
    if (parskip(type)) cycle
    if (ebegin(type) >= eend(type)) cycle
    Zix = Zindex(Zcomp, Ncomp, type)
    Nix = Nindex(Zcomp, Ncomp, type)
    Z = ZZ(Zcomp, Ncomp, type)
    A = AA(Zcomp, Ncomp, type)
    massstring='   '
    write(massstring,'(i3)') A
    finalnuclide=trim(nuc(Z))//adjustl(massstring)
    un=''
    col(1)='L'
    if (type /= 3 .and. type /= 6) then
      col(2)='T(L-1/2,L)'
      col(3)='T(L+1/2,L)'
      col(4)='T(L))'
      Ncol=4
    endif
    if (type == 3) then
      col(2)='T(L-1,L)'
      col(3)='T(L,L)'
      col(4)='T(L+1,L)'
      col(5)='T(L))'
      Ncol=5
    endif
    if (type == 6) then
      col(2)='T(L))'
      Ncol=2
    endif
    quantity=parname(type)//' transmission coefficients'
    topline=trim(finalnuclide)//' '//trim(quantity)
    tjlfile='transm.'//parsym(type)
    open (unit=1, file=tjlfile, status='replace')
    call write_header(topline,source,user,date,oformat)
    call write_residual(Z,A,finalnuclide)
    write(1,'("# parameters:")')
    call write_char(2,'particle',parname(type))
    if (flagtransen) then
      Nen =  eend(type) - ebegin(type) + 1
      call write_integer(2,'number of energies',Nen)
      do nen = ebegin(type), eend(type)
        e = real(egrid(nen) / specmass(Zix, Nix, type))
        write(1,'("# parameters:")')
        call write_real(2,'energy [MeV]',e)
        call write_datablock(quantity,Ncol,lmax(type, nen),col,un)
        write(*, '(/" Transmission coefficients for incident ", a8, " at ", f10.5, " MeV"/)') parname(type), e
!
! 1. Spin 1/2 particles: Neutrons, protons, tritons and Helium-3
!
        if (type /= 3 .and. type /= 6 .and. lmax(type, nen) >= 0) then
          write(*, '("   L   T(L-1/2,L)   T(L+1/2,L)    Tav(L)"/)')
          do l = 0, lmax(type, nen)
            write(1, '(i6, 9x, 3es15.6)') l, Tjl(type, nen, -1, l), Tjl(type, nen, 1, l), Tl(type, nen, l)
            write(*, '(1x, i3, 3es13.5)') l, Tjl(type, nen, -1, l), Tjl(type, nen, 1, l), Tl(type, nen, l)
          enddo
        endif
!
! 2. Spin 1 particles: Deuterons
!
        if (type == 3 .and. lmax(type, nen) >= 0) then
          write(*, '("   L    T(L-1,L)     T(L,L)       ", "T(L+1,L)     Tav(L)"/)')
          do l = 0, lmax(type, nen)
            write(1, '(i6, 9x, 4es15.6)') l, Tjl(type, nen, -1, l), Tjl(type, nen, 0, l), Tjl(type, nen, 1, l), Tl(type, nen, l)
            write(*, '(1x, i3, 4es13.5)') l, Tjl(type, nen, -1, l), Tjl(type, nen, 0, l), Tjl(type, nen, 1, l), Tl(type, nen, l)
          enddo
        endif
!
! 3. Spin 0 particles: Alpha-particles
!
        if (type == 6 .and. lmax(type, nen) >= 0) then
          write(*, '("   L     T(L)"/)')
          do l = 0, lmax(type, nen)
            write(1, '(i9, 6x, es15.6)') l, Tjl(type, nen, 0, l)
            write(*, '(1x, i3, es13.5)') l, Tjl(type, nen, 0, l)
          enddo
        endif
      enddo
    else
!
! ************ Transmission coefficients per angular momentum **********
!
!
! For each l-value, the whole set of transmission coefficients is given as a function of the energy and spin value.
!
      Nen =  lmax(type, eend(type)) + 1
      call write_integer(2,'number of angular L-values',Nen)
      do l = 0, lmax(type, eend(type))
        write(*, '(/" Transmission coefficients for incident ", a8, " and l= ", i2/)') parname(type), l
        write(1,'("# parameters:")')
        call write_integer(2,'L-value',l)
        Nen =  eend(type) - ebegin(type) + 1
        call write_datablock(quantity,Ncol,Nen,col,un)
!
! 1. Spin 1/2 particles: Neutrons, protons, tritons and Helium-3
!
        if (type /= 3 .and. type /= 6) then
          write(*, '("    Energy    T(L-1/2,L)   T(L+1/2,L)     ", "Tav(L)"/)')
          do nen = ebegin(type), eend(type)
            e = real(egrid(nen) / specmass(Zix, Nix, type))
            write(1, '(4es15.6)') e, Tjl(type, nen, - 1, l), Tjl(type, nen, 1, l), Tl(type, nen, l)
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
            write(1, '(5es15.6)') e, Tjl(type, nen, - 1, l), Tjl(type, nen, 0, l), Tjl(type, nen, 1, l), Tl(type, nen, l)
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
            write(1, '(2es15.6)') e, Tjl(type, nen, 0, l)
            write(*, '(1x, f10.5, es13.5)') e, Tjl(type, nen, 0, l)
          enddo
        endif
      enddo
    endif
    close(1)
!
! **************** Cross sections for inverse channels *****************
!
    quantity=parname(type)//' inverse reaction cross sections'
    topline=trim(finalnuclide)//' '//trim(quantity)
    crossfile='cross.'//parsym(type)
    un='mb'
    col(1)='E'
    un(1)='MeV'
    if (type == 1) then
      col(2)='total'
      col(3)='elastic'
      col(4)='reaction'
      col(5)='OMP reaction'
      Ncol=5
    else
      col(2)='reaction'
      col(3)='OMP reaction'
      Ncol=3
    endif
    open (unit=1, file=crossfile, status='replace')
    call write_header(topline,source,user,date,oformat)
    call write_residual(Z,A,finalnuclide)
    write(1,'("# parameters:")')
    call write_char(2,'particle',parname(type))
    Nen =  eend(type) - ebegin(type) + 1
    call write_datablock(quantity,Ncol,Nen,col,un)
    if (type == 1) then
      write(*, '(/" Total cross sections for ", a8/)') parname(1)
      write(*, '("      E        total      reaction    elastic", "   OMP reaction"/)')
      do nen = ebegin(type), eend(type)
        e = real(egrid(nen) / specmass(Zix, Nix, type))
        write(1, '(5es15.6)') e, xstot(1, nen), xselas(1, nen), xsreac(1, nen), xsopt(1, nen)
        write(*, '(1x, f10.5, 4es12.4)') e, xstot(1, nen), xsreac(1, nen), xselas(1, nen), xsopt(1, nen)
      enddo
    else
      write(*, '(/" Total cross sections for ", a8/)') parname(type)
      write(*, '("      E       reaction  OMP reaction"/)')
      do nen = ebegin(type), eend(type)
        e = real(egrid(nen) / specmass(Zix, Nix, type))
        write(1, '(3es15.6)') e, xsreac(type, nen), xsopt(type, nen)
        write(*, '(1x, f10.5, 2es12.4)') e, xsreac(type, nen), xsopt(type, nen)
      enddo
    endif
    close(1)
  enddo
  return
end subroutine inverseout
! Copyright A.J. Koning 2021
