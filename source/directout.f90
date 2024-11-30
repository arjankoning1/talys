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
  character(len=18) :: reaction   ! reaction
  character(len=132) :: topline    ! topline
  character(len=21) :: dirfile 
  character(len=13) :: Estr
  character(len=15) :: col(numlev2)     ! header
  character(len=15) :: un(numlev2)     ! header
  character(len=80) :: quantity   ! quantity
  integer   :: i                                 ! counter
  integer   :: iang                              ! running variable for angle
  integer   :: ilev                              ! counter for discrete levels
  integer   :: nen                               ! energy counter
  integer   :: Ncol
  integer   :: Nk
  integer   :: Nix                               ! neutron number index for residual nucleus
  integer   :: ilevel(numlev2)                   
  integer   :: plev(numlev2)                     ! parity of level
  integer   :: Zix                               ! charge number index for residual nucleus
  real(sgl) :: dad(numlev2, 0:numang)            ! direct angular distribution
  real(sgl) :: elev(numlev2)                     ! energy of level
  real(sgl) :: eoutlev(numlev2)                  ! emission energy
  real(sgl) :: jlev(numlev2)                     ! spin of level
  real(sgl) :: dlev(numlev2)                     ! 
  real(sgl) :: xslev(numlev2)                    ! 
!
! *********************** Inelastic cross sections *********************
!
  write(*,'(/" ++++++++++ DIRECT CROSS SECTIONS AND GIANT RESONANCES ++++++++++"/)')
  Zix = Zindex(0, 0, k0)
  Nix = Nindex(0, 0, k0)
  ilev = 0
  do i = 1, numlev2
    if (xsdirdisc(k0, i) /= 0.) then
      ilev = ilev + 1
      ilevel(ilev) = i
      elev(ilev) = edis(Zix, Nix, i)
      eoutlev(ilev) = eoutdis(k0, i)
      jlev(ilev) = jdis(Zix, Nix, i)
      plev(ilev) = parlev(Zix, Nix, i)
      dlev(ilev) = deform(Zix, Nix, i)
      xslev(ilev) = xsdirdisc(k0, i)
      if (flagang) then
        do iang = 0, nangle
          dad(ilev, iang) = directad(k0, i, iang)
        enddo
      endif
    endif
  enddo
  Nk = ilev
  col = ''
  un = ''
  Estr=''
  write(Estr,'(es13.6)') Einc
  dirfile = 'directE0000.000.out'
  write(dirfile(8:15), '(f8.3)') Einc
  write(dirfile(8:11), '(i4.4)') int(Einc)
  quantity='direct inelastic cross sections'
  reaction='('//parsym(k0)//','//parsym(k0)//"'"//')'
  topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
  open (unit = 1, file = dirfile, status = 'replace')
  call write_header(topline,source,user,date,oformat)
  call write_target
  call write_reaction(reaction,0.D0,0.D0,0,0)
  call write_real(2,'E-incident [MeV]',Einc)
  call write_real(2,'total discrete direct inelastic cross section [mb]',xsdirdisctot(k0))
  call write_real(2,'collective continuum inelastic cross section [mb]',xscollconttot(k0))
  col = ''
  un = ''
  col(1) = 'level'
  col(2) = 'energy'
  un(2) = 'MeV'
  col(3) = 'E-out'
  un(3) = 'MeV'
  col(4) = 'J/P'
  col(5) = 'cross section'
  un(5) = 'mb'
  col(6) = 'def. type'
  col(7) = 'def. par.'
  Ncol = 7
  call write_datablock(quantity,Ncol,Nk,col,un)
  do i = 1, Nk
     write(1, '(3x, i6, 6x, 2es15.6, f10.1, 1x, a1, 4x, es15.6, 6x, a1, 8x, es15.6)') ilevel(i), elev(i), eoutlev(i), &
 &    jlev(i), cparity(plev(i)), xslev(i), deftype(Zix, Nix), dlev(i)
  enddo
  if (flagang) then
    quantity='direct inelastic angular distributions'
    un = 'MeV'
    col = 'E             '
    col(1) = 'Angle'
    un(1) = ''
    do i = 1, ilev
      write(col(i+1)(3:7),'(f5.3)') elev(i)
      write(col(i+1)(9:12),'(f4.1)') jlev(i)
      write(col(i+1)(13:13),'(a1)') cparity(plev(i))
    enddo
    Ncol = ilev + 1
    Nk = nangle + 1
    call write_datablock(quantity,Ncol,Nk,col,un)
    do iang = 0, nangle
      write(1, '(5x, f5.1, 5x, 200es15.6)') angle(iang), (dad(i, iang), i = 1, ilev)
    enddo
  endif
!
! *********************** Giant resonances *****************************
!
  if ( .not. flaggiant) return
  quantity='giant resonance cross sections'
  un=''
  col(1)='Type'
  col(2)='Cross section'
  un(2)='mb'
  col(3)='Exc. energy'
  un(2)='MeV'
  col(3)='Emis. energy'
  un(2)='MeV'
  col(3)='Width'
  un(2)='MeV'
  col(3)='Deform. par.'
  Ncol = 6
  Nk = 4
  call write_datablock(quantity,Ncol,Nk,col,un)
  write(1, '("     GMR       ", 5es15.6)') xsgrcoll(k0, 0, 1), Egrcoll(0, 1), eoutgr(k0, 0, 1), Ggrcoll(0, 1), betagr(0, 1)
  write(1, '("     GQR       ", 5es15.6)') xsgrcoll(k0, 2, 1), Egrcoll(2, 1), eoutgr(k0, 2, 1), Ggrcoll(2, 1), betagr(2, 1)
  write(1, '("     LEOR      ", 5es15.6)') xsgrcoll(k0, 3, 1), Egrcoll(3, 1), eoutgr(k0, 3, 1), Ggrcoll(3, 1), betagr(3, 1)
  write(1, '("     HEOR      ", 5es15.6)') xsgrcoll(k0, 3, 2), Egrcoll(3, 2), eoutgr(k0, 3, 2), Ggrcoll(3, 2), betagr(3, 2)
! write(1, '(/" Total", f12.5/)') xsgrtot(k0)-xscollconttot(k0)
  if (flagddx) then
    quantity='average angular distributions'
    un=''
    col(1)='Angle'
    un(1)='degrees'
    col(2)='GMR'
    col(3)='GQR'
    col(4)='LEOR'
    col(5)='HEOR'
    Ncol = 5
    Nk = nanglecont + 1
    call write_datablock(quantity,Ncol,Nk,col,un)
    do iang = 0, nanglecont
      write(1, '(5es15.6)') anglecont(iang), grcollad(k0, 0, 1, iang), grcollad(k0, 2, 1, iang), &
 &      grcollad(k0, 3, 1, iang), grcollad(k0, 3, 2, iang)
    enddo
  endif
  if (flagspec) then
    quantity='Giant resonance spectra'
    un='mb'
    col(1)='E-out'
    un(1)='MeV'
    col(2)='Total'
    col(3)='GMR'
    col(4)='GQR'
    col(5)='LEOR'
    col(6)='HEOR'
    col(7)='Collective'
    Ncol = 7
    Nk = eend(k0) - ebegin(k0) + 1
    call write_datablock(quantity,Ncol,Nk,col,un)
    do nen = ebegin(k0), eend(k0)
      write(1, '(7es15.6)') egrid(nen), xsgr(k0, nen), xsgrstate(k0, 0, 1, nen), xsgrstate(k0, 2, 1, nen), &
 &      xsgrstate(k0, 3, 1, nen), xsgrstate(k0, 3, 2, nen), xscollcont(k0, nen)
    enddo
  endif
  close (unit = 1)
  call write_outfile(dirfile,flagoutall)
  return
end subroutine directout
! Copyright A.J. Koning 2021
