subroutine discreteout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of cross sections for discrete states
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
! Variables for energies
!   Ninclow        ! number of incident energies below Elow
! Variables for output
!   filechannels     ! flag for exclusive channel cross sections on separate file
!   filediscrete     ! flag for discrete level cross sections on separate file
! Variables for basic reaction
!   flagchannels     ! flag for exclusive channels calculation
! Variables for input energies
!   eninc            ! incident energy in MeV
!   nin              ! counter for incident energy
!   Ninc           ! number of incident energies
! Variables for main input
!   Atarget          ! mass number of target nucleus
!   k0               ! index of incident particle
!   Ltarget          ! excited level of target
! Variables for total cross sections
!   xsexclcont       ! exclusive single channel cross section for continuum
!   xsexclusive      ! exclusive single channel cross section
! Variables for energy grid
!   Einc             ! incident energy in MeV
! Variables for energies
!   eoutdis          ! outgoing energy of discrete state reaction
!   Ethresh          ! threshold incident energy for residual nucleus
!   Qres             ! Q - value for residual nucleus
! Variables for multiple emission
!   xsngn            ! total (projectile, gamma - ejectile) cross section
! Variables for binary reactions
!   xscompdisc       ! compound cross section for discrete state
!   xscompdisctot    ! compound cross section summed over discrete states
!   xscompound       ! total compound cross section
!   xsconttot        ! total cross section for continuum
!   xsdircont        ! direct cross section for continuum
!   xsdirect         ! total direct cross section
!   xsdisc           ! total cross section for discrete state
!   xsdisctot        ! total cross section summed over discrete states
! Variables for incident channel
!   dorigin          ! origin of direct cross section (Direct or Preeq)
!   xsbinary         ! cross section from initial compound to residual nucleus
!   xscompcont       ! compound cross section for continuum
!   xsdirdisc        ! direct cross section for discrete state direct cross section
!   xsdirdisctot     ! direct cross section summed over discrete states
! Variables for nuclides
!   Nindex           ! neutron number index for residual nucleus
!   parskip          ! logical to skip outgoing particle
!   Zindex           ! charge number index for residual nucleus
! Constants
!   cparity          ! parity (character)
!   parsym           ! symbol of particle
! Variables for levels
!   edis             ! energy of level
!   jdis             ! spin of level
!   parlev           ! parity of level
! Variables for level density
!   Nlast            ! last discrete level
! Variables for thermal cross sections
!   fxscompdisc      ! compound cross section for discrete state
!   fxsdirdisc       ! direct cross section for discrete state
!   fxsdisc          ! total cross section for discrete state
!   fxsdisctot       ! total cross section summed over discrete states
!   fxsexclcont      ! exclusive single channel cross section for contin
!   fxsexclusive     ! exclusive single channel cross section
!   fxsngn           ! total (projectile, gamma - ejectile) cross section
!
! *** Declaration of local data
!
  implicit none
  character(len=2) :: levelstring             
  character(len=6) :: discfile                ! file with elastic scattering angular distribution
  character(len=6) :: contfile                ! file for continuum
  character(len=6) :: totfile                 ! file with total cross sections
  character(len=9) :: reactstring(0:6)        ! reaction string
  character(len=18) :: reaction   ! reaction
  character(len=132) :: topline    ! topline
  character(len=15) :: col(5)     ! header
  character(len=15) :: un(5)     ! header
  character(len=80) :: quantity   ! quantity
  integer           :: MF
  integer           :: MT
  integer           :: MT0(0:6)
  integer           :: indent
  integer           :: id2
  integer           :: id4
  integer          :: Ncomp                   ! neutron number index for compound nucleus
  integer          :: nen                     ! energy counter
  integer          :: Ncol                    ! number of columns
  integer          :: nex                     ! excitation energy bin of compound nucleus
  integer          :: Nix                     ! neutron number index for residual nucleus
  integer          :: NL                      ! last discrete level
  integer          :: type                    ! particle type
  integer          :: Zcomp                   ! proton number index for compound nucleus
  integer          :: Zix                     ! charge number index for residual nucleus
!
! ************************* Make reaction string ***********************
!
  indent = 0
  id2 = indent + 2
  id4 = indent + 4
  MT0(0) = 101
  MT0(1) = 50
  MT0(2) = 600
  MT0(3) = 650
  MT0(4) = 700
  MT0(5) = 750
  MT0(6) = 800
  MF = 3
  un = 'mb'
  col(1)='E'
  un(1)='MeV'
  col(2)='xs'
  Ncol=2
  col(3)='Direct'
  col(4)='Compound'
  quantity='cross section'
  do type = 0, 6
    if (type == k0) then
      reactstring(type) = 'Inelastic'
    else
      reactstring(type) = '  ( , )  '
      write(reactstring(type)(4:4), '(a1)') parsym(k0)
      write(reactstring(type)(6:6), '(a1)') parsym(type)
    endif
  enddo
!
! ****************** Cross sections for discrete states ****************
!
  Zcomp = 0
  Ncomp = 0
  write(*, '(/" 5. Binary reactions to discrete levels", " and continuum")')
  do type = 0, 6
    if (parskip(type)) cycle
    if (xsbinary(type) == 0.) cycle
    Zix = Zindex(Zcomp, Ncomp, type)
    Nix = Nindex(Zcomp, Ncomp, type)
    write(*, '(/1x, a9, " cross sections:"/)') reactstring(type)
    write(*, '(" Inclusive:"/)')
    write(*, '(" Level Energy    E-out     J/P       Direct    ", "Compound      Total     Origin"/)')
    do nex = 0, Nlast(Zix, Nix, 0)
      if (type == k0 .and. nex == Ltarget) cycle
      write(*, '(1x, i2, 2f10.5, f7.1, a1, 3f12.5, 4x, a6)') nex, edis(Zix, Nix, nex), eoutdis(type, nex), jdis(Zix, Nix, nex), &
 &      cparity(parlev(Zix, Nix, nex)), xsdirdisc(type, nex), xscompdisc(type, nex), xsdisc(type, nex), dorigin(type, nex)
    enddo
    write(*, '(31x, 3("   ---------"))')
    write(*, '(" Discrete  ", a9, ":", 10x, 3f12.5)') reactstring(type), xsdirdisctot(type), xscompdisctot(type), xsdisctot(type)
    write(*, '(" Continuum ", a9, ":", 10x, 3f12.5)') reactstring(type), xsdircont(type), xscompcont(type), xsconttot(type)
    write(*, '(31x, 3("   ---------"))')
    write(*, '(" Total     ", a9, ":", 10x, 3f12.5/)') reactstring(type), xsdirect(type), xscompound(type), xsbinary(type)
    if (type /= 0) write(*, '(" (", a1, ",g", a1, ") cross section:", f12.5)') parsym(k0), parsym(type), xsngn(type)
    if (flagchannels) then
      write(*, '(/" Exclusive"/)')
      write(*, '(" Discrete  ", a9, ":", f12.5)') reactstring(type), xsdisctot(type)
      write(*, '(" Continuum ", a9, ":", f12.5)') reactstring(type), xsexclcont(type)
      write(*, '(" Total     ", a9, ":", f12.5)') reactstring(type), xsexclusive(type)
    endif
  enddo
!
! Write results to separate file
!
  do type = 0, 6
    if (parskip(type)) cycle
    Zix = Zindex(Zcomp, Ncomp, type)
    Nix = Nindex(Zcomp, Ncomp, type)
    do nex = 0, Nlast(Zix, Nix, 0)
      if (type == k0 .and. nex == Ltarget) cycle
      if (filediscrete(nex)) then
        discfile = '  .L00'
        write(discfile(1:2), '(2a1)') parsym(k0), parsym(type)
        write(discfile(5:6), '(i2.2)') nex
        if (nin == Ninclow + 1) then
          open (unit = 1, file = discfile, status = 'replace')
          levelstring='  '
          write(levelstring,'(i2)') nex
          reaction='('//parsym(k0)//','//parsym(type)//'_'//trim(adjustl(levelstring))//')'
          topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
          Ncol=4
          if ( type > 0) then
            MT =  MT0(type) + nex
          else
            MT = 0
          endif
          call write_header(indent,topline,source,user,date,oformat)
          call write_target(indent)
          call write_reaction(indent,reaction,Qres(Zix, Nix, nex),Ethresh(Zix, Nix, nex),MF,MT)
          call write_level(id2,-1,nex,edis(Zix, Nix, nex),jdis(Zix, Nix, nex),parlev(Zix, Nix, nex),0.)
          call write_quantity(id2,quantity)
          call write_datablock(id2,Ncol,Ninc,col,un)
          do nen = 1, Ninclow
            write(1, '(4es15.6)') eninc(nen), fxsdisc(nen, type, nex), fxsdirdisc(nen, type, nex), fxscompdisc(nen, type, nex)
          enddo
        else
          open (unit = 1, file = discfile, status = 'old', position = 'append')
        endif
        write(1, '(4es15.6)') Einc, xsdisc(type, nex), xsdirdisc(type, nex), xscompdisc(type, nex)
        close (unit = 1)
      endif
    enddo
  enddo
!
! Write continuum cross sections to separate file
!
  if (filechannels) then
    do type = 0, 6
      if (parskip(type)) cycle
      Zix = Zindex(Zcomp, Ncomp, type)
      Nix = Nindex(Zcomp, Ncomp, type)
!     if (type == k0) then
!       NL = 0
!     else
        NL = Nlast(Zix, Nix, 0)
!     endif
      contfile = '  .con'
      write(contfile(1:2), '(2a1)') parsym(k0), parsym(type)
      reaction='('//parsym(k0)//','//parsym(type)//'_con)'
      if (nin == Ninclow + 1) then
        open (unit = 1, file = contfile, status = 'replace')
        topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
        Ncol=4
        if (type == 0) MT = 0
        if (type == 1) MT= 91
        if (type > 1) then
          MT= MT0(type) + 49
        endif
        col(3)='Continuum'
        col(4)='('//parsym(k0)//',g'//parsym(type)//')'
        call write_header(indent,topline,source,user,date,oformat)
        call write_target(indent)
        call write_reaction(indent,reaction,Qres(Zix, Nix, NL),Ethresh(Zix, Nix, NL),MF,MT)
        call write_level(id2,-1,NL,edis(Zix, Nix, NL),jdis(Zix, Nix, NL),parlev(Zix, Nix, NL),0.)
        call write_quantity(id2,quantity)
        call write_datablock(id2,Ncol,Ninc,col,un)
        do nen = 1, Ninclow
          write(1, '(4es15.6)') eninc(nen), fxsexclcont(nen, type) + fxsngn(nen, type), fxsexclcont(nen, type), fxsngn(nen, type)
        enddo
      else
        open (unit = 1, file = contfile, status = 'old', position = 'append')
      endif
      write(1, '(4es15.6)') Einc, xsexclcont(type) + xsngn(type), xsexclcont(type), xsngn(type)
      close (unit = 1)
    enddo
  endif
!
! Write cumulated binary cross sections to separate file
!
  if (filechannels) then
    do type = 0, 6
      if (parskip(type)) cycle
      Zix = Zindex(Zcomp, Ncomp, type)
      Nix = Nindex(Zcomp, Ncomp, type)
      totfile = '  .tot'
      write(totfile(1:2), '(2a1)') parsym(k0), parsym(type)
      reaction='('//parsym(k0)//','//parsym(type)//')'
      if (nin == Ninclow + 1) then
        open (unit = 1, file = totfile, status = 'replace')
        topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
        Ncol=5
        col(3)='Discrete'
        col(4)='Continuum'
        col(5)='('//parsym(k0)//',g'//parsym(type)//')'
        call write_header(indent,topline,source,user,date,oformat)
        call write_target(indent)
        call write_reaction(indent,reaction,Qres(Zix, Nix, 0),Ethresh(Zix, Nix, 0),0,0)
        call write_quantity(id4,quantity)
        call write_datablock(id2,Ncol,Ninc,col,un)
        do nen = 1, Ninclow
          write(1, '(5es15.6)') eninc(nen), fxsexclusive(nen, type), fxsdisctot(nen, type), fxsexclcont(nen, type), &
 &          fxsngn(nen, type)
        enddo
      else
        open (unit = 1, file = totfile, status = 'old', position = 'append')
      endif
      write(1, '(5es15.6)') Einc, xsexclusive(type), xsdisctot(type), xsexclcont(type), xsngn(type)
      close (unit = 1)
    enddo
  endif
  return
end subroutine discreteout
! Copyright A.J. Koning 2021
