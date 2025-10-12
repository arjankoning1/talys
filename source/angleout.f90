subroutine angleout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of discrete angular distributions
!
! Author    : Arjan Koning
!
! 2025-10-08: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Variables for output
!   fileangle       ! designator for angular distributions on separate file
!   flagblockangle  ! flag to block angle files
! Variables for existence libraries
!   angexist        ! flag for existence of angular distributions
!   legexist        ! flag for existence of Legendre coefficients
! Variables for direct reactions
!   flaglegendre    ! flag for output of Legendre coefficients
! Variables for numerics
!   nangle          ! number of angles
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
!   Ltarget         ! excited level of target
!   Starget         ! symbol of target nucleus
! Variables for discrete levels
!   nlev            ! number of levels for nucleus
! Variables for energy grid
!   angle           ! angle in degrees
!   Einc            ! incident energy in MeV
! Variables for binary reactions
!   xsdisc          ! total cross section for discrete state
!   xsdisctot       ! total cross section summed over discrete states
! Variables for incident channel
!   cleg            ! compound nucleus Legendre coefficient
!   directad        ! direct angular distribution
!   dleg            ! direct reaction Legendre coefficient
!   ruth            ! elastic / Rutherford ratio
!   elasni          ! nuclear+interference term
! Variables for angular distributions
!   cleg0           ! Legendre coefficient normalized to the first one
!   compad          ! compound angular distribution
!   discad          ! discrete state angular distribution
!   tleg            ! total Legendre coefficient
!   tlegnor         ! total Legendre coefficient normalized to 1
! Variables to normalize compound nucleus cross section
!   J2end           ! end of J summation
! Variables for nuclides
!   Nindex          ! neutron number index for residual nucleus
!   parskip         ! logical to skip outgoing particle
!   Zindex          ! charge number index for residual nucleus
! Constants
!   parsym          ! symbol of particle
!
! *** Declaration of local data
!
  implicit none
  character(len=2)  :: levelstring
  character(len=21) :: discfile    ! file with elastic scattering angular distribution
  character(len=21) :: legfile     ! file with Legendre coefficients
  character(len=12) :: Estr
  character(len=18) :: reaction   ! reaction
  character(len=132) :: topline    ! topline
  character(len=15) :: col(6)     ! header
  character(len=15) :: un(6)     ! units
  character(len=80) :: quantity   ! quantity
  integer           :: MF
  integer           :: MT
  integer           :: indent
  integer           :: id2
  integer           :: id4
  integer           :: MT0(0:6)
  integer           :: i           ! counter
  integer           :: iang        ! running variable for angle
  integer           :: LL          ! counter for l value
  integer           :: Ncol        ! 
  integer           :: Nix         ! neutron number index for residual nucleus
  integer           :: type        ! particle type
  integer           :: Zix         ! charge number index for residual nucleus
!
! **************** Elastic scattering angular distribution *************
!
!
  MT0(0) = 101
  MT0(1) = 50
  MT0(2) = 600
  MT0(3) = 650
  MT0(4) = 700
  MT0(5) = 750
  MT0(6) = 800
  MF = 4
  MT = 2
  indent = 0
  id2 = indent + 2
  id4 = indent + 4
!
! 1. Legendre coefficients
!
  Estr=''
  write(Estr,'(es12.6)') Einc
  write(*, '(/" 8. Discrete state angular distributions")')
  if (flaglegendre) then
    write(*, '(/" 8a1. Legendre coefficients for elastic scattering"/)')
!
! Write results to separate file
!
    if (flagblockangle) then
      legfile = '  leg.L00'
      write(legfile(1:2), '(2a1)') parsym(k0), parsym(k0)
      write(legfile(8:9), '(i2.2)') Ltarget
      if (.not. legexist(k0, k0, Ltarget)) then
        legexist(k0, k0, Ltarget) = .true.
        open (unit=1, file=legfile, status='unknown')
      else
        open (unit=1, file=legfile, status='unknown', position='append')
      endif
    else
      legfile = '          leg.L00'
      write(legfile(1:2), '(2a1)') parsym(k0), parsym(k0)
      write(legfile(3:10), '(f8.3)') Einc
      write(legfile(3:6), '(i4.4)') int(Einc)
      write(legfile(16:17), '(i2.2)') Ltarget
      open (unit=1, file=legfile, status='unknown')
    endif
    reaction='('//parsym(k0)//',el)'
    quantity='Legendre coefficients'
    topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
    un = ''
    col(1)='L'
    col(2)='Total'
    col(3)='Direct'
    col(4)='Compound'
    col(5)='Normalized'
    col(6)='ENDF-6'
    Ncol=6
    call write_header(indent,topline,source,user,date,oformat)
    call write_target(indent)
    call write_reaction(indent,reaction,0.D0,0.D0,MF,MT)
    call write_real(id2,'E-incident [MeV]',Einc)
    call write_quantity(id2,quantity)
    call write_datablock(id2,Ncol,J2end+1,col,un)
    do LL = 0, J2end
      write(1, '(i6, 9x, 5es15.6)') LL, tleg(k0, Ltarget, LL), dleg(k0, Ltarget, LL), cleg(k0, Ltarget, LL), &
 &      tlegnor(k0, Ltarget, LL), cleg0(k0, Ltarget, LL)
    enddo
    close (unit = 1)
    call write_outfile(legfile,flagoutall)
  endif
!
! 2. Angular distributions
!
  write(*, '(/" 8a2. Elastic scattering angular distribution"/)')
  if (k0 == 1) then
!
! Write results to separate file
!
    if (flagblockangle) then
      discfile = 'nnang.L00'
      write(discfile(1:2), '(2a1)') parsym(k0), parsym(k0)
      write(discfile(8:9), '(i2.2)') Ltarget
      if (.not. angexist(k0, k0, Ltarget)) then
        angexist(k0, k0, Ltarget) = .true.
        open (unit=1, file=discfile, status='unknown')
      else
        open (unit=1, file=discfile, status='unknown', position='append')
      endif
    else
      discfile = 'nn        ang.L00'
      write(discfile(3:10), '(f8.3)') Einc
      write(discfile(3:6), '(i4.4)') int(Einc)
      write(discfile(16:17), '(i2.2)') Ltarget
      open (unit=1, file=discfile, status='unknown')
    endif
    quantity='angular distribution'
    reaction='('//parsym(k0)//',el)'
    topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
    col(1)='Angle'
    un(1)='deg'
    col(2)='xs'
    un(2)='mb/sr'
    col(3)='Direct'
    un(3)='mb/sr'
    col(4)='Compound'
    un(4)='mb/sr'
    Ncol=4
    call write_header(indent,topline,source,user,date,oformat)
    call write_target(indent)
    call write_reaction(indent,reaction,0.D0,0.D0,MF,MT)
    call write_real(id2,'E-incident [MeV]',Einc)
    call write_quantity(id2,quantity)
    call write_datablock(id2,Ncol,nangle+1,col,un)
    do iang = 0, nangle
      write(1, '(4es15.6)') angle(iang), discad(k0, Ltarget, iang), directad(k0, Ltarget, iang), compad(k0, Ltarget, iang)
    enddo
    close (unit = 1)
    call write_outfile(discfile,flagoutall)
  else
!
! Write results to separate file
!
    if (flagblockangle) then
      discfile = '  ang.L00'
      write(discfile(1:2), '(2a1)') parsym(k0), parsym(k0)
      write(discfile(8:9), '(i2.2)') Ltarget
      if (.not. angexist(k0, k0, Ltarget)) then
        angexist(k0, k0, Ltarget) = .true.
        open (unit=1, file=discfile, status='unknown')
      else
        open (unit=1, file=discfile, status='unknown', position='append')
      endif
    else
      discfile = '          ang.L00'
      write(discfile(1:2), '(2a1)') parsym(k0), parsym(k0)
      write(discfile(3:10), '(f8.3)') Einc
      write(discfile(3:6), '(i4.4)') int(Einc)
      write(discfile(16:17), '(i2.2)') Ltarget
      open (unit=1, file=discfile, status='unknown')
    endif
    quantity='angular distribution'
    reaction='('//parsym(k0)//',el)'
    topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
    col(1)='Angle'
    un(1)='deg'
    col(2)='xs'
    un(2)='mb/sr'
    col(3)='Direct'
    un(3)='mb/sr'
    col(4)='Compound'
    un(4)='mb/sr'
    col(5)='xs/Rutherford'
    un(5)=''
    col(6)='Nuc+interfer.'
    un(6)=''
    Ncol=6
    call write_header(indent,topline,source,user,date,oformat)
    call write_target(indent)
    call write_reaction(indent,reaction,0.D0,0.D0,MF,MT)
    call write_real(id2,'E-incident [MeV]',Einc)
    call write_quantity(id2,quantity)
    call write_datablock(id2,Ncol,nangle+1,col,un)
    do iang = 0, nangle
      write(1, '(6es15.6)') angle(iang), max(discad(k0, Ltarget, iang), directad(k0, Ltarget, iang)), &
 &      directad(k0, Ltarget, iang), compad(k0, Ltarget, iang), ruth(iang), elasni(iang)
    enddo
    close (unit = 1)
    call write_outfile(discfile,flagoutall)
  endif
!
! ************** Inelastic scattering angular distributions ************
!
  Zix = Zindex(0, 0, k0)
  Nix = Nindex(0, 0, k0)
!
! 1. Legendre coefficients
!
  if (flaglegendre) then
    write(*, '(/" 8b1. Legendre coefficients for inelastic scattering")')
    do i = 0, nlev(Zix, Nix)
      if (i == Ltarget) cycle
      if (xsdisc(k0, i) == 0.) cycle
!
! Write results to separate file
!
      if (fileangle(i)) then
        if (flagblockangle) then
          legfile = '  leg.L00'
          write(legfile(1:2), '(2a1)') parsym(k0), parsym(k0)
          write(legfile(8:9), '(i2.2)') i
          if (.not. legexist(k0, k0, i)) then
            legexist(k0, k0, i) = .true.
            open (unit=1, file=legfile, status='unknown')
          else
            open (unit=1, file=legfile, status='unknown', position='append')
          endif
        else
          legfile = '          leg.L00'
          write(legfile(1:2), '(2a1)') parsym(k0), parsym(k0)
          write(legfile(3:10), '(f8.3)') Einc
          write(legfile(3:6), '(i4.4)') int(Einc)
          write(legfile(16:17), '(i2.2)') i
          open (unit=1, file=legfile, status='unknown')
        endif
        quantity='Legendre coefficients'
        levelstring='  '
        write(levelstring,'(i2)') i
        reaction='('//parsym(k0)//','//parsym(k0)//'_'//trim(adjustl(levelstring))//')'
        topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
        MT =  MT0(k0) + i
        un = ''
        col(1)='L'
        col(2)='Total'
        col(3)='Direct'
        col(4)='Compound'
        col(5)='Normalized'
        col(6)='ENDF-6'
        Ncol=6
        call write_header(indent,topline,source,user,date,oformat)
        call write_target(indent)
        call write_reaction(indent,reaction,0.D0,0.D0,MF,MT)
        call write_real(id2,'E-incident [MeV]',Einc)
        call write_level(id2,-1,i,edis(Zix, Nix, i),jdis(Zix, Nix, i),parlev(Zix, Nix, i),0.)
        call write_quantity(id2,quantity)
        call write_datablock(id2,Ncol,J2end+1,col,un)
        do LL = 0, J2end
          write(1, '(i6, 9x, 5es15.6)') LL, tleg(k0, i, LL), dleg(k0, i, LL), cleg(k0, i, LL), tlegnor(k0, i, LL), cleg0(k0, i, LL)
        enddo
        close (unit = 1)
        call write_outfile(legfile,flagoutall)
      endif
    enddo
  endif
!
! 2. Angular distributions
!
  write(*, '(/" 8b2. Inelastic angular distributions")')
  do i = 0, nlev(Zix, Nix)
    if (i == Ltarget) cycle
    if (xsdisc(k0, i) == 0.) cycle
!
! Write results to separate file
!
    if (fileangle(i)) then
      if (flagblockangle) then
        discfile = '  ang.L00'
        write(discfile(1:2), '(2a1)') parsym(k0), parsym(k0)
        write(discfile(8:9), '(i2.2)') i
        if (.not. angexist(k0, k0, i)) then
          angexist(k0, k0, i) = .true.
          open (unit=1, file=discfile, status='unknown')
        else
          open (unit=1, file=discfile, status='unknown', position='append')
        endif
      else
        discfile = '          ang.L00'
        write(discfile(1:2), '(2a1)') parsym(k0), parsym(k0)
        write(discfile(3:10), '(f8.3)') Einc
        write(discfile(3:6), '(i4.4)') int(Einc)
        write(discfile(16:17), '(i2.2)') i
        open (unit=1, file=discfile, status='unknown')
      endif
      quantity='angular distribution'
      levelstring='  '
      write(levelstring,'(i2)') i
      reaction='('//parsym(k0)//','//parsym(k0)//'_'//trim(adjustl(levelstring))//')'
      topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
      MT =  MT0(k0) + i
      Ncol=4
      col(1)='Angle'
      un(1)='deg'
      col(2)='xs'
      un(2)='mb/sr'
      col(3)='Direct'
      un(3)='mb/sr'
      col(4)='Compound'
      un(4)='mb/sr'
      call write_header(indent,topline,source,user,date,oformat)
      call write_target(indent)
      call write_reaction(indent,reaction,0.D0,0.D0,MF,MT)
      call write_real(id2,'E-incident [MeV]',Einc)
      call write_level(id2,-1,i,edis(Zix, Nix, i),jdis(Zix, Nix, i),parlev(Zix, Nix, i),0.)
      call write_quantity(id2,quantity)
      call write_datablock(id2,Ncol,nangle+1,col,un)
      do iang = 0, nangle
        write(1, '(4es15.6)') angle(iang), discad(k0, i, iang), directad(k0, i, iang), compad(k0, i, iang)
      enddo
      close (unit = 1)
      call write_outfile(discfile,flagoutall)
    endif
  enddo
!
! ********** Non-inelastic scattering angular distributions ************
!
  write(*, '(/" 8c. Angular distributions for other reactions")')
  do type = 0, 6
    if (parskip(type)) cycle
    if (type == k0) cycle
    if (xsdisctot(type) == 0.) cycle
    Zix = Zindex(0, 0, type)
    Nix = Nindex(0, 0, type)
!
! 1. Legendre coefficients
!
    if (flaglegendre) then
      write(*, '(/" 8c1. Legendre coefficients for (", a1, ",", a1, ")")') parsym(k0), parsym(type)
      do i = 0, nlev(Zix, Nix)
        if (xsdisc(type, i) == 0.) cycle
!
! Write results to separate file
!
        if (fileangle(i)) then
          if (flagblockangle) then
            legfile = '  leg.L00'
            write(legfile(1:2), '(2a1)') parsym(k0), parsym(type)
            write(legfile(8:9), '(i2.2)') i
            if (.not. legexist(k0, type, i)) then
              legexist(k0, type, i) = .true.
              open (unit=1, file=legfile, status='unknown')
            else
              open (unit=1, file=legfile, status='unknown', position='append')
            endif
          else
            legfile = '          leg.L00'
            write(legfile(1:2), '(2a1)') parsym(k0), parsym(type)
            write(legfile(3:10), '(f8.3)') Einc
            write(legfile(3:6), '(i4.4)') int(Einc)
            write(legfile(16:17), '(i2.2)') i
            open (unit=1, file=legfile, status='unknown')
          endif
          quantity='Legendre coefficients'
          levelstring='  '
          write(levelstring,'(i2)') i
          reaction='('//parsym(k0)//','//parsym(type)//'_'//trim(adjustl(levelstring))//')'
          topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
          MT =  MT0(type) + i
          un = ''
          col(1)='L'
          col(2)='Total'
          col(3)='Direct'
          col(4)='Compound'
          col(5)='Normalized'
          col(6)='ENDF-6'
          Ncol=6
          call write_header(indent,topline,source,user,date,oformat)
          call write_target(indent)
          call write_reaction(indent,reaction,0.D0,0.D0,MF,MT)
          call write_real(id2,'E-incident [MeV]',Einc)
          call write_level(id2,-1,i,edis(Zix, Nix, i),jdis(Zix, Nix, i),parlev(Zix, Nix, i),0.)
          call write_quantity(id2,quantity)
          call write_datablock(id2,Ncol,J2end+1,col,un)
          do LL = 0, J2end
            write(1, '(i6, 9x, 5es15.6)') LL, tleg(type, i, LL), dleg(type, i, LL), cleg(type, i, LL), tlegnor(type, i, LL), &
 &            cleg0(type, i, LL)
          enddo
        close (unit = 1)
        call write_outfile(legfile,flagoutall)
      endif
      enddo
    endif
!
! 2. Angular distributions
!
    write(*, '(/" 8c2. (", a1, ",", a1, ") angular distributions")') parsym(k0), parsym(type)
    do i = 0, nlev(Zix, Nix)
      if (xsdisc(type, i) == 0.) cycle
!
! Write results to separate file
!
      if (fileangle(i)) then
        if (flagblockangle) then
          discfile = '  ang.L00'
          write(discfile(1:2), '(2a1)') parsym(k0), parsym(type)
          write(discfile(8:9), '(i2.2)') i
          if (.not. angexist(k0, type, i)) then
            angexist(k0, type, i) = .true.
            open (unit=1, file=discfile, status='unknown')
          else
            open (unit=1, file=discfile, status='unknown', position='append')
          endif
        else
          discfile = '          ang.L00'
          write(discfile(1:2), '(2a1)') parsym(k0), parsym(type)
          write(discfile(3:10), '(f8.3)') Einc
          write(discfile(3:6), '(i4.4)') int(Einc)
          write(discfile(16:17), '(i2.2)') i
          open (unit=1, file=discfile, status='unknown')
        endif
        quantity='angular distribution'
        levelstring='  '
        write(levelstring,'(i2)') i
        reaction='('//parsym(k0)//','//parsym(type)//'_'//trim(adjustl(levelstring))//')'
        topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
        MT =  MT0(type) + i
        Ncol=4
        col(1)='Angle'
        un(1)='deg'
        col(2)='xs'
        un(2)='mb/sr'
        col(3)='Direct'
        un(3)='mb/sr'
        col(4)='Compound'
        un(4)='mb/sr'
        call write_header(indent,topline,source,user,date,oformat)
        call write_target(indent)
        call write_reaction(indent,reaction,0.D0,0.D0,MF,MT)
        call write_real(id2,'E-incident [MeV]',Einc)
        call write_level(id2,-1,i,edis(Zix, Nix, i),jdis(Zix, Nix, i),parlev(Zix, Nix, i),0.)
        call write_quantity(id2,quantity)
        call write_datablock(id2,Ncol,nangle+1,col,un)
        do iang = 0, nangle
          write(1, '(4es15.6)') angle(iang), discad(type, i, iang), directad(type, i, iang), compad(type, i, iang)
        enddo
        close (unit = 1)
        call write_outfile(discfile,flagoutall)
      endif
    enddo
  enddo
  return
end subroutine angleout
! Copyright A.J. Koning 2021
