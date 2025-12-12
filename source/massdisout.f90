subroutine massdisout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of fission fragment yields
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
! Variables for input energies
!   eninc           ! incident energy in MeV
!   nin0            ! counter for incident energy
!   Ninc          ! number of incident energies
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
!   Ninit           ! neutron number of initial compound nucleus
!   Ztarget         ! charge number of target nucleus
! Variables for energies
!   Einc0           ! incident energy in MeV
! Constants
!   iso             ! counter for isotope
!   natstring       ! string extension for file names
!   nuc             ! symbol of nucleus
!   parsym          ! symbol of particle
! Variables for existence libraries
!   fpaexist        ! flag for existence of fission product per mass unit
!   fpexist         ! flag for existence of fission product
! Variables for mass distribution
!   fpeps           ! ratio for limit for fission product cross section
!   fpratio         ! fission product isomeric ratio
!   xsApost         ! post - neutron emission corrected cross section
!   xsApre          ! pre - neutron emission cross section
!   xsfpex          ! excitation energy spectrum per fission fragment
!   xstotpost       ! post - neutron emission fission product cross section
!   xstotpre        ! pre - neutron emission fission product cross section
!   xsZApost        ! post - neutron emission corrected isotopic cross section
!   xsZApre         ! pre - neutron emission isotopic cross section
!   yieldApost      ! post - neutron emission corrected fission yield
!   yieldApre       ! pre - neutron emission fission yield
!   yieldfpex       ! fission yield per isomer
!   yieldtotpost    ! post - neutron emission fission product yield
!   yieldtotpre     ! pre - neutron emission fission product yield
!   yieldZApost     ! post - neutron emission corrected isotopic yield
!   yieldZApre      ! pre - neutron emission isotopic yield
!
! *** Declaration of local data
!
  implicit none
  character(len=3) :: massstring
  character(len=6) :: finalnuclide
  character(len=12) :: isostring(0:1)    ! string to designate target isomer
  character(len=12) :: Estr
  character(len=132):: fpfile            ! file with fission product
  character(len=132):: fpyieldfile       ! file with fission yields
  character(len=18) :: reaction   ! reaction
  character(len=132) :: topline    ! topline
  character(len=15) :: col(8)     ! header
  character(len=15) :: un(8)     ! header
  character(len=80) :: quantity   ! quantity
  integer           :: MF
  integer           :: MT
  integer           :: i
  integer           :: indent
  integer           :: id2
  integer           :: id4
  integer           :: ia                ! mass number from abundance table
  integer           :: in                ! counter for neutrons
  integer           :: iz                ! charge number of residual nucleus
  integer           :: Ncol              ! number of columns
  integer           :: Nfy               ! number of fission yields
  integer           :: nen               ! energy counter
  integer           :: nex               ! excitation energy bin of compound nucleus
!
! ****************** Output of fission yields **************************
!
  indent = 0
  id2 = indent + 2
  id4 = indent + 4
  MF = 8
  MT = 454
  Estr=''
  write(Estr,'(es12.6)') Einc0
  un = ''
  col(1)='A'
  col(2)='FP_yield'
  col(3)='FF_yield'
  col(4)='FP_xs'
  un(4)='mb'
  col(5)='FF_xs'
  un(5)='mb'
  Ncol=5
  isostring(0) = 'ground state'
  isostring(1) = 'isomer      '
  write(*, '(/" ++++++++++ FISSION YIELDS ++++++++++"/)')
!
! Write results to separate files
!
  fpyieldfile = 'yieldA0000.000.fis'//natstring(iso)
  if (Einc0 < 0.001) then
    write(fpyieldfile(7:14), '(es8.2)') Einc0
  else
    write(fpyieldfile(7:14), '(f8.3)') Einc0
    write(fpyieldfile(7:10), '(i4.4)') int(Einc0)
  endif
  open (unit = 1, file = fpyieldfile, status = 'replace')
  quantity='fission yield'
  reaction='('//parsym(k0)//',f)'
  topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' per A'//' at '//Estr//' MeV'
  call write_header(indent,topline,source,user,date,oformat)
  call write_target(indent)
  call write_reaction(indent,reaction,0.D0,0.D0,MF,MT)
  call write_real(id2,'E-incident [MeV]',Einc0)
  call write_quantity(id2,quantity)
  call write_datablock(id2,Ncol,Atarget,col,un)
  write(*, '(" Fission yields as function of A"/)')
  write(*, '("   A       FP yield       FF yield ", "         FP xs          FF xs")')
  do ia = 1, Atarget
    write(*, '(i4, 4es15.6)') ia, yieldApost(ia), yieldApre(ia), xsApost(ia), xsApre(ia)
    write(1, '(i6, 9x, 4es15.6)') ia, yieldApost(ia), yieldApre(ia), xsApost(ia), xsApre(ia)
  enddo
  write(*, '(/" Tot", 4es15.6)') yieldtotpost, yieldtotpre, xstotpost, xstotpre
  write(1, '(/"#Tot           ", 4es15.6)') yieldtotpost, yieldtotpre, xstotpost, xstotpre
  close (unit = 1)
!
! Write ff/fp production
!
  write(*, '(/" Fission yields as function of Z, A"/)')
  write(*, '("    Z    A iso     FP yield       FF yield", "         FP xs          FF xs    Isom. Ratio")')
  Nfy=0
  do ia = 1, Atarget
    do iz = 1, Ztarget
      in = ia - iz
      if (in < 1 .or. in > Ninit) cycle
      if (xsZApre(iz, in) <= fpeps .and. xsZApost(iz, in) <= fpeps .and. .not. fpexist(1, iz, in, -1)) cycle
      do i = 1,2
        fpfile = 'fp000000.tot'//natstring(iso)
        if (i == 2) fpfile(1:1)='r'
        write(fpfile(3:8), '(2i3.3)') iz, ia
        massstring='   '
        write(massstring,'(i3)') ia
        finalnuclide=trim(nuc(iz))//trim(adjustl(massstring))
        if ( .not. fpexist(i, iz, in, -1)) then
          open (unit = 1, file = fpfile, status = 'replace')
          quantity='fission yield'
          if (i == 2) quantity='cross section'
          topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' for '//finalnuclide
          col(1)='E'
          un(1)='MeV'
          if (i == 1) then
            col(2)='FP_yield'
            col(3)='FF_yield'
            col(4)='FP_xs'
            col(5)='FF_xs'
            Ncol=5
          else
            col(2)='FP_xs'
            un(2)='mb'
            Ncol=2
          endif
          call write_header(indent,topline,source,user,date,oformat)
          call write_target(indent)
          if (i == 1) then
            call write_reaction(indent,reaction,0.D0,0.D0,MF,MT)
          else
            call write_reaction(indent,reaction,0.D0,0.D0,6,5)
          endif
          call write_residual(id2,iz,ia,finalnuclide)
          call write_quantity(id2,quantity)
          call write_datablock(id2,Ncol,Ninc,col,un)
          if (i == 1) then
            do nen = 1, nin0 - 1
              write(1, '(5es15.6)') eninc(nen), 0., 0., 0., 0.
            enddo
          else
            do nen = 1, nin0 - 1
              write(1, '(2es15.6)') eninc(nen), 0.
            enddo
          endif
        else
          open (unit = 1, file = fpfile, status = 'old', position = 'append')
        endif
        if (i == 1 .and. xsZApre(iz, in) >= fpeps .and. xsZApost(iz, in) >= fpeps) then
          write(*, '(2i5, i3, 4es15.6)') iz, ia, -1, yieldZApost(iz, in), yieldZApre(iz, in), xsZApost(iz, in), xsZApre(iz, in)
          Nfy=Nfy+1
        endif
        if (i == 1) then
          write(1, '(5es15.6)') Einc0, yieldZApost(iz, in), yieldZApre(iz, in), xsZApost(iz, in), xsZApre(iz, in)
        else
          write(1, '(2es15.6)') Einc0, xsZApost(iz, in)
        endif
        close (unit = 1)
        if (xsfpex(iz, in, 1) > 0.) then
          do nex = 0, 1
            write(fpfile(10:12), '("i", i2.2)') nex
            finalnuclide=trim(nuc(iz))//trim(adjustl(massstring))//isochar(nex)
            if ( .not. fpexist(i, iz, in, nex)) then
              open (unit = 1, file = fpfile, status = 'replace')
              reaction='('//parsym(k0)//',f)'
              topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' for '//finalnuclide
              col(1)='E'
              un(1)='MeV'
              if (i == 1) then
                col(2)='FP_yield'
                un(2)=''
                col(3)='FP_xs'
                un(3)='mb'
                col(4)='Ratio'
                un(4)=''
                Ncol=4
              else
                col(2)='FP_xs'
                un(2)='mb'
                Ncol=2
              endif
              call write_header(indent,topline,source,user,date,oformat)
              call write_target(indent)
              if (i == 1) then
                call write_reaction(indent,reaction,0.D0,0.D0,0,0)
              else
                call write_reaction(indent,reaction,0.D0,0.D0,6,5)
              endif
              call write_residual(id2,iz,ia,finalnuclide)
              if (nex == 0) then
                write(1, '("#   isomer: 0")') 
              else
                write(1, '("#   isomer: 1")') 
              endif
              call write_quantity(id2,quantity)
              call write_datablock(id2,Ncol,Atarget,col,un)
              if (i == 1) then
                do nen = 1, nin0 - 1
                  write(1, '(4es15.6)') eninc(nen), 0., 0., 0.
                enddo
              else
                do nen = 1, nin0 - 1
                  write(1, '(2es15.6)') eninc(nen), 0.
                enddo
              endif
              fpexist(i, iz, in, nex) = .true.
            else
              open (unit = 1, file = fpfile, status = 'old', position = 'append')
            endif
            if (xsZApre(iz, in) >= fpeps .and. xsZApost(iz, in) >= fpeps) then
              write(*, '(2i5, i3, 2(es15.6, 15x), es15.6)') iz, ia, nex, yieldfpex(iz, in, nex), xsfpex(iz, in, nex), &
 &            fpratio(iz, in, nex)
              Nfy=Nfy+1
            endif
            if (i == 1) then
              write(1, '(4es15.6)') Einc0, yieldfpex(iz, in, nex), xsfpex(iz, in, nex), fpratio(iz, in, nex)
            else
              write(1, '(2es15.6)') Einc0, xsZApost(iz, in)
            endif
            close (unit = 1)
          enddo
        endif
        if ( .not. fpexist(i, iz, in, -1)) fpexist(i, iz, in, -1) = .true.
      enddo
    enddo
!
! Write cumulative ff/fp production
!
    un = ''
    if (xsApre(ia) < fpeps .and. xsApost(ia) < fpeps .and. .not. fpaexist(ia)) cycle
    iz = 0
    fpfile = 'fp000.tot'//natstring(iso)
    write(fpfile(3:5), '(i3.3)') ia
    if ( .not. fpaexist(ia)) then
      fpaexist(ia) = .true.
      finalnuclide=trim(adjustl(massstring))
      open (unit = 1, file = fpfile, status = 'replace')
      reaction='('//parsym(k0)//',f)'
      quantity='fission yield'
      topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' for  A='//fpfile(3:5)
      col(1)='E'
      un(1)='MeV'
      col(2)='FP_yield'
      col(3)='FF_yield'
      col(4)='FP_xs'
      un(4)='mb'
      col(5)='FF_xs'
      un(5)='mb'
      Ncol=5
      call write_header(indent,topline,source,user,date,oformat)
      call write_target(indent)
      call write_reaction(indent,reaction,0.D0,0.D0,0,0)
      call write_residual(id2,iz,ia,finalnuclide)
      call write_quantity(id2,quantity)
      call write_datablock(id2,Ncol,Ninc,col,un)
      do nen = 1, nin0 - 1
        write(1, '(5es15.6)') eninc(nen), 0., 0., 0., 0.
      enddo
    else
      open (unit = 1, file = fpfile, status = 'old', position = 'append')
    endif
    write(1, '(5es15.6)') Einc0, yieldApost(ia), yieldApre(ia), xsApost(ia), xsApre(ia)
    close (unit = 1)
  enddo
  fpyieldfile = 'yieldZA0000.000.fis'//natstring(iso)
  if (Einc0 < 0.001) then
    write(fpyieldfile(8:15), '(es8.2)') Einc0
  else
    write(fpyieldfile(8:15), '(f8.3)') Einc0
    write(fpyieldfile(8:11), '(i4.4)') int(Einc0)
  endif
  open (unit = 1, file = fpyieldfile, status = 'replace')
  quantity='fission yield'
  topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' per Z, A'//' at '//Estr//' MeV'
  un = ''
  col(1)='Z'
  col(2)='A'
  col(3)='isomer'
  col(4)='FP_yield'
  col(5)='FF_yield'
  col(6)='FP_xs'
  un(6)='mb'
  col(7)='FF_xs'
  un(7)='mb'
  col(8)='Isomeric_ratio'
  Ncol=8
  MT = 459
  call write_header(indent,topline,source,user,date,oformat)
  call write_target(indent)
  call write_reaction(indent,reaction,0.D0,0.D0,MF,MT)
  call write_real(id2,'E-incident [MeV]',Einc0)
  call write_quantity(id2,quantity)
  call write_datablock(id2,Ncol,Nfy,col,un)
  do ia = 1, Atarget
    do iz = 1, Ztarget
      in = ia - iz
      if (in < 1 .or. in > Ninit) cycle
      if (xsZApre(iz, in) < fpeps .and. xsZApost(iz, in) < fpeps .and. .not. fpexist(1, iz, in, -1)) cycle
      if (xsZApre(iz, in) >= fpeps .and. xsZApost(iz, in) >= fpeps) &
 &      write(1, '(3(i6, 9x), 4es15.6)') iz, ia, -1, yieldZApost(iz, in), yieldZApre(iz, in), xsZApost(iz, in), xsZApre(iz, in)
      if (xsfpex(iz, in, 1) > 0.) then
        do nex = 0, 1
          if (xsZApre(iz, in) >= fpeps .and. xsZApost(iz, in) >= fpeps) &
 &          write(1, '(3(i6, 9x), 2(es15.6, 15x), es15.6)') iz, ia, nex, yieldfpex(iz, in, nex), xsfpex(iz, in, nex), &
 &          fpratio(iz, in, nex)
        enddo
      endif
    enddo
  enddo
  close (unit=1)
  write(*, '(/"Total          ", 4es15.6)') yieldtotpost, yieldtotpre, xstotpost, xstotpre
  return
end subroutine massdisout
! Copyright A.J. Koning 2021
