subroutine nudisout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of number of fission neutrons
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
! 2022-02-28: Output of average emission energies per fragment to separate files
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! All global variables
!   numneu        ! number of neutrons
!   numnu         ! number of neutrons from fission
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
!   Ztarget         ! charge number of target nucleus
! Variables for energies
!   Einc0           ! incident energy in MeV
! Constants
!   iso             ! counter for isotope
!   natstring       ! string extension for file names
!   parname         ! name of particle
!   parsym          ! symbol of particle
! Variables for mass distribution
!   EaverageA       ! average emission energy per A
!   EaverageZA      ! average emission energy per (Z,A)
!   nuA             ! nu per A
!   nuZA            ! nu per Z,A
!   nubar           ! average nu
!   Pdisnu          ! prompt fission neutrons distribution
!   Pdisnuav        ! average prompt fission neutrons distribution
!
! *** Declaration of local data
!
  implicit none
  character(len=18) :: reaction   ! reaction
  character(len=132) :: topline    ! topline
  character(len=15) :: Estr
  character(len=15) :: col(4)     ! header
  character(len=15) :: un(4)     ! header
  character(len=80) :: quantity   ! quantity
  character(len=8)  :: Estring   ! string for incident energy
  character(len=132):: nufile    ! file for nubar
  character(len=132):: Eline     ! line for incident energy
  integer           :: i         ! counter
  integer           :: k         ! counter
  integer           :: ia        ! mass number
  integer           :: in        ! neutron number
  integer           :: iz        ! charge number
  integer           :: Ncol      ! number of columns
  integer           :: type      ! particle type
!
! Write results to separate files
!
! nu per number, P(nu) and nubar as function of Z and A
!
  Estring = '0000.000'
  if (Einc0 < 0.001) then
    write(Estring(1:8), '(es8.2)') Einc0
  else
    write(Estring(1:8), '(f8.3)') Einc0
    write(Estring(1:4), '(i4.4)') int(Einc0)
  endif
  Eline = '# E-incident =          MeV'
  if (Einc0 < 0.001) then
    write(Eline(16:23), '(es8.2)') Einc0
  else
    write(Eline(16:23), '(f8.3)') Einc0
  endif
  Estr=''
  write(Estr,'(es13.6)') Einc0
  write(*, '(/" +++ NUMBER OF PROMPT FISSION NEUTRONS AND GAMMAS +++")')
  un = ''
  do type = 0, 6
    if (parskip(type)) cycle
    if (nubar(type) == 0.) cycle
    write(*, '(/" nubar for ", a8, f10.5)') parname(type), nubar(type)
!
! P(nu)
!
    nufile = 'Pnux'//Estring//'.fis'//natstring(iso)
    nufile(4:4) = parsym(type)
    open (unit = 1, file = nufile, status = 'replace')
    quantity='multiplicity'
    reaction='('//parsym(k0)//',f)'
    topline=trim(targetnuclide)//trim(reaction)//' prompt '//trim(adjustl(parname(type)))//' '//trim(quantity)//' at '//Estr//' MeV'
    col(1)='number'
    col(2)='nu'
    Ncol=2
    call write_header(topline,source,user,date,oformat)
    call write_target
    call write_reaction(reaction,0.d0,0.d0,0,0)
    call write_real(2,'E-incident [MeV]',Einc0)
    call write_char(2,'ejectile',parname(type))
    call write_real(2,'nubar-prompt',nubar(type))
    call write_real(2,'average P(nu)',Pdisnuav(type))
    call write_datablock(quantity,Ncol,numnu+1,col,un)
    write(*, '(/" P(nu) for ", a8/)') parname(type)
    write(*, '(" number    nu ")')
    do i = 0, numnu
      if (Pdisnu(type, i) > 0. .or. i <= 4) then
        write(*, '(i3, es15.6)') i, Pdisnu(type, i)
        write(1, '(i6, 9x, es15.6)') i, Pdisnu(type, i)
      endif
    enddo
    close (unit = 1)
    write(*, '(/" Average P(nu) for ", a8, f10.5/)') parname(type), Pdisnuav(type)
!
! nu(A)
!
    nufile = 'nuxA'//Estring//'.fis'//natstring(iso)
    nufile(3:3) = parsym(type)
    open (unit = 1, file = nufile, status = 'replace')
    quantity='multiplicity'
    reaction='('//parsym(k0)//',f)'
    topline=trim(targetnuclide)//trim(reaction)//' prompt '//trim(adjustl(parname(type)))//' '//trim(quantity)// &
 &    ' as function of mass - nu(A)'//' at '//Estr//' MeV'
    col(1)='A'
    col(2)='nu'
    Ncol=2
    call write_header(topline,source,user,date,oformat)
    call write_target
    call write_reaction(reaction,0.d0,0.d0,0,0)
    call write_real(2,'E-incident [MeV]',Einc0)
    call write_char(2,'ejectile',parname(type))
    call write_real(2,'nubar-prompt',nubar(type))
    call write_datablock(quantity,Ncol,Atarget,col,un)
    write(*, '(/"  nu(A) for ", a8/)') parname(type)
    write(*, '("  A        nu   ")')
    do ia = 1, Atarget
      write(*, '(i3, es17.6)') ia, nuA(type, ia)
      write(1, '(i6, 9x, es15.6)') ia, nuA(type, ia)
    enddo
    close (unit = 1)
!
! nu(Z,A)
!
    nufile = 'nuxZA'//Estring//'.fis'//natstring(iso)
    nufile(3:3) = parsym(type)
    open (unit = 1, file = nufile, status = 'replace')
    quantity='multiplicity'
    reaction='('//parsym(k0)//',f)'
    topline=trim(targetnuclide)//trim(reaction)//' prompt '//trim(adjustl(parname(type)))//' '//trim(quantity)// &
 &    ' as function of nuclide - nu(Z,A)'//' at '//Estr//' MeV'
    col(1)='Z'
    col(2)='A'
    col(3)='nu'
    Ncol=3
    call write_header(topline,source,user,date,oformat)
    call write_target
    call write_reaction(reaction,0.d0,0.d0,0,0)
    call write_real(2,'E-incident [MeV]',Einc0)
    call write_char(2,'ejectile',parname(type))
    write(*, '(/"  nu(Z,A) for ", a8/)') parname(type)
    write(*, '("  Z   A       nu")')
    k = 0
    do iz = 1, Ztarget
      do ia = iz + 1, Atarget
        in = ia - iz
        if (in < numneu) then
          if (nuZA(type, iz, in) > 0.) k = k + 1
        endif
      enddo
    enddo
    call write_datablock(quantity,Ncol,k,col,un)
    do iz = 1, Ztarget
      do ia = iz + 1, Atarget
        in = ia - iz
        if (in < numneu) then
          if (nuZA(type, iz, in) > 0.) then
            write(*, '(i3, i4, es17.6)') iz, ia, nuZA(type, iz, in)
            write(1, '(2(i6, 9x), es15.6)') iz, ia, nuZA(type, iz, in)
          endif
        endif
      enddo
    enddo
    close (unit = 1)
  enddo
  if (fymodel <=2) return
!
! E-average(Z,A) and E-average(A)
!
  nufile = 'EavZA'//Estring//'.fis'//natstring(iso)
  open (unit=1, file=nufile, status='replace')
  quantity='emission energy'
  reaction='('//parsym(k0)//',f)'
  topline=trim(targetnuclide)//trim(reaction)//' average '//trim(quantity)//'  as function of nuclide - E-av(Z,A)'// &
 &  ' at '//Estr//' MeV'
  col(1)='Z'
  col(2)='A'
  col(3)='gamma'
  un(3)='MeV'
  col(4)='neutron'
  un(4)='MeV'
  Ncol=4
  call write_header(topline,source,user,date,oformat)
  call write_target
  call write_reaction(reaction,0.d0,0.d0,0,0)
  call write_real(2,'E-incident [MeV]',Einc0)
  write(*,'(/"  Average emission energy per (Z,A)")')
  write(*,'(/"  Z  A      gamma    neutron"/)')
  k = 0
  do iz = 1, Ztarget
    do ia = iz+1, Atarget
      in = ia - iz
      if (in < numneu) then
        if (EaverageZA(0,iz,in) > 0. .or. EaverageZA(1,iz,in) > 0.) k = k + 1
      endif
    enddo
  enddo
  call write_datablock(quantity,Ncol,k,col,un)
  do iz = 1, Ztarget
    do ia = iz+1, Atarget
      in = ia - iz
      if (in < numneu) then
        if (EaverageZA(0,iz,in) > 0. .or. EaverageZA(1,iz,in) > 0.) then
          write(*, '(i3,i4,2es12.5)') iz, ia, (EaverageZA(type,iz,in), type=0,1)
          write(1, '(2(i6,9x),2es15.6)') iz, ia, (EaverageZA(type,iz,in), type=0,1)
        endif
      endif
    enddo
  enddo
  close (unit = 1)
  nufile='EavA'//Estring//'.fis'//natstring(iso)
  open (unit=1, file=nufile, status='replace')
  quantity='emission energy'
  reaction='('//parsym(k0)//',f)'
  topline=trim(targetnuclide)//trim(reaction)//' average '//trim(quantity)//' as function of mass - E-av(A)'//' at '//Estr//' MeV'
  col(1)='A'
  col(2)='gamma'
  col(3)='neutron'
  Ncol=3
  call write_header(topline,source,user,date,oformat)
  call write_target
  call write_reaction(reaction,0.d0,0.d0,0,0)
  call write_real(2,'E-incident [MeV]',Einc0)
  call write_datablock(quantity,Ncol,Atarget,col,un)
  write(*,'(/"  Average emission energy per A")')
  write(*,'(/"  A     gamma    neutron"/)')
  do ia = 1, Atarget
    write(*, '(i3,2es12.5)') ia, (EaverageA(type,ia), type=0,1)
    write(1, '(i6,9x,2es15.6)') ia, (EaverageA(type,ia), type=0,1)
  enddo
  close (unit = 1)
  return
end subroutine nudisout
! Copyright A.J. Koning 2021
