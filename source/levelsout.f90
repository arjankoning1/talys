subroutine levelsout(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of basic nuclear structure: masses and discrete levels
!
! Author    : Arjan Koning
!
! 2024-02-11: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl             ! single precision kind
! All global variables
!   numlev          ! maximum number of discrete levels
! Variables for discrete levels
!   nlev            ! number of levels for nucleus
! Variables for nuclides
!   AA              ! mass number of residual nucleus
!   ZZ              ! charge number of residual nucleus
! Constants
!   nuc             ! symbol of nucleus
! Variables for levels
!   bassign         ! flag for assignment of branching ratio
!   branchlevel     ! level to which branching takes place
!   branchratio     ! gamma - ray branching ratio to level
!   conv            ! conversion coefficient
!   edis            ! energy of level
!   ENSDF           ! string from original ENSDF discrete level file
!   eassign         ! flag for assignment of energy
!   jassign         ! flag for assignment of spin
!   jdis            ! spin of level
!   nbranch         ! number of branching levels
!   parlev          ! parity of level
!   passign         ! flag for assignment of parity
!   tau             ! lifetime of state in seconds
!
! *** Declaration of local data
!
  implicit none
  character(len=3)  :: massstring !
  character(len=6)  :: resstring !
  character(len=6)  :: finalnuclide !
  character(len=132) :: topline    ! topline
  character(len=15) :: col(11)     ! header
  character(len=15) :: taustring
  character(len=15) :: un(11)     ! header
  character(len=80) :: quantity   ! quantity
  character(len=20)  :: structurefile                 ! file with nuclear structure data
  integer            :: A                             ! mass number of target nucleus
  integer            :: i                             ! counter
  integer            :: k                             ! designator for particle
  integer            :: Ncol                          ! number of columns
  integer            :: Nentries                      ! number of entries
  integer            :: Nix                           ! neutron number index for residual nucleus
  integer            :: Z                             ! charge number of target nucleus
  integer            :: Zix                           ! charge number index for residual nucleus
  integer            :: indent
  integer            :: id2
  integer            :: id4
!
! ********************* Nuclear structure file *************************
!
  indent = 0
  id2 = indent + 2
  id4 = indent + 4
  Z = ZZ(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  massstring='   '
  write(massstring,'(i3)') A
  finalnuclide=trim(nuc(Z))//adjustl(massstring)
  resstring='000000'
  write(resstring(1:3),'(i3.3)') Z
  write(resstring(4:6),'(i3.3)') A
  if (flagblocklevels) then
    structurefile='levels.out'//natstring(iso)
    if (nin == Ninclow + 1 .and. Zix == 0 .and. Nix == 0) then
      open (unit = 1, file = structurefile, status = 'unknown')
    else
      open (unit = 1, file = structurefile, status = 'unknown', position = 'append')
    endif
  else
    structurefile = 'levels'//resstring//'.out'
    open (unit = 1, file = structurefile, status = 'replace')
  endif
  quantity='Nuclear masses and levels'
  topline=trim(finalnuclide)//' '//trim(quantity)
!
! Mass and separation energies
!
  call write_header(indent,topline,source,user,date,oformat)
  call write_residual(indent,Z,A,finalnuclide)
  call write_char(id2,'parameters','')
  call write_double(id4,'nuclear mass [amu]',nucmass(Zix,Nix))
  call write_double(id4,'neutron separation energy [MeV]',S(Zix,Nix,1))
  call write_double(id4,'proton separation energy [MeV]',S(Zix,Nix,2))
  call write_double(id4,'alpha separation energy [MeV]',S(Zix,Nix,6))
  call write_integer(id4,'number of excited levels',nlev(Zix,Nix))
!
! Discrete levels
!
  un = ''
  col(1)='Level'
  col(2)='Energy'
  un(2)='MeV'
  col(3)='Spin'
  col(4)='Parity'
  col(5)='Num._Branches'
  col(6)='Branch._level'
  col(7)='Branch._ratio'
  un(7)='%'
  col(8)='Intern._conv'
  col(9)='Half-life'
  un(9)='sec'
  col(10)='Assignment'
  col(11)='ENSDF'
  Ncol=11
  Nentries=0
  do i = 0, nlev(Zix, Nix)
    Nentries = Nentries + 1 + nbranch(Zix, Nix, i)
  enddo
  call write_quantity(id2,quantity)
  call write_datablock(id2,Ncol,Nentries,col,un)
  do i = 0, nlev(Zix, Nix)
    taustring=''
    if (tauripl(Zix, Nix, i) /= 0.) write(taustring(1:15),'(es15.6)') tauripl(Zix, Nix, i)
    write(1, '(i6, 9x, es15.6, 7x, f4.1, 13x, a1, 9x, i6, 51x, a15, 6x, 3a1, 1x, a18)') levnum(Zix, Nix, i),  &
 &    edis(Zix, Nix, i), jdis(Zix, Nix, i), cparity(parlev(Zix, Nix, i)), nbranch(Zix, Nix, i), taustring, &
 &    eassign(Zix, Nix, i), jassign(Zix, Nix, i), passign(Zix, Nix, i), ENSDF(Zix, Nix, i)
    do k = 1, nbranch(Zix, Nix, i)
      write(1, '(75x, "--->", i6, 9x, f8.4, 4x, es15.6, 26x, a1)') branchlevel(Zix, Nix, i, k),  & 
 &      branchratio(Zix, Nix, i, k) * 100., conv(Zix, Nix, i, k), bassign(Zix, Nix, i, k)
    enddo
  enddo
  close (unit = 1)
  return
end subroutine levelsout
! Copyright A.J. Koning 2024
