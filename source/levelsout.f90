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
  character(len=15) :: col(10)     ! header
  character(len=15) :: un(10)     ! header
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
!
! ********************* Nuclear structure file *************************
!
  Z = ZZ(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  massstring='   '
  write(massstring,'(i3)') A
  finalnuclide=trim(nuc(Z))//adjustl(massstring)
  resstring='000000'
  write(resstring(1:3),'(i3.3)') Z
  write(resstring(4:6),'(i3.3)') A
  structurefile = 'levels'//resstring//'.tot'
  open (unit = 1, file = structurefile, status = 'replace')
  quantity='Nuclear masses and levels'
  topline=trim(finalnuclide)//' '//trim(quantity)
!
! Mass and separation energies
!
  call write_header(topline,source,user,date,oformat)
  call write_residual(Z,A,finalnuclide)
  write(1,'("# parameters:")')
  call write_double(2,'nuclear mass [amu]',nucmass(Zix,Nix))
  call write_real(2,'neutron separation energy [MeV]',S(Zix,Nix,1))
  call write_real(2,'proton separation energy [MeV]',S(Zix,Nix,2))
  call write_real(2,'alpha separation energy [MeV]',S(Zix,Nix,6))
  call write_integer(2,'number of excited levels',nlev(Zix,Nix))
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
  col(8)='Half-life'
  un(8)='sec'
  col(9)='Assignment'
  col(10)='ENSDF'
  Ncol=10
  Nentries=0
  do i = 0, nlev(Zix, Nix)
    Nentries = Nentries + 1 + nbranch(Zix, Nix, i)
  enddo
  call write_datablock(quantity,Ncol,Nentries,col,un)
  do i = 0, nlev(Zix, Nix)
    if (tau(Zix, Nix, i) /= 0.) then
      write(1, '(i6, 9x, es15.6, 7x, f4.1, 13x, a1, 9x, i6, 36x, es15.6, 6x, 3a1, 1x, a18)') levnum(Zix, Nix, i),  &
 &      edis(Zix, Nix, i), jdis(Zix, Nix, i), cparity(parlev(Zix, Nix, i)), nbranch(Zix, Nix, i), tau(Zix, Nix, i), &
 &      eassign(Zix, Nix, i), jassign(Zix, Nix, i), passign(Zix, Nix, i), ENSDF(Zix, Nix, i)
    else
      write(1, '(i6, 9x, es15.6, 7x, f4.1, 13x, a1, 9x, i6, 36x, 15x, 6x, 3a1, 1x, a18)') i, edis(Zix, Nix, i), jdis(Zix, Nix, i), &
 &      cparity(parlev(Zix, Nix, i)), nbranch(Zix, Nix, i), eassign(Zix, Nix, i), jassign(Zix, Nix, i), passign(Zix, Nix, i), &
 &      ENSDF(Zix, Nix, i)
    endif
    do k = 1, nbranch(Zix, Nix, i)
      write(1, '(75x, "--->", i6, 9x, f8.4, 26x, a1)') branchlevel(Zix, Nix, i, k), branchratio(Zix, Nix, i, k) * 100., &
 &      bassign(Zix, Nix, i, k)
    enddo
  enddo
  close (unit = 1)
  return
end subroutine levelsout
! Copyright A.J. Koning 2024
