subroutine gammadecay(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Scheme for discrete gamma decay
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
!   sgl             ! single precision kind
! All global variables
!   numlev          ! maximum number of discrete levels
! Variables for discrete levels
!   flagelectron    ! flag for application of electron conversion coefficient
!   nlev            ! number of levels for nucleus
! Variables for nuclides
!   AA              ! mass number of residual nucleus
!   ZZ              ! charge number of residual nucleus
! Constants
!   cparity         ! parity (character)
!   nuc             ! symbol of nucleus
!   parsym          ! symbol of particle
! Variables for levels
!   bassign         ! flag for assignment of branching ratio
!   branchlevel     ! level to which branching takes place
!   branchratio     ! gamma - ray branching ratio to level
!   conv            ! conversion coefficient
!   edis            ! energy of level
!   ENSDF           ! string from original ENSDF discrete level file
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
  character(len=6)  :: finalnuclide !
  character(len=18) :: reaction   ! reaction
  character(len=132) :: topline    ! topline
  character(len=15) :: col(8)     ! header
  character(len=15) :: un(8)     ! header
  character(len=80) :: quantity   ! quantity
  character(len=7)   :: decayfile                     ! decay data file
  character(len=10)  :: discretefile                  ! file with discrete level data
  integer            :: A                             ! mass number of target nucleus
  integer            :: i                             ! counter
  integer            :: j                             ! counter
  integer            :: k                             ! designator for particle
  integer            :: l                             ! counter
  integer            :: N                             ! neutron number of residual nucleus
  integer            :: Ncol                          ! number of columns
  integer            :: ng                            ! number of gamma lines
  integer            :: Ngam(numlev)                  ! counter for discrete gamma rays
  integer            :: Nix                           ! neutron number index for residual nucleus
  integer            :: type                          ! particle type
  integer            :: Z                             ! charge number of target nucleus
  integer            :: Zix                           ! charge number index for residual nucleus
  real(sgl)          :: br(numlev, numlev*numlev/2)   ! branching ratio multiplied by initial flux
  real(sgl)          :: brtmp                         ! help variable
  real(sgl)          :: egam(numlev, numlev*numlev/2) ! gamma energy
  real(sgl)          :: egamtmp                       ! help variable
  real(sgl)          :: flux(0:numlev)                ! flux in level (normalized to 1)
  real(sgl)          :: total                         ! help variable
  real(sgl)          :: totalg                        ! help variable
  real(sgl)          :: yieldg(numlev)                ! total discrete gamma yield per level
!
! ********************* Construct gamma decay scheme *******************
!
  do i = 1, nlev(Zix, Nix)
    ng = 0
    do j = 0, i - 1
      flux(j) = 0.
    enddo
    flux(i) = 1.
    total = 0.
    totalg = 0.
    do j = i, 1, - 1
      if (tau(Zix, Nix, j) /= 0.) cycle
      do k = 1, nbranch(Zix, Nix, j)
        l = branchlevel(Zix, Nix, j, k)
        if (flux(j) == 0.) cycle
        ng = ng + 1
        egam(i, ng) = edis(Zix, Nix, j) - edis(Zix, Nix, l)
        br(i, ng) = branchratio(Zix, Nix, j, k) * flux(j)
        flux(l) = flux(l) + br(i, ng)
        total = total + br(i, ng)
        if (flagelectron) then
          totalg = totalg + br(i, ng) / (1. + conv(Zix, Nix, j, k))
        else
          totalg = total
        endif
      enddo
    enddo
    Ngam(i) = ng
    yieldg(i) = totalg
    if (total == 0.) cycle
!
! Normalize intensities to 1.
!
    do j = 1, ng
      br(i, j) = br(i, j) / total
    enddo
!
! Sort discrete gamma ray lines in descending order
!
    do j = 1, ng
      do k = 1, j
        if (egam(i, j) > egam(i, k)) then
          egamtmp = egam(i, k)
          brtmp = br(i, k)
          egam(i, k) = egam(i, j)
          br(i, k) = br(i, j)
          egam(i, j) = egamtmp
          br(i, j) = brtmp
        endif
      enddo
    enddo
  enddo
!
! Write decay information to file
!
  un = ''
  col(1)='E'
  un(1)='MeV'
  col(2)='fraction'
  Ncol=2
  quantity='Discrete gamma decay'
  Z = ZZ(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  type = 2 * Zix + Nix
  reaction='('//parsym(k0)//','//parsym(type)//')'
  decayfile = 'decay.'//parsym(type)
  massstring='   '
  write(massstring,'(i3)') A
  finalnuclide=trim(nuc(Z))//adjustl(massstring)
  open (unit = 1, file = decayfile, status = 'replace')
  topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
  call write_header(topline,source,user,date,oformat)
  call write_target
  call write_reaction(reaction,0.D0,0.D0,0,0)
  call write_residual(Z,A,finalnuclide)
  call write_datablock(quantity,Ncol,nlev(Zix, Nix),col,un)
  do i = 1, nlev(Zix, Nix)
    write(1, '(2(i6, 9x), es15.6)') i, Ngam(i), yieldg(i)
    do j = 1, Ngam(i)
      write(1, '(2es15.6)') egam(i, j), br(i, j)
    enddo
  enddo
  close (unit = 1)
!
! Write discrete level file for gamma transition probabilities
!
  N = A - Z
  discretefile = 'discrete.'//parsym(type)
  open (unit = 1, file = discretefile, status = 'replace')
  quantity='Discrete levels'
  topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
  un = ''
  col(1)='Number'
  col(2)='Energy'
  un(2)='MeV'
  col(3)='Spin'
  col(4)='Parity'
  col(5)='Branchin ratio'
  un(5)='%'
  col(6)='Half-life'
  un(6)='sec'
  col(7)='Assignment'
  col(8)='ENSDF'
  Ncol=8
  call write_header(topline,source,user,date,oformat)
  call write_target
  call write_residual(Z,A,finalnuclide)
  call write_datablock(quantity,Ncol,nlev(Zix, Nix),col,un)
  do i = 0, nlev(Zix, Nix)
    if (tau(Zix, Nix, i) /= 0.) then
      write(1, '(1x, i3, 1x, es12.5, 1x, f4.1, 3x, a1, 9x, i2, 15x,  es10.3, 7x, 2a1, a18)') levnum(Zix, Nix, i),  &
 &      edis(Zix, Nix, i), jdis(Zix, Nix, i), cparity(parlev(Zix, Nix, i)), nbranch(Zix, Nix, i), tau(Zix, Nix, i), &
 &      jassign(Zix, Nix, i), passign(Zix, Nix, i), ENSDF(Zix, Nix, i)
    else
      write(1, '(1x, i3, 1x, es12.5, 1x, f4.1, 3x, a1, 9x, i2, 32x, 2a1, a18)') i, edis(Zix, Nix, i), jdis(Zix, Nix, i), &
 &      cparity(parlev(Zix, Nix, i)), nbranch(Zix, Nix, i), jassign(Zix, Nix, i), passign(Zix, Nix, i), ENSDF(Zix, Nix, i)
    endif
    do k = 1, nbranch(Zix, Nix, i)
      write(1, '(33x, "--->", i3, 2x, f8.4, 18x, a1)') branchlevel(Zix, Nix, i, k), branchratio(Zix, Nix, i, k) * 100., &
 &      bassign(Zix, Nix, i, k)
    enddo
  enddo
  close (unit = 1)
  return
end subroutine gammadecay
! Copyright A.J. Koning 2021
