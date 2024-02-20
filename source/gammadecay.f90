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
  character(len=6)  :: resstring !
  character(len=6)  :: finalnuclide !
  character(len=132) :: topline    ! topline
  character(len=15) :: col(10)     ! header
  character(len=15) :: un(10)     ! header
  character(len=80) :: quantity   ! quantity
  character(len=15)  :: gammafile                     ! file with gamma decay data
  integer            :: A                             ! mass number of target nucleus
  integer            :: i                             ! counter
  integer            :: j                             ! counter
  integer            :: k                             ! designator for particle
  integer            :: l                             ! counter
  integer            :: Nentries                      ! number of entries
  integer            :: Ncol                          ! number of columns
  integer            :: ng                            ! number of gamma lines
  integer            :: Ngam(numlev)                  ! counter for discrete gamma rays
  integer            :: Nix                           ! neutron number index for residual nucleus
  integer            :: Z                             ! charge number of target nucleus
  integer            :: Zix                           ! charge number index for residual nucleus
  integer            :: plev(numlev, numlev*numlev/2) ! level
  integer            :: plevtmp                       ! help variable
  integer            :: dlev(numlev, numlev*numlev/2) ! level
  integer            :: dlevtmp                       ! help variable
  real(sgl)          :: br(numlev, numlev*numlev/2)   ! branching ratio multiplied by initial flux
  real(sgl)          :: brtmp                         ! help variable
  real(sgl)          :: penergy(numlev, numlev*numlev/2) ! level energy
  real(sgl)          :: penergytmp                       ! help variable
  real(sgl)          :: denergy(numlev, numlev*numlev/2) ! level energy
  real(sgl)          :: denergytmp                       ! help variable
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
        plev(i, ng) = j
        penergy(i, ng) = edis(Zix, Nix, j)
        dlev(i, ng) = l
        denergy(i, ng) = edis(Zix, Nix, l)
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
          plevtmp = plev(i, k)
          dlevtmp = dlev(i, k)
          penergytmp = penergy(i, k)
          denergytmp = denergy(i, k)
          egamtmp = egam(i, k)
          brtmp = br(i, k)
          plev(i, k) = plev(i, j)
          dlev(i, k) = dlev(i, j)
          penergy(i, k) = penergy(i, j)
          denergy(i, k) = denergy(i, j)
          egam(i, k) = egam(i, j)
          br(i, k) = br(i, j)
          plev(i, j) = plevtmp
          dlev(i, j) = dlevtmp
          penergy(i, j) = penergytmp
          denergy(i, j) = denergytmp
          egam(i, j) = egamtmp
          br(i, j) = brtmp
        endif
      enddo
    enddo
  enddo
!
! Write decay information to file
!
  quantity='Cumulative discrete gamma decay'
  Z = ZZ(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  resstring='000000'
  write(resstring(1:3),'(i3.3)') Z
  write(resstring(4:6),'(i3.3)') A
  gammafile = 'gamma'//resstring//'.tot'
  open (unit = 1, file = gammafile, status = 'replace')
  topline=trim(targetnuclide)//' '//trim(quantity)
  call write_header(topline,source,user,date,oformat)
  call write_target
  Nentries=0
  do i = 1, nlev(Zix, Nix)
    Nentries = Nentries + 1 + Ngam(i)
  enddo
  un = ''
  col(1)='Level'
  col(2)='Energy'
  un(2)='MeV'
  col(3)='Num._branch'
  col(4)='Total_yield'
  col(5)='Parent_level'
  col(6)='Daughter_level'
  col(7)='Gamma_energy'
  un(7)='MeV'
  col(8)='Branch._ratio'
  col(9)='Parent_energy'
  un(9)='MeV'
  col(10)='Daughter_energy'
  un(10)='MeV'
  Ncol=10
  call write_datablock(quantity,Ncol,Nentries,col,un)
  do i = 1, nlev(Zix, Nix)
    write(1, '(2(i6, 9x, es15.6))') i, edis(Zix, Nix, i), Ngam(i), yieldg(i)
    do j = 1, Ngam(i)
      write(1, '(60x, 2(i6, 9x), 4es15.6)') plev(i, j), dlev(i, j), egam(i, j), br(i, j), penergy(i, j), denergy(i, j)
    enddo
  enddo
  close (unit = 1)
  return
end subroutine gammadecay
! Copyright A.J. Koning 2021
