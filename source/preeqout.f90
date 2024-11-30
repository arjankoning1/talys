subroutine preeqout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of pre-equilibrium cross sections
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
! Variables for preequilibrium
!   flag2comp       ! flag for two - component pre - equilibrium model
!   flaggshell      ! flag for energy dependence of single particle level den
!   g               ! single - particle level density parameter
!   gn              ! single - particle neutron level density parameter
!   gp              ! single - particle proton level density parameter
!   pairmodel       ! model for preequilibrium pairing energy
! Variables for level density
!   alev            ! level density parameter
! Variables for energy grid
!   ebegin          ! first energy point of energy grid
!   egrid           ! outgoing energy grid
! Variables for energies
!   eend            ! last energy point of energy grid
!   Etotal          ! total energy of compound system (target + projectile)
! Variables for nuclides
!   parskip         ! logical to skip outgoing particle
! Constants
!   parname         ! name of particle
! Variables for incident channel
!   xspreeq         ! preeq. cross section per particle typ and outgoing energye
!   xspreeqsum      ! total preequilibrium cross section summed over particles
!   xspreeqtot      ! preequilibrium cross section per particle type
! Variables for preequilibrium initialization
!   Efermi          ! depth of Fermi well
!   maxexc          ! maximal exciton number
!   maxpar          ! maximal particle number
!   RnJ             ! spin distribution for particle - hole stat
!   RnJsum          ! (2J + 1) * sum over spin distributions
! Variables for preequilibrium
!   Esurf           ! well depth for surface interaction
!   preeqnorm       ! preequilibrium normalizati
!   xspreeqbu       ! preequilibrium cross section per particle type and outgoing energy for brea
!   xspreeqki       ! preequilibrium cross section per particle type and outgoing energy for knoc
!   xspreeqps       ! preequilibrium cross section per particle type and outgoing energy for pick
!   xspreeqtotbu    ! preequilibrium cross section per particle type for breakup
!   xspreeqtotki    ! preequilibrium cross section per particle type for knockout and inelas
!   xspreeqtotps    ! preequilibrium cross section per particle type for pickup and strippin
!   xsstep          ! preeq. cross section per particle type, stage and outgoing E
!   xssteptot       ! preequilibrium cross section per particle type and stage
!
! *** Declaration of local data
!
  implicit none
  logical   :: surfwell  ! flag for surface effects in finite well
  character(len=3)  :: massstring !
  character(len=6)  :: finalnuclide !
  character(len=20) :: phdfile    ! file for output
  character(len=20) :: preeqfile  ! file for output
  character(len=18) :: reaction   ! reaction
  character(len=132):: topline    ! topline
  character(len=15) :: col(14)     ! header
  character(len=15) :: un(14)     ! units
  character(len=80) :: quantity   ! quantity
  integer           :: Z             ! charge number of target nucleus
  integer           :: A             ! mass number of target nucleus
  integer           :: Ncol       ! number of columns
  integer   :: h         ! help variable
  integer   :: J         ! spin of level
  integer   :: k         ! designator for particle
  integer   :: n         ! exciton number
  integer   :: nen       ! energy counter
  integer   :: Nix       ! neutron number index for residual nucleus
  integer   :: p         ! particle number
  integer   :: type      ! particle type
  integer   :: Zix       ! charge number index for residual nucleus
  real(sgl) :: damp      ! shell damping factor
  real(sgl) :: Eex       ! excitation energy
  real(sgl) :: gs        ! single-particle level density parameter
  real(sgl) :: gsn       ! single-particle neutron level density parameter
  real(sgl) :: gsp       ! single-particle proton level density parameter
  real(sgl) :: ignatyuk  ! function for energy dependent level density parameter a
  real(sgl) :: nonpski   ! preequilibrium cross section without pickup etc.
  real(sgl) :: phdens    ! function for particle-hole state density
  real(sgl) :: phdens2   ! function for two-component particle-hole state density
  real(sgl) :: preeqpair ! pre-equilibrium pairing energy
!
! ************************ Pre-equilibrium *****************************
!
! 1. Output of particle-hole state densities
!
! ignatyuk  : function for energy dependent level density parameter a
! phdens2   : function for two-component particle-hole state density
!
  Zix = 0
  Nix = 0
  surfwell = .false.
  if (nin == Ninc) then
    Z = ZZ(Zix, Nix, 0) 
    A = AA(Zix, Nix, 0) 
    massstring='   '  
    write(massstring,'(i3)') A
    finalnuclide=trim(nuc(Z))//adjustl(massstring)
    write(*, '(/," ++++++++++ PARTICLE-HOLE STATE DENSITIES ++++++++++",/)')
    phdfile = 'ph_density.out'
    quantity = 'particle-hole state density'
    topline=trim(finalnuclide)//' '//trim(quantity)
    open (unit = 1, file = phdfile, status = 'replace')
    call write_header(topline,source,user,date,oformat)
    call write_residual(Z,A,finalnuclide)
    un='MeV^-1'
    col=''
    un(1)='MeV'
    un(2)='MeV'
    col(1)='Ex'
    col(2)='P(3)'
    if ( .not. flag2comp) then
      col(3)='g'
      col(4)='1p1h'
      col(5)='2p1h'
      col(6)='2p2h'
      col(7)='3p2h'
      col(8)='3p3h'
      col(9)='4p3h'
      col(10)='4p4h'
      col(11)='5p4h'
      Ncol=11
      call write_datablock(quantity,Ncol,int(Etotal),col,un)
      do nen = 1, int(Etotal)
        Eex = real(nen)
        gs = g(0, 0)
        if (flaggshell) gs = gs * ignatyuk(Zix, Nix, Eex, 0) / alev(0, 0)
        write(1, '(11es15.6)') Eex, preeqpair(Zix, Nix, 3, Eex, pairmodel), gs, &
 &        ((phdens(Zix, Nix, h + k, h, gs, Eex, Efermi, surfwell), k = 0, 1), h = 1, 4)
      enddo
    else
      col(3)='g(p)'
      col(4)='g(n)'
      col(5)='1p1h0p0h'
      col(6)='0p0h1p1h'
      col(7)='1p1h1p0h'
      col(8)='1p0h1p1h'
      col(9)='2p1h0p0h'
      col(10)='0p0h2p1h'
      col(11)='2p2h0p0h'
      col(12)='0p0h2p2h'
      col(13)='1p1h1p1h'
      Ncol=13
      call write_datablock(quantity,Ncol,int(Etotal),col,un)
      do nen = 1, int(Etotal)
        Eex = real(nen)
        gsp = gp(0, 0)
        gsn = gn(0, 0)
        if (flaggshell) then
          damp = ignatyuk(Zix, Nix, Eex, 0) / alev(0, 0)
          gsp = gsp * damp
          gsn = gsn * damp
        endif
        write(1, '(13es15.6)') Eex, preeqpair(Zix, Nix, 3, Eex, pairmodel), gsp, gsn, &
 &        phdens2(Zix, Nix, 1, 1, 0, 0, gsp, gsn, Eex, Efermi, surfwell), &
 &        phdens2(Zix, Nix, 0, 0, 1, 1, gsp, gsn, Eex, Efermi, surfwell), &
          phdens2(Zix, Nix, 1, 1, 1, 0, gsp, gsn, Eex, Efermi, surfwell), &
          phdens2(Zix, Nix, 1, 0, 1, 1, gsp, gsn, Eex, Efermi, surfwell), &
          phdens2(Zix, Nix, 2, 1, 0, 0, gsp, gsn, Eex, Efermi, surfwell), &
          phdens2(Zix, Nix, 0, 0, 2, 1, gsp, gsn, Eex, Efermi, surfwell), &
          phdens2(Zix, Nix, 2, 2, 0, 0, gsp, gsn, Eex, Efermi, surfwell), &
          phdens2(Zix, Nix, 0, 0, 2, 2, gsp, gsn, Eex, Efermi, surfwell), &
          phdens2(Zix, Nix, 1, 1, 1, 1, gsp, gsn, Eex, Efermi, surfwell)
      enddo
    endif
    quantity = 'particle-hole spin distribution'
    un = ''
    col(1)='n'
    do J = 0, 8
      col(2+J)='J=  '
      write(col(2+J)(3:4),'(i2)') J
    enddo
    col(11)='Sum'
    Ncol=11
    call write_datablock(quantity,Ncol,maxexc,col,un)
    do n = 1, maxexc
      write(1, '(3x, i6, 6x, 10es15.6)') n, ((2* J + 1)*RnJ(n, J), J = 0, 8), RnJsum(n)
    enddo
    close(1)
    call write_outfile(phdfile,flagoutall)
  endif
!
! 2. Output of pre-equilibrium cross sections
!
  preeqfile = 'preeq.out'
  if (preeqfirst) then
    open (unit = 1, file = preeqfile, status = 'replace')
    preeqfirst = .false.
  else
    open (unit = 1, file = preeqfile, status = 'old', position = 'append')
  endif
  un='mb'
  col=''
  col(2)='Total'
  col(3)='p=1'
  col(4)='p=2'
  col(5)='p=3'
  col(6)='p=4'
  col(7)='p=5'
  col(8)='p=6'
  col(9)='Exciton model'
  col(10)='Pickup/strip'
  col(11)='Knockout'
  col(12)='Breakup'
  Ncol=12
  write(*, '(/" ++++++++++ TOTAL PRE-EQUILIBRIUM CROSS SECTIONS ++++++++++",/)')
  if (preeqnorm /= 0.) write(*, '(/" Pre-equilibrium normalization factor: ", f8.5/)') preeqnorm
  col(1)='particle'
  un(1)=''
  reaction='('//parsym(k0)//',x)'
  quantity = 'preequilibrium cross sections'
  topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
  call write_header(topline,source,user,date,oformat)
  call write_target
  call write_reaction(reaction,0.D0,0.D0,0,0)
  call write_real(2,'E-incident [MeV]',Einc)
  call write_real(2,'Pre-equilibrium cross section [mb]', xspreeqsum)
  call write_datablock(quantity,Ncol,7,col,un)
  do type = 0, 6
    nonpski = xspreeqtot(type) - xspreeqtotps(type) - xspreeqtotki(type) - xspreeqtotbu(type)
    write(1, '(4x, a8, 4x, 11es15.6)') parname(type),xspreeqtot(type), (xssteptot(type, p), p = 1, maxpar), nonpski, &
   &    xspreeqtotps(type), xspreeqtotki(type), xspreeqtotbu(type)
  enddo
  quantity = 'preequilibrium emission spectra'
  do type = 0, 6
    if (parskip(type)) cycle
    if (ebegin(type) >= eend(type)) cycle
    col(1)='E-out'
    un(1)='MeV'
    reaction='('//parsym(k0)//',x'//parsym(type)//')'
    topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
    call write_header(topline,source,user,date,oformat)
    call write_target
    call write_reaction(reaction,0.D0,0.D0,0,0)
    call write_real(2,'E-incident [MeV]',Einc)
    call write_real(2,'Pre-equilibrium cross section [mb]', xspreeqtot(type))
    call write_datablock(quantity,Ncol,eend(type)-ebegin(type)+1,col,un)
    do nen = ebegin(type), eend(type)
      nonpski = xspreeq(type, nen) - xspreeqps(type, nen) - xspreeqki(type, nen) - xspreeqbu(type, nen)
      write(1, '(12es15.6)') egrid(nen), xspreeq(type, nen), (xsstep(type, p, nen), p = 1, maxpar), nonpski, &
 &      xspreeqps(type, nen), xspreeqki(type, nen), xspreeqbu(type, nen)
    enddo
  enddo
  close (unit = 1)  
  call write_outfile(preeqfile,flagoutall)
  write(*, '(/," Total pre-equilibrium cross section:", f12.5)') xspreeqsum
  return
end subroutine preeqout
! Copyright A.J. Koning 2021
