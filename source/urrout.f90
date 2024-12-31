subroutine urrout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of unresolved resonance parameters in separate files
!
! Author    : Gilles Noguere and Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl            ! single precision kind
! All global variables
!   numJ           ! maximum J - value
!   numl           ! number of l values
! Variables for energies
!   Ninclow      ! number of incident energies below Elow
! Variables for compound reactions
!   flagurrnjoy    ! normalization of URR parameters with NJOY method
!   lurr           ! maximal orbital angular momentum for URR
!   xscaptherm     ! thermal capture cross section
! Variables for input energies
!   eninc          ! incident energy in MeV
!   nin            ! counter for incident energy
!   Ninc           ! number of incident energies
! Variables for main input
!   Atarget        ! mass number of target nucleus
!   k0             ! index of incident particle
!   Ltarget        ! excited level of target
!   Ztarget        ! charge number of target nucleus
! Variables for OMP
!   Rprime         ! potential scattering radius
! Variables for OMP
!   RprimeU        ! potential scattering radius
! Variables for energy grid
!   Einc           ! incident energy in MeV
! Variables for compound nucleus from target
!   JmaxU          ! maximal total angular momentum
!   JminU          ! minimal total angular momentum
!   lmaxU          ! maximal orbital angular momentum
!   lminU          ! minimal orbital angular momentum
!   nulj           ! (l, j) number of degrees of freedom for URR calculation
!   Purrlj         ! (l, j) parity for URR calculation
! Variables for nuclides
!   Q              ! Q - value
! Constants
!   cparity        ! parity (character)
!   parsym         ! symbol of particle
! Variables for levels
!   jdis           ! spin of level
! Variables for resonance parameters
!   Dl             ! mean resonance spacing per l value
!   Dlj            ! mean resonance spacing per J, l value
! Variables for masses
!   S              ! separation energy
! Variables for URR
!   flagurrendf    ! flag for URR info to ENDF
!   strengthl      ! l neutron strength function
!   strengthlj     ! (l, j) neutron strength function
!   urrwidth       ! channel width in URR
!   xsurrN         ! URR cross section
!   xsurrT         ! URR cross section
! Variables for existence libraries
!   urrexist       ! flag for existence of URR
!
! *** Declaration of local data
!
  implicit none
  character(len=3)  :: lstring              ! string for l value
  character(len=13) :: Estr
  character(len=20) :: urrfile              ! file with URR parameters
  character(len=250) :: urrline(1000)
  character(len=18) :: reaction   ! reaction
  character(len=15) :: col(100)    ! header
  character(len=15) :: un(100)    ! header
  character(len=80) :: quantity   ! quantity
  character(len=80) :: extension
  character(len=132) :: topline    ! topline
  integer           :: J                    ! spin of level
  integer           :: istat
  integer           :: i
  integer           :: k
  integer           :: Ne
  integer           :: Nurr
  integer           :: Ncol
  integer           :: l                    ! multipolarity
  integer           :: nen                  ! energy counter
  integer           :: odd                  ! odd (1) or even (0) nucleus
  integer           :: type                 ! particle type
  real(sgl)         :: x(0:numl, 0:numJ)    ! help variable
  real(sgl)         :: xx(4)                ! x value
  real(sgl)         :: y(0:numl, 0:numJ)    ! coordinates of the point to test
!
! General output file
!
  urrline = ''
  reaction = '('//parsym(k0)//',urr)'
  quantity='URR parameter'
  un = ''
  col = ''
  col(1)='l'
  col(2)='J'
  col(3)='P'
  col(4)='D(l)'
  un(4)='eV'
  col(5)='D(l,J)'
  un(5)='eV'
  col(6)='S(l)'
  col(7)='S(l,J)'
  col(8)='Gx(l,J)'
  un(8)='eV'
  col(9)='Gn(l,J)'
  un(9)='eV'
  col(10)='Gg(l,J)'
  un(10)='eV'
  col(11)='Gf(l,J)'
  un(11)='eV'
  col(12)='Gp(l,J)'
  un(12)='eV'
  col(13)='Ga(l,J)'
  un(13)='eV'
  Ncol=13
  if (nin == Ninclow+1) then
    flagurrendf = .true.
    open (unit = 21, file = 'urr.dat', status = 'unknown')
  endif
  Estr=''
  write(Estr,'(es13.6)') Einc
  topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
  open (unit = 1, file = 'urr.tmp', status = 'unknown')
  call write_header(topline,source,user,date,oformat)
  call write_target
  call write_reaction(reaction,0.d0,0.d0,0,0)
  call write_real(2,'E-incident [MeV]',Einc)
  write(1,'("# parameters:")')
  call write_double(2,'neutron separation energy [MeV]',S(0, 0, 1))
  write(1,'("# observables:")')
  call write_real(2,'thermal capture cross section [mb]',xscaptherm(-1))
  call write_real(2,'potential scattering radius [fm]',Rprime)
  Ne = 0
  do l = 0, lurr
    do J = JminU(l), JmaxU(l)
      Ne = Ne + 1
    enddo
  enddo
  call write_datablock(quantity,Ncol,Ne,col,un)
  odd = mod(Atarget + 1, 2)
  do l = 0, lurr
    do J = JminU(l), JmaxU(l)
      write(1, '(3x,i6,12x,f5.1,9x,i6,4x,10es15.6)') l, J+0.5*odd, Purrlj(l, J), Dl(l), Dlj(l, J), strengthl(l), &
 &      strengthlj(l, J), urrwidth(1, l, J), urrwidth(3, l, J), &
 &      urrwidth(0, l, J), urrwidth( - 1, l, J), urrwidth(2, l, J), urrwidth(6, l, J)
    enddo
  enddo
  close (unit = 1)
  open (unit = 1, file = 'urr.tmp', status = 'unknown')
  k = 0
  do
    k = k + 1
    read(1,'(a)', iostat = istat) urrline(k)
    if (istat == -1) exit
  enddo
  Nurr = k - 1
  close (unit = 1, status = 'delete')
  do k = 1, Nurr
    write(21,'(a)') trim(urrline(k))
  enddo
  if (nin == Ninc) close (unit = 21)
!
! Output of (l,J) dependent widths, spacings and strength functions in separate files
!
  do type = -1, 6
    if (type == -1 .and. .not.flagfission) cycle
    if (type == 5 .and. Q(2) <= 0.) cycle
    if (type == 6 .and. Q(6) <= 0.) cycle
    do l = 0, lurr
      do J = 0, numJ
        x(l, J) = 0.
        y(l, J) = urrwidth(type, l, J)
        if (type ==  - 1 .or. type == 1) x(l, J) = nulj(type, l, J)
        if (type == 2) y(l, J) = Dlj(l, J)
        if (type == 3) x(l, J) = nulj(0, l, J)
        if (type == 4) y(l, J) = strengthlj(l, J)
        if (type == 5) then
          x(l, J) = nulj(2, l, J)
          y(l, J) = urrwidth(2, l, J)
        endif
        if (type == 6) x(l, J) = nulj(6, l, J)
      enddo
    enddo
    do l = lminU, min(lmaxU, lurr)
      urrfile = '                    '
      lstring = 'l00'
      write(lstring(2:3), '(i2.2)') l
      extension =''
      if (type == -1) then
        urrfile = 'urrfiswidth.'//lstring
        extension='average fission width'
      endif
      if (type == 0) then
        urrfile = 'urrgamwidth.'//lstring
        extension='average radiation width'
      endif
      if (type == 1) then
        urrfile = 'urrcomwidth.'//lstring
        extension='average competitive width'
      endif
      if (type == 2) then
        urrfile = 'urrspacinglj.'//lstring
        extension='mean level spacing per l,J'
      endif
      if (type == 3) then
        urrfile = 'urrneuwidth.'//lstring
        extension='reduced neutron width'
      endif
      if (type == 4) then
        urrfile = 'urrneustrengthlj.'//lstring
        extension='neutron strength function'
      endif
      if (type == 5) then
        urrfile = 'urrprowidth.'//lstring
        extension='average proton width'
      endif
      if (type == 6) then
        urrfile = 'urralpwidth.'//lstring
        extension='average alpha width'
      endif
      k = 1 + JmaxU(l) - JminU(l)
      un = ''
      col(1)='E'
      un(1)='MeV'
      do i=1,k
        col(3*(i-1)+2)='J'
        col(3*(i-1)+3)='nu'
        col(3*(i-1)+4)='width'
        if (type == 2) col(3*(i-1)+4)='spacing'
        if (type == 4) col(3*(i-1)+4)='strength'
      enddo
      Ncol=1+3*k
      if ( .not. urrexist(type, l)) then
        urrexist(type, l) = .true.
        topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' - '//extension
        open (unit = 1, file = urrfile, status = 'replace')
        call write_header(topline,source,user,date,oformat)
        call write_target
        call write_reaction(reaction,0.d0,0.d0,0,0)
        call write_integer(2,'l-value',l)
        call write_datablock(quantity,Ncol,Ninc,col,un)
        do nen = 1, Ninclow
          write(1, '(100es15.6)') eninc(nen), (J + 0.5 * odd, x(l, J), y(l, J), J = JminU(l), JmaxU(l))
        enddo
        do nen = Ninclow + 1, nin - 1
          write(1, '(es15.6, " !!! not calculated")') eninc(nen)
        enddo
      else
        open (unit = 1, file = urrfile, status = 'old', position = 'append')
      endif
      write(1, '(100es15.6)') Einc, (J + 0.5 * odd, x(l, J), y(l, J), J = JminU(l), JmaxU(l))
      close (unit = 1)
    enddo
  enddo
!
! Output of l-dependent mean level spacing and neutron strength function in separate file
!
  do type = 7, 8
    do l = lminU, min(lmaxU, lurr)
      urrfile = '                    '
      lstring = 'l00'
      write(lstring(2:3), '(i2.2)') l
      col(1)='E'
      un(1)='MeV'
      un(2)=''
      Ncol=2
      if (type == 7) then
        urrfile = 'urrspacingl.'//lstring
        xx(1) = Dl(l)
        extension='mean level spacing per l'
        col(2)='spacing'
      else
        urrfile = 'urrneustrengthl.'//lstring
        xx(1) = strengthl(l)
        extension='neutron strength function per l'
        col(2)='strength'
      endif
      if ( .not. urrexist(type, l)) then
        urrexist(type, l) = .true.
        topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' - '//extension
        open (unit = 1, file = urrfile, status = 'replace')
        call write_header(topline,source,user,date,oformat)
        call write_target
        call write_reaction(reaction,0.d0,0.d0,0,0)
        call write_integer(2,'l-value',l)
        call write_datablock(quantity,Ncol,Ninc,col,un)
        do nen = 1, Ninclow
          write(1, '(2es15.6)') eninc(nen), xx(1)
        enddo
        do nen = Ninclow + 1, nin - 1
          write(1, '(es15.6, " !!! not calculated")') eninc(nen)
        enddo
      else
        open (unit = 1, file = urrfile, status = 'old', position = 'append')
      endif
      write(1, '(2es15.6)') Einc, xx(1)
      close (unit = 1)
    enddo
  enddo
!
! Output of URR cross section from NJOY method and TALYS
!
  un = 'mb'
  col(1)='E'
  un(1)='MeV'
  col(2)='xs'
  col(3)='elastic'
  col(4)='capture'
  col(5)='fission'
  Ncol=5
  do type = 9, 11
    if (type /= 10 .and. .not. flagurrnjoy) cycle
    quantity='cross section'
    if (type == 9) then
      topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' with formalism from NJOY'
      urrfile = 'urrnjoy.tot         '
      do j = 1, 4
        xx(j) = xsurrN(j)
      enddo
    endif
    if (type == 10) then
      topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' with formalism from TALYS'
      urrfile = 'urrtalys.tot        '
      do j = 1, 4
        xx(j) = xsurrT(j)
      enddo
    endif
    if (type == 11) then
      quantity='cross section ratio'
      topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' TALYS:NJOY'
      urrfile = 'urrratio.tot        '
      do j = 1, 4
        if (xsurrN(j) > 0.) then
          xx(j) = xsurrT(j) / xsurrN(j)
        else
          xx(j) = 1.
        endif
      enddo
    endif
    if ( .not. urrexist(type, 0)) then
      urrexist(type, 0) = .true.
      open (unit = 1, file = urrfile, status = 'replace')
      call write_header(topline,source,user,date,oformat)
      call write_target
      call write_reaction(reaction,0.d0,0.d0,0,0)
      write(1,'("# observables:")')
      if (type == 9) call write_real(2,'potential scattering radius [fm]',RprimeU)
      if (type == 10) call write_real(2,'potential scattering radius [fm]',Rprime)
      if (type == 11) call write_real(2,'ratio of potential scattering radius ',Rprime/RprimeU)
      call write_datablock(quantity,Ncol,Ninc,col,un)
      do nen = 1, Ninclow
        write(1, '(5es15.6)') eninc(nen), xx(1), xx(2), xx(4), xx(3)
      enddo
      do nen = Ninclow + 1, nin - 1
        write(1, '(es15.6, " !!! not calculated")') eninc(nen)
      enddo
    else
      open (unit = 1, file = urrfile, status = 'old', position = 'append')
    endif
    write(1, '(5es15.6)') Einc, xx(1), xx(2), xx(4), xx(3)
    close (unit = 1)
  enddo
  return
end subroutine urrout
! Copyright A.J. Koning 2023
