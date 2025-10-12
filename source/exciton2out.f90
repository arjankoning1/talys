subroutine exciton2out
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of two-component exciton model parameters
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
!   sgl           ! single precision kind
! Variables for preequilibrium
!   M2constant    ! constant for matrix element in exciton model
!   Rnunu         ! ratio for two - component matrix element
!   Rnupi         ! ratio for two - component matrix element
!   Rpinu         ! ratio for two - component matrix element
!   Rpipi         ! ratio for two - component matrix element
! Variables for main input
!   Ainit         ! mass number of initial compound nucleus
! Constants
!   hbar          ! Planck's constant / 2.pi in MeV.s
!   parname       ! name of particle
! Variables for preequilibrium initialization
!   maxpar        ! maximal particle number
! Variables for exciton model
!   M2nunu        ! square of neutron - neutron matrix element
!   M2nupi        ! square of neutron - proton matrix element
!   M2pinu        ! square of proton - neutron matrix element
!   M2pipi        ! square of proton - proton matrix element
!   wemispart2    ! two - component emission rate per par
!   wemistot2     ! total two - component emission rate p
! Variables for preequilibrium
!   Ecomp         ! total energy of composite system
!   p0            ! initial particle number
!   pnu0          ! initial neutron number
!   ppi0          ! initial proton number
!   Spre          ! time - integrated strength of two - component exciton state
!
! *** Declaration of local data
!
  implicit none
  character(len=30) :: excfile
  character(len=18) :: reaction   ! reaction
  character(len=15) :: col(12)    ! header
  character(len=15) :: un(12)    ! units
  character(len=80) :: quantity   ! quantity
  character(len=132) :: topline    ! topline
  integer            :: h            ! help variable
  integer            :: hnu          ! neutron hole number
  integer            :: hpi          ! proton hole number
  integer            :: n            ! exciton number
  integer            :: p            ! particle number
  integer            :: pnu          ! neutron particle number
  integer            :: ppi          ! proton particle number
  integer            :: type         ! particle type
  integer            :: Ncol
  integer            :: Nk
  integer            :: indent
  integer            :: id2
  integer            :: id4
  real(sgl)          :: lambdanupi   ! neutron-proton transition rate for n --> n
  real(sgl)          :: lambdanuplus ! neutron transition rate for n --> n+2
  real(sgl)          :: lambdapinu   ! proton-neutron transition rate for n --> n
  real(sgl)          :: lambdapiplus ! proton transition rate for n --> n+2
!
! ************************ Exciton model *******************************
!
  indent = 0
  id2 = indent + 2
  id4 = indent + 4
  write(*, '(/" ++++++++++ TWO-COMPONENT EXCITON MODEL ++++++++++",/)')
!
! 1. Output of matrix element
!
! matrix    : subroutine for matrix element for exciton model
!
  quantity='two-component exciton model'
  reaction= '('//parsym(k0)//',x)'
  topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
  if (flagblockpreeq) then
    excfile='exciton.out'//natstring(iso)
    if (nin == Ninclow + 1) then
      open (unit = 1, file = excfile, status = 'unknown')
    else
      open (unit = 1, file = excfile, status = 'unknown', position = 'append')
    endif
  else
    excfile='exciton0000.000.out'//natstring(iso)
    write(excfile(8:15), '(f8.3)') Einc
    write(excfile(8:11), '(i4.4)') int(Einc)
    open (unit = 1, file = excfile, status = 'replace')
  endif
  call write_header(indent,topline,source,user,date,oformat)
  call write_target(indent)
  quantity='matrix element'
  call write_char(indent,'parameters','')
  call write_real(id2,'E-incident [MeV]',Einc)
  call write_real(id2,'E-compound [MeV]',Ecomp)
  call write_real(id2,'Surface effective well depth [MeV]',Esurf)
  call write_real(id2,'Constant for matrix element',M2constant)
  call write_real(id2,'p-p ratio for matrix element',Rpipi)
  call write_real(id2,'n-n ratio for matrix element',Rnunu)
  call write_real(id2,'p-n ratio for matrix element',Rpinu)
  call write_real(id2,'n-p ratio for matrix element',Rnupi)
  un = ''
  col(1) = 'p(p)'
  col(2) = 'h(p)'
  col(3) = 'p(n)'
  col(4) = 'h(n)'
  col(5) = 'M2pipi'
  col(6) = 'M2nunu'
  col(7) = 'M2pinu'
  col(8) = 'M2nupi'
  Ncol = 8
  Nk = 0
  do p = p0, maxpar
    do ppi = ppi0, maxpar
      hpi = ppi - ppi0
      do pnu = pnu0, maxpar
        if (ppi + pnu == p) Nk = Nk + 1
      enddo
    enddo
  enddo
  call write_quantity(id2,quantity)
  call write_datablock(id2,Ncol,Nk,col,un)
  do p = p0, maxpar
    do ppi = ppi0, maxpar
      hpi = ppi - ppi0
      do pnu = pnu0, maxpar
        if (ppi + pnu == p) then
          hnu = pnu - pnu0
          h = hpi + hnu
          n = p + h
          call matrix(Ainit, n)
          write(1, '(4(3x,i6,6x), 4es15.6)') ppi, hpi, pnu, hnu, M2pipi, M2nunu, M2pinu, M2nupi
        endif
      enddo
    enddo
  enddo
!
! 2. Output of emission rates or escape widths
!
  quantity='emission rate'
  do type = 0,6
    col(5 + type) = parname(type)
    un(5+type) = 'sec^-1'
  enddo
  col(12) = 'Total'
  un(12) = 'sec^-1'
  Ncol = 12
  call write_quantity(id2,quantity)
  call write_datablock(id2,Ncol,Nk,col,un)
  do p = p0, maxpar
    do ppi = ppi0, maxpar
      hpi = ppi - ppi0
      do pnu = pnu0, maxpar
        if (ppi + pnu == p) then
          hnu = pnu - pnu0
          write(1, '(4(3x, i6, 6x), 8es15.6)') ppi, hpi, pnu, hnu, (wemispart2(type, ppi, hpi, pnu, hnu), type = 0, 6), &
 &          wemistot2(ppi, hpi, pnu, hnu)
        endif
      enddo
    enddo
  enddo
  quantity='escape width'
  do type = 0,6
    un(5+type) = 'MeV'
  enddo
  col(12) = 'Total'
  un(12) = 'MeV'
  call write_quantity(id2,quantity)
  call write_datablock(id2,Ncol,Nk,col,un)
  do p = p0, maxpar
    do ppi = ppi0, maxpar
      hpi = ppi - ppi0
      do pnu = pnu0, maxpar
        if (ppi + pnu == p) then
          hnu = pnu - pnu0
          write(1, '(4(3x, i6, 6x), 8es15.6)') ppi, hpi, pnu, hnu, &
 &          (wemispart2(type, ppi, hpi, pnu, hnu) * hbar, type = 0, 6), wemistot2(ppi, hpi, pnu, hnu) * hbar
        endif
      enddo
    enddo
  enddo
!
! 3. Output of transition rates or damping widths and total widths
!
  quantity='internal transition rate'
  col(5) = 'lambdapiplus'
  col(6) = 'lambdanuplus'
  col(7) = 'lambdapinu'
  col(8) = 'lambdanupi'
  do type = 1,4
    un(4+type) = 'sec^-1'
  enddo
  Ncol = 8
  call write_quantity(id2,quantity)
  call write_datablock(id2,Ncol,Nk,col,un)
  do p = p0, maxpar
    do ppi = ppi0, maxpar
      hpi = ppi - ppi0
      do pnu = pnu0, maxpar
        if (ppi + pnu == p) then
          hnu = pnu - pnu0
          write(1, '(4(3x, i6, 6x), 4es15.6)') ppi, hpi, pnu, hnu, lambdapiplus(0, 0, ppi, hpi, pnu, hnu), &
 &          lambdanuplus(0, 0, ppi, hpi, pnu, hnu), lambdapinu(0, 0, ppi, hpi, pnu, hnu), lambdanupi(0, 0, ppi, hpi, pnu, hnu)
        endif
      enddo
    enddo
  enddo
  quantity='damping width'
  col(5) = 'gammapiplus'
  col(6) = 'gammanuplus'
  col(7) = 'gammapinu'
  col(8) = 'gammanupi'
  do type = 1,4
    un(4+type) = 'MeV'
  enddo
  call write_quantity(id2,quantity)
  call write_datablock(id2,Ncol,Nk,col,un)
  do p = p0, maxpar
    do ppi = ppi0, maxpar
      hpi = ppi - ppi0
      do pnu = pnu0, maxpar
        if (ppi + pnu == p) then
          hnu = pnu - pnu0
          write(1, '(4(3x, i6, 6x), 4es15.6)') ppi, hpi, pnu, hnu, lambdapiplus(0, 0, ppi, hpi, pnu, hnu) * hbar, &
 &          lambdanuplus(0, 0, ppi, hpi, pnu, hnu) * hbar, lambdapinu(0, 0, ppi, hpi, pnu, hnu) * hbar, &
 &          lambdanupi(0, 0, ppi, hpi, pnu, hnu) * hbar
        endif
      enddo
    enddo
  enddo
  quantity='total width'
  col(5) = 'gammatot'
  Ncol = 5
  call write_quantity(id2,quantity)
  call write_datablock(id2,Ncol,Nk,col,un)
  do p = p0, maxpar
    do ppi = ppi0, maxpar
      hpi = ppi - ppi0
      do pnu = pnu0, maxpar
        if (ppi + pnu == p) then
          hnu = pnu - pnu0
          write(1, '(4(3x, i6, 6x), es15.6)') ppi, hpi, pnu, hnu, hbar * (lambdapiplus(0, 0, ppi, hpi, pnu, hnu) + &
 &          lambdanuplus(0, 0, ppi, hpi, pnu, hnu) + lambdapinu(0, 0, ppi, hpi, pnu, hnu) + &
 &          lambdanupi(0, 0, ppi, hpi, pnu, hnu) + wemistot2(ppi, hpi, pnu, hnu))
        endif
      enddo
    enddo
  enddo
!
! 4. Output of lifetimes of exciton states
!
  quantity='lifetime'
  col(5) = 'Strength'
  un(5) = 'sec'
  call write_quantity(id2,quantity)
  call write_datablock(id2,Ncol,Nk,col,un)
  do p = p0, maxpar
    do ppi = ppi0, maxpar
      hpi = ppi - ppi0
      do pnu = pnu0, maxpar
        if (ppi + pnu == p) then
          hnu = pnu - pnu0
          write(1, '(4(3x, i6, 6x), es15.6)') ppi, hpi, pnu, hnu, Spre(ppi, hpi, pnu, hnu)
        endif
      enddo
    enddo
  enddo
  close(1)
  call write_outfile(excfile,flagoutall)
  return
end subroutine exciton2out
! Copyright A.J. Koning 2021
