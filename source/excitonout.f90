subroutine excitonout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of exciton model parameters
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
! Variables for main input
!   Ainit         ! mass number of initial compound nucleus
! Constants
!   hbar          ! Planck's constant / 2.pi in MeV.s
!   parname       ! name of particle
! Variables for exciton model initialization
!   Qfactor       ! Q - factor for neutron / proton distinction
! Variables for preequilibrium
!   Ecomp         ! total energy of composite system
!   p0            ! initial particle number
! Variables for preequilibrium initialization
!   maxpar        ! maximal particle number
! Variables for exciton model
!   depletion     ! depletion factor at each stage
!   M2            ! square of matrix element
!   tauexc        ! lifetime of exciton state
!   wemispart     ! emission rate per particle and exciton number
!   wemistot      ! total emission rate per exciton number
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
  integer           :: h          ! help variable
  integer           :: n          ! exciton number
  integer           :: p          ! particle number
  integer           :: type       ! particle type
  integer            :: Ncol
  integer            :: Nk
  real(sgl)         :: lambdaplus ! transition rate for n --> n+2
!
! ************************ Exciton model *******************************
!
  write(*, '(/" ++++++++++ EXCITON MODEL ++++++++++",/)')
!
! 1. Output of matrix element
!
! matrix    : subroutine for matrix element for exciton model
!
  quantity='one-component exciton model'
  reaction= '('//parsym(k0)//',x)'
  topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
  excfile='exciton0000.000.tot'//natstring(iso)
  write(excfile(8:15), '(f8.3)') Einc
  write(excfile(8:11), '(i4.4)') int(Einc)
  open (unit = 1, file = excfile, status = 'replace')
  call write_header(topline,source,user,date,oformat)
  call write_target
  write(1,'("# parameters:")')
  call write_real(2,'E-incident [MeV]',Einc)
  call write_real(2,'E-compound [MeV]',Ecomp)
  call write_real(2,'Surface effective well depth [MeV]',Esurf)
  call write_real(2,'Constant for matrix element',M2constant)
  un = ''
  col(1) = 'p'
  col(2) = 'h'
  col(3) = 'M2'     
  Ncol = 3
  Nk = maxpar - p0 + 1
  call write_quantity(quantity)
  call write_datablock(Ncol,Nk,col,un)
  do p = p0, maxpar
    h = p - p0
    n = p + h
    call matrix(Ainit, n)
    write(1, '(2(3x, i6, 6x), es15.6)') p, h, M2
  enddo
!
! 2. Output of Q-factors
!
  quantity='Q-factors'
  do type = 1,6
    col(2 + type) = parname(type)
  enddo
  Ncol = 8
  call write_quantity(quantity)
  call write_datablock(Ncol,Nk,col,un)
  do p = p0, maxpar
    h = p - p0
    write(1, '(2(3x, i6, 6x), 6es15.6)') p, h, (Qfactor(type, p), type = 1, 6)
  enddo
!
! 3. Output of emission rates or escape widths
!
  quantity='emission rates'
  do type = 0,6
    col(3 + type) = parname(type)
    un(3+type) = 'sec^-1'
  enddo
  col(10) = 'Total'
  un(10) = 'sec^-1'
  Ncol = 10
  call write_quantity(quantity)
  call write_datablock(Ncol,Nk,col,un)
  do p = p0, maxpar
    h = p - p0
    write(1, '(2(3x, i6, 6x), 8es15.6)') p, h, (wemispart(type, p, h), type = 0, 6), wemistot(p, h)
  enddo
  quantity='escape widths'
  do type = 0,6
    un(3+type) = 'MeV'
  enddo
  col(10) = 'Total'
  un(10) = 'MeV'
  call write_quantity(quantity)
  call write_datablock(Ncol,Nk,col,un)
  do p = p0, maxpar
    h = p - p0
    write(1, '(2(3x, i6, 6x), 8es15.6)') p, h, (wemispart(type, p, h) * hbar, type = 0, 6), wemistot(p, h) * hbar
  enddo
!
! 4. Output of transition rates or damping widths and total widths
!
  quantity='internal transition rates'
  col(3) = 'lambdaplus'
  un(3) = 'sec^-1'
  Ncol = 3
  call write_quantity(quantity)
  call write_datablock(Ncol,Nk,col,un)
  do p = p0, maxpar
    h = p - p0
    write(1, '(2(3x, i6, 6x), es15.6)') p, h, lambdaplus(0, 0, p, h)
  enddo
  quantity='damping widths'
  col(3) = 'gammaplus'
  un(3) = 'MeV'
  call write_quantity(quantity)
  call write_datablock(Ncol,Nk,col,un)
  do p = p0, maxpar
    h = p - p0
    write(1, '(2(3x, i6, 6x), es15.6)') p, h, lambdaplus(0, 0, p, h)*hbar
  enddo
  quantity='total widths'
  call write_quantity(quantity)
  call write_datablock(Ncol,Nk,col,un)
  do p = p0, maxpar
    h = p - p0
    write(1, '(2(3x, i6, 6x), es15.6)') p, h, (lambdaplus(0, 0, p, h) + wemistot(p, h)) * hbar
  enddo
!
! 5. Output of depletion factors
!
  quantity='depletion factors'
  col(3) = 'depletion'
  un(3) = ''
  call write_quantity(quantity)
  call write_datablock(Ncol,Nk,col,un)
  do p = p0, maxpar
    h = p - p0
    write(1, '(2(3x, i6, 6x), es15.6)') p, h, depletion(p, h)
  enddo
!
! 6. Output of lifetimes of exciton states
!
  quantity='lifetimes'
  col(3) = 'mean lifetime'
  un(3) = 'sec'   
  call write_quantity(quantity)
  call write_datablock(Ncol,Nk,col,un)
  do p = p0, maxpar
    h = p - p0
    write(1, '(2(3x, i6, 6x), es15.6)') p, h, tauexc(p, h)
  enddo
  close(1)
  call write_outfile(excfile,flagoutall)
  return
end subroutine excitonout
! Copyright A.J. Koning 2021
