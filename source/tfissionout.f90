subroutine tfissionout(Zcomp, Ncomp, nex)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of fission transmission coefficients
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
! Variables for compound nucleus from target
!   Exinc    ! excitation energy of entrance bin
! Variables for excitation energy grid
!   maxJ     ! maximal J - value
! Variables for fission transmission coefficients
!   denfis    ! fission level density
!   gamfis    ! fission width
!   taufis    ! fission lifetime
!   tfis      ! fission transmission coefficient for Hill - Wheeler magnitude
! Variables for nuclides
!   AA       ! mass number of residual nucleus
!   NN       ! neutron number of residual nucleus
!   ZZ       ! charge number of residual nucleus
! Constants
!   nuc      ! symbol of nucleus
!
! *** Declaration of local data
!
  implicit none
  character(len=3)  :: massstring
  character(len=6)  :: finalnuclide ! 
  character(len=80)  :: fisfile               ! fission file
  character(len=132) :: topline    ! topline
  character(len=15) :: col(9)     ! header
  character(len=15) :: un(9)     ! header
  character(len=80) :: quantity   ! quantity
  integer :: A     ! mass number of target nucleus
  integer :: J     ! spin of level
  integer :: J2    ! 2 * J
  integer :: Ncol  ! 
  integer :: N     ! neutron number of residual nucleus
  integer :: Ncomp ! neutron number index for compound nucleus
  integer :: nex   ! excitation energy bin of compound nucleus
  integer :: odd   ! odd (1) or even (0) nucleus
  integer :: Z     ! charge number of target nucleus
  integer :: Zcomp ! proton number index for compound nucleus
  integer :: indent
  integer :: id2
  integer :: id4
!
! *************** Output of fission transmission coefficients **********
!
  indent = 0
  id2 = indent + 2
  id4 = indent + 4
  Z = ZZ(Zcomp, Ncomp, 0)
  N = NN(Zcomp, Ncomp, 0)
  A = AA(Zcomp, Ncomp, 0)
  fisfile = 'fis000000.trans'
  write(fisfile(4:6), '(i3.3)') Z
  write(fisfile(7:9), '(i3.3)') A
  massstring='   '
  write(massstring,'(i3)') A            
  finalnuclide=trim(nuc(Z))//adjustl(massstring)
  quantity='transmission coefficients'
  if (tfisexist(Zcomp, Ncomp)) then
    open (unit = 1, file = fisfile, status = 'old', position = 'append')
  else
    open (unit = 1, file = fisfile, status = 'replace')
    tfisexist(Zcomp, Ncomp) = .true.
    topline=trim(finalnuclide)//' fission '//trim(quantity)
    call write_header(indent,topline,source,user,date,oformat)
    call write_residual(indent,Z,A,finalnuclide) 
  endif
  un = ''
  un(4)= 'eV'
  un(5)= 'eV'
  un(6)= 'sec'
  un(7)= 'sec'
  un(8)= 'MeV^-1'
  un(9)= 'MeV^-1'
  col(1)= 'J'
  col(2)= 'T(J,-)'
  col(3)= 'T(J,+)'
  col(4)= 'Gamma(J,-)'
  col(5)= 'Gamma(J,+)'
  col(6)= 'tau(J,-)'
  col(7)= 'tau(J,+)'
  col(8)= 'density(J,-)'
  col(9)= 'density(J,+)'
  Ncol = 9
  call write_quantity(id2,quantity)
  call write_real(id4,'Excitation energy [MeV]',Exinc)
  call write_datablock(id2,Ncol,maxJ(Zcomp, Ncomp, nex) + 1,col,un)
  write(*, '(/" Fission transmission coefficients for Z=", i3, &
 &  " N=", i3, " (",a,") and an excitation energy of ", f8.3, " MeV"/)') Z, N, trim(finalnuclide), Exinc
  odd = mod(A, 2)
  do J = 0, maxJ(Zcomp, Ncomp, nex)
    J2 = 2 * J + odd
    write(1, '(2x, f4.1, 9x, 8es15.6)') 0.5*J2, tfis(J, -1), tfis(J, 1), gamfis(J, -1), gamfis(J, 1), taufis(J, -1), &
 &  taufis(J, 1), denfis(J, -1), denfis(J, 1)
  enddo
  close (1)
  call write_outfile(fisfile,flagoutall)
  return
end subroutine tfissionout
! Copyright A.J. Koning 2021
