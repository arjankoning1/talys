subroutine productionout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of particle production cross sections
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
! Variables for output
!   filefission     ! flag for fission cross sections on separate file
!   filetotal       ! flag for total cross sections on separate file
! Variables for fission
!   flagfission     ! flag for fission
! Variables for input energies
!   eninc           ! incident energy in MeV
!   nin             ! counter for incident energy
!   Ninc          ! number of incident energies
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
! Variables for energies
!   Ninclow       ! number of incident energies below Elow
! Variables for total cross sections
!   xsfistot        ! total fission cross section
! Variables for incident channel
!   multiplicity    ! particle multiplicity
!   xsparticle      ! total particle production cross section
! Variables for energy grid
!   Einc            ! incident energy in MeV
! Variables for nuclides
!   parskip         ! logical to skip outgoing particle
!   Q               ! Q - value
! Constants
!   iso             ! counter for isotope
!   natstring       ! string extension for file names
!   parname         ! name of particle
!   parsym          ! symbol of particle
!
! *** Declaration of local data
!
  implicit none
  character(len=13) :: totfile    ! file with total cross sections
  character(len=15) :: fisfile    ! fission file
  character(len=18) :: reaction   ! reaction
  character(len=15) :: col(3)    ! header
  character(len=15) :: un(3)    ! header
  character(len=80) :: quantity   ! quantity
  character(len=132) :: topline    ! topline
  integer           :: MF
  integer           :: MT
  integer           :: Ncol       ! counter
  integer           :: nen        ! energy counter
  integer           :: type       ! particle type
!
! ************** Total particle production cross sections **************
!
  MF = 3
  quantity='cross section'
  col(1)='E'
  un(1)='MeV'
  col(2)='xs'
  un(2)='mb'
  col(3)='Multiplicity'
  un(3)=''
  Ncol=3
  write(*, '(/" 3. Total particle production cross sections"/)')
  do type = 0, 6
    if (parskip(type)) cycle
    write(*, '(1x, a8, "=", es12.5, "    Multiplicity=", es12.5)') parname(type), xsparticle(type), multiplicity(type)
!
! Write results to separate file
!
    if (filetotal) then
      totfile = ' prod.tot'//natstring(iso)
      write(totfile(1:1), '(a1)') parsym(type)
      reaction='('//parsym(k0)//',x'//parsym(type)//')'
      if (nin == Ninclow + 1) then
        open (unit = 1, file = totfile, status = 'replace')
        topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
        if (type == 0) MT = 202
        if (type == 1) MT = 201
        if (type > 1) MT = 201 + type
        call write_header(topline,source,user,date,oformat)
        call write_target
        call write_reaction(reaction,0.d0,0.d0,MF,MT)
        call write_quantity(quantity)
        call write_datablock(Ncol,Ninc,col,un)
        do nen = 1, Ninclow
          write(1, '(3es15.6)') eninc(nen), 0., 0.
        enddo
      else
        open (unit = 1, file = totfile, status = 'old', position = 'append')
      endif
      write(1, '(3es15.6)') Einc, xsparticle(type), multiplicity(type)
      close (unit = 1)
    endif
  enddo
!
! Total fission cross section
!
  if (flagfission) then
    write(*, '(" fission =", es12.5)') xsfistot
!
! Write results to separate file
!
    if (filefission) then
      fisfile = 'fission.tot'//natstring(iso)
      reaction='('//parsym(k0)//',f)'
      Ncol=2
      if (nin == Ninclow + 1) then
        open (unit = 1, file = fisfile, status = 'replace')
        topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
        MT = 18
        call write_header(topline,source,user,date,oformat)
        call write_target
        call write_reaction(reaction,0.d0,0.d0,MF,MT)
        call write_quantity(quantity)
        call write_datablock(Ncol,Ninc,col,un)
        do nen = 1, Ninclow
          write(1, '(2es15.6)') eninc(nen), 0.
        enddo
      else
        open (unit = 1, file = fisfile, status = 'old', position = 'append')
      endif
      write(1, '(2es15.6)') Einc, xsfistot
      close (unit = 1)
    endif
  endif
  return
end subroutine productionout
! Copyright A.J. Koning 2021
