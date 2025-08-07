subroutine pfnsout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of prompt fission neutron spectrum
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
! Variables for nuclides
!   parskip         ! logical to skip outgoing particle
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
! Variables for fission
!   fymodel         ! fission yield model, 1: Brosa 2: GEF
! Variables for energies
!   Einc0           ! incident energy in MeV
! Constants
!   iso             ! counter for isotope
!   natstring       ! string extension for file names
!   parname         ! name of particle
!   parsym          ! symbol of particle
! Variables for mass distribution
!   Eavpfns         ! average energy of prompt fission neutrons spectrum
!   Epfns           ! energy of PFNS
!   Epfnsaverage    ! average energy of PFN
!   maxpfns         ! maximum energy of prompt fission neutrons spectrum
!   NEpfns          ! number of energies for PFNS
!   pfns            ! prompt fission neutron spectrum
!   pfnscm          ! prompt fission neutron spectrum in CM
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: pfnsfile  ! file for nubar
  character(len=12) :: Estr
  character(len=18) :: reaction   ! reaction
  character(len=132) :: topline    ! topline
  character(len=15) :: col(4)     ! header
  character(len=15) :: un(4)     ! header
  character(len=80) :: quantity   ! quantity
  integer           :: MF
  integer           :: MT
  integer           :: indent
  integer           :: id2
  integer           :: id4
  integer           :: Ncol      ! number of columns
  integer           :: nen       ! energy counter
  integer           :: type      ! particle type
!
! Write PFNS, PFGS, etc.
!
  indent = 0
  id2 = indent + 2
  id4 = indent + 4
  MF = 5
  Estr=''
  write(Estr,'(es12.6)') Einc
  col(1)='E-out'
  un(1)='MeV'
  col(2)='spectrum'
  un(2)='mb/MeV'
  col(3)='Maxwell_ratio'
  un(3)=''
  col(4)='spectrum_CM'
  un(4)='mb/MeV'
  Ncol=4
  quantity='emission spectrum'
  do type = 0, 6
    if (parskip(type)) cycle
    if (pfnsmodel == 1 .and. type /= 1) cycle
    write(*, '(/" +++ Prompt fission ", a8, " spectrum +++")') parname(type)
    pfnsfile = 'pfxs0000.000.fis'//natstring(iso)
    pfnsfile(3:3) = parsym(type)
    if (Einc0 < 0.001) then
      write(pfnsfile(5:12), '(es8.2)') Einc0
    else
      write(pfnsfile(5:12), '(f8.3)') Einc0
      write(pfnsfile(5:8), '(i4.4)') int(Einc0)
    endif
    open (unit = 1, file = pfnsfile, status = 'replace')
    reaction='('//parsym(k0)//',f)'
    topline=trim(targetnuclide)//trim(reaction)//' prompt fission '//trim(adjustl(parname(type)))//' '//trim(quantity)// &
 &    ' at '//Estr//' MeV'
    call write_header(indent,topline,source,user,date,oformat)
    call write_target(indent)
    if (type ==1) then
      MT = 18
    else
      MT = 0
    endif
    call write_reaction(indent,reaction,0.D0,0.D0,MF,MT)
    call write_char(id2,'ejectile',parname(type))
    call write_real(id2,'E-incident [MeV]',Einc)
    call write_real(id2,'E-average [MeV]',Eavpfns(type))
    call write_quantity(id2,quantity)
    call write_datablock(id2,Ncol,NEpfns,col,un)
    write(*, '(/" E-average         = ", f8.3, " MeV")') Eavpfns(type)
    write(*, '(" Weighted E-average= ", f8.3, " MeV")') Epfnsaverage(type)
    write(*, '(/"       E-out         spectrum    Maxwell ratio   spectrum_CM"/)')
    do nen = 1, NEpfns
      write(*, '(4es15.4)') Epfns(nen), pfns(type, nen), maxpfns(type, nen), pfnscm(type, nen)
      write(1, '(4es15.6)') Epfns(nen), pfns(type, nen), maxpfns(type, nen), pfnscm(type, nen)
    enddo
    close (unit = 1)
    call write_outfile(pfnsfile,flagoutall)
  enddo
  return
end subroutine pfnsout
! Copyright A.J. Koning 2021
