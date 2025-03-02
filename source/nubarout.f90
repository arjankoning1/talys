subroutine nubarout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of average number of fission neutrons
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
! Variables for input energies
!   eninc           ! incident energy in MeV
!   nin0            ! counter for incident energy
!   Ninc            ! number of incident energies
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
! Variables for energies
!   Einc0           ! incident energy in MeV
!   Ninclow       ! number of incident energies below Elow
! Constants
!   iso             ! counter for isotope
!   natstring       ! string extension for file names
!   parname         ! name of particle
!   parsym          ! symbol of particle
! Variables for thermal cross sections
!   fnubar          ! nubar
! Variables for mass distribution
!   nubar           ! average nu
! Variables for existence libraries
!   nubarexist      ! flag for existence of nubar file
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: nufile    ! file for nubar
  character(len=18) :: reaction   ! reaction
  character(len=132) :: topline    ! topline
  character(len=15) :: col(4)     ! header
  character(len=15) :: un(4)     ! header
  character(len=80) :: quantity   ! quantity
  integer           :: MF
  integer           :: MT
  integer           :: Ncol      ! number of columns
  integer           :: nen       ! energy counter
  integer           :: type      ! particle type
!
! Write results to separate files
!
! nu per number, P(nu) and nubar as function of Z and A
!
  MF = 1
  quantity='multiplicity'
  un = ''
  col(1)='E'
  un(1)='MeV'
  col(2)='nubar'
  Ncol=2
  if (nin0 == Ninclow+1 .and. Ninclow > 0) then
    do nen = 1, Ninclow
      do type = 0, 6
        fnubar(nen, type) = nubar(type)
      enddo
    enddo
  endif
  write(*, '(/" +++ AVERAGE NUMBER OF PROMPT FISSION NEUTRONS +++"/)')
  do type = 0, 6
    if (nubar(type) == 0.) cycle
!
! Write nubar
!
    nufile = 'nubarx.tot'//natstring(iso)//'       '
    nufile(6:6) = parsym(type)
    if ( .not. nubarexist(type)) then
      nubarexist(type) = .true.
      open (unit = 1, file = nufile, status = 'replace')
      reaction='('//parsym(k0)//',f)'
      topline=trim(targetnuclide)//trim(reaction)//' average '//trim(adjustl(parname(type)))//' '//trim(quantity)
      if (type == 1) then
        MT = 456
      else
        MT = 0
      endif
      call write_header(topline,source,user,date,oformat)
      call write_target
      call write_reaction(reaction,0.d0,0.d0,MF,MT)
      call write_char(2,'ejectile',parname(type))
      call write_quantity(quantity)
      call write_datablock(Ncol,Ninc,col,un)
      do nen = 1, Ninclow
        write(1, '(2es15.6)') eninc(nen), fnubar(nen, type)
      enddo
      do nen = Ninclow + 1, nin0 - 1
        write(1, '(2es15.6)') eninc(nen), 0.
      enddo
    else
      open (unit = 1, file = nufile, status = 'old', position = 'append')
    endif
    write(1, '(2es15.6)') Einc0, nubar(type)
    write(*, '(a8, es12.5)') parname(type), nubar(type)
    close (unit = 1)
  enddo
  return
end subroutine nubarout
! Copyright A.J. Koning 2021
