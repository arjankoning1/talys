subroutine binaryout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of binary cross sections
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
! Variables for energies
!   Ninclow       ! number of incident energies below Elow
! Variables for output
!   filetotal       ! flag for total cross sections on separate file
!   flagcompo       ! flag for output of cross section components
! Variables for input energies
!   eninc           ! incident energy in MeV
!   nin             ! counter for incident energy
!   Ninc          ! number of incident energies
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
! Variables for energy grid
!   Einc            ! incident energy in MeV
! Variables for incident channel
!   xsbinary        ! cross section from initial compound to residual nucleus
!   xsdirdisctot    ! direct cross section summed over discrete states
!   xsgrtot         ! total smoothed giant resonance cross section
!   xspreeqtot      ! preequilibrium cross section per particle type
! Variables for nuclides
!   parskip         ! logical to skip outgoing particle
! Constants
!   parname         ! name of particle
!   parsym          ! symbol of particle
! Variables for direct capture
!   xsracape        ! direct radiative capture cross section
! Variables for thermal cross sections
!   fxsbinary       ! cross section from initial compound to residual n
!
! *** Declaration of local data
!
  implicit none
  character(len=20) :: Qstring    ! Q-value
  character(len=10) :: binfile    ! file for binary output
  character(len=18) :: reaction   ! reaction
  character(len=132) :: topline    ! topline
  character(len=15) :: col(8)     ! header
  character(len=15) :: un(8)     ! units
  character(len=80) :: quantity   ! quantity
  integer           :: Ncol       ! number of columns
  integer           :: nen        ! energy counter
  integer           :: type       ! particle type
  real(sgl)         :: xsc        ! interpolated cross section
  real(sgl)         :: xsrac      ! direct radiative capture cross section
!
! ******************* Binary non-elastic channels **********************
!
  quantity='cross section'
  write(*, '(/" 2. Binary non-elastic cross sections ", "(non-exclusive)")')
  if (flagcompo) then
    write(*, '(36x, " Direct  Preequilibrium Compound  Dir. Capt.")')
  else
    write(*, '()')
  endif
  do type = - 1, 6
    if (parskip(type)) cycle
    if (type == 0) then
      xsrac = xsracape
    else
      xsrac = 0.
    endif
    if (flagcompo .and. type >= 0) then
      xsc = max(xsbinary(type) - xsdirdisctot(type) - xspreeqtot(type) - xsgrtot(type) - xsrac, 0.)
      write(*, '(1x, a8, "=", es12.5, 12x, 4es12.5)') parname(type), xsbinary(type), xsdirdisctot(type), &
 &      xspreeqtot(type) + xsgrtot(type), xsc, xsrac
    else
      write(*, '(1x, a8, "=", es12.5)') parname(type), xsbinary(type)
    endif
  enddo
!
! Write results to separate file
!
  if (filetotal) then
    binfile = 'binary.tot'
    if (nin == Ninclow + 1) then
      open (unit = 1, file = binfile, status = 'replace')
      reaction='('//parsym(k0)//',bin)'
      topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' - binary'
      col(1)='E'
      un(1)='MeV'
      do type=0,6
        col(2+type)=parname(type)
        un(2+type)= 'mb'
      enddo
      Ncol=8
      call write_header(topline,source,user,date,oformat)
      call write_target
      call write_reaction(reaction,0.D0,0.D0,0,0)
      do type = 0, 6
        Qstring = 'Q('//parsym(k0)//','//parsym(type)//') [MeV]'
        call write_real(2,Qstring,Q(type))
      enddo
      call write_datablock(quantity,Ncol,Ninc,col,un)
      do nen = 1, Ninclow
        write(1, '(8es15.6)') eninc(nen), (fxsbinary(nen, type), type = 0, 6)
      enddo
    else
      open (unit = 1, file = binfile, status = 'old', position = 'append')
    endif
    write(1, '(8es15.6)') Einc, (xsbinary(type), type = 0, 6)
    close (unit = 1)
  endif
  return
end subroutine binaryout
! Copyright A.J. Koning 2021
