subroutine totalout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of total cross sections
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
!   filetotal        ! flag for total cross sections on separate file
! Variables for input energies
!   eninc            ! incident energy in MeV
!   flaginitpop      ! flag for initial population distribution
!   nin              ! counter for incident energy
!   Ninc             ! number of incident energies
! Variables for main input
!   Atarget          ! mass number of target nucleus
!   k0               ! index of incident particle
! Variables for gamma rays
!   flagracap        ! flag for radiative capture model
! Variables for multiple emission
!   xsinitpop        ! initial population cross section
! Variables for energy grid
!   Einc             ! incident energy in MeV
! Variables for energies
!   eninccm          ! center - of - mass incident energy in MeV
!   flaggiant        ! flag for collective contribution from giant resonances
!   Ninclow        ! number of incident energies below Elow
! Variables for binary reactions
!   xscompel         ! compound elastic cross section
!   xscompnonel      ! total compound non - elastic cross section
!   xselastot        ! total elastic cross section (shape + compound)
!   xsnonel          ! non - elastic cross
! Variables for incident channel
!   xsdirdiscsum     ! total direct cross section
!   xselasinc        ! total elastic cross section (neutrons only) for inc. channel
!   xsgrsum          ! sum over giant resonance cross sections
!   xspreeqsum       ! total preequilibrium cross section summed over particles
!   xsreacinc        ! reaction cross section for incident channel
!   xstotinc         ! total cross section (neutrons only) for incident channel
! Constants
!   iso              ! counter for isotope
!   natstring        ! string extension for file names
!   parsym           ! symbol of particle
! Variables for direct capture
!   xsracape         ! direct radiative capture cross section
! Variables for thermal cross sections
!   fxscompel        ! compound elastic cross section
!   fxscompnonel     ! total compound non - elastic cross section
!   fxsdirdiscsum    ! total direct cross section
!   fxselasinc       ! total elastic cross section (neutrons only) for i
!   fxselastot       ! total elastic cross section (neutrons only) for i
!   fxsnonel         ! non - elastic cross section for incident channel
!   fxspreeqsum      ! total preequilibrium cross section summed over pa
!   fxsracape        ! direct capture cross section
!   fxsreacinc       ! reaction cross section for incident channel
!   fxstotinc        ! total cross section (neutrons only) for incident
!
! *** Declaration of local data
!
  implicit none
  character(len=18) :: reaction   ! reaction
  character(len=15) :: col(11)    ! header
  character(len=15) :: un(11)     ! header
  character(len=18) :: totfile    ! file with total cross sections
  character(len=80) :: quantity   ! quantity
  character(len=132) :: topline   ! topline
  integer           :: nen        ! energy counter
  integer           :: Ncol       ! counter
!
! ********************* Total cross sections ***************************
!
  if (Einc >= 0.001) then
    write(*, '(/" ########### REACTION SUMMARY FOR E=", f10.5, " ###########"/)') Einc
  else
    write(*, '(/" ########### REACTION SUMMARY FOR E=", es12.5, " ###########"/)') Einc
  endif
  if (flaginitpop) then
    write(*, '(" 1. Initial population cross section =", es12.5/)') xsinitpop
    return
  endif
  write(*, '(" Center-of-mass energy: ", f8.3/)') eninccm
  write(*, '(" 1. Total (binary) cross sections"/)')
  if (k0 <= 1) write(*, '(" Total           =", es12.5)') xstotinc
  if (k0 <= 1) write(*, '("   Shape elastic   =", es12.5)') xselasinc
  write(*, '("   Reaction        =", es12.5)') xsreacinc
  write(*, '("     Compound elastic=", es12.5)') xscompel
  write(*, '("     Non-elastic     =", es12.5)') xsnonel
  write(*, '("       Direct          =", es12.5)') xsdirdiscsum
  write(*, '("       Pre-equilibrium =", es12.5)') xspreeqsum
  if (flagracap) write(*, '("       Direct Capture  =", es12.5)') xsracape
  if (flaggiant) write(*, '("       Giant resonance =", es12.5)') xsgrsum
  write(*, '("       Compound non-el =", es12.5)') xscompnonel
  if (k0 <= 1) write(*, '("     Total elastic   =", es12.5)') xselastot
!
! Write results to separate file
!
  if (filetotal) then
    quantity='cross section'
    un = 'mb'
    col(1)='E'
    un(1)='MeV'
    col(2)='xs'
    col(2)='Non-elastic'
    col(3)='Elastic'
    col(4)='Total'
    col(5)='Compound_elast.'
    col(6)='Shape_elastic'
    col(7)='Reaction'
    col(8)='Compound_nonel.'
    col(9)='Direct'
    col(10)='Preequilibrium'
    col(11)='Direct_capture'
    Ncol=11
    totfile = 'all.tot'//natstring(iso)
    if (nin == Ninclow + 1) then
      open (unit = 1, file = totfile, status = 'replace')
      reaction='('//parsym(k0)//',all)'
      topline=trim(targetnuclide)//trim(reaction)//' general '//trim(quantity)//'s [mb]'
      call write_header(topline,source,user,date,oformat)
      call write_target
      call write_reaction(reaction,0.d0,0.d0,0,0)
      call write_datablock(quantity,Ncol,Ninc,col,un)
      do nen = 1, Ninclow
        write(1, '(11es15.6)') eninc(nen), fxsnonel(nen), fxselastot(nen), fxstotinc(nen), fxscompel(nen), &
 &        fxselasinc(nen), fxsreacinc(nen), fxscompnonel(nen), fxsdirdiscsum(nen), fxspreeqsum(nen),fxsracape(nen)
      enddo
    else
      open (unit = 1, file = totfile, status = 'old', position = 'append')
    endif
    write(1, '(11es15.6)') Einc, xsnonel, xselastot, xstotinc, xscompel, xselasinc, xsreacinc, xscompnonel, &
 &    xsdirdiscsum, xspreeqsum, xsracape
    close (unit = 1)
!
! Total cross sections (i.e. from OMP) only
!
    quantity='cross section'
    col(1)='E'
    un(1)='MeV'
    col(2)='xs'
    un(2)='mb'
    Ncol=2
    totfile = 'total.tot'//natstring(iso)
    reaction='('//parsym(k0)//',tot)'
    if (nin == Ninclow + 1) then
      open (unit = 1, file = totfile, status = 'replace')
      topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
      call write_header(topline,source,user,date,oformat)
      call write_target
      call write_reaction(reaction,0.d0,0.d0,3,1)
      call write_datablock(quantity,Ncol,Ninc,col,un)
      do nen = 1, Ninclow
        write(1, '(2es15.6)') eninc(nen), fxstotinc(nen)
      enddo
    else
      open (unit = 1, file = totfile, status = 'old', position = 'append')
    endif
    write(1, '(2es15.6)') Einc, xstotinc
    close (unit = 1)
!
! Elastic cross sections only
!
    totfile = 'elastic.tot'//natstring(iso)
    reaction='('//parsym(k0)//',el)'
    if (nin == Ninclow + 1) then
      open (unit = 1, file = totfile, status = 'replace')
      topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
      call write_header(topline,source,user,date,oformat)
      call write_target
      call write_reaction(reaction,0.d0,0.d0,3,2)
      call write_datablock(quantity,Ncol,Ninc,col,un)
      do nen = 1, Ninclow
        write(1, '(2es15.6)') eninc(nen), fxselastot(nen)
      enddo
    else
      open (unit = 1, file = totfile, status = 'old', position = 'append')
    endif
    write(1, '(2es15.6)') Einc, xselastot
    close (unit = 1)
!
! Nonelastic cross sections only
!
    totfile = 'nonelastic.tot'//natstring(iso)
    reaction='('//parsym(k0)//',non)'
    if (nin == Ninclow + 1) then
      open (unit = 1, file = totfile, status = 'replace')
      topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
      call write_header(topline,source,user,date,oformat)
      call write_target
      call write_reaction(reaction,0.d0,0.d0,3,3)
      call write_datablock(quantity,Ncol,Ninc,col,un)
      do nen = 1, Ninclow
        write(1, '(2es15.6)') eninc(nen), fxsnonel(nen)
      enddo
    else
      open (unit = 1, file = totfile, status = 'old', position = 'append')
    endif
    write(1, '(2es15.6)') Einc, xsnonel
    close (unit = 1)
!
! Reaction cross sections only
!
    totfile = 'reaction.tot'//natstring(iso)
    reaction='('//parsym(k0)//',reac)'
    if (nin == Ninclow + 1) then
      open (unit = 1, file = totfile, status = 'replace')
      topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
      call write_header(topline,source,user,date,oformat)
      call write_target
      call write_reaction(reaction,0.d0,0.d0,0,0)
      call write_datablock(quantity,Ncol,Ninc,col,un)
      do nen = 1, Ninclow
        write(1, '(2es15.6)') eninc(nen), fxsreacinc(nen)
      enddo
    else
      open (unit = 1, file = totfile, status = 'old', position = 'append')
    endif
    write(1, '(2es15.6)') Einc, xsreacinc
    close (unit = 1)
  endif
  return
end subroutine totalout
! Copyright A.J. Koning 2021
