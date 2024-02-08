subroutine spectraout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of particle spectra
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
!   filespectrum    ! designator for spectrum on separate file
!   flagblock       ! flag to block spectra, angle and gamma files
! Variables for existence libraries
!   spexist1        ! flag for existence of spectra
!   spexist2        ! flag for existence of spectra
! Variables for basic reaction
!   flaglabddx      ! flag for calculation of DDX in LAB system
!   flagrecoil      ! flag for calculation of recoils
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
! Variables for incident channel
!   xsparticle      ! total particle production cross section
! Variables for energies
!   Ninclow        ! number of incident energies below Elow
! Variables for energy grid
!   ebegin          ! first energy point of energy grid
!   Einc            ! incident energy in MeV
! Variables for input energies
!   eninc            ! incident energy in MeV
!   nin              ! counter for incident energy
!   Ninc             ! number of incident energies
! Variables for spectra
!   buratio         ! break - up ratio
!   Eaverage        ! average outgoing energy
!   eendout         ! last energy point of energy grid
!   espec           ! outgoing energy grid
!   preeqratio      ! pre - equilibrium ratio
!   xscompout       ! compound emission angular distribution
!   xsdiscout       ! smoothed angular distribution for discrete state
!   xsmpreeqout     ! multiple preequilibrium angular distribution
!   xspreeqbuout    ! preequilibrium cross section for breakup
!   xspreeqkiout    ! preequilibrium cross section for knockout and inelastic
!   xspreeqout      ! preequilibrium angular distribution per particle type an
!   xspreeqpsout    ! preequilibrium cross section for pickup and stripping
!   xssumout        ! cross section summed over mechanisms
! Variables for nuclides
!   parskip         ! logical to skip outgoing particle
!   Q               ! Q - value
! Constants
!   iso             ! counter for isotope
!   natstring       ! string extension for file names
!   parname         ! name of particle
!   parsym          ! symbol of particle
! Variables for recoil
!   Eejlab          ! center of ejectile lab bin
!   iejlab          ! number of ejectile lab bins
!   xsejlab         ! LAB ejectile spectrum
!   xsejlabint      ! LAB energy - integrated spectrum
!
! *** Declaration of local data
!
  implicit none
  character(len=10) :: Efile       ! file with average energies
  character(len=21) :: specfile    ! file with composite particle spectra
  character(len=13) :: Estr
  character(len=18) :: reaction   ! reaction
  character(len=132) :: topline    ! topline
  character(len=15) :: col(11)     ! header
  character(len=15) :: un(11)     ! header
  character(len=80) :: quantity   ! quantity
  integer           :: MF
  integer           :: MT
  integer           :: nen         ! energy counter
  integer           :: type        ! particle type
  integer           :: Ncol        ! number of columns
!
! ***************************** Spectra ********************************
!
  MF = 6
  MT = 5
  Estr=''
  write(Estr,'(es13.6)') Einc
  un = 'mb/MeV'
  col(1)='E-out'
  un(1)='MeV'
  col(2)='xs'
  col(3)='Direct'
  col(4)='Preequilibrium'
  col(5)='Multiple_preeq'
  col(6)='Compound'
  col(7)='Preeq_ratio'
  un(7)=''
  col(8)='Break-up_ratio'
  un(8)=''
  col(9)='Stripping'
  col(10)='Knock-out'
  col(11)='Break-up'
  Ncol=7
  quantity='emission spectrum'
  write(*, '(/" 7. Composite particle spectra")')
  do type = 0, 6
    if (parskip(type)) cycle
    if (xsparticle(type) == 0.) cycle
    write(*, '(/" Spectra for outgoing ", a8/)') parname(type)
    if (k0 <= 2 .and. type <= 2) then
      write(*, '("  Energy   Total       Direct    Pre-equil.  Mult. preeq  Compound"/)')
      do nen = ebegin(type), eendout(type)
        write(*, '(f8.3, 5es12.5)') espec(type, nen), xssumout(type, nen), xsdiscout(type, nen), &
 &        xspreeqout(type, nen), xsmpreeqout(type, nen), xscompout(type, nen)
      enddo
    else
      write(*, '("  Energy   Total       Direct    Pre-equil.  Mult. preeq  Compound    Stripping   Knock-out   Break-up"/)')
      do nen = ebegin(type), eendout(type)
        write(*, '(f8.3, 8es12.5)') espec(type, nen), xssumout(type, nen), xsdiscout(type, nen), &
 &        xspreeqout(type, nen), xsmpreeqout(type, nen), xscompout(type, nen), xspreeqpsout(type, nen), &
 &        xspreeqkiout(type, nen), xspreeqbuout(type, nen)
      enddo
    endif
    if (flagrecoil .and. flaglabddx) then
      write(*, '(/" LAB spectra for outgoing ", a8/)') parname(type)
      write(*, '("  Energy   Cross section"/)')
      do nen = 1, iejlab(type)
        write(*, '(f8.3, es12.5)') Eejlab(type, nen), xsejlab(type, nen)
      enddo
      write(*, '(/" Energy-integrated cross section:", es12.5/)') xsejlabint(type)
    endif
!
! Write results to separate file
!
    if (filespectrum(type)) then
      if (flagblock) then
        specfile = ' spec.tot'//natstring(iso)
        write(specfile(1:1), '(a1)') parsym(type)
        if (.not. spexist1(type)) then
          spexist1(type) = .true.
          open (unit=1, file=specfile, status='unknown')
        else
          open (unit=1, file=specfile, status='unknown', position='append')
        endif
      else
        specfile=' spec0000.000.tot'//natstring(iso)
        write(specfile(1:1), '(a1)') parsym(type)
        write(specfile(6:13), '(f8.3)') Einc
        write(specfile(6:9), '(i4.4)') int(Einc)
        open (unit=1, file=specfile, status='unknown')
      endif
      reaction='('//parsym(k0)//',x'//parsym(type)//')'
      topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
      call write_header(topline,source,user,date,oformat)
      call write_target
      call write_reaction(reaction,0.D0,0.D0,MF,MT)
      call write_real(2,'E-incident [MeV]',Einc)
      call write_real(2,'E-average [MeV]',Eaverage(type))
      if (k0 <= 2 .and. type <= 2) then
        Ncol=7
        call write_datablock(quantity,Ncol,eendout(type)-ebegin(type)+1,col,un)
        do nen = ebegin(type), eendout(type)
          write(1, '(7es15.6)') espec(type, nen), xssumout(type, nen), xsdiscout(type, nen), &
 &          xspreeqout(type, nen), xsmpreeqout(type, nen), xscompout(type, nen), preeqratio(type, nen)
        enddo
      else
        Ncol=11
        call write_datablock(quantity,Ncol,eendout(type)-ebegin(type)+1,col,un)
        do nen = ebegin(type), eendout(type)
          write(1, '(11es15.6)') espec(type, nen), xssumout(type, nen), xsdiscout(type, nen), &
 &          xspreeqout(type, nen), xsmpreeqout(type, nen), xscompout(type, nen), preeqratio(type, nen), &
 &          buratio(type, nen), xspreeqpsout(type, nen), xspreeqkiout(type, nen), xspreeqbuout(type, nen)
        enddo
      endif
      close (unit = 1)
      if (flagrecoil .and. flaglabddx) then
        if (flagblock) then
          specfile = ' spec.lab'//natstring(iso)
          write(specfile(1:1), '(a1)') parsym(type)
          if (.not. spexist2(type)) then
            spexist2(type) = .true.
            open (unit=1, file=specfile, status='unknown')
          else
            open (unit=1, file=specfile, status='unknown', position='append')
          endif
        else
          specfile=' spec0000.000.lab'//natstring(iso)
          write(specfile(1:1), '(a1)') parsym(type)
          write(specfile(6:13), '(f8.3)') Einc
          write(specfile(6:9), '(i4.4)') int(Einc)
          open (unit=1, file=specfile, status='unknown')
        endif
        quantity='emission spectrum in LAB frame'
        topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
        Ncol=2
        call write_header(topline,source,user,date,oformat)
        call write_target
        call write_reaction(reaction,0.D0,0.D0,MF,MT)
        call write_real(2,'E-incident [MeV]',Einc)
        call write_datablock(quantity,Ncol,iejlab(type),col,un)
        do nen = 1, iejlab(type)
          write(1, '(2es15.6)') Eejlab(type, nen), xsejlab(type, nen)
        enddo
        close (unit = 1)
      endif
    endif
  enddo
  write(*, '(/" Average emission energies"/)')
  do type = 0, 6
    if (parskip(type)) cycle
    write(*, '(1x, a8, 4x, f8.3)') parname(type), Eaverage(type)
    if (filespectrum(type)) then
      Efile = 'Eaverage.'//parsym(type)
      if (nin == Ninclow + 1) then
        open (unit = 1, file = Efile, status = 'replace')
        quantity='average emission energy'
        reaction='('//parsym(k0)//',x'//parsym(type)//')'
        topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
        Ncol=2
        col(1)='E'
        col(2)='E-average'
        un(2)='MeV'
        call write_header(topline,source,user,date,oformat)
        call write_target
        call write_reaction(reaction,0.D0,0.D0,0,0)
        call write_datablock(quantity,Ncol,Ninclow,col,un)
        do nen = 1, Ninclow
          write(1, '(2es15.6)') eninc(nen), 0.
        enddo
      else
        open (unit = 1, file = Efile, status = 'old', position = 'append')
      endif
      write(1, '(2es15.6)') Einc, Eaverage(type)
      close (unit = 1)
    endif
  enddo
  return
end subroutine spectraout
! Copyright A.J. Koning 2021
