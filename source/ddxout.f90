subroutine ddxout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of double-differential cross sections
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
!   sgl              ! single precision kind
! All global variables
!   numen2           ! maximum number of outgoing energies
! Variables for existence libraries
!   ddxexist1        ! flag for existence of DDX
!   ddxexist2        ! flag for existence of DDX
!   ddxexist3        ! flag for existence of DDX
!   ddxexist4        ! flag for existence of DDX
! Variables for output
!   ddxacount        ! counter for double - differential cross section files
!   ddxecount        ! counter for double - differential cross section files
!   ddxmode            ! mode for DDX: 0: None, 1: Angular distributions, 2: Spectra per angle, 3: Both
!   fileddxa         ! designator for double - differential cross sections on separate file
!   fileddxe         ! designator for double - differential cross sections on separate file
!   flagblock        ! flag to block spectra, angle and gamma files
! Variables for basic reaction
!   flaglabddx       ! flag for calculation of DDX in LAB system
!   flagrecoil       ! flag for calculation of recoils
! Variables for numerics
!   nanglecont       ! number of angles for continuum
! Variables for main input
!   Atarget          ! mass number of target nucleus
!   k0               ! index of incident particle
! Variables for incident channel
!   xsparticle       ! total particle production cross section
! Variables for energy grid
!   anglecont        ! angle in degrees for continuum
!   deltaE           ! energy bin around outgoing energies
! Variables for energy grid
!   ebegin           ! first energy point of energy grid
!   Einc             ! incident energy in MeV
! Variables for spectra
!   eendout          ! last energy point of energy grid
!   espec            ! outgoing energy grid
!   xscompoutad      ! compound emission angular distribution
!   xsdiscoutad      ! smoothed angular distribution for discrete state
!   xsmpreeqoutad    ! multiple preequilibrium angular distribution
!   xspreeqoutad     ! preequilibrium angular distribution per particle type
!   xssumout         ! cross section summed over mechanisms
!   xssumoutad       ! angular distribution summed over mechanisms
! Variables for nuclides
!   parskip          ! logical to skip outgoing particle
! Constants
!   iso              ! counter for isotope
!   natstring        ! string extension for file names
!   parname          ! name of particle
!   parsym           ! symbol of particle
! Variables for ECIS
!   anginc           ! angle increment
! Variables for recoil
!   ddxejlab         ! array containing the ddx spectrum of light part
!   Eejlab           ! center of ejectile lab bin
!   iejlab           ! number of ejectile lab bins
!
! *** Declaration of local data
!
  implicit none
  character(len=28) :: ddxfile         ! file for DDX
  character(len=13) :: Estr
  character(len=18) :: reaction   ! reaction
  character(len=132) :: topline    ! topline
  character(len=15) :: col(6)     ! header
  character(len=15) :: un(6)     ! header
  character(len=80) :: quantity   ! quantity
  integer           :: i               ! level
  integer           :: MF
  integer           :: MT
  integer           :: iang            ! running variable for angle
  integer           :: Ncol            ! number of columns
  integer           :: nen             ! energy counter
  integer           :: type            ! particle type
  real(sgl)         :: angf            ! help variable
  real(sgl)         :: enf             ! help variable
  real(sgl)         :: Eo(0:numen2)    ! outgoing energy grid based on excitation energy
  real(sgl)         :: Eout            ! outgoing energy
  real(sgl)         :: fac             ! factor
  real(sgl)         :: xs1             ! help variable
  real(sgl)         :: xs2             ! help variable
  real(sgl)         :: xs3             ! help variable
  real(sgl)         :: xs4             ! help variable
  real(sgl)         :: xs5             ! help variable
  real(sgl)         :: xsa             ! help variable
  real(sgl)         :: xsb             ! help variable
!
! ******************* Double-differential cross sections ***************
!
! locate       : subroutine to find value in ordered table
!
! 1. Angular distributions per outgoing energy
!
  MF = 6
  MT = 5
  Estr=''
  write(Estr,'(es13.6)') Einc
  un='mb/MeV/sr'
  col(1)='Angle'
  un(1)='deg'
  col(2)='xs'
  col(3)='Direct'
  col(4)='Preequilibrium'
  col(5)='Multiple_preeq'
  col(6)='Compound'
  Ncol=6
  quantity='double-differential cross section'
  if (ddxmode == 1 .or. ddxmode == 3) then
    write(*, '(/" 9. Double-differential cross sections per outgoing energy")')
    do type = 1, 6
      if (parskip(type)) cycle
      if (xsparticle(type) == 0.) cycle
      do nen = ebegin(type), eendout(type)
        Eo(nen) = espec(type, nen)
        if (xssumout(type, nen) == 0.) cycle
        write(*, '(/" DDX for outgoing ", a8, " at ", f8.3, " MeV"/)') parname(type), Eo(nen)
        write(*, '(" Angle   Total      Direct     Pre-equil.  Mult. preeq   Compound"/)')
        do iang = 0, nanglecont
          write(*, '(1x, f5.1, 5es12.5)') anglecont(iang), xssumoutad(type, nen, iang), xsdiscoutad(type, nen, iang), &
 &          xspreeqoutad(type, nen, iang), xsmpreeqoutad(type, nen, iang), xscompoutad(type, nen, iang)
        enddo
      enddo
      do i = 1, ddxecount(type)
        enf = fileddxe(type, i)
        call locate(Eo, ebegin(type), eendout(type), enf, nen)
        fac = (enf - Eo(nen)) / deltaE(nen)
        if (flagblock) then
          ddxfile = ' ddx.MeV'//natstring(iso)
          write(ddxfile(1:1), '(a1)') parsym(type)
          if (.not. ddxexist1(type)) then
            ddxexist1(type) = .true.
            open (unit=1, file=ddxfile, status='unknown')
          else
            open (unit=1, file=ddxfile, status='unknown', position='append')
          endif
        else
          ddxfile = ' ddxE0000.000E0000.0.MeV'//natstring(iso)
          write(ddxfile(1:1), '(a1)') parsym(type)
          write(ddxfile(6:13), '(f8.3)') Einc
          write(ddxfile(6:9), '(i4.4)') int(Einc)
          write(ddxfile(15:20), '(f6.1)') enf
          write(ddxfile(15:18), '(i4.4)') int(enf)
          open (unit=1, file=ddxfile, status='unknown')
        endif
        write(ddxfile(1:1), '(a1)') parsym(type)
        write(ddxfile(6:13), '(f8.3)') Einc
        write(ddxfile(6:9), '(i4.4)') int(Einc)
        write(ddxfile(15:20), '(f6.1)') enf
        write(ddxfile(15:18), '(i4.4)') int(enf)
        reaction='('//parsym(k0)//',x'//parsym(type)//')'
        topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
        call write_header(topline,source,user,date,oformat)
        call write_target
        call write_reaction(reaction,0.D0,0.D0,MF,MT)
        call write_real(2,'E-incident [MeV]',Einc)
        call write_real(2,'E-emission [MeV]',enf)
        call write_datablock(quantity,Ncol,nanglecont+1,col,un)
        do iang = 0, nanglecont
          xsa = xssumoutad(type, nen, iang)
          xsb = xssumoutad(type, nen + 1, iang)
          xs1 = xsa + fac * (xsb - xsa)
          xsa = xsdiscoutad(type, nen, iang)
          xsb = xsdiscoutad(type, nen + 1, iang)
          xs2 = xsa + fac * (xsb - xsa)
          xsa = xspreeqoutad(type, nen, iang)
          xsb = xspreeqoutad(type, nen + 1, iang)
          xs3 = xsa + fac * (xsb - xsa)
          xsa = xsmpreeqoutad(type, nen, iang)
          xsb = xsmpreeqoutad(type, nen + 1, iang)
          xs4 = xsa + fac * (xsb - xsa)
          xsa = xscompoutad(type, nen, iang)
          xsb = xscompoutad(type, nen + 1, iang)
          xs5 = xsa + fac * (xsb - xsa)
          write(1, '(6es15.6)') anglecont(iang), xs1, xs2, xs3, xs4, xs5
        enddo
        close (unit = 1)
      enddo
!
! Results in LAB frame
!
      if (flagrecoil .and. flaglabddx) then
        do nen = 1, iejlab(type)
          Eo(nen) = Eejlab(type, nen)
          write(*, '(/" DDX for outgoing ", a8, " at ", f8.3, " MeV in LAB frame")') parname(type), Eo(nen)
          write(*, '(" Angle   Cross section"/)')
          do iang = 0, nanglecont
            write(*, '(1x, f5.1, es12.5)') anglecont(iang), ddxejlab(type, nen, iang)
          enddo
        enddo
        do i = 1, ddxecount(type)
          enf = fileddxe(type, i)
          call locate(Eo, 1, iejlab(type), enf, nen)
          fac = (enf - Eo(nen)) / deltaE(nen)
          if (flagblock) then
            ddxfile = ' ddx.lab'//natstring(iso)
            write(ddxfile(1:1), '(a1)') parsym(type)
            if (.not. ddxexist2(type)) then
              ddxexist2(type) = .true.
              open (unit=1, file=ddxfile, status='unknown')
            else
              open (unit=1, file=ddxfile, status='unknown', position='append')
            endif
          else
            ddxfile = ' ddxE0000.000E0000.0.lab'//natstring(iso)
            write(ddxfile(1:1), '(a1)') parsym(type)
            write(ddxfile(6:13), '(f8.3)') Einc
            write(ddxfile(6:9), '(i4.4)') int(Einc)
            write(ddxfile(15:20), '(f6.1)') enf
            write(ddxfile(15:18), '(i4.4)') int(enf)
            open (unit=1, file=ddxfile, status='unknown')
          endif
          write(ddxfile(1:1), '(a1)') parsym(type)
          write(ddxfile(6:13), '(f8.3)') Einc
          write(ddxfile(6:9), '(i4.4)') int(Einc)
          write(ddxfile(15:20), '(f6.1)') enf
          write(ddxfile(15:18), '(i4.4)') int(enf)
          quantity='double-differential cross section in LAB frame'
          reaction='('//parsym(k0)//',x'//parsym(type)//')'
          topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
          call write_header(topline,source,user,date,oformat)
          call write_target
          call write_reaction(reaction,0.D0,0.D0,MF,MT)
          call write_real(2,'E-incident [MeV]',Einc)
          call write_real(2,'E-emission [MeV]',enf)
          call write_datablock(quantity,Ncol,nanglecont+1,col,un)
          do iang = 0, nanglecont
            xsa = ddxejlab(type, nen, iang)
            xsb = ddxejlab(type, nen + 1, iang)
            xs1 = xsa + fac * (xsb - xsa)
            write(1, '(2es15.6)') anglecont(iang), xs1
          enddo
          close (unit = 1)
        enddo
      endif
    enddo
  endif
!
! 2. Emission spectra per outgoing angle
!
  if (ddxmode == 2 .or. ddxmode == 3) then
    col(1)='E-out'
    un(1)='MeV'
    anginc = 180. / nanglecont
    write(*, '(/" 9. Double-differential cross sections per outgoing angle")')
    do type = 1, 6
      if (parskip(type)) cycle
      if (xsparticle(type) == 0.) cycle
      do iang = 0, nanglecont
        write(*, '(/" DDX for outgoing ", a8, " at ", f7.3, " degrees")') parname(type), anglecont(iang)
        write(*, '(/"    E-out    Total      Direct     Pre-equil. Mult. preeq   Compound"/)')
        do nen = ebegin(type), eendout(type)
          Eout = espec(type, nen)
          write(*, '(1x, f8.3, 5es12.5)') Eout, xssumoutad(type, nen, iang), xsdiscoutad(type, nen, iang), &
 &          xspreeqoutad(type, nen, iang), xsmpreeqoutad(type, nen, iang), xscompoutad(type, nen, iang)
        enddo
      enddo
      do i = 1, ddxacount(type)
        angf = fileddxa(type, i)
        call locate(anglecont, 0, nanglecont, angf, iang)
        fac = (angf - anglecont(iang)) / anginc
        if (flagblock) then
          ddxfile = ' ddx.deg'//natstring(iso)
          write(ddxfile(1:1), '(a1)') parsym(type)
          if (.not. ddxexist3(type)) then
            ddxexist3(type) = .true.
            open (unit=1, file=ddxfile, status='unknown')
          else
            open (unit=1, file=ddxfile, status='unknown', position='append')
          endif
        else
          ddxfile = ' ddxE0000.000A000.0.deg'//natstring(iso)
          write(ddxfile(1:1), '(a1)') parsym(type)
          write(ddxfile(6:13), '(f8.3)') Einc
          write(ddxfile(6:9), '(i4.4)') int(Einc)
          write(ddxfile(15:19), '(f5.1)') angf
          write(ddxfile(15:17), '(i3.3)') int(angf)
          open (unit=1, file=ddxfile, status='unknown')
        endif
        write(ddxfile(1:1), '(a1)') parsym(type)
        write(ddxfile(6:13), '(f8.3)') Einc
        write(ddxfile(6:9), '(i4.4)') int(Einc)
        write(ddxfile(15:19), '(f5.1)') angf
        write(ddxfile(15:17), '(i3.3)') int(angf)
        quantity='double-differential cross section'
        reaction='('//parsym(k0)//',x'//parsym(type)//')'
        topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
        call write_header(topline,source,user,date,oformat)
        call write_target
        call write_reaction(reaction,0.D0,0.D0,MF,MT)
        call write_real(2,'E-incident [MeV]',Einc)
        call write_real(2,'Angle [deg]',angf)
        call write_datablock(quantity,Ncol,eendout(type)-ebegin(type)+1,col,un)
        do nen = ebegin(type), eendout(type)
          Eout = espec(type, nen)
          xsa = xssumoutad(type, nen, iang)
          xsb = xssumoutad(type, nen, iang + 1)
          xs1 = xsa + fac * (xsb - xsa)
          xsa = xsdiscoutad(type, nen, iang)
          xsb = xsdiscoutad(type, nen, iang + 1)
          xs2 = xsa + fac * (xsb - xsa)
          xsa = xspreeqoutad(type, nen, iang)
          xsb = xspreeqoutad(type, nen, iang + 1)
          xs3 = xsa + fac * (xsb - xsa)
          xsa = xsmpreeqoutad(type, nen, iang)
          xsb = xsmpreeqoutad(type, nen, iang + 1)
          xs4 = xsa + fac * (xsb - xsa)
          xsa = xscompoutad(type, nen, iang)
          xsb = xscompoutad(type, nen, iang + 1)
          xs5 = xsa + fac * (xsb - xsa)
          write(1, '(6es15.6)') Eout, xs1, xs2, xs3, xs4, xs5
        enddo
        close (unit = 1)
      enddo
!
! Results in LAB frame
!
      if (flagrecoil .and. flaglabddx) then
        do iang = 0, nanglecont
          write(*, '(/" DDX for outgoing ", a8, " at ", f8.3, " degrees in LAB frame"/)') parname(type), anglecont(iang)
          write(*, '("  Energy  Cross section"/)')
          do nen = 1, iejlab(type)
            write(*, '(1x, f8.3, es12.5)') Eejlab(type, nen), ddxejlab(type, nen, iang)
          enddo
        enddo
        do i = 1, ddxacount(type)
          angf = fileddxa(type, i)
          call locate(anglecont, 0, nanglecont, angf, iang)
          fac = (angf - anglecont(iang)) / anginc
          if (flagblock) then
            ddxfile = ' ddx.lab'//natstring(iso)
            write(ddxfile(1:1), '(a1)') parsym(type)
            if (.not. ddxexist4(type)) then
              ddxexist4(type) = .true.
              open (unit=1, file=ddxfile, status='unknown')
            else
              open (unit=1, file=ddxfile, status='unknown', position='append')
            endif
          else
            ddxfile = ' ddxE0000.000A000.0.lab'//natstring(iso)
            write(ddxfile(1:1), '(a1)') parsym(type)
            write(ddxfile(6:13), '(f8.3)') Einc
            write(ddxfile(6:9), '(i4.4)') int(Einc)
            write(ddxfile(15:19), '(f5.1)') angf
            write(ddxfile(15:17), '(i3.3)') int(angf)
            open (unit=1, file=ddxfile, status='unknown')
          endif
          ddxfile = ' ddxE0000.000A000.0.lab'//natstring(iso)
          write(ddxfile(1:1), '(a1)') parsym(type)
          write(ddxfile(6:13), '(f8.3)') Einc
          write(ddxfile(6:9), '(i4.4)') int(Einc)
          write(ddxfile(15:19), '(f5.1)') angf
          write(ddxfile(15:17), '(i3.3)') int(angf)
          quantity='double-differential cross section in LAB frame'
          reaction='('//parsym(k0)//',x'//parsym(type)//')'
          topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
          call write_header(topline,source,user,date,oformat)
          call write_target
          call write_reaction(reaction,0.D0,0.D0,MF,MT)
          call write_real(2,'E-incident [MeV]',Einc)
          call write_real(2,'Angle [deg]',angf)
          call write_datablock(quantity,Ncol,eendout(type)-ebegin(type)+1,col,un)
          do nen = 1, iejlab(type)
            xsa = ddxejlab(type, nen, iang)
            xsb = ddxejlab(type, nen, iang + 1)
            xs1 = xsa + fac * (xsb - xsa)
            write(1, '(6es15.6)') Eejlab(type, nen), xs1
          enddo
        enddo
        close (unit = 1)
      endif
    enddo
  endif
  return
end subroutine ddxout
! Copyright A.J. Koning 2021
