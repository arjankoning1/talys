subroutine recoilout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of recoils
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
!   dbl           ! double precision kind
! Variables for output
!   filerecoil    ! flag for recoil spectra on separate file
!   flagblock     ! flag to block spectra, angle and gamma files
! Variables for existence libraries
!   recexist      ! flag for existence of recoils
! Variables for numerics
!   maxenrec      ! number of recoil energies
!   maxN          ! maximal number of neutrons away from initial compound nucleus
!   maxZ          ! maximal number of protons away from initial compound nucleus
!   xseps         ! limit for cross sections
! Variables for main input
!   Atarget       ! mass number of target nucleus
!   k0            ! index of incident particle
! Variables for energy grid
!   Einc          ! incident energy in MeV
! Variables for incident channel
!   xselasinc     ! total elastic cross section (neutrons only) for inc. channel
!   xspopnuc      ! population cross section per nucleus
! Variables for nuclides
!   AA            ! mass number of residual nucleus
!   ZZ            ! charge number of residual nucleus
! Constants
!   iso           ! counter for isotope
!   natstring     ! string extension for file names
!   nuc           ! symbol of nucleus
!   parN          ! neutron number of particle
!   parsym        ! symbol of particle
!   parZ          ! charge number of particle
! Variables for recoil
!   Erec          ! recoil energy
!   recoilint     ! total recoil integrated over spectrum
!   specrecoil    ! recoil spectrum
!
! *** Declaration of local data
!
  implicit none
  character(len=3)  :: massstring !
  character(len=6)  :: finalnuclide !
  character(len=13) :: Estr
  character(len=18) :: reaction   ! reaction
  character(len=132) :: topline    ! topline
  character(len=15) :: col(2)     ! header
  character(len=15) :: un(2)     ! header
  character(len=80) :: quantity   ! quantity
  character(len=9)  :: recstring  ! string
  character(len=29) :: recfile    ! file with recoil spectra
  integer           :: A          ! mass number of target nucleus
  integer           :: MF
  integer           :: MT
  integer           :: Ncol       ! number of columns
  integer           :: Ncomp      ! neutron number index for compound nucleus
  integer           :: nen        ! energy counter
  integer           :: Z          ! charge number of target nucleus
  integer           :: Zcomp      ! proton number index for compound nucleus
  real(dbl)         :: sumcm      ! total residual production in the CM frame
!
! ***************************** Spectra ********************************
!
  MF = 6
  MT = 5
  Estr=''
  write(Estr,'(es13.6)') Einc
  col(1)='E-out'
  un(1)='MeV'
  col(2)='xs'
  un(2)='mb/MeV'
  Ncol=2
  quantity='emission spectrum'
  reaction='('//parsym(k0)//',x)'
  write(*, '(/" 8. Recoil spectra")')
  do Zcomp = 0, maxZ
    do Ncomp = 0, maxN
      if (xspopnuc(Zcomp, Ncomp) < xseps) cycle
      Z = ZZ(Zcomp, Ncomp, 0)
      A = AA(Zcomp, Ncomp, 0)
      massstring='   '
      write(massstring,'(i3)') A
      finalnuclide=trim(nuc(Z))//adjustl(massstring)
      write(*, '(/" Recoil Spectrum for ", i3, a2/)') A, nuc(Z)
      write(*, '("   Energy   Cross section"/)')
      do nen = 0, maxenrec
        write(*, '(1x, f8.3, es12.5)') Erec(Zcomp, Ncomp, nen), specrecoil(Zcomp, Ncomp, nen)
      enddo
      write(*, '(/" Integrated recoil spectrum       : ", es12.5)') recoilint(Zcomp, Ncomp)
      if (Zcomp == parZ(k0) .and. Ncomp == parN(k0)) then
        sumcm = xspopnuc(Zcomp, Ncomp) + xselasinc
      else
        sumcm = xspopnuc(Zcomp, Ncomp)
      endif
      write(*, '(" Residual production cross section: ", es12.5)') sumcm
!
! Write results to separate file
!
      if (filerecoil) then
        recstring = 'rec000000'
        write(recstring(4:9), '(2i3.3)') Z, A
        if (flagblock) then
          recfile = recstring//'.tot'//natstring(iso)
          if (.not. recexist(Zcomp,Ncomp)) then
            recexist(Zcomp, Ncomp) = .true.
            open (unit=1, file=recfile, status='unknown')
          else
            open (unit=1, file=recfile, status='unknown', position='append')
          endif
        else
          recfile = recstring//'spec0000.000.tot'//natstring(iso)
          write(recfile(14:21), '(f8.3)') Einc
          write(recfile(14:17), '(i4.4)') int(Einc)
          open (unit=1,file=recfile,status='unknown')
        endif
        topline=trim(targetnuclide)//trim(reaction)//' recoil '//trim(quantity)//' at '//Estr//' MeV'
        call write_header(topline,source,user,date,oformat)
        call write_target
        call write_reaction(reaction,0.D0,0.D0,MF,MT)
        call write_real(2,'E-incident [MeV]',Einc)
        call write_residual(Z,A,finalnuclide)
        call write_datablock(quantity,Ncol,maxenrec+1,col,un)
        do nen = 0, maxenrec
          write(1, '(2es15.6)') Erec(Zcomp, Ncomp, nen), specrecoil(Zcomp, Ncomp, nen)
        enddo
        close (unit = 1)
      endif
    enddo
  enddo
  return
end subroutine recoilout
! Copyright A.J. Koning 2021
