subroutine gamdisout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of discrete gamma-ray intensities
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
!   sgl            ! single precision kind
! All global variables
!   numlev         ! maximum number of discrete levels
! Variables for existence libraries
!   gamexist       ! flag for existence of gamma production cross section
! Variables for numerics
!   maxN           ! maximal number of neutrons away from initial compound nucleus
!   maxZ           ! maximal number of protons away from initial compound nucleus
!   popeps         ! limit for population cross sections
! Variables for input energies
!   eninc          ! incident energy in MeV
!   nin            ! counter for incident energy
!   Ninc         ! number of incident energies
! Variables for main input
!   Atarget        ! mass number of target nucleus
!   k0             ! index of incident particle
! Variables for output
!   filegamdis     ! flag for gamma-ray intensities on separate file
! Variables for energies
!   Ethresh        ! threshold incident energy for residual nucleus
!   Ninclow        ! number of incident energies below Elow
! Variables for multiple emission
!   xsgamdis       ! discrete gamma-ray cross section
!   xsgamdistot    ! total discrete gamma-ray cross section
! Variables for energy grid
!   Einc           ! incident energy in MeV
! Variables for incident channel
!   xspopnuc       ! population cross section per nucleus
! Variables for nuclides
!   AA             ! mass number of residual nucleus
!   ZZ             ! charge number of residual nucleus
! Constants
!   cparity        ! parity (character)
!   nuc            ! symbol of nucleus
!   parsym         ! symbol of particle
! Variables for levels
!   edis           ! energy of level
!   jdis           ! spin of level
!   parlev         ! parity of level
!
! *** Declaration of local data
!
  implicit none
  character(len=2) :: levelstring1   ! level string
  character(len=2) :: levelstring2   ! level string
  character(len=3)  :: massstring !
  character(len=6)  :: finalnuclide !
  character(len=18) :: reaction   ! reaction
  character(len=132) :: topline    ! topline
  character(len=15) :: col(5)     ! header
  character(len=15) :: un(5)     ! header
  character(len=80) :: quantity   ! quantity
  character(len=19) :: gamfile    ! giant resonance parameter file
  integer           :: A          ! mass number of target nucleus
  integer           :: i1         ! value
  integer           :: i2         ! value
  integer           :: indent
  integer           :: id2
  integer           :: id4
  integer           :: id6
  integer           :: Ncol       ! number of columns
  integer           :: Ncomp      ! neutron number index for compound nucleus
  integer           :: nen        ! energy counter
  integer           :: Z          ! charge number of target nucleus
  integer           :: Zcomp      ! proton number index for compound nucleus
  real(sgl)         :: Egam       ! gamma energy
!
! **************** Discrete gamma-ray intensities **********************
!
  indent = 0
  id2 = indent + 2
  id4 = indent + 4
  id6 = indent + 6
  col(1)='E'
  un(1)='MeV'
  col(2)='xs'
  un(2)='mb'
  Ncol=2
  quantity='cross section'
  reaction='('//parsym(k0)//',x)'
  write(*, '(/" 10. Gamma-ray intensities")')
  do Zcomp = 0, maxZ
    do Ncomp = 0, maxN
      if (xspopnuc(Zcomp, Ncomp) < popeps) cycle
      Z = ZZ(Zcomp, Ncomp, 0)
      A = AA(Zcomp, Ncomp, 0)
      write(*, '(/" Nuclide: ", i3, a2/)') A, nuc(Z)
      write(*, '("     Initial level          Final level", "     Gamma Energy  Cross section "/)')
      write(*, '("  no.  J/Pi    Ex         no.  J/Pi    Ex"/)')
      do i1 = 0, numlev
        do i2 = 0, i1
          if (xsgamdis(Zcomp, Ncomp, i1, i2) == 0.) cycle
          Egam = edis(Zcomp, Ncomp, i1) - edis(Zcomp, Ncomp, i2)
          write(*, '(1x, i3, 2x, f4.1, a1, f8.4, "  --->", i3, 2x, f4.1, a1, f8.4, f11.5, es15.6)') i1, jdis(Zcomp, Ncomp, i1), &
 &          cparity(parlev(Zcomp, Ncomp, i1)), edis(Zcomp, Ncomp, i1), i2, jdis(Zcomp, Ncomp, i2), &
 &          cparity(parlev(Zcomp, Ncomp, i2)), edis(Zcomp, Ncomp, i2), Egam, xsgamdis(Zcomp, Ncomp, i1, i2)
        enddo
      enddo
      write(*, '(/"  Total", 47x, es15.6)') xsgamdistot(Zcomp, Ncomp)
    enddo
  enddo
!
! Write results on separate files
!
  if (filegamdis) then
    do Zcomp = 0, maxZ
      do Ncomp = 0, maxN
        Z = ZZ(Zcomp, Ncomp, 0)
        A = AA(Zcomp, Ncomp, 0)
        massstring='   '
        write(massstring,'(i3)') A
        finalnuclide=trim(nuc(Z))//trim(adjustl(massstring))
        do i1 = 0, numlev
          do i2 = 0, i1
            if (xsgamdis(Zcomp, Ncomp, i1, i2) == 0. .and. .not. gamexist(Zcomp, Ncomp, i1, i2)) cycle
            Egam = edis(Zcomp, Ncomp, i1) - edis(Zcomp, Ncomp, i2)
            gamfile = 'gam000000L00L00.tot'
            write(gamfile(4:9), '(2i3.3)') Z,A
            write(gamfile(11:12), '(i2.2)') i1
            write(gamfile(14:15), '(i2.2)') i2
            if ( .not. gamexist(Zcomp, Ncomp, i1, i2)) then
              gamexist(Zcomp, Ncomp, i1, i2) = .true.
              levelstring1='  '
              write(levelstring1,'(i2)') i1
              levelstring2='  '
              write(levelstring2,'(i2)') i2
              reaction='('//parsym(k0)//',xg_'//trim(adjustl(levelstring1))//'-'//trim(adjustl(levelstring2))//')'
              open (unit = 1, file = gamfile, status = 'unknown')
              topline=trim(targetnuclide)//trim(reaction)//trim(finalnuclide)//' gamma-ray production '//trim(quantity)
              call write_header(indent,topline,source,user,date,oformat)
              call write_target(indent)
              call write_reaction(indent,reaction,Qres(Zcomp, Ncomp, i1),Ethresh(Zcomp, Ncomp, i1),0,0)
              call write_level(id2,-1,levnum(Zcomp, Ncomp, i1),edis(Zcomp, Ncomp, i1), &
 &              jdis(Zcomp, Ncomp, i1),parlev(Zcomp, Ncomp, i1),0.)
              call write_level(id4,-1,levnum(Zcomp, Ncomp, i2),edis(Zcomp, Ncomp, i2), &
 &              jdis(Zcomp, Ncomp, i2),parlev(Zcomp, Ncomp, i2),0.)
              call write_real(id2,'gamma energy [MeV]',Egam)
              call write_residual(id2,Z,A,finalnuclide)
              call write_quantity(id2,quantity)
              call write_datablock(id2,Ncol,Ninc,col,un)
              do nen = 1, Ninclow
                write(1, '(2es15.6)') eninc(nen), 0.
              enddo
              do nen = Ninclow + 1, nin - 1
                write(1, '(2es15.6)') eninc(nen), 0.
              enddo
            else
              open (unit = 1, file = gamfile, status = 'old', position = 'append')
            endif
            write(1, '(2es15.6)') Einc, xsgamdis(Zcomp, Ncomp, i1, i2)
            close (unit = 1)
          enddo
        enddo
      enddo
    enddo
  endif
  return
end subroutine gamdisout
! Copyright A.J. Koning 2023
