subroutine residualout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of residual production cross sections
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
! Variables for energies
!   Ninclow       ! number of incident energies below Elow
! Variables for numerics
!   maxN            ! maximal number of neutrons away from initial compound nucleus
!   maxZ            ! maximal number of protons away from initial compound nucleus
!   xseps           ! limit for cross sections
! Variables for output
!   fileresidual    ! flag for residual production cross sections on separate file
!   flagcompo       ! flag for output of cross section components
! Variables for fission
!   flagfission     ! flag for fission
! Variables for input energies
!   eninc           ! incident energy in MeV
!   flaginitpop     ! flag for initial population distribution
!   nin             ! counter for incident energy
!   Ninc            ! number of incident energies
! Variables for main input
!   Ainit           ! mass number of initial compound nucleus
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
! Variables for existence libraries
!   rpexist         ! flag for existence of residual production cross section
!   rpisoexist      ! flag for existence of isomeric residual production cross section
! Variables for total cross sections
!   xsfistot        ! total fission cross section
! Variables for energy grid
!   Einc            ! incident energy in MeV
! Variables for energies
!   Ethresh         ! threshold incident energy for residual nucleus
!   Qres            ! Q - value for residual nucleus
! Variables for multiple emission
!   xsinitpop       ! initial population cross section
!   xspopcomp       ! compound population cross section per nucleus
!   xspoppreeq      ! preequilibrium population cross section per nucleus
! Variables for binary reactions
!   xsnonel         ! non - elastic cross
!   xspopdir        ! direct population cross section per nucleus
! Variables for incident channel
!   maxA            ! maximal number of nucleons away from initial compound nucleus
!   xsbranch        ! branching ratio for isomeric cross section
!   xsmassprod      ! residual production cross section per mass unit
!   xspopex         ! population cross section summed over spin and parity
!   xspopnuc        ! population cross section per nucleus
!   xsresprod       ! total residual production ( = reaction) cross section
! Variables for nuclides
!   AA              ! mass number of residual nucleus
!   ZZ              ! charge number of residual nucleus
! Constants
!   iso             ! counter for isotope
!   natstring       ! string extension for file names
!   nuc             ! symbol of nucleus
!   parN            ! neutron number of particle
!   parsym          ! symbol of particle
!   parZ            ! charge number of particle
! Variables for levels
!   edis            ! energy of level
!   tau             ! lifetime of state in seconds
! Variables for level density
!   Nlast           ! last discrete level
! Variables for thermal cross sections
!   fxsbranch       ! branching ratio for isomeric cross section
!   fxspopex        ! population cross section summed over spin and par
!   fxspopnuc       ! population cross section per nucleus
! Variables for masses
!   nucmass         ! mass of nucleus
!
! *** Declaration of local data
!
  implicit none
  logical           :: tauflag
  character(len=3)  :: massstring ! 
  character(len=6)  :: finalnuclide ! 
  character(len=18) :: reaction   ! reaction
  character(len=16) :: isofile    ! file with isomeric cross section
  character(len=16) :: rpfile     ! file with residual production cross sections
  character(len=132) :: topline    ! topline
  character(len=15) :: col(5)     ! header
  character(len=15) :: un(5)     ! header
  character(len=80) :: quantity   ! quantity
  integer           :: A          ! mass number of target nucleus
  integer           :: Acomp      ! mass number index for compound nucleus
  integer           :: Ares       ! mass number of residual nucleus
  integer           :: kiso       ! counter for isomer
  integer           :: Ncomp      ! neutron number index for compound nucleus
  integer           :: Ncol       ! number of columns
  integer           :: nen        ! energy counter
  integer           :: nex        ! excitation energy bin of compound nucleus
  integer           :: Z          ! charge number of target nucleus
  integer           :: Zcomp      ! proton number index for compound nucleus
!
! *************** Residual production cross sections *******************
!
  un = 'mb'
  col(1)='E'
  un(1)='MeV'
  col(2)='xs'
  Ncol=2
  col(3)='Direct'
  col(4)='Preequilibrium'
  col(5)='Compound'
  quantity='cross section'
  reaction='('//parsym(k0)//',x)'
  write(*, '(/" 4. Residual production cross sections"/)')
  write(*, '("   a. Per isotope"/)')
  write(*, '("   Z   A  nuclide    total     level   ", "isomeric    isomeric    lifetime")')
  write(*, '(17x, "cross section", 7x, "cross section  ratio"/)')
  do Acomp = 0, maxA
    do Zcomp = 0, maxZ
      Ncomp = Acomp - Zcomp
      if (Ncomp < 0 .or. Ncomp > maxN) cycle
      if (xspopnuc(Zcomp, Ncomp) < xseps) cycle
!
! A. Total and ground state
!
      Z = ZZ(Zcomp, Ncomp, 0)
      A = AA(Zcomp, Ncomp, 0)
      write(*, '(1x, 2i4, " (", i3, a2, ")", es12.5, "    0   ", es12.5, f9.5)') Z, A, A, nuc(Z), xspopnuc(Zcomp, Ncomp), &
 &      xspopex(Zcomp, Ncomp, 0), xsbranch(Zcomp, Ncomp, 0)
!
! B. Per isomer
!
      do nex = 1, Nlast(Zcomp, Ncomp, 0)
        if (tau(Zcomp, Ncomp, nex) /= 0.) then
          write(*, '(31x, i3, 3x, es12.5, f9.5, 2x, es12.5, " sec. ")')  levnum(Zcomp, Ncomp, nex), xspopex(Zcomp, Ncomp, nex), &
 &          xsbranch(Zcomp, Ncomp, nex), tau(Zcomp, Ncomp, nex)
        endif
      enddo
    enddo
  enddo
  write(*, '(/"   b. Per mass"/)')
  write(*, '("   A  cross section"/)')
  do Acomp = 0, maxA
    if (xsmassprod(Acomp) > xseps) then
      Ares = Ainit - Acomp
      write(*, '(1x, i4, es12.5)') Ares, xsmassprod(Acomp)
    endif
  enddo
!
! ************* Check of residual production cross section *************
!
  write(*, '(/" Total residual production cross section:", f14.7)') xsresprod
  if (flagfission) then
    write(*, '(" Total fission cross section            :", f14.7)') xsfistot
    write(*, '(" Fission + res. production cross section:", f14.7)') xsresprod + xsfistot
  endif
  if (flaginitpop) then
    write(*, '(" Initial population cross section       :", f14.7)') xsinitpop
  else
    write(*, '(" Non-elastic cross section              :", f14.7)') xsnonel
  endif
!
! Write results to separate file
!
  if (fileresidual) then
    do Acomp = 0, maxA
      do Zcomp = 0, maxZ
        Ncomp = Acomp - Zcomp
        if (Ncomp < 0 .or. Ncomp > maxN) cycle
        if (xspopnuc(Zcomp, Ncomp) < xseps .and. .not. rpexist(Zcomp, Ncomp)) cycle
!
! A. Total
!
        if (Zcomp == parZ(k0) .and. Ncomp == parN(k0)) then
          nex = 1
        else
          nex = 0
        endif
        Z = ZZ(Zcomp, Ncomp, 0)
        A = AA(Zcomp, Ncomp, 0)
        rpfile = 'rp000000.tot'//natstring(iso)
        write(rpfile(3:8), '(2i3.3)') Z, A
        massstring='   '
        write(massstring,'(i3)') A
        finalnuclide=trim(nuc(Z))//adjustl(massstring)
        if ( .not. rpexist(Zcomp, Ncomp)) then
          rpexist(Zcomp, Ncomp) = .true.
          open (unit = 1, file = rpfile, status = 'unknown')
          topline=trim(targetnuclide)//trim(reaction)//trim(finalnuclide)//' '//trim(quantity)
          if (flagcompo) then
            col(3)='Direct'
            Ncol=5
          else
            Ncol=2
          endif
          call write_header(topline,source,user,date,oformat)
          call write_target
          call write_reaction(reaction,Qres(Zcomp, Ncomp, nex),Ethresh(Zcomp, Ncomp, nex),6,5)
          call write_residual(Z,A,finalnuclide)
          call write_double(2,'mass [amu]',nucmass(Zcomp, Ncomp))
          call write_datablock(quantity,Ncol,Ninc,col,un)
          if (flagcompo) then
            do nen = 1, Ninclow
              write(1, '(5es15.6)') eninc(nen), fxspopnuc(nen, Zcomp, Ncomp), 0., 0., fxspopnuc(nen, Zcomp, Ncomp)
            enddo
            do nen = Ninclow + 1, nin - 1
              write(1, '(5es15.6)') eninc(nen), 0., 0., 0., 0.
            enddo
          else
            do nen = 1, Ninclow
              write(1, '(2es15.6)') eninc(nen), fxspopnuc(nen, Zcomp, Ncomp)
            enddo
            do nen = Ninclow + 1, nin - 1
              write(1, '(2es15.6)') eninc(nen), 0.
            enddo
          endif
        else
          open (unit = 1, file = rpfile, status = 'old', position = 'append')
        endif
        if (flagcompo) then
          write(1, '(5es15.6)') Einc, xspopnuc(Zcomp, Ncomp), xspopdir(Zcomp, Ncomp), &
 &          xspoppreeq(Zcomp, Ncomp), xspopcomp(Zcomp, Ncomp)
        else
          write(1, '(2es15.6)') Einc, xspopnuc(Zcomp, Ncomp)
        endif
        close (unit = 1)
!
! B. Per ground state and isomer
!
        tauflag = .false.
        do nex=1,Nlast(Zcomp,Ncomp,0)
          if (tau(Zcomp,Ncomp,nex).ne.0.)  tauflag = .true.
        enddo
        if (tauflag) then
          do nex = 0, Nlast(Zcomp, Ncomp, 0)
            if (nex == 0 .or. tau(Zcomp, Ncomp, nex) /= 0.) then
              isofile = 'rp000000.L00'//natstring(iso)
              write(isofile(3:8), '(2i3.3)') Z, A
              write(isofile(11:12), '(i2.2)') levnum(Zcomp, Ncomp, nex)
              if ( .not. rpisoexist(Zcomp, Ncomp, nex)) then
                rpisoexist(Zcomp, Ncomp, nex) = .true.
                open (unit = 1, file = isofile, status = 'replace')
                if (nex == 0) then
                  kiso = 0
!                 write(1, '("# ", a1, " + ", a, ": Production of ", a, " - Ground state")') &
!&                  parsym(k0), trim(targetnuclide), trim(finalnuclide)
                else
                  kiso = kiso + 1
!                 write(1, '("# ", a1, " + ", a, ": Production of ", a, " - Level", i3, f12.5, " MeV")') &
!&                  parsym(k0), trim(targetnuclide), trim(finalnuclide), levnum(Zcomp, Ncomp, nex), edis(Zcomp, Ncomp, nex)
                endif
                finalnuclide=trim(nuc(Z))//trim(adjustl(massstring))//isochar(kiso)
                topline=trim(targetnuclide)//trim(reaction)//trim(finalnuclide)//' '//trim(quantity)
                col(3)='Isomeric_ratio'
                un(3)=''
                Ncol=3
                call write_header(topline,source,user,date,oformat)
                call write_target
                call write_reaction(reaction,Qres(Zcomp, Ncomp, nex),Ethresh(Zcomp, Ncomp, nex),6,5)
                call write_residual(Z,A,finalnuclide)
                call write_double(2,'mass [amu]',nucmass(Zcomp, Ncomp))
                call write_level(2,kiso,levnum(Zcomp, Ncomp, nex),edis(Zcomp, Ncomp, nex), &
 &                jdis(Zcomp, Ncomp, nex),parlev(Zcomp, Ncomp, nex),tau(Zcomp, Ncomp, nex))
                call write_datablock(quantity,Ncol,Ninc,col,un)
                do nen = 1, Ninclow
                  write(1, '(3es15.6)') eninc(nen), fxspopex(nen, Zcomp, Ncomp, nex), fxsbranch(nen, Zcomp, Ncomp, nex)
                enddo
                do nen = Ninclow + 1, nin - 1
                  write(1, '(3es15.6)') eninc(nen), 0., fxsbranch(max(Ninclow, 1), Zcomp, Ncomp, nex)
                enddo
              else
                open (unit = 1, file = isofile, status = 'old', position = 'append')
              endif
              write(1, '(3es15.6)') Einc, xspopex(Zcomp, Ncomp, nex), xsbranch(Zcomp, Ncomp, nex)
              close (unit = 1)
            endif
          enddo
        endif
      enddo
    enddo
  endif
  return
end subroutine residualout
! Copyright A.J. Koning 2021
